import numpy as np
import scipy as sp
import os
import matplotlib.pyplot as plt
import scipy.signal as sig
import pickle as pkl
import copy
import scipy.optimize as opt
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

def genLorentzian(pulse_width,order):
    #returns a lorentzian pulse that integrates to 1 and then is made zero mean
    
    gam = np.sqrt(pulse_width/2.*(0.5**(1./order))/(1-0.5**(1./order)));
    pulse = np.linspace(-pulse_width*0.5,pulse_width*0.5,int(pulse_width*1));
    pulse = 1/np.pi/gam*(gam**2/(pulse**2+gam**2))
    
    return (pulse**order-np.mean(pulse**order))/sum(pulse**order)
    
def optimize_sum(center,wavelet):
    wavelet_copy = wavelet;
    wavelet_copy[len(wavelet_copy)/2] = center
    return (sum(wavelet_copy)**2)

def ricker_wavelet(pulse_width,order):
    f = 1./pulse_width
    t = np.linspace(1*pulse_width,-1*pulse_width,int(pulse_width*2+1));
    A = (1-2*np.pi**2*f**2*t**2)*np.exp(-np.pi**2*f**2*t**2)
    A[len(A)/2] = opt.minimize(optimize_sum,A[len(A)/2],args=A)['x'][0]
    return A/sum(A**2)
    

#def customWavelet(pulse_width,relative_width):

def lowFatButter(cutoff,fs,order=5):
    nyq = 0.5*fs
    normalCutoff = cutoff/nyq
    b,a = sig.butter(order,normalCutoff,btype='low',analog=False);
    return b,a

def highFatButter(cutoff,fs,order=5):
    nyq = 0.5*fs
    normalCutoff = cutoff/nyq
    b,a = sig.butter(order,normalCutoff,btype='high',analog=False);
    return b,a
    
def closest(num,arr):
    return arr[np.argmin((arr-num)**2)]

def chi_sq(line1,line2):
    return sum((line1-line2)**2/len(line1))

def get_cali_curve(lpdf):
    cali_curve = np.zeros(len(lpdf[0]));
    std_curve = np.zeros(len(lpdf[0]));
    for i in range(len(lpdf[0])):
        cali_curve[i] = np.mean(lpdf[:,i])
        std_curve[i] = np.std(lpdf[:,i])
    return cali_curve/np.max(cali_curve), std_curve/np.max(cali_curve)    

def discriminator(peak_pdf, cali_curve, amp_f, curve_f):
    #normalize cali_curve
    cali_curve = cali_curve/max(cali_curve)

    metric = np.zeros(len(peak_pdf))

    if max(peak_pdf)==0:
        return 0,0,0,0
    return [max(peak_pdf)**amp_f,sum(peak_pdf/max(peak_pdf)*cali_curve**curve_f)]#,sum((cali_curve-peak_pdf/sum(peak_pdf))**2),sum(peak_pdf);

def discriminator_complex(peak_pdf, cali_curve, amp_f, curve_f):
    cali_curve = cali_curve/max(cali_curve)
    metric = np.zeros(len(peak_pdf))

    if max(peak_pdf)==0:
        return 0,0,0,0
    return max(peak_pdf)**amp_f,sum(peak_pdf/max(peak_pdf)*cali_curve**curve_f),sum((cali_curve-peak_pdf/sum(peak_pdf))**2),sum(peak_pdf);
    #used to be sum(peak_pdf)

def simple_rule(amp,trend):
    if amp>80 and trend>0.7:
        return True
    else:
        return False

# purges peaks based on threshold set at a pulse width
def purge(pws,pdfs,spike_loc,pw,thresh):
    for i,val in enumerate(pws):
        if val == pw:
            j = 0;
            while j<len(pdfs):
                if pdfs[j,i] < thresh:
                    pdfs = np.delete(pdfs,j,axis=0);
                    spike_loc = np.delete(spike_loc,j);
                else:
                    j+=1;
    return pdfs,spike_loc

def std_thresh(discs,thresh=[1.4,1]):
    mean = np.mean(discs,axis=0);
    std = np.std(discs,axis=0);
    spike = np.zeros(len(discs));
    for i,val in enumerate(discs):
        if val[0]>mean[0]+std[0]*thresh[0] and val[1]>mean[1]+std[1]*thresh[1]:
            spike[i] = 1;
        else:
            spike[i] = 0
    return spike



def stat(dist1,dist2):
    return (np.mean(dist1)-np.mean(dist2))/(np.sqrt(np.std(dist1)**2/len(dist1)+np.sqrt(np.std(dist2)**2/len(dist2))))


def get_opt(lpdf,tpdf,amp_f,curve_f):
    ccurve = get_cali_curve(lpdf)
    lmet = []
    for i in lpdf:
        lmet.append(discriminator(i,ccurve,amp_f,curve_f))
    tmet = []
    for i in tpdf:
        tmet.append(discriminator(i,ccurve,amp_f,curve_f))
    return stat(lmet,tmet),lmet,tmet


def pdf(signal,pulse_widths,pw0,prelabeled_peaks=[],wavelet_func = ricker_wavelet):
    
    lpeak = prelabeled_peaks

    #signal2 = np.convolve(signal,genLorentzian(pw0,1),mode='same');
    signal2 = np.convolve(signal,wavelet_func(pw0,1),mode='same');
    
    
    pdf = np.zeros([len(prelabeled_peaks),len(pulse_widths)])    

    tpeak = sig.find_peaks(-signal2,prominence=40.)[0]#chagned from 10 to 50
    tpeakcopy = copy.deepcopy(tpeak)
    totalpdf = np.zeros([len(tpeak),len(pulse_widths)]);
    
    for i in range(len(lpeak)):
        np.delete(tpeak,np.argmin(np.abs(prelabeled_peaks[i]-tpeak)));
    
    labeled_peaks = np.zeros(len(lpeak));
    ptl = 10
    for i in range(len(lpeak)):
        loc = lpeak[i]
        proms = sig.peak_prominences(-signal2,sig.find_peaks(-signal2[loc-ptl:loc+ptl+1])[0]+(loc-ptl),wlen=int(1.5*pw0))[0];
        labeled_peaks[i] =sig.find_peaks(-signal2[loc-ptl:loc+ptl+1])[0][np.argmax(proms)]+loc-ptl
        lpeak[i] = labeled_peaks[i]
    
    not_det = []
    
    
    for p in range(len(pulse_widths)):
        
        pulse_width = pulse_widths[p] 
        
        #signal2 = np.convolve(signal,genLorentzian(pulse_width,1),mode='same');
        signal2 = np.convolve(signal,wavelet_func(pulse_width,1),mode='same');

        labeled_peaks = np.zeros(len(lpeak));
        
        ptl = 3; #peak tracking tolerance
        
        #updating peak locations
        for i in range(len(lpeak)):
            loc = lpeak[i]
            try:
                labeled_peaks[i] = loc-ptl+closest(ptl,sig.find_peaks(-signal2[loc-ptl:loc+ptl+1])[0])
                lpeak[i] = labeled_peaks[i]
            except:
                not_det.append(loc)
           
        for i in range(len(lpeak)):
            try:
                pdf[i,p] = sig.peak_prominences(-signal2,[int(lpeak[i])],wlen=int(pulse_width*1.5))[0][0]
            except:
                pdf[i,p] = 0
                
        total_peaks = np.zeros(len(tpeak))
        
        #updating peak locations
        for i in range(len(tpeak)):
            loc = tpeak[i]
            if not loc in not_det:
                if(len(sig.find_peaks(-signal2[loc-ptl:loc+ptl+1])[0])==0):
                    #if peak no longer detected within range
                    not_det.append(tpeak[i]);
                    totalpdf[i,p] = 0
                elif loc-ptl+closest(ptl,sig.find_peaks(-signal2[loc-ptl:loc+ptl+1])[0]) in lpeak:
                    #if peak is part of labeled peaks
                    not_det.append(tpeak[i]);
                    totalpdf[i,p] = 0
                else: 
                    total_peaks[i] = loc-ptl+closest(ptl,sig.find_peaks(-signal2[loc-ptl:loc+ptl+1])[0])
                    #udpate peak locations
                    tpeak[i] = total_peaks[i];
                    try:
                        totalpdf[i,p] = sig.peak_prominences(-signal2,[int(tpeak[i])],wlen=int(pulse_width*1.5))[0][0]
                    except:
                        totalpdf[i,p] = 0;
        
    return pulse_widths,pdf,totalpdf,lpeak,tpeakcopy
    
def pdf_complex(signal,pulse_widths,pw0,prelabeled_peaks=[]):
    
    lpeak = prelabeled_peaks
    
    freq = 1/(pw0*40e-6)
    i,q = sig.gausspulse(np.linspace(40e-6*-max(pulse_widths),40e-6*max(pulse_widths),2*max(pulse_widths)),bw=1,fc=freq,retquad=True)
    wavelet = i+1j*q
    signal2 = np.convolve(signal,wavelet,mode='same');

    pdf = np.zeros([len(prelabeled_peaks),len(pulse_widths)])    

    tpeak = sig.find_peaks(-signal2,prominence=10.)[0]
    tpeakcopy = copy.deepcopy(tpeak)
    totalpdf = np.zeros([len(tpeak),len(pulse_widths)]);
    
    for i in range(len(lpeak)):
        np.delete(tpeak,np.argmin(np.abs(prelabeled_peaks[i]-tpeak)));
    
    labeled_peaks = np.zeros(len(lpeak));
    ptl = 10
    for i in range(len(lpeak)):
        loc = lpeak[i]
        proms = sig.peak_prominences(-signal2,sig.find_peaks(-signal2[loc-ptl:loc+ptl+1])[0]+(loc-ptl),wlen=int(1.5*pw0))[0];
        labeled_peaks[i] =sig.find_peaks(-signal2[loc-ptl:loc+ptl+1])[0][np.argmax(proms)]+loc-ptl
        lpeak[i] = labeled_peaks[i]
    
    not_det = []
    
    
    for p in range(len(pulse_widths)):
        
        pulse_width = pulse_widths[p] 
        
        i,q = sig.gausspulse(np.linspace(40e-6*-pulse_width,40e-6*pulse_width,2*pulse_width),bw=1,fc=freq,retquad=True)
        wavelet = i+1j*q
        signal2 = np.convolve(signal,wavelet,mode='same');
        
        labeled_peaks = np.zeros(len(lpeak));
        
        ptl = 3; #peak tracking tolerance
        
        #updating peak locations
        for i in range(len(lpeak)):
            loc = lpeak[i]
            try:
                labeled_peaks[i] = loc-ptl+closest(ptl,sig.find_peaks(-signal2[loc-ptl:loc+ptl+1])[0])
                lpeak[i] = labeled_peaks[i]
            except:
                not_det.append(loc)
           
        print lpeak
        for i in range(len(lpeak)):
            try:
                pdf[i,p] = sig.peak_prominences(np.absolute(signal2),[int(lpeak[i])],wlen=int(pulse_width*1.5))[0][0]
            except:
                pdf[i,p] = 0
                
        total_peaks = np.zeros(len(tpeak))
        
        #updating peak locations
        for i in range(len(tpeak)):
            loc = tpeak[i]
            if not loc in not_det:
                if(len(sig.find_peaks(-signal2[loc-ptl:loc+ptl+1])[0])==0):
                    #if peak no longer detected within range
                    not_det.append(tpeak[i]);
                    totalpdf[i,p] = 0
                elif loc-ptl+closest(ptl,sig.find_peaks(-signal2[loc-ptl:loc+ptl+1])[0]) in lpeak:
                    #if peak is part of labeled peaks
                    not_det.append(tpeak[i]);
                    totalpdf[i,p] = 0
                else: 
                    total_peaks[i] = loc-ptl+closest(ptl,sig.find_peaks(-signal2[loc-ptl:loc+ptl+1])[0])
                    #udpate peak locations
                    tpeak[i] = total_peaks[i];
                    try:
                        totalpdf[i,p] = sig.peak_prominences(np.absolute(signal2),[int(tpeak[i])],wlen=int(pulse_width*1.5))[0][0]
                    except:
                        totalpdf[i,p] = 0;
        
    return pulse_widths,pdf,totalpdf,lpeak,tpeakcopy