import numpy as np
import scipy as sp
import os
import matplotlib.pyplot as plt
import scipy.signal as sig
import pickle as pkl
import copy

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

def genLorentzian(pulse_width,order):
    #returns a lorentzian pulse that integrates to 1 and then is made zero mean
    
    gam = np.sqrt(pulse_width/2.*(0.5**(1./order))/(1-0.5**(1./order)));
    pulse = np.linspace(-pulse_width*0.5,pulse_width*0.5,int(pulse_width*1));
    pulse = 1/np.pi/gam*(gam**2/(pulse**2+gam**2))
    
    return (pulse**order-np.mean(pulse**order))/sum(pulse**order)

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

def nBiasedConvolve(data, conv):
    #takes data and convolutes with conv, while making sure the mean of both data and conv in the window is zero mean
    convData = []
    for i in range(0,len(data)-len(conv)):
        temp_data = data[i:i+len(conv)]
        std_data = (temp_data-np.mean(temp_data))/np.std(temp_data)
        temp_data= temp_data-np.mean(temp_data)
        unbiased_data_neg = temp_data#-max(temp_data)
        convData.append(sum(unbiased_data_neg*conv));
    return np.array(convData)
    
def closest(num,arr):
    return arr[np.argmin(arr-num)]

def peak_prom(data,peaks):
    prom = np.zeros(len(peaks));
    for i in range(len(peaks)):
        try:
            prom[i] = sig.peak_prominences(data,[int(peaks[i])])[0][0]
        except:
            prom[i] = 0
    return prom
    
def mean_lines(arr):
    mean_line = []
    for i in range(len(arr[0])):
        mean_line.append(np.mean(arr[:,i]))
    return np.array(mean_line)

def chi_sq(line1,line2):
    return sum((line1-line2)**2/len(line1))

def get_cali_curve(lpdf):
    cali_curve = np.zeros(len(lpdf[0]));
    for i in range(len(lpdf[0])):
        cali_curve[i] = np.mean(lpdf[:,i])
    return cali_curve/np.max(cali_curve)


def discriminator(peak_pdf, cali_curve, amp_f, curve_f):
    #normalize cali_curve
    cali_curve = cali_curve/max(cali_curve)

    metric = np.zeros(len(peak_pdf))

    if max(peak_pdf)==0:
        return 0,0
    return max(peak_pdf)**amp_f,sum(peak_pdf/sum(peak_pdf)*cali_curve**curve_f);

def simple_rule(amp,trend):
    if amp>80 and trend>0.7:
        return True
    else:
        return False

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


def pdf(signal,pulse_widths,pw0,prelabeled_peaks=[]):
    
    lpeak = prelabeled_peaks
    
    signal2 = nBiasedConvolve(signal,genLorentzian(pw0,1)-np.mean(genLorentzian(pw0,1)))
    signal2 = np.append(np.zeros(int(pw0*0.5)+1),signal2);

    pdf = np.zeros([len(prelabeled_peaks),len(pulse_widths)])    

    tpeak = sig.find_peaks(-signal2,prominence=10.)[0]
    tpeakcopy = copy.deepcopy(tpeak)
    totalpdf = np.zeros([len(tpeak),len(pulse_widths)]);
    
    for i in range(len(lpeak)):
        np.delete(tpeak,np.argmin(np.abs(prelabeled_peaks[i]-tpeak)));
    
    labeled_peaks = np.zeros(len(lpeak));
    for i in range(len(lpeak)):
        loc = lpeak[i]
        labeled_peaks[i] = loc-10+closest(10,sig.argrelextrema(-signal2[loc-10:loc+11],np.greater)[0])
        #update peak locations
        lpeak[i] = labeled_peaks[i]
    
    not_det = []
    
    
    for p in range(len(pulse_widths)):
    #pulse_width = 41.;
        pulse_width = pulse_widths[p] 
        
        signal2 = nBiasedConvolve(signal,genLorentzian(pulse_width,1)-np.mean(genLorentzian(pulse_width,1)))
        signal2=np.append(np.zeros(int(pulse_width*0.5)+1),signal2);
        
        labeled_peaks = np.zeros(len(lpeak));
        for i in range(len(lpeak)):
            loc = lpeak[i]
            try:
                labeled_peaks[i] = loc-5+closest(5,sig.argrelextrema(-signal2[loc-5:loc+6],np.greater)[0])
                lpeak[i] = labeled_peaks[i]
            except:
                not_det.append(loc)
            #update peak locations
           
        
        lprom = peak_prom(-signal2,labeled_peaks);
        print lprom
        for i in range(len(lprom)):
            pdf[i,p] = lprom[i];

        total_peaks = np.zeros(len(tpeak))
        
        for i in range(len(tpeak)):
            loc = tpeak[i]
            if not loc in not_det:
                if(len(sig.argrelextrema(-signal2[loc-5:loc+6],np.greater)[0])==0):
                    not_det.append(tpeak[i]);
                elif loc-5+closest(5,sig.argrelextrema(-signal2[loc-5:loc+6],np.greater)[0]) in lpeak:
                    not_det.append(tpeak[i]);
                else: 
                    total_peaks[i] = loc-5+closest(5,sig.argrelextrema(-signal2[loc-5:loc+6],np.greater)[0])
                    #udpate peak locations
                    tpeak[i] = total_peaks[i];
                  
        tprom = peak_prom(-signal2,total_peaks)
        
        for i in range(len(tprom)):
            totalpdf[i,p] = tprom[i];

        #pdf.append(metric2)
    return pulse_widths,pdf,totalpdf,lpeak,tpeakcopy