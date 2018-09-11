
# coding: utf-8

# In[1]:


from astropy.io import ascii
import numpy as np
import pandas as pd
import scipy
import matplotlib.pyplot as plt


# In[2]:


#reads/displays ascii formatted txt file and names the data "spectra"
spectra=ascii.read("../lab0_spectral_data.txt")
#print(spectra)


# In[3]:


plt.figure()
plt.semilogy(spectra['Am-241'])
plt.semilogy(spectra['Cs-137'])
plt.title('$^{241}Am$ and $^{137}Cs$ Spectra by Channel')
plt.xlabel('Channel')
plt.ylabel('Counts')
plt.legend()
plt.savefig('../images/SpecUncal')
#plt.show()


# In[4]:


#Peak finding script - from https://gist.github.com/endolith/250860
import sys
from numpy import NaN, Inf, arange, isscalar, asarray, array

def peakdet(v, delta, x = None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html

    Returns two arrays

    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.

    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.

    """
    maxtab = []
    mintab = []

    if x is None:
        x = arange(len(v))

    v = asarray(v)

    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')

    if not isscalar(delta):
        sys.exit('Input argument delta must be a scalar')

    if delta <= 0:
        sys.exit('Input argument delta must be positive')

    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN

    lookformax = True

    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]

        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return array(maxtab), array(mintab)


# In[5]:


#Finds and plots peaks for Cs-137 and Am-241
AmPeak, mintab = peakdet(spectra['Am-241'],39000)
CsPeak, mintab = peakdet(spectra['Cs-137'],39000)
plt.figure()
plt.semilogy(spectra['Am-241'])
plt.semilogy(spectra['Cs-137'])
plt.scatter(array(AmPeak)[:,0], array(AmPeak)[:,1], color='blue')
plt.scatter(array(CsPeak)[:,0], array(CsPeak)[:,1], color='orange')
plt.title('$^{241}Am$ and $^{137}Cs$ Identified Peaks by Channel')
plt.xlabel('Channel')
plt.ylabel('Counts')
plt.legend()
plt.savefig('../images/Peaks')
#plt.show()
#print(AmPeak)
#print(CsPeak)
x1 = float(AmPeak[:,0])
x2 = float(CsPeak[:,0])


# In[6]:


# 2 Point Linear Energy Calibration using 59.541 keV peak from Am-241 and 661.657 keV peak from Cs-137
peakE = [59.541, 661.657]
channel = [x1, x2]
slope, intercept = np.polyfit(channel, peakE, 1)
#print(slope, intercept)



# In[7]:


#Calibrated Energy Range
Energy = range(8192)*slope + intercept


#Plots Linear fit of Energy
plt.figure()
plt.plot(range(8192),Energy)
plt.title('Linear Fit')
plt.xlabel('Channel')
plt.ylabel('Energy (keV)')
plt.savefig('../images/LinearFit')
# In[8]:


#Plots Ba-133 Spectra vs Calibrated Energy
plt.figure()
plt.semilogy(Energy,spectra['Ba-133'])
plt.title('$^{133}Ba$ Calibrated Energy Spectrum')
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.legend()
plt.savefig('../images/BaSpecCal')
#plt.show()


# In[9]:


#Finds Ba Peaks
BaPeak, mintab = peakdet(spectra['Ba-133'],5000)
plt.figure()
plt.semilogy(spectra['Ba-133'])
plt.scatter(array(BaPeak)[:,0], array(BaPeak)[:,1], color='blue')
plt.title('$^{133}Ba$ Identified Peaks by Channel')
plt.xlabel('Channel')
plt.ylabel('Counts')
plt.legend()
plt.savefig('../images/BaPeaks')
#plt.show()
#print(BaPeak)


# In[10]:


E1 = (float(BaPeak[2,0]))*slope + intercept
E2 = (float(BaPeak[3,0]))*slope + intercept
E3 = (float(BaPeak[4,0]))*slope + intercept
E4 = (float(BaPeak[5,0]))*slope + intercept
E5 = (float(BaPeak[6,0]))*slope + intercept
#print(E1, E2, E3, E4, E5)


# In[11]:


BaCalibrated = np.asarray([E1, E2, E3, E4, E5])
BaExpected = np.asarray([80.9979, 276.3989, 302.8508, 356.0129, 383.8485])
#print(BaCalibrated, BaExpected)


CalVsEx = {'Expected Energy (keV)': BaExpected,
            'Calibrated Energy (keV)': BaCalibrated}

ascii.write(CalVsEx, output='../text/Comparison.tex', overwrite=True, Writer=ascii.Latex, names=['Expected Energy (keV)','Calibrated Energy (keV)'], col_align='|lr|',
            latexdict={'preamble': r'\begin{center}',
                       'tablefoot': r'\end{center}',
                       'tabletype': 'table*'})

# In[12]:


#Calculates Absolute, Relative, and Percent Errors
abs_error = abs(BaExpected-BaCalibrated)
rel_error = abs_error/abs(BaExpected)
percent_error = 100*rel_error
#print(abs_error)
#print(rel_error)
#print(percent_error)


# In[13]:


dat = [abs_error, rel_error, percent_error]

columns = ('80.9979 keV', '276.3989 keV', '302.8508 keV', '356.0129 keV', '383.8485 keV')
rows = ('Absolute Error', 'Relative Error', 'Percent Error')

fig, ax = plt.subplots()

# Hide axes
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
ax.set_frame_on(False)


ax.table(cellText=dat, rowLabels=rows, colLabels=columns, loc='center')
plt.savefig('../images/ErrorAnalysis')
#plt.show()


# In[14]:


data = {'\gamma Energy (keV)': [80.9979, 276.3989, 302.8508, 356.0129, 383.8485],
            '\epsilon': abs_error,
            '\eta': rel_error,
            '\delta': percent_error}

ascii.write(data, output='../text/ErrorAnalysis.tex', overwrite=True, Writer=ascii.Latex, names=['\gamma Energy (keV)','\epsilon', '\eta', '\delta'], col_align='|lr|',
            latexdict={'preamble': r'\begin{center}',
                       'tablefoot': r'\end{center}',
                       'tabletype': 'table*'})
