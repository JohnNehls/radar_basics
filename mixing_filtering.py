#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl
from scipy import fft
from scipy import  signal

# simulation Parameters
COMPLEX = True
fmax = 2e9
fsample = 2*fmax
dt_sample = 1/fsample
Tmax = 200e-9 # length of the pulse
Nsample = int(fsample*Tmax) # Nyquist

# time vector
t = np.linspace(0,Tmax,Nsample) ## ISSUE: real dt and real sample freq differ due to this!

# baseband signal
f0 = 300e6
phi = np.pi/3

# result below work with real and complex baseband signal
if COMPLEX:
    base = np.exp(1j*(2*np.pi*f0*t+phi)) # complex
else:
    base = np.cos(2*np.pi*f0*t+phi) # real

# local oscillator (LO)
f_lo = 1e9
phi = np.pi
LO = np.cos(2*np.pi*f_lo*t+phi)

# rf/RF
rf = LO*base
RF = fft.fftshift(fft.fft(rf))/Nsample

# band-pass filter
f1, f2 = f_lo+100e6, f_lo+500e6
numtaps = 21
bp_coef = signal.firwin(numtaps, [f1, f2], pass_zero=False, fs=fsample)

# bandpassed rf/RF
# Fastest filter application. https://scipy.github.io/old-wiki/pages/Cookbook/ApplyFIRFilter.html
rf_bp = signal.fftconvolve(rf, bp_coef, mode="same")
RF_bp = fft.fftshift(fft.fft(rf_bp))/Nsample


# plot the rf vs rf_bp
fig, ax = pl.subplots(2,1)
ax[0].set_title("filter RF")
ax[0].plot(t,np.real(rf),'-k',label='rf')
ax[0].plot(t,np.real(rf_bp),'--r',label='rf_bp')
ax[0].set_xlabel("time [s]")
ax[0].set_ylabel("voltage")
ax[0].legend()

k = np.linspace(-1/(2*dt_sample), 1/(2*dt_sample), Nsample)  # frequency
ax[1].plot(k,np.real(RF),'-k', label='rf')
ax[1].plot(k,np.real(RF_bp),'--r', label='rf_bp')
ax[1].set_xlabel("frequency [Hz]")
ax[1].set_ylabel("voltage")
ax[1].legend()


# Convert back to baseband
rf_bplo = LO*rf_bp
RF_bplo = fft.fftshift(fft.fft(rf_bplo))/Nsample

# filter
f1, f2 = 100e6, 500e6
numtaps = 21
bp_coef = signal.firwin(numtaps, [f1, f2], pass_zero=False, fs=fsample)
rf_bplobp = signal.fftconvolve(rf_bplo, bp_coef, mode="same")
RF_bplobp = fft.fftshift(fft.fft(rf_bplobp))/Nsample

# compare to original baseband
fig, ax = pl.subplots(2,1)
ax[0].set_title("filter baseband")
ax[0].plot(t,np.real(base),'-k',label='base')
ax[0].plot(t,np.real(rf_bplobp),'--r',label='rf_bplobp')
ax[0].set_xlabel("time [s]")
ax[0].set_ylabel("voltage")
ax[0].legend()

k = np.linspace(-1/(2*dt_sample), 1/(2*dt_sample), Nsample)  # frequency
ax[1].plot(k,np.real(RF_bplo),'-k', label='rf_bplo')
ax[1].plot(k,np.real(RF_bplobp),'--r', label='rf_bplobp')
ax[1].set_xlabel("frequency [Hz]")
ax[1].set_ylabel("voltage")
ax[1].legend()
