#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as pl
from scipy.signal import correlate

C = 3e8 # m/s
PI = np.pi
COMPLEX = True # likely will stay in complex in the future

# simulation Parameters
fmax = 5000e6 # max baseband signal # default: 500e6 higher values  used for cleaner plots
if COMPLEX:
    fsample = fmax # Complex is 2 samples, real and imag
else:
    fsample = 2*fmax
dt_sample = 1/fsample
Tmax = 50e-9 # length of the pulse
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

# distance
# set it so it is pi out of phase: there and back takes T/2 seconds
distance_pi = 1/f0/4*C
dist_factor = 1.8
distance = dist_factor*distance_pi
print(f"For {dist_factor}pi phase shift at {f0=}, {distance=}")
# with I/Q we should be able to allow phase [0, 2*pi)

# phase calc functions
def distance2phase(distance, c, f0):
    return 2*np.pi*distance/c*f0

def there_and_back_phase(distance, c, f0):
    return distance2phase(2*distance, c, f0)

# propagation loss multiplication
def k_there(distance):
    # solve at which distance the factor is 1: np.sqrt(1/(4*PI)) = 0.28209479177387814
    return (4*PI*(distance+0.28209479177387814)**2)**(-1)

def k_there_and_back(distance):
    return k_there(distance)**2

# plot the loss vs distance
ds = np.linspace(0,10,1000)
pl.figure()
pl.title("There and back 4piR^2 loss")
pl.plot(ds, k_there_and_back(ds))
plt.xlabel("target distance [m]")
plt.ylabel("multiplcation factor")
plt.yscale('log')
pl.grid()

# plot the two signals
recievedSignal=base*k_there_and_back(distance)*np.exp(1j*there_and_back_phase(distance, C, f0))
pl.figure()
pl.plot(t, np.real(base), '-k')
pl.plot(t, np.real(recievedSignal), '--r')
pl.grid()

# find the phase difference between two signals
def calc_phase_difference(signal1, signal2, f0, t):
    '''Find a better way to do this without passing full time array'''
    xcorr = correlate(signal1, signal2)
    am = np.argmax(xcorr)
    return np.remainder(2*np.pi*t[am]*f0,2*np.pi)

# test signal with complex baseband: base2 = base*np.exp(1j*np.pi/2)
print(f"Estimated phase difference = {calc_phase_difference(base, recievedSignal, f0, t)}")
