#!/usr/bin/env python

import numpy as np

#TODO improve names (make consistent)
from propagation import phase_thereback, propagation_phaseloss, calc_phase_difference
from propagation import unambigous_distance

C = 3e8 # m/s

def unambigous_speed(f0, c, PRF):
    '''Unambigeous speed for complex (I/Q) signal'''
    return unambigous_distance(f0, c)*PRF

def basebandPulse(t, phi, complex):
    ''''''
    if complex:
        base = np.exp(1j*(2*np.pi*f0*t+phi)) # complex
    else:
        base = np.cos(2*np.pi*f0*t+phi) # real
    return base

def phase2distance(phase, c, f0):
    ''''''
    return c*phase/(2*np.pi*f0)

def phase2distance_thereback(phase, c, f0):
    ''''''
    return phase2distance(phase, c, f0)/2

def min_speed(PRF, Npulses, c, f0, phase_threshold):
    '''Not sure how to fill this out'''
    CPI = Npulses/PRF
    dist = phase2distance(phase_threshold, c , f0)
    return dist/CPI


# simulation Parameters
COMPLEX = True
fmax = 100e9
fsample = 2*fmax
dt_sample = 1/fsample
Tmax = 50e-9 # length of the pulse
Nsample = int(fsample*Tmax) # Nyquist

# time vector
t = np.linspace(0,Tmax,Nsample) ## ISSUE: real dt and real sample freq differ due to this!

# baseband signal
f0 = 300e6
phi = np.pi/3

PRF = 1.34333e3 # Hz
Npulses = 3
t_start = 0 # seconds
pulse_starts = []

print(f"Unambigeous distance {unambigous_speed(f0, C, PRF)} m/s")

PRI = 1/PRF
for i in range(Npulses):
    pulse_starts.append(i*PRI+t_start)

pulses = [basebandPulse(t+offset, phi, COMPLEX) for offset in pulse_starts]

def targetDist_t(x0, v, t):
    return x0 + v*t

# target
x0 = 17.2452
v = 671.665
v = 17.1
target_locations = [targetDist_t(x0,v,t) for t in pulse_starts]

return_signals = [propagation_phaseloss(target_locations[i], C, f0)*pulses[i] for i in range(Npulses)]

import matplotlib.pyplot as pl
fig, ax = pl.subplots(2,1)
for sig in pulses:
    ax[0].plot(t,np.real(sig))

for sig in return_signals:
    ax[1].plot(t,np.real(sig))

# plot each transmitted and normalized recieved pulse
fig, ax = pl.subplots(Npulses,1)
for i in range(Npulses):
    ax[i].plot(t,np.real(pulses[i]))
    ax[i].plot(t,np.real(return_signals[i])/np.real(return_signals[i]).max(),'--')

return_phases = [phase_thereback(target_locations[i], C, f0)%(2*np.pi) for i in range(Npulses)]
est_return_phases = [calc_phase_difference(pulses[i], return_signals[i], f0, t[1]-t[0])
                     for i in range(Npulses)]

print(f"{return_phases=}")
print(f"{est_return_phases=}")

print("")
print(f"ret diff{[(j-i)%(2*np.pi) for i, j in zip(return_phases[:-1], return_phases[1:])] }")
print(f"est diff{[(j-i)%(2*np.pi) for i, j in zip(est_return_phases[:-1], est_return_phases[1:])] }")

# I thought there would be a way to get the phase difference from first to last pulse to get the minium detectable speed, but that does not seem to be the case from above.
#min_speed(PRF, Npulses, C, f0, 0.24)
