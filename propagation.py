#!/usr/bin/env python
import numpy as np
from scipy.signal import correlate

# phase calc functions
def distance2phase(distance, c, f0):
    ''''''
    return 2*np.pi*distance/c*f0

def phase_thereback(distance, c, f0):
    ''''''
    return distance2phase(2*distance, c, f0)

# propagation loss multiplication
def k_there(distance):
    ''''''
    # solve at which distance the factor is 1: np.sqrt(1/(4*PI)) = 0.28209479177387814
    return (4*np.pi*(distance+0.28209479177387814)**2)**(-1)

def k_thereback(distance):
    ''''''
    return k_there(distance)**2

def propagation_phaseloss(distance, c, f0):
    ''''''
    return k_thereback(distance)*np.exp(1j*phase_thereback(distance, c, f0))

def unambigous_distance(f0, C):
    '''Unambigeous distance for complex (I/Q) signal'''
    return 1/f0/2*C

# find the phase difference between two signals
def calc_phase_difference(signal1, signal2, f0, dt):
    '''better way?'''
    xcorr = correlate(signal1, signal2, mode='full', method='fft')
    am = np.argmax(xcorr)
    T = dt*am
    return np.remainder(2*np.pi*T*f0,2*np.pi)

def main():
    import matplotlib.pyplot as plt

    C = 3e8 # m/s
    COMPLEX = True # likely will stay in complex in the future
    # simulation Parameters
    fmax = 10e9 # max baseband signal # default: 500e6 higher values  used for cleaner plots
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
    dist_factor = 1.0
    distance = dist_factor*distance_pi
    print(f"For {dist_factor}pi phase shift at {f0=}, {distance=}")
    # with I/Q we should be able to allow phase [0, 2*pi)

    # plot the loss vs distance
    ds = np.linspace(0,10,1000)
    plt.figure()
    plt.title("There and back 4piR^2 loss")
    plt.plot(ds, k_thereback(ds))
    plt.xlabel("target distance [m]")
    plt.ylabel("multiplcation factor")
    plt.yscale('log')
    plt.grid()

    # plot the two signals
    # recievedSignal=base*k_thereback(distance)*np.exp(1j*phase_thereback(distance, C, f0))
    recievedSignal=base*propagation_phaseloss(distance, C, f0)
    plt.figure()
    plt.plot(t, np.real(base), '-k')
    plt.plot(t, np.real(recievedSignal), '--r')
    plt.grid()

    # test signal with complex baseband: base2 = base*np.exp(1j*np.pi/2)
    print(f"Estimated phase difference = {calc_phase_difference(base, recievedSignal, f0, t[1]-t[0])}")


if __name__ == "__main__":
    main()




