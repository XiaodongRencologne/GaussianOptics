import numpy as np


def Gaussian_propagation2d(w0,z,Lambda):
    # confocal distance, in grasp it is denoted by 'b'
    z_c= np.pi*w0**2/Lambda
    # beam radius at distance of z
    w_z=w0*np.sqrt(1+(z/z_c)**2)
    # radius of wavefront curvature at z
    R_z=z+z_c**2/z
    # beam parameter
    q=z+1j*z_c
    # phasor
    np.arctan(z/z_c)
