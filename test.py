# %%
import numpy as np

def Gaussian_propagation2d(w0,z,Lambda):
    # confocal distance, in grasp it is denoted by 'b'
    z_c= np.pi*w0**2/Lambda
    # beam radius at distance of z
    w_z=w0*np.sqrt(1+(z/z_c)**2)
    # radius of wavefront curvature at z
    if z==0:
        R_z=999999999999999999999999999999
    else:
        R_z=z+z_c**2/z
    # beam parameter
    q=z+1j*z_c
    # phasor
    phi=np.arctan(z/z_c)

    # beamsize
    theta=Lambda/(np.pi*w0)

    def Gbeam(x,y):
        k=2*np.pi/Lambda
        E=np.sqrt(2/np.pi/w_z**2)*np.exp(-(x**2+y**2)/w_z**2-1j*k*z-1j*k/2*(x**2+y**2)/R_z+1j*phi)
        return E

    return z_c, w_z, R_z, q, phi, Gbeam

def ThinLens(win, din, Lambda, para='f', value=10):
    # para:  'f' or 'dout'
    z_c, w_z, Rin, q1, phi1, b1 =Gaussian_propagation2d(win,din,Lambda)
    if para=='f':
        f=value       
        M_max=f/z_c
        factor=(din/f-1)**2+(z_c/f)**2
        dout=f*(1+(din/f-1)/factor)
        wout=win/np.sqrt(factor)
        Rout=1/(1/Rin-1/f)
        M=wout/win
        print(M)
        print(M_max)
        print(1/np.sqrt(factor))
        print(dout)
        print(Rout)

    elif para=='dout':
        dout=value
        d=din+dout
        pass




# %%
ThinLens(1.4, 0, 1.0162456203, para='f', value=33.73731665)
# %%
