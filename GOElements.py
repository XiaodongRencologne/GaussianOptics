# %%
import numpy as np
import GaussianOptics
from GaussianOptics import ThinLens, DrawBeamCountour
import transforms3d

from matplotlib import cm
import matplotlib.pyplot as plt
import pyvista as pv
pv.set_jupyter_backend('trame')#('static')#

# %%
def Mirror(GO_para,Angle,Scale=5):
    Angle=Angle/180*np.pi
    Rin=GO_para['Rin']
    Rout=GO_para['Rout']
    a=np.abs(Rin+Rout)/2
    c=np.sqrt(Rin**2+Rout**2-2*Rin*Rout*np.cos(Angle))/2
    #eccentricity
    e=c/a
    b=np.sqrt(np.abs(a**2-c**2))
    p=b**2/a
    
    if e>1.0:
        GO_para['Mtype']='Hyperb'
        phi0=np.arccos(((2*c)**2+Rin**2-Rout**2)/(2*Rin*2*c))
        print(phi0*180/np.pi)
        size=GO_para['Lambda']/np.pi/GO_para['win']*3
        phi_max=phi0+size
        phi_min=phi0-size
        P_max=np.pi-np.arctan(b/a)-1/180*np.pi
        if np.abs(phi_max)>=P_max:
            phi_max=P_max
        elif np.abs(phi_min)>=P_max:
            phi_min=-P_max
    else:
        GO_para['Mtype']='Ellip'
        phi0=np.pi-np.arccos(((2*c)**2+Rin**2-Rout**2)/(2*Rin*2*c))
        size=GO_para['Lambda']/np.pi/GO_para['win']*3
        phi_max=phi0+size
        phi_min=phi0-size

    t=np.linspace(phi_min,phi_max,101)
    r=p/(1+e*np.cos(t))
    plt.plot(r*np.cos(t),r*np.sin(t),'k-')
    plt.axis('equal')
    plt.show()
# %%
'''CHAMP informations'''
Lambda1=1.0162456203
win2=1.4
#win2=2.2
din2=33.73731665
f2=67.80690232
angle2=90
P2=ThinLens(win2, din2, f2, Lambda1)
Mirror(P2, angle2)

# %%
