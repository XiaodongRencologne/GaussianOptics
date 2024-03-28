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
def Ellipsoid(GO_para,Angle, Scale=5):
    Angle=Angle/180*np.pi
    Rin=GO_para['Rin']
    Rout=GO_para['Rout']
    a=np.abs(Rin+Rout)/2
    c=np.sqrt(Rin**2+Rout**2-2*Rin*Rout*np.cos(Angle))/2
    #eccentricity
    e=c/a
    b=np.sqrt(np.abs(a**2-c**2))
    p=b**2/a
    
    phi0=np.arccos(((2*c)**2+Rin**2-Rout**2)/(2*Rin*2*c))
    size=GO_para['Lambda']/np.pi/GO_para['win']*3
    phi_max=phi0+size
    phi_min=phi0-size

    print(phi_max*180/np.pi)
    print(phi_min*180/np.pi)

    t=np.linspace(phi_min,phi_max,101)
    r=p/(1+e*np.cos(t))
    plt.plot(r*np.cos(t),r*np.sin(t),'k--')
    plt.axis('equal')




# %%
'''CHAMP informations'''
Lambda1=1.0162456203
win2=1.4
#win2=2.2
din2=33.73731665
f2=67.80690232
angle2=90
P2=ThinLens(win2, din2, f2, Lambda1)
Ellipsoid(P2, angle2)

# %%

def Mirror(focal_length,Rin,Rout,angle,Theta_beam,Size=5,Mirror_type='concave'):
    angle=angle*np.pi/180;
    if Mirror_type=='concave':
        a=np.abs(Rin+Rout)/2;
        b=np.sqrt(Rin*np.abs(Rout)/2*(1+np.cos(angle)));
        c=1/2*np.sqrt(Rin**2+Rout**2+2*Rin*np.abs(Rout)*np.cos(angle));
        p=b**2/a;
        e=c/a;
        phi_0=np.arccos((Rin**2+4*c**2-Rout**2)/(4*c*Rin))
        phi_max=phi_0+Size/2*Theta_beam;
        phi_min=phi_0-Size/2*Theta_beam;
        t=np.linspace(phi_min,phi_max,101);
        
        dx=np.linspace(0,1,101)
        rho_0=p/(1+e*np.cos(phi_0))*dx;  
        print('phi_0:',phi_0/np.pi*180);
        plt.plot(rho_0*np.cos(phi_0),rho_0*np.sin(phi_0));
        
        
        # plot the mirror size;
        Curve=p/(1+e*np.cos(t))
        plt.plot(Curve*np.cos(t),Curve*np.sin(t),'r-')
        #plot the elliptic 
        t=np.linspace(0,2*np.pi,101);
        print('p',p)
        Curve=p/(1+e*np.cos(t))
        plt.plot(Curve*np.cos(t),Curve*np.sin(t),'k--');
        
    else:
        a=np.abs(Rin+Rout)/2;
        b=np.sqrt(Rin*np.abs(Rout)/2*(1-np.cos(angle)));

        c=1/2*np.sqrt(Rin**2+Rout**2+2*Rin*np.abs(Rout)*np.cos(angle));
        p=b**2/a;
        e=c/a;
        phi_0=np.pi-np.arccos((Rin**2+4*c**2-Rout**2)/(4*c*Rin))
        phi_max=phi_0+Size/2*Theta_beam;
        phi_min=phi_0-Size/2*Theta_beam;
        t=np.linspace(phi_min,phi_max,101);
        
        dx=np.linspace(0,1,101)
        rho_0=p/(1+e*np.cos(phi_0))*dx;  
        print('phi_0:',phi_0/np.pi*180);
        plt.plot(rho_0*np.cos(phi_0),rho_0*np.sin(phi_0));
        
        
        # plot the mirror size;
        Curve=p/(1+e*np.cos(t))
        plt.plot(Curve*np.cos(t),Curve*np.sin(t),'r-')
        #plot the elliptic 
        t=np.linspace(0,2*np.pi,101);
        print('p',p)
        Curve=p/(1+e*np.cos(t))
        plt.plot(Curve*np.cos(t),Curve*np.sin(t),'k--');
        
    print('e:',e,'\n','a:',a,'\n','b:',b,'\n','c:',c)
    return phi_0,Size/2*Theta_beam,e,a,b,c;

def Ellips(focal_length,Rin,Rout,angle,Theta_beam,Size=5):
    angle=angle*np.pi/180;

    a=np.abs(Rin+Rout)/2;
    b=np.sqrt(Rin*np.abs(Rout)/2*(1+np.cos(angle)));
    c=1/2*np.sqrt(Rin**2+Rout**2-2*Rin*np.abs(Rout)*np.cos(angle));
    p=b**2/a;
    e=c/a;
    phi_0=np.arccos((Rin**2+4*c**2-Rout**2)/(4*c*Rin))
    phi_max=phi_0+Size/2*Theta_beam;
    phi_min=phi_0-Size/2*Theta_beam;
    t=np.linspace(phi_min,phi_max,101);
        
    dx=np.linspace(0,1,101)
    rho_0=p/(1+e*np.cos(phi_0))*dx;  
    print('phi_0:',phi_0/np.pi*180);
    plt.plot(rho_0*np.cos(phi_0),rho_0*np.sin(phi_0));
        
        
    # plot the mirror size;
    Curve=p/(1+e*np.cos(t))
    plt.plot(Curve*np.cos(t),Curve*np.sin(t),'r-')
    #plot the elliptic 
    t=np.linspace(0,2*np.pi,101);
    print('p',p)
    Curve=p/(1+e*np.cos(t))
    plt.plot(Curve*np.cos(t),Curve*np.sin(t),'k--');
        
    print('e:',e,'\n','a:',a,'\n','b:',b,'\n','c:',c)
    return phi_0,Size/2*Theta_beam,e,a,b,c;

    
def center_projector(beta,alpha,e,a,c):
    y0=(1+e)*(e*np.cos(beta)+np.cos(alpha))*np.sin(beta)*(a-c)
    Y=(1-e**2)*np.sin(beta)**2+(np.cos(beta)+e*np.cos(alpha))**2
    y0=y0/Y;
    
    a2=((1+e)*np.sin(alpha)*(a-c))**2/Y
    
    b2=((1+e)*(np.cos(beta)+e*np.cos(alpha))*np.sin(alpha)*(a-c))**2/Y**2;
    
    return y0,np.sqrt(a2),np.sqrt(b2);

size=5