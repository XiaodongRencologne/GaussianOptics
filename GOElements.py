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
def Mirror(GO_para,Angle, Scale=5):
    Angle=Angle/180*np.pi
    a=np.abs(GO_para['Rin']+GO_para['Rout'])/2
    
    Rout



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