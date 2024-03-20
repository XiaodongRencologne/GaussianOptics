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

def ThinLens(win, din, f, Lambda):
    # para:  'f' or 'din'
    # para:  
    z_c, w_z, Rin, q1, phi1, b1 =Gaussian_propagation2d(win,din,Lambda)
    Mirror_para={'Lambda':Lambda,
                 'f':f,
                 'win':win,'din':din,'Rin':Rin,'q1':q1, 'phi1': phi1,
                 'wout':None,'dout':None,'Rout':0,'q2': 0, 'phi2': None,
                 'M':0}
  
    M_max=f/z_c
    factor=(din/f-1)**2+(z_c/f)**2
    dout=f*(1+(din/f-1)/factor)
    wout=win/np.sqrt(factor)
    Rout=1/(1/f-1/Rin)
    M=wout/win
    z_out_c=np.pi*wout**2/Lambda
    q2=dout+1j*z_out_c

    Mirror_para['wout']=wout
    Mirror_para['dout']=dout
    Mirror_para['Rout']=Rout
    Mirror_para['q2']=q2
    Mirror_para['phi2']=np.arctan(dout/z_out_c)
    Mirror_para['M']=M
    for key in Mirror_para.keys():
        print(key,':',Mirror_para[key])
    return Mirror_para



def ThinLens_f_M(win, wout, f, Lambda):
        M=wout/win
        z_c_in = np.pi*win**2/Lambda
        z_c_out= np.pi*wout**2/Lambda
        Mirror_para1={'Lambda':Lambda,
                        'f':f,
                        'win':win,'din':None,'Rin':None,'q1':None, 'phi1':None,
                        'wout':wout,'dout':None,'Rout':None,'q2': None, 'phi2': None,
                        'M':M}
        Mirror_para2={'Lambda':Lambda,
                        'f':f,
                        'win':win,'din':None,'Rin':None,'q1':None, 'phi1':None,
                        'wout':wout,'dout':None,'Rout':None,'q2': None, 'phi2': None,
                        'M':M}
        f0=np.pi*win*wout/Lambda
        # pos
        din=f+1/M*np.sqrt(f**2-f0**2);Rin=din*(1+(np.pi*win**2/Lambda/din)**2)
        dout=f+M*np.sqrt(f**2-f0**2);Rout=1/(1/f-1/Rin)
        q1=din+1j*z_c_in
        q2=dout+1j*z_c_out
        phi1=np.arctan(din/z_c_in)
        phi2=np.arctan(dout/z_c_out)
        Mirror_para1['din']=din
        Mirror_para1['Rin']=Rin
        Mirror_para1['dout']=dout
        Mirror_para1['Rout']=Rout
        Mirror_para1['q1']=q1
        Mirror_para1['q2']=q2
        Mirror_para1['phi1']=phi1
        Mirror_para1['phi2']=phi2

        # neg
        din=f-1/M*np.sqrt(f**2-f0**2);Rin=din*(1+(np.pi*win**2/Lambda/din)**2)
        dout=f-M*np.sqrt(f**2-f0**2);Rout=1/(1/f-1/Rin)
        q1=din+1j*z_c_in
        q2=dout+1j*z_c_out
        phi1=np.arctan(din/z_c_in)
        phi2=np.arctan(dout/z_c_out)
        Mirror_para2['din']=din
        Mirror_para2['Rin']=Rin
        Mirror_para2['dout']=dout
        Mirror_para2['Rout']=Rout
        Mirror_para2['q1']=q1
        Mirror_para2['q2']=q2
        Mirror_para2['phi1']=phi1
        Mirror_para2['phi2']=phi2
        Mirror_para2['din']=din
        Mirror_para2['dout']=dout
        for key in Mirror_para1.keys():
             print(key,':',Mirror_para1[key])
        for key in Mirror_para2.keys():
             print(key,':',Mirror_para2[key])
        return Mirror_para1, Mirror_para2

def ThinLens_d_M(win, wout, d, Lambda):
        M=wout/win
        z_c_in = np.pi*win**2/Lambda
        z_c_out= np.pi*wout**2/Lambda
        f0=np.pi*win*wout/Lambda
        Mirror_para1={'Lambda':Lambda,
                        'f':None,
                        'win':win,'din':None,'Rin':None,'q1':None, 'phi1':None,
                        'wout':wout,'dout':None,'Rout':None,'q2': None, 'phi2': None,
                        'M':M}
        Mirror_para2={'Lambda':Lambda,
                        'f':None,
                        'win':win,'din':None,'Rin':None,'q1':None, 'phi1':None,
                        'wout':wout,'dout':None,'Rout':None,'q2': None, 'phi2': None,
                        'M':M}
        if M!=1.0:
            f=(np.sqrt((M-1/M)**2*f0**2+d**2)*(M+1/M)-2*d)/(M-1/M)**2
            din=(d-f*(1-M**2))/(1+M**2);Rin=din*(1+(np.pi*win**2/Lambda/din)**2)
            dout=(M**2*d+f*(1-M**2))/(1+M**2);Rout=1/(1/f-1/Rin)
            q1=din+1j*z_c_in
            q2=dout+1j*z_c_out
            phi1=np.arctan(din/z_c_in)
            phi2=np.arctan(dout/z_c_out)
            Mirror_para1['f']=f
            Mirror_para1['din']=din
            Mirror_para1['Rin']=Rin
            Mirror_para1['dout']=dout
            Mirror_para1['Rout']=Rout
            Mirror_para1['q1']=q1
            Mirror_para1['q2']=q2
            Mirror_para1['phi1']=phi1
            Mirror_para1['phi2']=phi2

            f=(-np.sqrt((M-1/M)**2*f0**2+d**2)*(M+1/M)-2*d)/(M-1/M)**2
            din=(d-f*(1-M**2))/(1+M**2);Rin=din*(1+(np.pi*win**2/Lambda/din)**2)
            dout=(M**2*d+f*(1-M**2))/(1+M**2);Rout=1/(1/f-1/Rin)
            q1=din+1j*z_c_in
            q2=dout+1j*z_c_out
            phi1=np.arctan(din/z_c_in)
            phi2=np.arctan(dout/z_c_out)
            Mirror_para2['f']=f
            Mirror_para2['din']=din
            Mirror_para2['Rin']=Rin
            Mirror_para2['dout']=dout
            Mirror_para2['Rout']=Rout
            Mirror_para2['q1']=q1
            Mirror_para2['q2']=q2
            Mirror_para2['phi1']=phi1
            Mirror_para2['phi2']=phi2

        elif M==1.0:
            f=d/4+f0**2/d
            din=(d-f*(1-M**2))/(1+M**2);Rin=din*(1+(np.pi*win**2/Lambda/din)**2)
            dout=(M**2*d+f*(1-M**2))/(1+M**2);Rout=1/(1/f-1/Rin)
            q1=din+1j*z_c_in
            q2=dout+1j*z_c_out
            phi1=np.arctan(din/z_c_in)
            phi2=np.arctan(dout/z_c_out)
            Mirror_para1['f']=f
            Mirror_para1['din']=din
            Mirror_para1['Rin']=Rin
            Mirror_para1['dout']=dout
            Mirror_para1['Rout']=Rout
            Mirror_para1['q1']=q1
            Mirror_para1['q2']=q2
            Mirror_para1['phi1']=phi1
            Mirror_para1['phi2']=phi2
        for key in Mirror_para1.keys():
             print(key,':',Mirror_para1[key])
        for key in Mirror_para2.keys():
             print(key,':',Mirror_para2[key])
        return Mirror_para1, Mirror_para2
# %%
p=ThinLens(1.4, 33.73731665, 67.80690232, 1.0162456203)
# %%
p1,p2=ThinLens_f_M(1.4, 1.4, 67.80690232, 1.0162456203)
# %%
p=ThinLens(1.4, 135.34254813273228, 67.80690232, 1.0162456203)
# %%
p=ThinLens(1.4, 0.2712565072677364, 67.80690232, 1.0162456203)
# %%
p1,p2=ThinLens_d_M(1.4, 1.4, 0.2712565072677364+0.2712565072677176, 1.0162456203)
# %%
