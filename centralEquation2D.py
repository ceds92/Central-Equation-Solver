# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 09:06:35 2022

@author: jced0001
"""
import numpy as np
import numpy.linalg as npl
import matplotlib.pyplot as plt
import potentials as pot

# Parameters
rs = 101                                                                        # Real-space sampling
ks = 15                                                                         # K-space sampling
N  = 10                                                                         # Number of terms to consider in Fourier space
v0 = -3/4                                                                       # Characteristic potential strength (eV)
m_eff = 0.400                                                                   # Effective electron mass me*/me
L = 1.0                                                                         # Nearest neighbor spacing (nm)
r = 0.3                                                                         # Radii of potential wells

# Physical constants
pi = np.pi
i  = complex(0,1)
h_bar = 0.1973269804e3                                                          # Planck's constant (eV.nm)
m_e   = 0.51099895e6                                                            # Electron rest mass (eV)

# Build real-space and UC of periodic potential
# a,X,V   = pot.muffin_top(L, r, rs)
# a,X,V   = pot.hexagonal_primitive(L, r, rs)
a,X,V   = pot.kagome_primitive(L, r, rs)
V      *= v0
a1,a2   = a
norm_a  = npl.norm(a,axis=1)
x1,x2   = X
X1,X2   = np.meshgrid(x1,x2)
dx1,dx2 = x1[1]-x1[0], x2[1]-x2[0]                                              # dx,dy used for integrating

# Reciprocal lattice vectors
c  = 2*pi/(a1[0]*a2[1] - a1[1]*a2[0])
b1 = c*np.array([a2[1],-a2[0]])
b2 = c*np.array([-a1[1],a1[0]])

vg = np.zeros(((2*N+1),(2*N+1)),dtype = 'complex_')                                         
C  = np.zeros_like(X1,dtype = 'complex_')

KEG = []
Gn = np.arange(-N,N+1)[...,None]*b1
Gm = np.arange(-N,N+1)[...,None]*b2
print("Calculating Fourier Coefficients...")
for n in np.arange(0,2*N+1):
    for m in np.arange(0,2*N+1):
        vg[n,m] = np.sum(V*np.exp(-i*(Gn[n][0]*X1 + Gm[m][1]*X2 + Gn[n][1]*X2 + Gm[m][0]*X1)))*dx1*dx2
        KEG.append(Gn[n] + Gm[m])
        if(n%int(0.2*(2*N+1)) == 0):
            if(m == 0):
                print(int(100*n/(2*N+1)),'%')
        
# Rebuild potential from coefficients to check...
x = 3*x1
y = 3*x2
X,Y = np.meshgrid(x,y)
print("Rebuilding Potential...")
for n in np.arange(0,2*N+1):
    for m in np.arange(0,2*N+1):
        C += vg[n,m]*np.exp(i*(Gn[n][0] + Gm[m][0])*X)*np.exp(i*(Gn[n][1] + Gm[m][1])*Y)
        if(n%int(0.2*(2*N+1)) == 0):
            if(m == 0):
                print(int(100*n/(2*N+1)),'%')

mid = int((1+(2*N+1)**2)/2)
bound = int(mid/2)
KEG = np.array(KEG)[mid-bound:mid+bound+1]

vv = vg.reshape((len(vg)*len(vg[0])))
UG = np.zeros((mid,mid),dtype='complex_')
print("Constructing Potential Matrix")
for nv, v in enumerate(vv):
    if(nv < mid):
        diag = np.diag([v]*(nv+1),mid-1-nv)
    else:
        diag = np.diag([v]*((2*N+1)*(2*N+1) - nv),mid-1-nv)
    UG += diag
    
    if((nv%int(0.1*len(vv))) == 0): print(int(100*nv/len(vv))+1,'%')
    
# Set up k space
maxkx = abs(b1[0]/2) + abs(b2[0]/2)
maxky = abs(b1[1]/2) + abs(b2[1]/2)
maxkk = np.array((max([maxkx,maxky]),max([maxkx,maxky])))
kk = np.linspace(-maxkk,maxkk,ks)
Ek = np.zeros((mid,ks,ks))
print("Calculating Kinetic Terms and Eigen Values...")
Coeff = np.array(np.zeros(mid**2*len(kk)**2),dtype = 'complex_').reshape((mid,mid,len(kk),len(kk))) # Fourier coefficients of our wavefn
for kn,kx in enumerate(kk[:,0]):
    for km,ky in enumerate(kk[:,1]):
        kg = np.array([kx,ky]) - KEG
        diag = ((h_bar**2)/(2*m_e*m_eff)) * (npl.norm(kg,axis=1)**2)
        
        HE = np.diag(diag,0)
        
        H = HE + UG
        vals,vecs = npl.eigh(H)
        Ek[:,kn,km] = vals
        
        Coeff[:,:,kn,km] = vecs
        
        k = kn*len(kk)+km
        if((k%int(0.1*len(kk)**2)) == 0): print(int(100*k/len(kk)**2)+1,'%')

# %%
# Wavefunction
E    = -0.1                                                                       # Energy of the wfn to plot (eV)
E = int(100*np.min(Ek[0]))/100
# E = 2
emin = E - 0.025
emax = E + 0.025
B = ((Ek >= emin) & (Ek <= emax))
psi_total = np.zeros_like(X1)
R = np.sqrt(X1**2 + X2**2)
for row in range(mid):
    for kn,kx in enumerate(kk[:,0]):
        for km,ky in enumerate(kk[:,1]):
            if(not B[row,kn,km]): continue
            k  = np.array([kx,ky])
            kg = k+KEG
            kg1 = kg[:,0].reshape((mid,1))
            kg2 = kg[:,1].reshape((mid,1))
            # exp = np.exp(i*(np.matmul(kg1,x1.reshape((1,rs))) + np.matmul(kg2,x2.reshape((1,rs)))))
            exp = np.exp(i*(np.matmul(kg1,x1.reshape((1,rs))) + np.matmul(kg2,x2.reshape((1,rs)))))
            psi_k = np.sum(Coeff[:,row,kn,km].reshape((mid,1))*exp,0)
            psi_total += abs(psi_k)
    
    if((row%int(0.05*mid)) == 0): print(int(100*row/mid)+1,'%')

# %%
# Plotting
fig = plt.figure()
ax1 = fig.add_subplot(2,2,1)
ax2 = fig.add_subplot(2,2,3)
ax3 = fig.add_subplot(1,2,2,projection='3d')
# ax4 = fig.add_subplot(2,2,4)

extent=np.array([min(x1),max(x1),min(x2),max(x2)])

# ax1.imshow(V,extent=extent)
ax1.imshow(np.real(C) + np.imag(C),extent=extent*3)

mid = int(len(x1)/2)
u = []
oo = np.array(x1[mid],x2[mid])
oo = oo -a1/2 -a2/2
u.append(oo + a1)
u.append(u[0] + a2)
u.append(u[1] - a1)
u.append(u[2] - a2)
u.append(u[0])

for p in range(5):
    xx = u[p][0],u[p+1][0]
    yy = u[p][1],u[p+1][1]
    ax1.plot(xx,yy,c='red',linewidth=1)
    if(p==3): break


numbands = 3
K1,K2 = np.meshgrid(kk[:,0],kk[:,1])
for ne,e in enumerate(Ek):
    ax3.plot_surface(K1, K2, e, linewidth=0, antialiased=True)
    # ax2.plot(kk[:,0]*a1[0]/pi,e[:,int(ks/2)])
    ax2.plot(kk[:,1]*a2[1]/pi,e[int(ks/2),:])
    # daxis = np.sqrt(kk[:,1]**2 + kk[:,0]**2)
    # daxis[0:15] *= -1
    # ax2.plot(daxis,np.diag(e[:,:]))
    if(ne == numbands - 1):
        break