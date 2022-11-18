# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 09:27:40 2022

@author: jced0001
"""
import numpy as np
import numpy.linalg as npl
import matplotlib.pyplot as plt

## Constants
h_bar = 0.1973269804e3                                                          # eV nm
m_e   = 0.51099895e6                                                            # eV
pi    = np.pi                                                                   # 3.14
i = complex(0,1)                                                                # Make i a complex i

# Physical parameters
potStrength = 1
a   = 2                                                                         # Period of the varying potential (nm)
L = a/2                                                                         # Half period of the varying potential (nm)
m_eff = 1
s = 101                                                                         # Real-space sampling
x = np.linspace(-L,L,s)                                                         # Real-space x coordinate

# Define the structure of our periodic potential
V  = (np.heaviside(x+a/10,1)+np.heaviside(a/10 - x,1))                          # Single Wall
# V += 0.5*(np.heaviside(x+a/10 - L,1)+np.heaviside(a/10 - x - L,1))              # Add a second wall
V -= np.min(V)
V /= np.max(V)
V -= V[0]
V *= potStrength
# Set up figure
plt.close('all')
fig = plt.figure()
ax1 = fig.add_subplot(2,2,1)
ax2 = fig.add_subplot(2,2,2)
ax3 = fig.add_subplot(2,2,3)
ax4 = fig.add_subplot(2,2,4)

# When solving the central equation numerically like this, it's helpful to
# expamd the periodic potential, V, using Fourier.
# See http://faculty.otterbein.edu/DRobertson/compsci/tise-stud.pdf or
# http://solidstate.mines.edu/videonotes/VN_9_1.pdf. or
# https://silo.tips/download/the-central-equation (best one)
#  
# I used an and bn to calculate vg... just helped me understand what's going on
# See https://www.mathsisfun.com/calculus/fourier-series.html
# 
# vg are the coefficients we need.
# 
N  = 50                                                                         # Number of terms to consider in Fourier expansion. Goes from -N to N
dx = x[1] - x[0]                                                                # dx used for integrating
an = np.zeros(2*N+1)
bn = np.zeros(2*N+1)
for n in np.arange(-N,N+1):
    an[n+N] = (1/L)*np.sum(V*np.cos(pi*x*n/L))*dx
    bn[n+N] = (1/L)*np.sum(V*np.sin(pi*x*n/L))*dx

# A  = np.zeros_like(x)
# B  = np.zeros_like(x)
vg = np.zeros(2*N+1,dtype = 'complex_')                                         
C  = np.zeros_like(x,dtype = 'complex_')
for n in np.arange(-N,N+1):
    # A += an[n+N]*np.cos(n*pi*x/L)
    # B += bn[n+N]*np.sin(n*pi*x/L)
    
    Gn = n*pi/L                                                                 # Convert n to G since we're summing over reciprocal lattice vectors
    vg[n+N] = potStrength*(an[n+N] + i*bn[n+N])/2
    C += vg[n+N]*np.exp(i*Gn*x)

# F = A + B                                                                     # Equivalent to C

ax1.plot(x,np.real(C) + np.imag(C)); ax1.set_title(r'$V(x) = \sum{v_{G}e^{iGx}}$') # Plot the potential from the coefficients we just calculated
ax1.set_xlabel("x (nm)")
ax1.set_ylabel("Energy (eV)")

# 
ks = 200
K  = np.linspace(-pi/a,pi/a,2*ks +1)                                            # k-sampling
Gn = np.arange(-N,N+1)*pi/L                                                     # Reciprocal lattice vectors
HV = np.zeros((int(len(Gn)/2)+1,int(len(Gn)/2)+1),dtype = 'complex_')           # Potential component of the Hamiltonian
for n,v in enumerate(vg):
    if(n < N+1):                                                                # Create the upper triangle of the matrix
        diag = np.diag(np.array([v]*(n+1)),N-n)
    else:                                                                       # Create the lower triangle of the matrix
        diag = np.diag(np.array([v]*(2*N+1-n)),N-n)
    HV += diag

Ek  = np.zeros((N+1,2*ks +1))                                                   # Kinetic terms
Coeff = np.array(np.zeros((N+1)*(N+1)*len(K)),dtype = 'complex_').reshape((N+1,N+1,len(K))) # Fourier coefficients of our wavefn
mid = int(N)
diag_g = Gn[mid-int(N/2):mid+int(N/2)+1]                                        # Reciprocal vectors used on the main diagonal
for n,k in enumerate(K):                                                        # Now we're calculating the coefficients of our wavefn for each momentum
    diag = ((h_bar**2)/(2*m_e*m_eff)) * ((k + diag_g))**2                       # See equation: x in y
    HE   = np.diag(diag,0)                                                      # 
    
    H = HE + HV
    vals,vecs = npl.eigh(H)
    Ek[:,n] = vals
    Coeff[:,:,n] = vecs

# Band Structure
numBands = 4
for b,e in enumerate(Ek):
    ax3.plot(K,e,label='n = ' + str(b)); ax3.set_title(r'Band Structure')
    if(b + 1 == numBands): break
ax3.legend()
ax3.set_xlabel("Wave Vector k")
ax3.set_ylabel("Energy (eV)")

# Wavefunction
E    = -0.1                                                                       # Energy of the wfn to plot (eV)
E = int(100*np.min(Ek[0]))/100
# E = 2
emin = E - 0.025
emax = E + 0.025
B = ((Ek >= emin) & (Ek <= emax))
psi_total = np.zeros(s)
for row in range(len(Coeff[0])):
    for col in range(len(Coeff[0][0])):
        if(not B[row,col]): continue
        kg = K[col] + diag_g.reshape((N+1,1))
        psi_k = sum(np.exp(i*np.matmul(kg,x.reshape((1,s))))*Coeff[:,row,col].reshape((N+1,1)),0)
        # psi_total += np.real(psi_k)
        psi_total += abs(psi_k)
            
P = psi_total**2
# P = np.real(psi_total)
if(np.sum(P) > 0):
    P /= (np.sum(P)*dx)

ax2.plot(x,(V - np.min(V))/np.max((V - np.min(V))),linestyle='dashed', label=r'$V(x)$'); ax2.set_title(r'E = ' + str(E) + " eV")
ax2.plot(x,P,label=r'$|\psi|^{2}$')
# ax2.legend()
ax2.set_xlabel("x (nm)")
ax2.set_ylabel(r'$|\psi|^{2}$')

# LDOS
x0  = 0                                                                         # Location of the LDOS spectrum (nm)
x0idx = np.argmin(abs(x - x0))
Exmax = np.max(Ek[numBands-1])
Exmin = np.min(Ek[0])
w = 0.025*(Exmax - Exmin)
ovrfl = 3*w
Ex = np.linspace(Exmin-ovrfl,Exmax+ovrfl,200)
DOS  = np.zeros_like(Ex)
LDOS = np.zeros_like(Ex)
for n,e in enumerate(Ex):
    emin = e - w
    emax = e + w
    B = ((Ek >= emin) & (Ek <= emax))
    psi_total = np.zeros(s)
    dos  = 0
    ldos = 0
    for row in range(len(Coeff[0])):
        for col in range(len(Coeff[0][0])):
            if(not B[row,col]): continue
            kg = K[col] + diag_g.reshape((N+1,1))
            psi_k = sum(np.exp(i*np.matmul(kg,x.reshape((1,s))))*Coeff[:,row,col].reshape((N+1,1)),0)
            dos  += np.sum(abs(psi_k))
            ldos += abs(psi_k[x0idx])
    DOS[n]  = dos
    LDOS[n] = ldos

LDOS *= 100/np.max(DOS)
DOS  /= np.max(DOS)
ax4.plot(DOS ,Ex,label='DOS')
ax4.plot(LDOS,Ex,label='100xLDOS (x=' + str(x0) + "nm)")
ax4.set_title("DOS and LDOS at x = "   + str(x0) + " nm")
ax4.set_xlabel("LDOS (arb)")
ax4.set_ylabel("Energy (eV)")
ax4.legend()