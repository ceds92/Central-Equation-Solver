# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 17:33:43 2022

@author: jced0001
"""
import numpy as np
import numpy.linalg as npl
import pickle
import global_

i     = complex(0,1)                                                            # Complex i
pi    = np.pi                                                                   # 3.14
h_bar = 0.1973269804e3                                                          # Planck's constant (eV.nm)
m_e   = 0.51099895e6                                                            # Electron rest mass (eV)
class EPWE():
    def __init__(self,ks=15,N=5,m_eff=1):
        self.N      = N                                                         # Number of terms to consider in Fourier space
        self.ks     = ks                                                        # K-space sampling
        self.m_eff  = m_eff                                                     # Effective electron mass me*/me
        
        self.Ek     = []                                                        # This is where the kinetic terms go
        self.Coeff  = []                                                        # This is where the coefficients of the wavefunction go
        self.psi    = {}                                                        # This is where calculated wavefunctions get stored
        self.C      = {}                                                        # This is where the reconstructed potential will go
        
        self.LDOS    = []                                                       # Computed LDOS curves go here
        self.LDOSEx  = []
        self.LDOSx0s = []
        
        self.valid = False
        self.running = {"main":False,
                        "LDOS":False,
                        "map": False,
                        "potential": False}
    
    def stop(self,name):
        exec("global_." + name + "_running.clear()")
        self.running[name] = False
        
    def run(self,a,X,V,initiator=None):
        """
        Function to execute EPWE for a given potential with lattice vectors
        a=[a1,a2] in the space defined by X=[x1,x2]

        Parameters
        ----------
        a : Lattice vectors (np.array([a1,a2]))
        X : Real space ([x1,x2])
        V : 2D potential map

        """
        name = "main"
        self.running[name] = True
        self.valid = False
        
        G,K,bz   = self.brillouinZone(a,initiator)                              # Determine the first Brillouin zne
        if(self.checkEventFlags(name)): self.stop(name); return                 # Sim has been stopped. Clear the running flag and return
        
        Gnm,vg   = self.fourierCoefs(a,X,V,G,initiator)                         # Determine the Fourier coefficients of the potential
        if(self.checkEventFlags(name)): self.stop(name); return                 # Sim has been stopped. Clear the running flag and return
        
        U        = self.constructPotentialMatrix(vg,initiator)                  # Construct the potential matrix used to solve the TISE for a periodic potential
        if(self.checkEventFlags(name)): self.stop(name); return                 # Sim has been stopped. Clear the running flag and return
        
        Ek,Coeff = self.solveTISE(Gnm,K,U,bz,initiator)                         # Solve the TISE
        if(self.checkEventFlags(name)): self.stop(name); return                 # Sim has been stopped. Clear the running flag and return
        
        self.a = a                                                              # Primitive lattice vectors
        self.X = X                                                              # Real space [x1,x2]
        self.V = V                                                              # Potential map
        
        self.K     = K                                                          # k space defined by [k1,k2]
        self.bz    = bz                                                         # Frist Brillouin zone
        self.G     = G                                                          # Reciprocal lattice vectors [Gn,Gm]
        self.Gnm   = Gnm                                                        # All combinations of reciprocal vectors
        self.vg    = vg                                                         # Fourier coefficients of the periodic potentil
        self.U     = U                                                          # The potential matrix
        self.Ek    = Ek                                                         # Kinetic terms
        self.Coeff = Coeff                                                      # Wavefunction coefficients
        
        self.valid = True
        
        if(initiator): initiator.simFinished(True)
        
        self.stop("main"); return                                               # Sim has finished. Clear the running flag and return
        
    def save(self,name="EPWESIM"):
        """
        This function saves the simulation results
        
        Parameters
        ----------
        name : epwe filename
        
        """
        if(not self.valid):
            print("Error, simulation not valid. Need to rerun before saving.")
            return
        
        if(not name.endswith('.epwe')): name += '.epwe'
        
        pkldict = {"N"      : self.N,
                   "ks"     : self.ks,
                   "m_eff"  : self.m_eff,
                   "a"      : self.a,
                   "X"      : self.X,
                   "V"      : self.V,
                   "K"      : self.K,
                   "bz"     : self.bz,
                   "G"      : self.G,
                   "Gnm"    : self.Gnm,
                   "vg"     : self.vg,
                   "U"      : self.U,
                   "Ek"     : self.Ek,
                   "Coeff"  : self.Coeff,
                   "psi"    : self.psi,
                   "C"      : self.C,
                   "LDOS"   : self.LDOS,
                   "LDOSEx" : self.LDOSEx,
                   "LDOSx0s": self.LDOSx0s}
        
        pickle.dump(pkldict,open(name,'wb'))
    
    def load(self,name):
        """
        This function loads previous simulation results

        Parameters
        ----------
        name : epwe filename
        
        """
        pkldict = pickle.load(open(name,'rb'))
        
        self.N      = pkldict["N"]
        self.ks     = pkldict["ks"]
        self.m_eff  = pkldict["m_eff"]
        self.a      = pkldict["a"]
        self.X      = pkldict["X"]
        self.V      = pkldict["V"]
        self.K      = pkldict["K"]
        self.bz     = pkldict["bz"]
        self.G      = pkldict["G"]
        self.Gnm    = pkldict["Gnm"]
        self.vg     = pkldict["vg"]
        self.U      = pkldict["U"]
        self.Ek     = pkldict["Ek"]
        self.Coeff  = pkldict["Coeff"]
        self.psi    = pkldict["psi"]
        self.C      = pkldict["C"]
        self.LDOS   = pkldict["LDOS"]
        self.LDOSEx = pkldict["LDOSEx"]
        self.LDOSx0s= pkldict["LDOSx0s"]
        
        self.valid = True
    
    def brillouinZone(self,a,initiator=None):
        """
        This function calculates the first Brillouin zone from the primitive 
        lattice vectors

        Parameters
        ----------
        a : Primitive lattice vectors [a1,a2]

        Returns
        -------
        G  : Reciprocal lattice vectors [Gn,Gm]
        K  : K space defined by [k1,k2]
        bz : 2D map of the first Brillouin zone

        """
        if(initiator and not self.checkEventFlags("main")):
            initiator.progress("Mapping the first Brillouin zone")
            
        print("Mapping the first Brillouin zone")
        N  = self.N
        ks = self.ks
        
        # Reciprocal lattice vectors
        c  = 2*pi/(a[0][0]*a[1][1] - a[0][1]*a[1][0])
        b1 = c*np.array([a[1][1],-a[1][0]])                                     # Primitive reciprocal lattice vector b1
        b2 = c*np.array([-a[0][1],a[0][0]])                                     # Primitive reciprocal lattice vector b2
        
        Gn = np.arange(-N,N+1)[...,None]*b1                                     # Array of 2N+1 integer multiples of the reciprocal lattice vector b1
        Gm = np.arange(-N,N+1)[...,None]*b2                                     # Array of 2N+1 integer multiples of the reciprocal lattice vector b2
        
        # Set up k space
        na = npl.norm(a,axis=1)
        k1,k2 = np.transpose(np.linspace(-2*pi/na,2*pi/na,ks))                  # K sampling
        K1,K2 = np.meshgrid(k1,k2)
        gamma = int(ks/2)
        
        bz = np.zeros_like(K1)
        for q in [-1,0,1]:
            for w in [-1,0,1]:
                if(q==0 and w==0): continue
                G = q*b1 + w*b2                                                 # Reciprocal lattice vector
                p = G/2
                if(G[1] == 0):
                    if(p[0] < 0):
                        bz += K1 > p[0]
                    if(p[0] > 0):
                        bz += K1 < p[0]
                else:
                    m = -G[0]/G[1]
                    c = p[1] - m*p[0]
                    bz += K2 > (m*K1 + c)
                keepValue = bz[gamma,gamma]
                mask = bz == keepValue
                bz[~mask] = -1
                
        bz = np.array(bz > 0)
        
        print("Done")
        return np.array([Gn,Gm]), np.array([k1,k2]), bz
    
    def fourierCoefs(self,a,X,V,G,initiator=None):
        """
        This function calculates 2N+1 Fourier coefficients for a periodic 2D
        function.

        Parameters
        ----------
        N : Number of coefficients to calculate (2*N+1)
        V : 2D np.array containing the unit cell of the periodic function
        a : Lattice vectors [a1,a2]
        X : Real space [x1,x2]

        Returns
        -------
        G  : Every combination of the reciprocal vectors used
        vg : Matrix containing all of the Fourier coefficients

        """
        if(initiator and not self.checkEventFlags("main")):
            initiator.progress("Constructing Fourier terms...0 %")
        print("Calculating Fourier Coefficients...")
        N = self.N
        
        x1,x2 = X                                                               # Extract the real-space axis
        X1,X2 = np.meshgrid(x1,x2)                                              # The mesh
        dx1,dx2 = x1[1]-x1[0], x2[1]-x2[0]                                      # dx,dy used for integrating
        
        Gnm   = []                                                              # Initialse array containing all combinations of the reciprocal vectors
        Gn,Gm = G                                                               # Reciprocal lattice vectors
        vg = np.zeros(((2*N+1),(2*N+1)),dtype = 'complex_')                     # Initialise the matrix which will contain the fourier coefficients of the potential, V                                     
        
        for n in np.arange(0,2*N+1):                                            # Summing over all combinations of Gn and Gm...
            exp_n = np.exp(-i*(Gn[n][0]*X1 + Gn[n][1]*X2))                      # Exponential term for this Gn
            
            for m in np.arange(0,2*N+1):
                exp_m = np.exp(-i*(Gm[m][0]*X1 + Gm[m][1]*X2))                  # Exponential term for this Gm
                vg[n,m] = np.sum(V*exp_n*exp_m)*dx1*dx2                         # Integrate to get the n,m coefficient
                vg[n,m] = vg[n,m]/np.cross(a[0],a[1])                           # Normalise by dividing by the area of the unit cell
                
                Gnm.append(Gn[n] + Gm[m])
                if(n%int(0.2*(2*N+1)) == 0):
                    if(m == 0):
                        progress = int(100*n/(2*N+1))
                        print(progress,'%')
                        if(initiator and not self.checkEventFlags("main")):
                            progress = "Calculating Fourier terms..." + str(progress) + '%'
                            initiator.progress(progress)
            
                if(self.checkEventFlags("main")): return 0,0                    # Sim has been stopped. Clear the running flag and return
        
        print("Done")
        return np.array(Gnm),vg
        
    def constructPotentialMatrix(self,vg,initiator=None):
        """
        This function constructs the potential matrix used in solving the TISE

        Parameters
        ----------
        vg : Fourier coefficients of the periodic potential

        Returns
        -------
        U : The potential matrix used in solving the TISE

        """
        if(initiator and not self.checkEventFlags("main")):
            initiator.progress("Constructing potential matrix...0 %")
        print("Constructing Potential Matrix")
        N = self.N
        
        n = int((1+(2*N+1)**2)/2)                                               # The dimension of the potential matrix is
        U = np.zeros((n,n),dtype='complex_')
        vv = vg.reshape((len(vg)*len(vg[0])))
        for nv, v in enumerate(vv):
            if(nv < n):
                diag = np.diag([v]*(nv+1),n-1-nv)
            else:
                diag = np.diag([v]*((2*N+1)*(2*N+1) - nv),n-1-nv)
            U += diag
            
            if((nv%int(0.1*len(vv))) == 0):
                progress = int(100*nv/len(vv))+1
                print(progress,'%')
                if(initiator and not self.checkEventFlags("main")):
                    progress = "Constructing potential Matrix..." + str(progress) + '%'
                    initiator.progress(progress)

            if(self.checkEventFlags("main")): return 0                          # Sim has been stopped. Clear the running flag and return
        
        print("Done")
        return U
    
    def solveTISE(self,Gnm,K,U,bz,initiator):
        """
        This function solves the TISE for the system with periodic potential, V

        Parameters
        ----------
        Gnm : Combinations of the reciprocal lattice vectors
        K   : k space, defined as [k1,k2]
        U   : Potential matrix
        bz  : Map of the first Brillouin zone

        Returns
        -------
        Ek    : Kinetic matrix containing eigen values/energies
        Coeff : Fourier coefficients of the eigenstates/wavefunctions

        """
        if(initiator and not self.checkEventFlags("main")):
            initiator.progress("Solving TISE...0 %")
        print("Calculating Kinetic Terms and Eigen Values...")
        N     = self.N
        ks    = self.ks
        m_eff = self.m_eff
        n     = int((1+(2*N+1)**2)/2)                                           # The dimension of the hamiltonian
        
        k1,k2 = K
        
        Ek = np.zeros((n,ks,ks))*np.nan
        Coeff = np.array(np.zeros(n**2*len(k1)**2),dtype = 'complex_')          # Fourier coefficients of our wavefn will be stored here
        Coeff = Coeff.reshape((n,n,len(k1),len(k2)))                            # Can't make an np.array this shape unless we do it like this for some reason
        for kn,kx in enumerate(k1):
            for km,ky in enumerate(k2):
                if(bz[kn,km] == 0): continue
                kg = np.array([kx,ky]) - Gnm[n-int(n/2):n+int(n/2)+1]
                diag = ((h_bar**2)/(2*m_e*m_eff))*(npl.norm(kg,axis=1)**2)

                H0 = np.diag(diag,0)
                H  = H0 + U
                
                vals,vecs        = npl.eigh(H)
                Ek[:,kn,km]      = vals
                Coeff[:,:,kn,km] = vecs
                
                k = kn*ks+km
                if(k > 10):
                    if((k%int(0.1*ks**2)) == 0):
                        progress = int(100*k/ks**2)+1
                        print(progress,'%')
                        if(initiator and not self.checkEventFlags("main")):
                            progress = "Solving TISE..." + str(progress) + '%'
                            initiator.progress(progress)
            
            if(self.checkEventFlags("main")): return 0,0                        # Sim has been stopped. Clear the running flag and return
                
        print("Done")
        return Ek,Coeff
    
    def rebuildPotential(self,scale,initiator=None):
        """
        This function reconstructs the potential, V, from the Fourier 
        coefficients calculated in 'fourierCoefs'

        Parameters
        ----------
        scale : Factor to plot the rebuilt potential over an extended real 
                space. (i.e. X => X*scale)

        Returns
        -------
        C       : The reconstructed potential
        extent  : Real space extent (used for plotting, e.g. imshow)

        """
        self.running["potential"] = True
        
        print("Rebuilding Potential...")
        N     = self.N
        vg    = self.vg
        Gn,Gm = self.G
        
        x1,x2 = self.X*scale
        X1,X2 = np.meshgrid(x1,x2)
        C  = np.zeros_like(X1,dtype = 'complex_')
        
        for n in np.arange(0,2*N+1):
            exp_n = np.exp(i*(Gn[n][0]*X1 + Gn[n][1]*X2))
            
            for m in np.arange(0,2*N+1):
                exp_m = np.exp(i*(Gm[m][0]*X1 + Gm[m][1]*X2))
                C += vg[n,m]*exp_n*exp_m
                if(n%int(0.2*(2*N+1)) == 0):
                    if(m == 0):
                        progress = int(100*n/(2*N+1))
                        print(progress,'%')
                        if(initiator and not self.checkEventFlags("potential")):
                            initiator.progress(progress)
            
                if(initiator and self.checkEventFlags("potential")):
                    self.stop("potential")
                    return                                                      # Sim has been stopped. Clear the running flag and return
        
        print("Done")
        
        C = np.real(C) + np.imag(C)
        extent = np.array([min(x1),max(x1),min(x2),max(x2)])
        
        self.C = {"C"      : C,
                  "X"      : np.array([x1,x2]),
                  "extent" : extent}
        
        if(initiator):
            initiator.finish(True,C,extent,np.array([x1,x2]))
        
        self.stop("potential")
        return C,extent,np.array([x1,x2])
    
    def getWavefunction(self,E,scale=1,initiator=None):
        """
        This function computes the real-space density distribution within a
        given energy range, E.

        Parameters
        ----------
        E       : Can be either a single valued energy (float32) or two-valued
                  np.array([emin,emax])
        scale   : Factor to plot the distribution over an extended real space.
                  (i.e. number of unit cells)

        Returns
        -------
        psi    : Electron density distribution
        extent : Extent used for plotting (e.g. imshow)
        X      : Scaled real space

        """
        self.running["map"] = True
        
        print("Calculating Wavefunction...")
        N     = self.N
        bz    = self.bz
        k1,k2 = self.K
        Ek    = self.Ek
        Gnm   = self.Gnm
        Coeff = self.Coeff
        
        n     = int((1+(2*N+1)**2)/2)                                           # The dimension of the hamiltonian
        bound = int(n/2)
        
        B = Ek == E
        energy = E
        if(np.array(E).shape):
            emin = E[0]; emax = E[1]
            B = ((Ek >= emin) & (Ek <= emax))
            energy = int(1000*np.sum(E)/2)/1000
        
        B = B*bz
        
        x1,x2 = self.X*scale                                                    # Extract the real-space axis
        X1,X2 = np.meshgrid(x1,x2)                                              # The mesh
        
        count = 0
        psi   = np.zeros_like(X1)                                               # Sum of the states within the energy window emin to emax wil go here
        for row in range(n):
            for kn,kx in enumerate(k1):
                for km,ky in enumerate(k2):
                    if(not B[row,kn,km]): continue
                    k  = np.array([kx,ky])
                    kg = k + Gnm[n-bound:n+bound+1]
                    kg1 = kg[:,0].reshape((n,1))
                    kg2 = kg[:,1].reshape((n,1))
                    
                    exp = np.exp(i*(kg1[...,None]*X1 + kg2[...,None]*X2))
                    psi_k = np.sum(exp*Coeff[:,row,kn,km][...,None,None],0)
                    psi  += abs(psi_k)**2
                    
                    count += 1
                    if(np.sum(B) > 10):
                        if((count%int(0.1*np.sum(B))) == 0):
                            progress = int(100*count/np.sum(B))+1
                            print(progress,'%')
                            if(initiator and not self.checkEventFlags("map")):
                                initiator.progress(progress)
        
                    if(initiator and self.checkEventFlags("map")):
                        self.stop("map")
                        return                                                  # Sim has been stopped. Clear the running flag and return
                
        extent = np.array([min(x1),max(x1),min(x2),max(x2)])
        
        self.psi[energy] = {"E"      : E,
                            "X"      : np.array([x1,x2]),
                            "scale"  : scale,
                            "extent" : extent,
                            "psi"    : psi}
        
        print("Done")
        
        if(initiator):
            initiator.finish(True,psi.copy(),extent.copy(),np.array([x1,x2]).copy())
        
        self.stop("map")
        return psi.copy(),extent.copy(),np.array([x1,x2]).copy()
    
    def getLDOS(self,Exmin,Exmax,dE,pts,x0s,initiator=""):
        """
        This function computes the local density of states at given x0 values

        Parameters
        ----------
        Exmin : Lower bound of the energy range (eV)
        Exmax : Upper bound of the energy range (eV)
        dE    : Energy broadening term
        pts   : Number of points within the energy range
        x0s   : List of locations [np.array([x1,x2],np.array([x1n,x2n]),...] to 
                generate LDOS

        Returns
        -------
        LDOSs : Array of LDOS curves where each element is a curve
                corresponding to an x0
        Ex    : Energy axis to plot LDOS

        """
        self.running["LDOS"] = True
        
        N     = self.N
        X     = self.X
        bz    = self.bz
        k1,k2 = self.K
        Ek    = self.Ek
        Gnm   = self.Gnm
        Coeff = self.Coeff
        
        n     = int((1+(2*N+1)**2)/2)                                           # The dimension of the hamiltonian
        bound = int(n/2)
        
        x1,x2 = X
        X1,X2 = np.meshgrid(x1,x2)
        
        ovrfl = 0.03*(Exmax-Exmin)
        Ex = np.linspace(Exmin-ovrfl,Exmax+ovrfl,pts)
        LDOSs = np.zeros((len(x0s),len(Ex)))
        for nx0,x0 in enumerate(x0s):
            x0idx = np.argmin(abs(X - x0[...,None]),axis=1)
            print("Calculating LDOS at position",x1[x0idx[0]],',',x2[x0idx[1]])
        
            for ne,e in enumerate(Ex):
                emin = e - dE/2
                emax = e + dE/2
                B = ((Ek >= emin) & (Ek <= emax)) * bz
                ldos = 0
                for row in range(n):
                    for kn,kx in enumerate(k1):
                        for km,ky in enumerate(k2):
                            if(not B[row,kn,km]): continue
                            k  = np.array([kx,ky])
                            kg = k + Gnm[n-bound:n+bound+1]
                            kg1 = kg[:,0].reshape((n,1))
                            kg2 = kg[:,1].reshape((n,1))
                            
                            exp = np.exp(i*(kg1[...,None]*x1[x0idx[0]] + -kg2[...,None]*x2[x0idx[1]]))
                            psi_k = np.sum(exp*Coeff[:,row,kn,km][...,None,None],0)
                            
                            ldos += abs(psi_k)**2
                            
                        if(initiator and self.checkEventFlags("LDOS")):
                            self.stop("LDOS")
                            return                                              # Sim has been stopped. Clear the running flag and return
                    
                LDOSs[nx0,ne] = ldos
                if(ne%int(0.1*len(Ex)) == 0):
                    progress = int(100*ne/len(Ex))
                    print(progress,'%')
                    progressStr = "Point "+str(nx0+1)+'/'+str(len(x0s))+': '+str(progress)+' %'
                    if(initiator and not self.checkEventFlags("LDOS")):
                        initiator.progress(progressStr)
        
        LDOSs /= np.max(LDOSs)
        self.LDOS    = LDOSs
        self.LDOSEx  = Ex
        self.LDOSx0s = x0s
        
        print("Done")
        if(initiator):
            initiator.finish(True,LDOSs,Ex)
        
        self.stop("LDOS")
        return LDOSs,Ex
        
    def checkEventFlags(self,name):
        ldict = {}
        exec("running_" + name + " = global_." + name + "_running.is_set()",globals(),ldict)
        if(not ldict['running_'+name]): return True