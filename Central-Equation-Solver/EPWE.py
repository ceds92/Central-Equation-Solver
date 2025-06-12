# %%
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
                        "Sweep":False,
                        "map": False,
                        "potential": False}
    
    def stop(self,name,initiator=None):
        if(initiator): exec("global_." + name + "_running.clear()")
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
        if(initiator and self.checkEventFlags(name)): self.stop(name,initiator); return # Sim has been stopped. Clear the running flag and return
        
        Gnm,vg   = self.fourierCoefs(a,X,V,G,initiator)                         # Determine the Fourier coefficients of the potential
        if(initiator and self.checkEventFlags(name)): self.stop(name,initiator); return # Sim has been stopped. Clear the running flag and return
        
        U        = self.constructPotentialMatrix(vg,initiator)                  # Construct the potential matrix used to solve the TISE for a periodic potential
        if(initiator and self.checkEventFlags(name)): self.stop(name,initiator); return # Sim has been stopped. Clear the running flag and return
        
        Ek,Coeff = self.solveTISE(Gnm,K,U,bz,initiator)                         # Solve the TISE
        if(initiator and self.checkEventFlags(name)): self.stop(name,initiator); return # Sim has been stopped. Clear the running flag and return
        
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
        
        self.stop("main",initiator)                                             # Sim has finished. Clear the running flag and return
        
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
        
        bz = np.ones_like(K1, dtype=bool)
        for q in [-1, 0, 1]:
            for w in [-1, 0, 1]:
                if (q == 0 and w == 0): continue
                G = q * b1 + w * b2
                condition = (K1 * G[0] + K2 * G[1]) < (G[0]**2 + G[1]**2) / 2
                bz &= condition

        bz = np.array(bz > 0)
        
        print("Done2")
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
        vg = np.zeros(((2*N+1),(2*N+1)),dtype = 'complex')                     # Initialise the matrix which will contain the fourier coefficients of the potential, V                                     
        
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
            
                if(initiator and self.checkEventFlags("main")): return 0,0      # Sim has been stopped. Clear the running flag and return
        
        print("Done")
        return np.array(Gnm),vg
        
    def constructPotentialMatrix(self, vg, initiator=False):
        """
        Constructs the potential matrix for solving the TISE using plane wave expansion in 2D.

        Parameters
        ----------
        vg : ndarray
            Fourier coefficients of the periodic potential, shape (2*N + 1, 2*N + 1).
            vg[n, m] corresponds to G = (n - N)*b1 + (m - N)*b2, where n, m = 0 to 2*N.

        Returns
        -------
        U : ndarray
            Potential matrix, shape (num_G, num_G), where num_G = (2*N + 1)^2.
            U[i,j] = V_{G_i - G_j}.
        """
        print("Constructing Potential Matrix")
        N = self.N  # Assuming N is an attribute of the class
        num_G = (2*N + 1)**2  # Total number of reciprocal lattice vectors
        U = np.zeros((num_G, num_G), dtype=complex)  # Initialize U with complex type

        # Loop over all possible reciprocal lattice vector pairs
        for n_i in range(-N, N + 1):
            for m_i in range(-N, N + 1):
                # Map (n_i, m_i) to a single index i
                i = (n_i + N) * (2*N + 1) + (m_i + N)
                for n_j in range(-N, N + 1):
                    for m_j in range(-N, N + 1):
                        # Map (n_j, m_j) to a single index j
                        j = (n_j + N) * (2*N + 1) + (m_j + N)
                        # Compute G_diff components
                        k = n_i - n_j  # Difference in b1 direction
                        l = m_i - m_j  # Difference in b2 direction
                        # Check if G_diff is within the computed vg range
                        if -N <= k <= N and -N <= l <= N:
                            # Map k, l to vg indices: (n - N) = k => n = k + N
                            U[i, j] = vg[k + N, l + N]
                        else:
                            # Set to 0 for differences outside vg's range
                            U[i, j] = 0

        print("Done")
        return U

    def solveTISE(self, Gnm, K, U, bz, initiator=False):
        """
        Solves the TISE for a periodic potential in 2D using the plane wave expansion method.

        Parameters
        ----------
        Gnm : ndarray, shape ((2N+1)^2, 2)
            Combinations of reciprocal lattice vectors.
        K : list of ndarray
            k-space points [k1, k2], each of length ks.
        U : ndarray, shape ((2N+1)^2, (2N+1)^2)
            Potential matrix in the plane wave basis.
        bz : ndarray, shape (ks, ks)
            Boolean map of the first Brillouin zone (1 = inside, 0 = outside).

        Returns
        -------
        Ek : ndarray, shape ((2N+1)^2, ks, ks)
            Band energies (eigenvalues) for each k-point.
        Coeff : ndarray, shape ((2N+1)^2, (2N+1)^2, ks, ks)
            Wavefunction coefficients (eigenvectors).
        """
        if(initiator and not self.checkEventFlags("main")):
            initiator.progress("Solving TISE...0 %")
        print("Calculating Kinetic Terms and Eigen Values...")
        N = self.N
        ks = self.ks
        m_eff = self.m_eff
        num_G = (2 * N + 1) ** 2  # Correct basis size

        k1, k2 = K
        Ek = np.full((num_G, ks, ks), np.nan)  # Energies, NaN outside BZ
        Coeff = np.zeros((num_G, num_G, ks, ks), dtype=complex)  # Coefficients

        cnt = 0
        total = len(k1) * len(k2)
        for kn, kx in enumerate(k1):
            for km, ky in enumerate(k2):
                cnt += 1
                if((cnt%int(0.1*total)) == 0):
                    progress = int(100*cnt/total)
                    print(progress,'%')
                    if(initiator and not self.checkEventFlags("main")):
                        progress = "Solving TISE..." + str(progress) + '%'
                        initiator.progress(progress)

                if not bz[km, kn]:
                    continue  # Skip points outside BZ

                k = np.array([kx, ky])
                kg = k + Gnm  # Shape: (num_G, 2), all basis vectors
                diag = (h_bar**2 / (2 * m_e * m_eff)) * np.sum(kg**2, axis=1)
                H0 = np.diag(diag)
                H = H0 + U  # U must be (num_G, num_G)
                vals, vecs = np.linalg.eigh(H)
                Ek[:, kn, km] = vals
                Coeff[:, :, kn, km] = vecs
                
            if(initiator and self.checkEventFlags("main")): return 0,0          # Sim has been stopped. Clear the running flag and return

        print("Done")
        return Ek, Coeff
    
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
        C  = np.zeros_like(X1,dtype = 'complex')
        
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
                    self.stop("potential",initiator)
                    return                                                      # Sim has been stopped. Clear the running flag and return
        
        print("Done")
        
        C = np.real(C) + np.imag(C)
        extent = np.array([min(x1),max(x1),min(x2),max(x2)])
        
        self.C = {"C"      : C,
                  "X"      : np.array([x1,x2]),
                  "extent" : extent}
        
        if(initiator):
            initiator.finish(True,C,extent,np.array([x1,x2]))
        
        self.stop("potential",initiator)
        return C,extent,np.array([x1,x2])
    def getWavefunction(self, E, scale=1, initiator=None):
        """
        Computes the spatial distribution of electron density integrated over an energy window.

        Parameters
        ----------
        E : np.ndarray
            Energy window [Emin, Emax] in eV, where Emin and Emax define the range to integrate over.

        Returns
        -------
        rho : np.ndarray
            2D array of electron density (shape: (len(x1), len(x2))) in arbitrary units.
        X : np.ndarray
            Real-space grid coordinates [x1, x2] for plotting.
        """
        self.running["map"] = True
        print("Calculating Wavefunction...")

        i = 1j  # Imaginary unit

        # Extract simulation data
        N = self.N
        X = self.X
        bz = self.bz
        k1, k2 = self.K
        Ek = self.Ek  # Shape: (num_bands, len(k1), len(k2))
        Gnm = self.Gnm  # Shape: ((2*N + 1)**2, 2)
        Coeff = self.Coeff  # Shape: ((2*N + 1)**2, (2*N + 1)**2, len(k1), len(k2))

        n = (2 * N + 1) ** 2  # Number of plane waves (bands)
        Emin, Emax = E  # Unpack energy window
        energy = (Emax - Emin)/2

        # Real-space grid
        x1, x2 = X*scale
        X1, X2 = np.meshgrid(x1, x2)  # Shape: (len(x2), len(x1))
        rho = np.zeros(X1.shape, dtype=float)  # Electron density map

        print(f"Calculating wavefunction density for E = [{Emin}, {Emax}] eV")

        # Mask for states within energy window and BZ
        B = ((Ek >= Emin) & (Ek <= Emax)) & bz[np.newaxis, :, :]  # Shape: (n, len(k1), len(k2))

        # Loop over k-points and bands
        for kn, kx in enumerate(k1):
            for km, ky in enumerate(k2):
                if not np.any(B[:, kn, km]):  # Skip if no states in window at this k
                    continue
                k = np.array([kx, ky])
                kg = k + Gnm  # Shape: (n, 2)
                # Compute phase factor for all positions at once
                phase = np.exp(i * (kg[:, 0, None, None] * X1 + kg[:, 1, None, None] * X2))  # Shape: (n, len(x2), len(x1))
                for row in range(n):
                    if not B[row, kn, km]:
                        continue
                    # Wavefunction: psi = sum_G c_G * exp(i (k+G) · r)
                    psi = np.sum(Coeff[:, row, kn, km, None, None] * phase, axis=0)  # Shape: (len(x2), len(x1))
                    rho += np.abs(psi) ** 2  # Add density contribution
                    
                    if(initiator and self.checkEventFlags("map")):
                        self.stop("map",initiator)
                        return                                                  # Sim has been stopped. Clear the running flag and return
                    
        extent = np.array([min(x1),max(x1),min(x2),max(x2)])
        self.psi[energy] = {"E"      : E,
                            "X"      : np.array([x1,x2]),
                            "scale"  : scale,
                            "extent" : extent,
                            "psi"    : rho}
        
        print("Done")
        
        if(initiator):
            initiator.finish(True,rho.copy(),extent.copy(),np.array([x1,x2]).copy())
        
        self.stop("map",initiator)
        
        return rho.copy(),extent.copy(),np.array([x1,x2]).copy()
    
    def getLDOS(self, Exmin, Exmax, dE, pts, x0s, initiator=False):
        """
        Computes the local density of states at given x0 values.

        Parameters
        ----------
        Exmin : float
            Lower bound of the energy range (eV)
        Exmax : float
            Upper bound of the energy range (eV)
        dE    : float
            Energy broadening term
        pts   : int
            Number of points within the energy range
        x0s   : list of np.ndarray
            List of positions [np.array([x1,x2]), ...] to generate LDOS

        Returns
        -------
        LDOSs : np.ndarray
            Array of LDOS curves corresponding to each x0
        Ex    : np.ndarray
            Energy axis for plotting LDOS
        """
        if(not initiator): self.running["LDOS"] = True
        else: self.running[initiator.name] = True

        i = 1j  # Define imaginary unit

        N = self.N
        X = self.X
        bz = self.bz
        k1, k2 = self.K
        Ek = self.Ek
        Gnm = self.Gnm  # Shape: ((2*N + 1)**2, 2)
        Coeff = self.Coeff  # Shape: ((2*N + 1)**2, (2*N + 1)**2, len(k1), len(k2))

        n = (2*N + 1)**2  # Correct number of plane waves

        # Real-space coordinates
        x1, x2 = X
        X1, X2 = np.meshgrid(x1, x2)
        
        ovrfl = 0.03 * (Exmax - Exmin)
        Ex = np.linspace(Exmin - ovrfl, Exmax + ovrfl, pts)
        LDOSs = np.zeros((len(x0s), len(Ex)))
        
        for nx0, x0 in enumerate(x0s):
            x0idx = np.argmin(np.abs(X - x0[..., None]), axis=1)
            print("Calculating LDOS at position", x1[x0idx[0]], ',', x2[x0idx[1]])
            
            for ne, e in enumerate(Ex):
                emin = e - dE / 2
                emax = e + dE / 2
                B = ((Ek >= emin) & (Ek <= emax)) & bz[np.newaxis, :, :]
                ldos = 0
                for row in range(n):
                    for kn, kx in enumerate(k1):
                        for km, ky in enumerate(k2):
                            if not B[row, kn, km]:
                                continue
                            k = np.array([kx, ky])
                            kg = k + Gnm  # Shape: (n, 2)
                            exp = np.exp(i * (kg[:, 0] * x1[x0idx[0]] + kg[:, 1] * x2[x0idx[1]]))
                            psi_k = np.sum(Coeff[:, row, kn, km] * exp)
                            ldos += np.abs(psi_k)**2
                LDOSs[nx0, ne] = ldos
        
        LDOSs /= np.max(LDOSs)
        self.LDOS = LDOSs
        self.LDOSEx = Ex
        self.LDOSx0s = x0s
        
        print("Done")
        if(initiator):
            initiator.finish(True,LDOSs,Ex)
            self.stop(initiator.name,initiator)
        else:
            self.stop("LDOS",initiator)
        return LDOSs, Ex

    def checkEventFlags(self,name):
        ldict = {}
        exec("running_" + name + " = global_." + name + "_running.is_set()",globals(),ldict)
        if(not ldict['running_'+name]): return True
