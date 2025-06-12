# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 06:50:05 2022

@author: jced0001
"""

from Panel import Panel
import customtkinter as ctk
import numpy as np
import numpy.linalg as npl
import matplotlib.pyplot as plt
from   tkinter import filedialog
import pickle

class BandStructure2D(Panel):
    numBands  = 1
    bandOrder = False
    # Inset
    showBZ = False
    insetPos = np.array([0.65, 0.65, 0.3, 0.3]);
    insetColours = ['black','white'] + plt.rcParams['axes.prop_cycle'].by_key()['color']
    ###########################################################################
    # Constructor
    ###########################################################################
    def __init__(self, master, width, height, dpi, mainPanel):
        super().__init__(master, width, height, dpi, mainPanel=mainPanel)
        self.buttons()
        
    ###########################################################################
    # Panel
    ###########################################################################
    def buttons(self):
        self.btn = {
            "Add":    ctk.CTkButton(self.master, text="Add band",       command=self.addBand),
            "Remove": ctk.CTkButton(self.master, text="Remove band",    command=self.removeBand),
            "BZ":     ctk.CTkButton(self.master, text="Inset BZ",       command=self.toggleBZ),
            "BOrder": ctk.CTkComboBox(self.master,values=["Band order"],command=self.orderBands),        # Dropdown to change order of bands
            "Export": ctk.CTkButton(self.master, text="Exp .pk",        command=self.exportPickle),
            "PNG":    ctk.CTkButton(self.master, text="Exp PNG",        command=super().exportPNG), # Export the canvas to png
            "Close":  ctk.CTkButton(self.master, text="Close",          command=self.destroy)
            }
        
    def buttonHelp(self):
        overlayValues = ["Band Order","Eigen energy","Tracked bands"]
        self.btn['BOrder'].configure(values=overlayValues,fg_color=['#3B8ED0', '#1F6AA5'])

        helpStr = "Close this panel"
        self.btn['Close'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Add a band to the plot"
        self.btn['Add'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Remove a band from the plot"
        self.btn['Remove'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Export the main panel plot as a png"
        self.btn['PNG'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
    ###########################################################################
    # Update and Plotting
    ###########################################################################
    def update(self,fromThread=False):
        if(not self.active): return
        
        self.ax.cla()                                                           # Clear the axis
        self.ax.set_position([0.09, 0.07, 0.87, 0.9])                           # Make it take up the whole canvas
        self.showBandStructure()
        self.addInset(fromThread)
        self.canvas.figure = self.fig                                           # Assign the figure to the canvas
        self.canvas.draw()   
    
    def showBandStructure(self):
        if(self.mainPanel.sim and self.mainPanel.sim.valid):
            self.path_indices_in_bz,self.M,self.K,self.Gamma,self.k,self.path_kx,self.path_ky = self.getHighSymPath()
            self.highSymPoints = [0,self.M,self.K,self.Gamma]
            self.highSymLabels = [r'$\Gamma$','M','K',r'$\Gamma$']

            if(not self.bandOrder):
                Ek = np.array([self.mainPanel.sim.Ek[:, j, i] for i, j in self.path_indices_in_bz])
                for b in range(self.numBands):
                    Ekb = Ek[:,b]
                    self.ax.plot(self.k,Ekb,'o',linewidth=3,markersize=3)
                    self.ax.set_xticks(self.highSymPoints)
                    self.ax.set_xticklabels(self.highSymLabels)
            else:
                trackedEnergies = self.getTrackedBands()
                for b in range(self.numBands):
                    self.ax.plot(self.k, trackedEnergies[:, b],linewidth=3)
                self.ax.set_xticks(self.highSymPoints)
                self.ax.set_xticklabels(self.highSymLabels)
            
    def addInset(self,fromThread):
        self.btn["BZ"].configure(fg_color=['#3B8ED0', '#1F6AA5'])
        if(fromThread):                                                         # Workaround for bug where a non-main thread tries to update the BZ inset.
            if(self.showBZ):                                                    # If the BZ inset is currently showing and the update was casued by a non-main thread
                self.btn["BZ"].configure(fg_color='Red')                        # Make the button red to warn the user that the BZ inset may not reflect the latest changes
            return
        if(len(self.fig.axes) > 1):                                             # Remove the previous inset from the sxm figure if there was one
            self.fig.axes[1].cla()
            self.fig.axes[1].remove()
        if(not self.showBZ): return
        print("ShowingBZ")
        insetFig = plt.figure(figsize=(self.width/self.dpi,self.height/self.dpi),dpi=self.dpi)
        insetAx = insetFig.add_subplot(111)                                     # Take the axes
        k1,k2 = self.mainPanel.sim.K
        insetAx.imshow(self.mainPanel.sim.bz, origin='lower', extent=[k1[0], k1[-1], k2[0], k2[-1]])
        insetAx.axis("off")
        insetAx.plot(self.path_kx,self.path_ky)
        insetAx.remove()                                                        # Remove from the temporary figure
        insetAx.figure = self.fig                                               # Point the axes to the main figure
        
        insetAx.set_title("")                                                   # Get rid of the title
        
        self.fig.axes.append(insetAx)                                           # Append the axis to the sxm figure
        self.fig.add_axes(insetAx)
        
        plt.close(insetFig)
            
        self.fig.axes[1].set_position(self.insetPos)                            # Adjust the size and position of the inset
        
    ###########################################################################
    # Misc
    ###########################################################################
    def load(self):
        pass
    
    def addBand(self):
        if(self.mainPanel.sim and self.mainPanel.sim.valid):
            self.numBands += 1
            self.update()
            
    def removeBand(self):
        if(self.numBands == 1): return
        if(self.mainPanel.sim and self.mainPanel.sim.valid):
            self.numBands -= 1
            self.update()
    
    def toggleBZ(self):
        self.showBZ = not self.showBZ
        self.update()
        
    def orderBands(self,option):
        if(option == "Eigen energy"):  self.bandOrder = False
        if(option == "Tracked bands"): self.bandOrder = True
        
        self.btn['BOrder'].set("Band order")

        self.update()

    def find_closest_index(self, k_point, k1, k2):
        K1, K2 = np.meshgrid(k1, k2)
        dist = (K1 - k_point[0])**2 + (K2 - k_point[1])**2
        return np.unravel_index(np.argmin(dist), K1.shape)

    def get_path_indices(self,start_idx, end_idx, num_points=40):
        start = np.array(start_idx)
        end = np.array(end_idx)
        path = []
        for t in np.linspace(0, 1, num_points, endpoint=True):
            idx = np.round(start + t * (end - start)).astype(int)
            # Ensure indices stay within bounds
            idx[0] = np.clip(idx[0], 0, self.mainPanel.sim.ks-1)
            idx[1] = np.clip(idx[1], 0, self.mainPanel.sim.ks-1)
            path.append(tuple(idx))
        return list(dict.fromkeys(path))  # Remove duplicates, preserve order

    def getHighSymPath(self):
        N     = self.mainPanel.sim.N
        G1,G2 = self.mainPanel.sim.G
        k1,k2 = self.mainPanel.sim.K
        bz    = self.mainPanel.sim.bz
        b1,b2 = G1[N+1], G2[N+1]

        # For a perfect hexagonal lattice
        gamma_point = np.array([0,0])
        m_point     = b1/2
        k_point     = (2*b1 + b2)/3

        # Convert coordinates to indexes
        gamma_idx = self.find_closest_index(gamma_point, k1, k2)
        m_idx = self.find_closest_index(m_point, k1, k2)
        k_idx = self.find_closest_index(k_point, k1, k2)

        path_gamma_m = self.get_path_indices(gamma_idx, m_idx)
        path_m_k = self.get_path_indices(m_idx, k_idx)
        path_k_gamma = self.get_path_indices(k_idx, gamma_idx)

        path_gamma_m_in_bz  = [(i, j) for i, j in path_gamma_m if bz[i, j]]
        path_m_k_in_bz      = [(i, j) for i, j in path_m_k if bz[i, j]]
        path_k_gamma_in_bz  = [(i, j) for i, j in path_k_gamma if bz[i, j]]
        path_indices_in_bz = path_gamma_m_in_bz + path_m_k_in_bz + path_k_gamma_in_bz

        path_kx = [k1[j] for i, j in path_indices_in_bz]
        path_ky = [k2[i] for i, j in path_indices_in_bz]

        k = [0.0]
        path_k = np.array([path_kx, path_ky]).T
        for n in range(1, len(path_k)):
            dk = np.linalg.norm(path_k[n] - path_k[n-1])
            k.append(k[-1] + dk)
        k = np.array(k)

        M = len(path_gamma_m_in_bz) - 1
        K = len(path_m_k_in_bz) + M
        Gamma = len(path_indices_in_bz) - 1
        
        return path_indices_in_bz, k[M], k[K], k[Gamma], k, path_kx, path_ky
    
    def getTrackedBands(self):
        # # coeffs_along_path: Shape (n_points, num_G, num_G), eigenvectors in columns
        coeffs_along_path = np.array([self.mainPanel.sim.Coeff[:, :, j, i] for i, j in self.path_indices_in_bz])
        energies_along_path = np.array([self.mainPanel.sim.Ek[:, j, i] for i, j in self.path_indices_in_bz])
        # Number of bands
        num_bands = (2 * self.mainPanel.sim.N + 1) ** 2  # e.g. 529 for N=11
        n_points = len(self.path_indices_in_bz)

        # Initialize tracked energies
        tracked_energies = np.zeros((n_points, num_bands))
        tracked_energies[0] = energies_along_path[0]  # Start with first point unsorted

        # Track bands using overlap
        band_indices = np.arange(num_bands)  # Initial order at k=0
        for n in range(1, n_points):
            prev_coeffs = coeffs_along_path[n-1, :, :]  # Shape (num_G, num_G)
            curr_coeffs = coeffs_along_path[n, :, :]    # Shape (num_G, num_G)
            
            # Compute overlap matrix between eigenvectors
            overlaps = np.abs(np.dot(prev_coeffs.T.conj(), curr_coeffs))  # Shape (num_G, num_G)
            
            # Match bands by maximum overlap
            new_indices = np.zeros(num_bands, dtype=int)
            used = set()
            for i in range(num_bands):
                row = overlaps[band_indices[i]]
                j = np.argmax(row)
                while j in used:  # Avoid reusing indices
                    row[j] = -1
                    j = np.argmax(row)
                new_indices[i] = j
                used.add(j)
            
            band_indices = new_indices
            tracked_energies[n] = energies_along_path[n][band_indices]

        return tracked_energies
    
    def exportPickle(self):
        if(not len(self.highSymLabels)): return
        
        default = 'BS.pk'
        path = filedialog.asksaveasfilename(title="Save as",initialfile=default)
        if(not path.endswith('.pk')): path += '.pk'
        
        BS = []
        Ek = np.array([self.mainPanel.sim.Ek[:, j, i] for i, j in self.path_indices_in_bz])
        for b in range(self.numBands):
            Ekb = Ek[:,b]
            BS.append(Ekb)

        trackedEnergies = self.getTrackedBands()
        BS_tracked = []
        for b in range(self.numBands):
            BS_tracked.append(trackedEnergies[:, b])

        pkldict = {'highSymLabels' : self.highSymLabels,
                   'highSymPoints' : self.highSymPoints,
                   'BS': BS,
                   'BS_tracked' : BS_tracked,
                   'k' : self.k,
                   'Ek_2D': Ek,
                   'Ek': self.mainPanel.sim.Ek}
        
        pickle.dump(pkldict,open(path,'wb'))