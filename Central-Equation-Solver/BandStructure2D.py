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
    numBands = 1
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
            "Add":    ctk.CTkButton(self.master, text="Add band",    command=self.addBand),
            "Remove": ctk.CTkButton(self.master, text="Remove band", command=self.removeBand),
            "BZ":     ctk.CTkButton(self.master, text="Inset BZ",    command=self.toggleBZ),
            "Export": ctk.CTkButton(self.master, text="Exp .pk",     command=self.exportPickle),
            "PNG":    ctk.CTkButton(self.master, text="Exp PNG",     command=super().exportPNG), # Export the canvas to png
            "Close":  ctk.CTkButton(self.master, text="Close",       command=self.destroy)
            }
        
    def buttonHelp(self):
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
    def update(self):
        if(not self.active): return
        
        self.ax.cla()                                                           # Clear the axis
        self.ax.set_position([0.09, 0.07, 0.87, 0.9])                           # Make it take up the whole canvas
        self.showBandStructure()
        if(len(self.fig.axes) > 1):                                             # Remove the previous inset from the sxm figure if there was one
            self.fig.axes[1].cla()
            self.fig.axes[1].remove()
        if(self.showBZ): self.addInset()
        self.canvas.figure = self.fig                                           # Assign the figure to the canvas
        self.canvas.draw()   
    
    def showBandStructure(self):
        if(self.mainPanel.sim and self.mainPanel.sim.valid):
            # a   = self.mainPanel.sim.a
            # lim = np.pi/npl.norm(a,axis=1)
            # k1,k2 = self.mainPanel.sim.K
            # ny = int(len(k2)/2)
            self.gamma_m_k_gamma,M,K,G = self.getHighSymPath()
            self.highSymPoints = [0,M,K,G]
            self.highSymLabels = [r'$\Gamma$','M','K',r'$\Gamma$']
            for b in range(self.numBands):
                Ekb = self.mainPanel.sim.Ek[b]
                # self.ax.plot(k1/lim[0], Ekb[:,ny])
                # self.ax.set_xlim([-1.5,1.5])
                self.ax.plot(Ekb[self.gamma_m_k_gamma[1],self.gamma_m_k_gamma[0]])
                self.ax.set_xticks(self.highSymPoints)
                self.ax.set_xticklabels(self.highSymLabels)
            
    def addInset(self):
        insetFig = plt.figure(figsize=(self.width/self.dpi,self.height/self.dpi),dpi=self.dpi)
        insetAx = insetFig.add_subplot(111)                                     # Take the axes
        insetAx.imshow(self.mainPanel.sim.bz)
        insetAx.axis("off")
        insetAx.plot(*self.gamma_m_k_gamma)
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
        
    def getHighSymPath(self):
        N = self.mainPanel.sim.N
        ks = self.mainPanel.sim.ks
        b1 = self.mainPanel.sim.G[0][N+1]
        b2 = self.mainPanel.sim.G[1][N+1]
        k1,k2 = self.mainPanel.sim.K
        
        gammaPoint = 0*b1
        mPoint = b1/2
        mpPoint = (b1+b2)/2
        
        m_m = -mPoint[0]/mPoint[1]
        c_m = mPoint[1] - m_m*mPoint[0]
        # y_m = m_m*k1 + c_m
        
        m_mp = -mpPoint[0]/mpPoint[1]
        c_mp = mpPoint[1] - m_mp*mpPoint[0]
        # y_mp = m_mp*k1 + c_mp
        
        x = (c_m - c_mp)/(m_mp - m_m)
        y = m_m*x + c_m
        
        if(abs(m_mp) == np.inf):
            x = mpPoint[0]
            y = m_m*x + c_m
        
        if(abs(m_m) == np.inf):
            x = mPoint[0]
            y = m_mp*x + c_mp
            
        kPoint = [x,y]
        
        gammaIndex  = np.array([np.argmin(abs(k1 - gammaPoint[0])),np.argmin(abs(k2 - gammaPoint[1]))])
        mIndex      = np.array([np.argmin(abs(k1 - mPoint[0])),np.argmin(abs(k2 - mPoint[1]))])
        mpIndex     = np.array([np.argmin(abs(k1 - mpPoint[0])),np.argmin(abs(k2 - mpPoint[1]))])
        kIndex      = np.array([np.argmin(abs(k1 - kPoint[0])),np.argmin(abs(k2 - kPoint[1]))])
        
        mIndex[0] -= 1
        mIndex[1] += 1
        mpIndex[0] -= 1
        kIndex[0] -= 1
        kIndex[1] += 1
        
        xx = np.arange(ks)
        m = (mIndex[1] - gammaIndex[1])/(mIndex[0] - gammaIndex[0])
        c = mIndex[1] - m*mIndex[0]
        gamma_to_m = (m*xx + c).astype(int)
        mask = (xx >= gammaIndex[0]) & (xx <= mIndex[0])
        gamma_to_m = np.array([xx[mask],gamma_to_m[mask]])
        M = len(gamma_to_m[0]) - 1
        
        m = (mIndex[1] - kIndex[1])/(mIndex[0] - kIndex[0])
        c = mIndex[1] - m*mIndex[0]
        m_to_k = (m*xx + c).astype(int)
        mask = (xx >= mIndex[0]) & (xx <= kIndex[0])
        m_to_k = np.array([xx[mask],m_to_k[mask]])
        K = M + len(m_to_k[0])
        
        m = (gammaIndex[1] - kIndex[1])/(gammaIndex[0] - kIndex[0])
        c = gammaIndex[1] - m*gammaIndex[0]
        k_to_gamma = (m*xx + c).astype(int)
        mask = (xx >= gammaIndex[0]) & (xx <= kIndex[0])
        k_to_gamma = np.array([xx[mask],k_to_gamma[mask]])
        k_to_gamma = np.fliplr(k_to_gamma)
        G = K + len(k_to_gamma[0])
        
        gamma_m_k_gamma = np.concatenate((gamma_to_m,m_to_k,k_to_gamma),axis=1)
        
        return gamma_m_k_gamma,M,K,G
    
    def exportPickle(self):
        if(not len(self.highSymLabels)): return
        
        default = 'BS.pk'
        path = filedialog.asksaveasfilename(title="Save as",initialfile=default)
        if(not path.endswith('.pk')): path += '.pk'
        
        BS = []
        for b in range(self.numBands):
            Ekb = self.mainPanel.sim.Ek[b]
            BS.append(Ekb[self.gamma_m_k_gamma[1],self.gamma_m_k_gamma[0]])
            
        pkldict = {'highSymLabels' : self.highSymLabels,
                   'highSymPoints' : self.highSymPoints,
                   'BS'  : BS}
        
        pickle.dump(pkldict,open(path,'wb'))