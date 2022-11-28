# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 06:50:05 2022

@author: jced0001
"""

from Panel import Panel
import customtkinter as ctk
import numpy as np
import numpy.linalg as npl

class BandStructure2D(Panel):
    numBands = 1
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
            "Close":  ctk.CTkButton(self.master, text="Close",        command=self.destroy)
            }
        
    def buttonHelp(self):
        helpStr = "Close this panel"
        self.btn['Close'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Add a band to the plot"
        self.btn['Add'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Remove a band from the plot"
        self.btn['Remove'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
    ###########################################################################
    # Update and Plotting
    ###########################################################################
    def update(self):
        if(not self.active): return
        
        self.ax.cla()                                                           # Clear the axis
        self.ax.set_position([0.09, 0.07, 0.87, 0.9])                           # Make it take up the whole canvas
        self.showBandStructure()
        self.canvas.figure = self.fig                                           # Assign the figure to the canvas
        self.canvas.draw()   
    
    def showBandStructure(self):
        if(self.mainPanel.sim and self.mainPanel.sim.valid):
            a   = self.mainPanel.sim.a
            lim = np.pi/npl.norm(a,axis=1)
            k1,k2 = self.mainPanel.sim.K
            ny = int(len(k2)/2)
            for b in range(self.numBands):
                Ekb = self.mainPanel.sim.Ek[b]
                self.ax.plot(k1/lim[0], Ekb[:,ny])
                self.ax.set_xlim([-1.5,1.5])
    
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
        