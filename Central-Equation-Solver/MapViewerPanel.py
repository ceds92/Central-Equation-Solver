# -*- coding: utf-8 -*-
"""
Created on Sat Nov 26 11:29:07 2022

@author: jced0001
"""

from Panel import Panel
import customtkinter as ctk
import numpy as np

class MapViewerPanel(Panel):
    scaleBar = True
    plotCaption = True
    
    ###########################################################################
    # Constructor
    ###########################################################################
    def __init__(self, master, width, height, dpi, mainPanel):
        super().__init__(master, width, height, dpi, mainPanel=mainPanel,length=8,btnSize=2)
        self.buttons()
        self.map = 0
        self.im  = []
        self.init = True

    ###########################################################################
    # Panel
    ###########################################################################
    def buttons(self):
        self.btn = {
            "Prev":     ctk.CTkButton(self.master, text="Prev",         command=lambda : self.cycleMap(-1)),
            "Next":     ctk.CTkButton(self.master, text="Next",         command=lambda : self.cycleMap(1)),
            "cmap":     ctk.CTkButton(self.master, text="viridis",      command=super()._cmap), # Button to cycle through colour maps
            "Overlay":  ctk.CTkComboBox(self.master,values=["Overlay"], command=self.overlay),  # Dropdown to change overlay display
            "PNG":      ctk.CTkButton(self.master, text="Exp PNG",      command=super().exportPNG),   # Export the canvas to png
            "Close":    ctk.CTkButton(self.master, text="Close",        command=self.destroy)
            }
    
        overlayValues = ["Overlay","Caption","Scale Bar"]
        self.btn['Overlay'].configure(values=overlayValues,fg_color=['#3B8ED0', '#1F6AA5'])
        
    def buttonHelp(self):
        helpStr = "Close this panel"
        self.btn['Close'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Change the colour map"
        self.btn['cmap'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Show previous map"
        self.btn['Prev'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Show next map"
        self.btn['Next'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Show/hide overlay features"
        self.btn['Overlay'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Export the panel figure as a png"
        self.btn['PNG'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))

    ###########################################################################
    # Update and Plotting
    ###########################################################################
    def update(self):
        if(not self.active): return
        
        self.ax.cla()                                                           # Clear the axis
        self.ax.set_position([0, 0, 1, 1])                                      # Make it take up the whole canvas
        self.ax.axis('off')                                                     # Hide the axis
        self.showMap()
        self.updateOverlay()
        self.canvas.figure = self.fig                                           # Assign the figure to the canvas
        self.canvas.draw()
    
    def showMap(self):
        if(not (self.mainPanel.sim and self.mainPanel.sim.valid)): return
        if(not len(self.mainPanel.sim.psi)): return
        if(not len(self.im)):
            if(not self.cycleMap(0)): return
            
        cmap = self.cmaps[self.cmap][1]
        self.ax.imshow(self.im,extent=self.extent,cmap=cmap())
        
    def updateOverlay(self):
        if(not len(self.im)): return
        left, right = self.ax.get_xlim();                                       # Remember plot extent
        bottom, top = self.ax.get_ylim();                                       # Remember plot extent
        
        UC = self.mainPanel.getUC(self.X)
        for line in UC:
            self.ax.plot(line[0],line[1],c='red',linewidth=1)
        
        if(self.scaleBar): super().addPlotScalebar()                            # Add a scale bar to the plot
        if(self.plotCaption):                                                   # Caption the image with Vbias and Iset
            E  = int(1000*np.sum(self.E/2))/1000
            dE = int(1000*(self.E[1] - self.E[0]))/1000
            plotCaption  = r'$E = $'  + str(E)    + ' eV'                       # Show bias in a box in the top left of the image
            plotCaption += r'; $dE = $' + str(dE) + ' eV'                       # Show setpoint current in top left
            rx = abs(self.extent[1] - self.extent[0])
            ry = abs(self.extent[3] - self.extent[2])
            pos = [0.025*rx + self.extent[0],self.extent[3] - 0.1*ry]           # Put the caption in the top left
            super().addPlotCaption(plotCaption, pos)
        
        if(self.mainPanel.ldosPanel.active):
            for x0 in self.mainPanel.ldosPanel.x0s:
                self.ax.plot(x0[0],x0[1],'x',markersize=12)
                
        self.ax.set_xlim((left,right)); self.ax.set_ylim((bottom,top))          # Put back extent
    
    def cycleMap(self,m):
        if(not (self.mainPanel.sim and self.mainPanel.sim.valid)):
            self.updateHelpLabel("Error: Run a simulation then generate Maps to view them here.")
            return False
        if(not len(self.mainPanel.sim.psi)):
            self.updateHelpLabel("Error: Generate maps in the Maps panel first.")
            return False
        
        numMaps = len(list(self.mainPanel.sim.psi.keys()))
        
        if((self.map + m >= 0) and (self.map + m < numMaps)):
            self.map += m
        
        if(self.map > numMaps): self.map = 0
        
        energies = list(self.mainPanel.sim.psi.keys())
        energies.sort()
        
        self.E  = self.mainPanel.sim.psi[energies[self.map]]['E']
        self.im = self.mainPanel.sim.psi[energies[self.map]]['psi']
        self.X  = self.mainPanel.sim.psi[energies[self.map]]['X']
        self.extent = self.mainPanel.sim.psi[energies[self.map]]['extent']
        
        if(m == 0): return True
        
        self.update()
        
    ###########################################################################
    # Misc
    ###########################################################################
    def load(self):
        pass
    
    def overlay(self,option):
        if(not self.init):
            self.btn['Overlay'].set("Overlay")
            return
        if(option == "Caption"):    self.toggleCaption()
        if(option == "Scale Bar"):  self.toggleScaleBar()
        self.btn['Overlay'].set("Overlay")
        
    def toggleCaption(self):
        self.plotCaption = not self.plotCaption
        self.update()
    
    def toggleScaleBar(self):
        self.scaleBar = not self.scaleBar
        self.update()