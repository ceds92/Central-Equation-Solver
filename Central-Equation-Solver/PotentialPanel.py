# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 08:00:45 2022

@author: jced0001
"""

from Panel import Panel
import customtkinter as ctk

class PotentialPanel(Panel):
    scaleBar = True
    ###########################################################################
    # Constructor
    ###########################################################################
    def __init__(self, master, width, height, dpi, mainPanel):
        super().__init__(master, width, height, dpi, mainPanel=mainPanel)
        self.buttons()
        self.rebuiltV = []
        
    ###########################################################################
    # Panel
    ###########################################################################
    def buttons(self):
        self.btn = {
            "Refresh": ctk.CTkButton(self.master, text="Refresh",       command=self.refresh),
            "cmap":    ctk.CTkButton(self.master, text="viridis",       command=super()._cmap),     # Button to cycle through colour maps
            "PNG":     ctk.CTkButton(self.master, text="Exp PNG",       command=super().exportPNG), # Export the canvas to png
            "Overlay": ctk.CTkComboBox(self.master,values=["Overlay"],  command=self.overlay),      # Dropdown to change overlay display
            "Close":   ctk.CTkButton(self.master, text="Close",         command=self.destroy)
            }
    
        overlayValues = ["Overlay","Scale Bar"]
        self.btn['Overlay'].configure(values=overlayValues,fg_color=['#3B8ED0', '#1F6AA5'])
        
    def buttonHelp(self):
        helpStr = "Reconstruct the potential"
        self.btn['Refresh'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Change the colour map"
        self.btn['cmap'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Export the main panel plot as a png"
        self.btn['PNG'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Close this panel"
        self.btn['Close'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
    
    def refresh(self):
        if(self.mainPanel.sim and self.mainPanel.sim.valid):
            self.rebuiltV,self.extent,self.X = self.mainPanel.sim.rebuildPotential(scale=2)
            self.update()
    
    ###########################################################################
    # Update and Plotting
    ###########################################################################
    def update(self):
        if(not self.active): return
        
        self.ax.cla()                                                           # Clear the axis
        self.ax.set_position([0, 0, 1, 1])                                      # Make it take up the whole canvas
        self.showV()
        self.updateOverlay()
        self.ax.axis('off')                                                     # Hide the axis
        self.canvas.figure = self.fig                                           # Assign the figure to the canvas
        self.canvas.draw()
        
    def showV(self):
        if(not len(self.rebuiltV)): return
        cmap = self.cmaps[self.cmap][1]
        self.ax.imshow(self.rebuiltV,extent=self.extent,cmap=cmap())
    
    def updateOverlay(self):
        if(not len(self.rebuiltV)): return
        UC = self.mainPanel.getUC(self.X)
        for line in UC:
            self.ax.plot(line[0],line[1],c='red',linewidth=1)
            
        if(self.scaleBar): super().addPlotScalebar()                            # Add a scale bar to the plot
            
    ###########################################################################
    # Misc
    ###########################################################################
    def load(self):
        if(not len(self.mainPanel.sim.C.keys())): return
        self.rebuiltV = self.mainPanel.sim.C["C"]
        self.extent   = self.mainPanel.sim.C["extent"]
        self.X        = self.mainPanel.sim.C["X"]
        self.update()
        
    def overlay(self,option):
        if(option == "Scale Bar"):  self.toggleScaleBar()
        self.btn['Overlay'].set("Overlay")
    
    def toggleScaleBar(self):
        self.scaleBar = not self.scaleBar
        self.update()