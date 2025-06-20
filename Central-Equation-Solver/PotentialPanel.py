# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 08:00:45 2022

@author: jced0001
"""

from Panel import Panel
import customtkinter as ctk

class PotentialPanel(Panel):
    scaleBar = True
    init = False
    ###########################################################################
    # Constructor
    ###########################################################################
    def __init__(self, master, width, height, dpi, mainPanel):
        super().__init__(master, width, height, dpi, mainPanel=mainPanel)
        super().initGlobs(name="potential")
        self.buttons()
        self.rebuiltV = []
        self.init = 1
        
    ###########################################################################
    # Panel
    ###########################################################################
    def buttons(self):
        self.btn = {
            "Rebuild": ctk.CTkButton(self.master, text="Rebuild",       command=self.refresh),
            "cmap":    ctk.CTkButton(self.master, text="viridis",       command=super()._cmap),     # Button to cycle through colour maps
            "Export":  ctk.CTkComboBox(self.master,values=["Export"],   command=self.export),       # Dropdown to export the figure
            "Overlay": ctk.CTkComboBox(self.master,values=["Overlay"],  command=self.overlay),      # Dropdown to change overlay display
            "Close":   ctk.CTkButton(self.master, text="Close",         command=self.destroy)
            }
    
        exportValues = ["Export","PNG","Pickle"]
        self.btn['Export'].configure(values=exportValues,fg_color=['#3B8ED0', '#1F6AA5'])
        
        overlayValues = ["Overlay","Scale Bar"]
        self.btn['Overlay'].configure(values=overlayValues,fg_color=['#3B8ED0', '#1F6AA5'])
        
    def buttonHelp(self):
        helpStr = "Rebuild the potential from the N coefficients"
        self.btn['Rebuild'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Change the colour map"
        self.btn['cmap'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Export figure"
        self.btn['Export'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Show/hide overlay features"
        self.btn['Overlay'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Close this panel"
        self.btn['Close'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
    
    def refresh(self):
        if(not(self.mainPanel.sim and self.mainPanel.sim.valid)):
            self.updateHelpLabel("Error: simulation not valid. Need to rerun before rebuilding potential.")
            return
        
        if(self.mainPanel.sim.running["main"]):
            self.updateHelpLabel("Error: Wait for the simulation to finish running.")
            return
        
        
        if(self.mainPanel.sim.running["potential"]):
            super().stop()
            while(self.mainPanel.sim.running["potential"]):
                print("waiting to stop")
            self.btn['Rebuild'].configure(text="Rebuild")
            self.updateHelpLabel("Stopped!")
            return
        
        self.btn['Rebuild'].configure(text="STOP")
        self.updateHelpLabel("Rebuilding...")
        
        func = lambda : self.mainPanel.sim.rebuildPotential(scale=2,initiator=self)
        super().threadTask(func)
        
    def finish(self,success,rebuiltV="",extent="",X=""):
        if(not success):
            self.updateHelpLabel("Error: Cannot run while another process is running.")
            return
        
        self.btn['Rebuild'].configure(text="Rebuild")
        
        self.X        = X
        self.extent   = extent
        self.rebuiltV = rebuiltV
        self.updateHelpLabel("Done!")
        self.update()
    
    def progress(self,progress):
        self.updateHelpLabel("Rebuilding... " + str(progress) + " %")
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
    
    def export(self,option):
        if(option == "PNG"): super().exportPNG()
        if(option == "Pickle"): self.exportPickle()
        self.btn['Export'].set("Export")
    
    def exportPickle(self):
        if(not len(self.rebuiltV)): return
        
        pklDict = {"rebuiltV": self.rebuiltV,
                   "extent"  : self.extent}
        
        super().exportPickle(pklDict=pklDict,initialfile="rebuiltV")
        
    def toggleScaleBar(self):
        self.scaleBar = not self.scaleBar
        self.update()