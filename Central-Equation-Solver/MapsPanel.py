# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 20:19:38 2022

@author: jced0001
"""

from Panel import Panel
import customtkinter as ctk
import numpy as np

class MapsPanel(Panel):
    scaleBar = True
    plotCaption = True
    band = 0
    ###########################################################################
    # Constructor
    ###########################################################################
    def __init__(self, master, width, height, dpi, mainPanel):
        super().__init__(master, width, height, dpi, mainPanel=mainPanel,length=8,btnSize=2)
        super().initGlobs(name="map")
        self.buttons()
        self.forms = {}
        self.psi = []
        
    ###########################################################################
    # Panel
    ###########################################################################
    def buttons(self):
        self.btn = {
            "Prev":     ctk.CTkButton(self.master, text="Prev band",    command=lambda : self.cycleBand(-1)),
            "Reset":    ctk.CTkButton(self.master, text="Reset band",   command=lambda : self.cycleBand(0)),
            "Next":     ctk.CTkButton(self.master, text="Next band",    command=lambda : self.cycleBand(1)),
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
        
        helpStr = "Adjust preset energies to the previous band"
        self.btn['Prev'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Adjust preset energies to the next band"
        self.btn['Next'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Adjust preset energies to the current band"
        self.btn['Reset'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Show/hide overlay features"
        self.btn['Overlay'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Export the panel figure as a png"
        self.btn['PNG'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
    
    def special(self):
        params = []
        params.append(['Emin','Emax','Scale'])
        self.buildForm(name="Reference", params=params)
        self.showForm(name="Reference")
        
    def removeSpecial(self):
        self.hideForm()
    
    def buildForm(self,name,params,row=7):
        self.forms[name] = {"labels"  : [],
                            "entries" : [],
                            "buttons" : []}
        idr = 0;
        for idr,p in enumerate(params):
            for idp,param in enumerate(p):
                self.forms[name]['labels'].append([ctk.CTkLabel(self.master, text=param),row+idr,self.pos + 2*idp])
                self.forms[name]['entries'].append([ctk.CTkEntry(self.master),row+idr,self.pos+2*idp+1])
        
        idr += 1
        self.forms[name]['buttons'] = []
        self.forms[name]['buttons'].append([ctk.CTkButton(self.master, text="RUN", command=lambda n=name: self.submitForm(n)),row+idr,self.pos])
        
    ###########################################################################
    # Update and Plotting
    ###########################################################################
    def update(self):
        if(not self.active): return
        
        self.ax.cla()                                                           # Clear the axis
        self.ax.set_position([0, 0, 1, 1])                                      # Make it take up the whole canvas
        self.ax.axis('off')                                                     # Hide the axis
        self.showPsi()
        self.updateOverlay()
        self.canvas.figure = self.fig                                           # Assign the figure to the canvas
        self.canvas.draw()
    
    def showPsi(self):
        if(not len(self.psi)): return
        cmap = self.cmaps[self.cmap][1]
        self.ax.imshow(self.psi,extent=self.extent,cmap=cmap())
            
    def updateOverlay(self):
        if(not len(self.psi)): return
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
    ###########################################################################
    # Form Actions
    ###########################################################################
    def showForm(self,name):
        for l in self.forms[name]['labels']:
            l[0].grid(row=l[1],column=l[2],columnspan=1)
            l[0].configure(width=int(self.width/self.length),height=27)
            
        for idx,e in enumerate(self.forms[name]['entries']):
            e[0].grid(row=e[1],column=e[2],columnspan=1)
            e[0].configure(width=int(self.width/self.length),height=27)
            
            e[0].delete(0,ctk.END)
            e[0].insert(0,[0,0.1,2][idx])
        
        for idx,b in enumerate(self.forms[name]['buttons']):
            b[0].grid(row=b[1],column=b[2] + 2*idx,columnspan=2)
            b[0].configure(width=int(self.width/self.length),height=27)
            
    def hideForm(self,name=""):
        names = [name]
        if(not names[0]): names = self.forms.keys()
        
        for name in names:
            for l in self.forms[name]['labels']:
                l[0].grid_forget()
                
            for e in self.forms[name]['entries']:
                e[0].grid_forget()
            
            for b in self.forms[name]['buttons']:
                b[0].grid_forget()
                
    def submitForm(self,name):
        self.updateHelpLabel("Calculating...")
        params = []
        for e in self.forms[name]['entries']:
            try:
                params.append(np.float32(e[0].get()))
            except:
                self.updateHelpLabel("Error in form: All values must be numeric.")
                return
            
        if(not(self.mainPanel.sim and self.mainPanel.sim.valid)):
            self.updateHelpLabel("Error: simulation not valid. Need to rerun before calculting map.")
            return
        
        if(self.mainPanel.sim.running["main"]):
            self.updateHelpLabel("Error: Wait for the simulation to finish running.")
            return
        
        if(self.mainPanel.sim.running["map"]):
            super().stop()
            while(self.mainPanel.sim.running["map"]):
                print("waiting to stop")
            self.forms["Reference"]['buttons'][0][0].configure(text="RUN")
            self.updateHelpLabel("Stopped!")
            return
        
        self.forms["Reference"]['buttons'][0][0].configure(text="STOP")
        self.updateHelpLabel("Generating map...")
        self.E = np.array([params[0],params[1]])
        scale = params[2]
        func = lambda : self.mainPanel.sim.getWavefunction(self.E,scale,initiator=self)
        super().threadTask(func)
        
    def finish(self,success,psi="",extent="",X=""):
        if(not success):
            self.updateHelpLabel("Error: Cannot run while another process is running.")
            return
        
        self.forms["Reference"]['buttons'][0][0].configure(text="RUN")
        
        self.X      = X
        self.psi    = psi
        self.extent = extent
        self.updateHelpLabel("Done!")
        self.update()
    
    def progress(self,progress):
        self.updateHelpLabel("Generating map... " + str(progress) + " %")
        
    def overlay(self,option):
        if(option == "Caption"):        self.toggleCaption()
        if(option == "Scale Bar"):      self.toggleScaleBar()
        self.btn['Overlay'].set("Overlay")
        
    def toggleCaption(self):
        self.plotCaption = not self.plotCaption
        self.update()
    
    def toggleScaleBar(self):
        self.scaleBar = not self.scaleBar
        self.update()
        
    def cycleBand(self,b):
        if(self.mainPanel.sim and self.mainPanel.sim.valid and self.band + b >= 0):
            self.band += b
            
            emin = np.nanmin(self.mainPanel.sim.Ek[self.band])
            emax = np.nanmax(self.mainPanel.sim.Ek[self.band])
            for idx,e in enumerate(self.forms["Reference"]['entries']):
                e[0].grid(row=e[1],column=e[2],columnspan=1)
                e[0].configure(width=int(self.width/self.length),height=27)
                e[0].delete(0,ctk.END)
                e[0].insert(0,[emin,emax][idx])
                if(idx==1): break
        
        self.updateHelpLabel("Band:" + str(self.band))
    
    def load(self):
        pass