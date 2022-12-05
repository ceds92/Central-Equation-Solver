# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 07:59:08 2022

@author: jced0001
"""

from Panel import Panel
import customtkinter as ctk
import numpy as np
from scipy.signal import savgol_filter as savgol

class LDOSPanel(Panel):
    ###########################################################################
    # Constructor
    ###########################################################################
    def __init__(self, master, width, height, dpi, mainPanel):
        super().__init__(master, width, height, dpi, mainPanel=mainPanel,length=8,btnSize=2)
        super().initGlobs(name="LDOS")
        self.buttons()
        self.forms = {}
        self.x0s   = []
        self.sg_pts  = 1
        self.sg_poly = 1
        self.tunnellingFactor = 0
        self.gridLines = True
        self.LDOS  = []
        self.smootheLDOS = []
        self.init = 1
    
    ###########################################################################
    # Panel
    ###########################################################################
    def buttons(self):
        self.btn = {
            "Addx0": ctk.CTkButton(self.master, text="Add LDOS",command=self.addx0),
            "Undo":  ctk.CTkButton(self.master, text="Undo",    command=self.undo),
            "Grid":  ctk.CTkButton(self.master, text="Grid",    command=self.toggleGrid),
            "PNG":   ctk.CTkButton(self.master, text="Exp PNG", command=super().exportPNG), # Export the canvas to png
            "Close": ctk.CTkButton(self.master, text="Close",   command=self.destroy)
            }
        
        helpStr = "Place an LDOS marker on the potential in the left figure"
        self.btn['Addx0'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
    
        helpStr = "Undo the last LDOS marker"
        self.btn['Undo'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Close this panel"
        self.btn['Close'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Toggle grid lines on the plot"
        self.btn['Grid'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Export the panel figure as a png"
        self.btn['PNG'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
    def special(self):
        params = []
        params.append(['Emin','Emax','dE'])
        params.append(['Pts'])
        self.buildForm(name="LDOSForm", params=params)
        self.showForm(name="LDOSForm")
        
        self.smoothSlider = ctk.CTkSlider(self.master, orientation=ctk.HORIZONTAL, from_=0, to=15, width=420, command=self.smoothing) # Slider to select which bias/sweep signal slice to look show
        self.smoothSlider.grid(row=10,column=self.pos,columnspan=8,rowspan=1)   # Make it take up the entire length of the panel
        self.smoothSlider.set(0)
        
        self.expSlider = ctk.CTkSlider(self.master, orientation=ctk.HORIZONTAL, from_=0, to=15, width=420, command=self.exponential) # Slider to select which bias/sweep signal slice to look show
        self.expSlider.grid(row=12,column=self.pos,columnspan=8,rowspan=1)      # Make it take up the entire length of the panel
        self.expSlider.set(0)

    def removeSpecial(self):
        self.hideForm()
        self.smoothSlider.grid_forget()
        self.expSlider.grid_forget()
    
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
        self.ax.set_position([0.09, 0.07, 0.87, 0.9])                           # Make it take up the whole canvas
        self.plotLDOS()
        if(self.gridLines):
            self.ax.grid(b=True, which='major', linestyle='-')
            self.ax.grid(b=True, which='minor', linestyle='--')
            self.ax.minorticks_on()
        self.canvas.figure = self.fig                                           # Assign the figure to the canvas
        self.canvas.draw()    
    
    def plotLDOS(self):
        if(not len(self.LDOS)): return
        for LDOS in self.smoothedLDOS:
            exponential = 0.1*np.exp(self.Ex*self.tunnellingFactor/5) - 0.1
            self.ax.plot(self.Ex,LDOS + exponential)
        self.ax.plot(self.Ex,exponential,linestyle = 'dashed',linewidth=1,alpha=0.5)
            
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
            e[0].insert(0,[0,1,0.03,200][idx])
        
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
                
        if(not (self.mainPanel.sim and self.mainPanel.sim.valid)):
            self.updateHelpLabel("Error: simulation not valid. Need to rerun before calculting LDOS.")
            return
        
        if(self.mainPanel.sim.running["main"]):
            self.updateHelpLabel("Error: Wait for the simulation to finish running.")
            return
        
        if(self.mainPanel.sim.running["LDOS"]):
            super().stop()
            while(self.mainPanel.sim.running["LDOS"]):
                print("waiting to stop")
            self.forms["LDOSForm"]['buttons'][0][0].configure(text="RUN")
            self.updateHelpLabel("Stopped!")
            return
        
        self.forms["LDOSForm"]['buttons'][0][0].configure(text="STOP")
        self.updateHelpLabel("Calculating LDOS...")
        
        emin,emax,dE,pts = params
        x0s = np.array(self.x0s)
        
        if(not len(x0s)):
            self.updateHelpLabel("Add at least one marker on the main figure before running")
            return
        
        func = lambda : self.mainPanel.sim.getLDOS(emin,emax,dE,int(pts),x0s,initiator=self)
        super().threadTask(func)
        
    def finish(self,success,LDOS="",Ex=""):
        if(not success):
            self.updateHelpLabel("Error: Cannot run while another process is running.")
            return
        
        self.forms["LDOSForm"]['buttons'][0][0].configure(text="RUN")
        
        self.Ex   = Ex
        self.LDOS = LDOS
        self.smootheLDOS = self.smoothing(event=-1)
        
        self.update()
        self.mainPanel.update(upd=[0])
        
        self.updateHelpLabel("Done!")
    
    def progress(self,progressStr):
        self.updateHelpLabel("Calculating LDOS... " + str(progressStr))
    ###########################################################################
    # Adding x0
    ###########################################################################
    def addx0(self):
        self.updateHelpLabel("Place the LDOS marker on the potential map to the left\n"
                             + "Left click to place\nRight click to cancel")
        self.mainPanel.addx0Bind()
        
    def setx0(self,x0):
        self.x0s.append(np.array(x0))
    
    def undo(self):
        if(not len(self.x0s)): return
        del self.x0s[-1]
        self.mainPanel.update(upd=[0])
    
    ###########################################################################
    # Smoothing
    ###########################################################################
    def smoothing(self,event):
        if(event >=0):
            self.sg_pts = 2*int(event) + 1                                      # Change the bias on a slider event
        self.smoothedLDOS = self.LDOS.copy()
        if(self.sg_pts > self.sg_poly):
            for l,LDOS in enumerate(self.LDOS):
                self.smoothedLDOS[l] = savgol(LDOS,self.sg_pts,self.sg_poly,deriv=0)
        self.update()
    ###########################################################################
    # Misc
    ###########################################################################
    def toggleGrid(self):
        self.gridLines = not self.gridLines
        self.update()
    
    def exponential(self,event):
        self.tunnellingFactor = event
        self.update()
        
    def load(self):
        if(len(self.mainPanel.sim.LDOS)):
            self.LDOS = self.mainPanel.sim.LDOS
            self.Ex   = self.mainPanel.sim.LDOSEx
            self.x0s  = list(self.mainPanel.sim.LDOSx0s)
            
            self.smootheLDOS = self.smoothing(event=-1)
            self.update()