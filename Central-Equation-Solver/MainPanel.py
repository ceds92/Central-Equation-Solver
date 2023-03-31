# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 09:27:19 2022

@author: jced0001
"""

from Panel import Panel
from MapsPanel import MapsPanel
from BandStructure3D import BandStructure3D
from BandStructure2D import BandStructure2D
from PotentialPanel import PotentialPanel
from MapViewerPanel import MapViewerPanel
from LDOSPanel import LDOSPanel
from SweepPanel import SweepPanel
from FitPanel import FitPanel
import customtkinter as ctk
from   tkinter import filedialog
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import potentials as pot
from EPWE import EPWE
import pickle
import global_

matplotlib.use("TkAgg")


class MainPanel(Panel):
    mplibColours = plt.rcParams['axes.prop_cycle'].by_key()['color'] + ['white'] + ['black']
    bound = False                                                               # Used so only one function can be bound at a time
    
    ###########################################################################
    # Constructor
    ###########################################################################
    def __init__(self, master, width=512, height=512, dpi=96, scaleFactor=1):
        self.init()
        super().__init__(master, width=width, height=height, dpi=dpi, length=8, btnSize=2, scaleFactor=scaleFactor)
        super().initGlobs(name="main")
        self.buildSubPanels()
        self.buttons()
        super().create()
        self.canvas.figure = self.fig
        self.canvas.draw()
        
    def init(self):
        self.init = True
        self.scaleBar = True
        self.plotCaption = True
        self.sim = EPWE()
        self.currentstsPos = []
        self.forms = {}
        self.potentialType  = "Hexagonal"
        self.potentialTypes = {"Hexagonal" : [pot.hexagonal_primitive,  [1,0.3,101,-2]],
                               "Kagome"    : [pot.kagome_primitive,     [1,0.4,101,-2]],
                               "Muffin Tin": [pot.muffin_primitive,     [1,0.25,101,-2.5]],
                               "Big Star"  : [pot.HATCu_bigstar,        [0.93,0.9,201,1,0.55,0.2,0.4,-0.07]],
                               "Small Star": [pot.HATCu_smallstar,      [0.68,0.9,201,1,0.55,0.2,0.4,-0.07]],
                               "SimParams" : [self.run,                 [31,5,0.4]]}
        
    ###########################################################################
    # Panel
    ###########################################################################
    def buildSubPanels(self):                                                   # Sub panels are build in here. add new ones here first
        commonParams = [self.master, self.width, self.height, self.dpi, self]
        
        self.mapsPanel = MapsPanel(*commonParams)
        self.bs3DPanel = BandStructure3D(*commonParams)
        self.bs2DPanel = BandStructure2D(*commonParams)
        self.potentialPanel = PotentialPanel(*commonParams)
        self.ldosPanel = LDOSPanel(*commonParams)
        self.mapViewerPanel = MapViewerPanel(*commonParams)
        self.sweepPanel = SweepPanel(*commonParams)
        self.fitPanel = FitPanel(*commonParams)
        
        self.panels = [self.mapsPanel,self.bs3DPanel,self.bs2DPanel,self.potentialPanel,self.ldosPanel,self.mapViewerPanel,self.sweepPanel,self.fitPanel]
        
    def buttons(self):
        self.btn = {
            "Potential":ctk.CTkComboBox(self.master,values=["Potential"],   command=self.potselect),      # Dropdown to select potential map type
            "Overlay":  ctk.CTkComboBox(self.master,values=["Overlay"],     command=self.overlay),        # Dropdown to change overlay display
            "cmap":     ctk.CTkButton(self.master, text="viridis",          command=super()._cmap),       # Button to cycle through colour maps
            "Export":   ctk.CTkComboBox(self.master,values=["Export"],      command=self.export),         # Dropdown to export the figure
            "OpenPanel":ctk.CTkComboBox(self.master,values=["Open Panel"],  command=self.openPanel),      # Dropdown to open other panels
            "Save":     ctk.CTkButton(self.master, text="Save",             command=self.save),           # Save all session to a .epwe file
            "Load":     ctk.CTkButton(self.master, text="Load",             command=self.load),           # Load a .epwe file
            "Quit":     ctk.CTkButton(self.master, text="Quit",             command=self.quit)            # Button to quit the program
            }
        
        potentialValues = ["Potential Type"] + list(self.potentialTypes.keys())
        potentialValues.remove("SimParams")
        self.btn['Potential'].configure(values=potentialValues,fg_color=['#3B8ED0', '#1F6AA5'])
        self.btn['Potential'].set(self.potentialType)
        
        openPanelValues = ["Open Panel","Rebuilt Potential","2D Bands","3D Bands","Map Generator","Map Viewer","LDOS","Sweep","Fitting"]
        self.btn['OpenPanel'].configure(values=openPanelValues,fg_color=['#3B8ED0', '#1F6AA5'])
        
        overlayValues = ["Overlay","Caption","Scale Bar"]
        self.btn['Overlay'].configure(values=overlayValues,fg_color=['#3B8ED0', '#1F6AA5'])

        exportValues = ["Export","PNG","Pickle"]
        self.btn['Export'].configure(values=exportValues,fg_color=['#3B8ED0', '#1F6AA5'])
    
    def buttonHelp(self):
        helpStr = "Select a potntial type"
        self.btn['Potential'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Change the colour map"
        self.btn['cmap'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Export figure"
        self.btn['Export'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Show/hide overlay features"
        self.btn['Overlay'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Open another processing window"
        self.btn['OpenPanel'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Save this session to a .epwe file"
        self.btn['Save'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Load a previous session from a .epwe file"
        self.btn['Load'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Quit the program"
        self.btn['Quit'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
    def reorderPanels(self,destroyPos,destroyLength):
        for panel in self.panels:
            if(panel.active):
                if(panel.pos > destroyPos):
                    panel.destroy()
                    panel.pos -= destroyLength
                    panel.create()
    
    def getPos(self):
        pos = self.length
        for panel in self.panels:
            if(panel.active):
                pos += panel.length
        return pos
    
    def adjustWindowSize(self):
        windowWidth = self.width
        for panel in self.panels:
            if(panel.active):
                windowWidth += panel.width
        self.master.geometry("%dx%d" % (windowWidth,850))
    
    def special(self):
        for name,func in self.potentialTypes.items():
            params = func[0]()
            row = 9
            btnName = "Update Potential"
            if(name == "SimParams"):
                row = 7
                btnName = "Run Simulation"
            self.buildForm(name=name, params=params,row=row,btnName=btnName)
        
        self.showForm(name="SimParams")
        self.showForm(name=self.potentialType)
        self.forms["SimParams"]['buttons'][0][0].configure(fg_color='Red')
        
    def removeSpecial(self):
        self.hideForm()
    
    def buildForm(self,name,params,row,btnName):
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
        self.forms[name]['buttons'].append([ctk.CTkButton(self.master, text=btnName, command=lambda n=name: self.submitForm(n)),row+idr,self.pos])
        
    ###########################################################################
    # Main Panel Updates
    ###########################################################################
    def update(self,upd=[-1]):
        if(-1 in upd or 0 in upd):                                              # upd=0 selects mainPanel
            self.ax.cla()                                                       # Clear the axis
            self.updatePotential()                                              # Show the processed (tilt-corrected, filtered, etc.) sxm image
            self.ax.set_position([0, 0, 1, 1])                                  # Make it take up the whole canvas
            self.ax.axis('off')                                                 # Hide the axis
            
        if((-1 in upd or 1 in upd) and self.bs2DPanel.active):                  # upd=1 selects lineProfile panel
            self.bs2DPanel.update()
            
        if((-1 in upd or 2 in upd) and self.bs3DPanel.active):                  # upd=1 selects lineProfile panel
            self.bs3DPanel.update()
            
        if(-1 in upd or 0 in upd):                                              # upd=0 selects mainPanel. Update the overlay at the end, since other panels may contribute to the overlay
            self.updateOverlay()                                                # Add things to the foreground (e.g. scalebar and plot label)
        
        self.canvas.figure = self.fig                                           # Assign the figure to the canvas
        self.canvas.draw()                                                      # Redraw the canvas with the updated figure
    
    def updateOverlay(self):
        left, right = self.ax.get_xlim();                                       # Remember plot extent
        bottom, top = self.ax.get_ylim();                                       # Remember plot extent
        
        if(self.ldosPanel.active):
            for c,x0 in enumerate(self.ldosPanel.x0s):
                self.ax.plot(x0[0],x0[1],'.',markersize=14,c=self.mplibColours[c])
            if(len(self.currentstsPos)):
                self.ax.plot(self.currentstsPos[0],self.currentstsPos[1],'x',markersize=14)
                
        if(self.sweepPanel.active):
            for c,x0 in enumerate(self.sweepPanel.x0s):
                self.ax.plot(x0[0],x0[1],'*',markersize=14,c=self.mplibColours[c])
            if(len(self.currentstsPos)):
                self.ax.plot(self.currentstsPos[0],self.currentstsPos[1],'x',markersize=14)
        
        UC = self.getUC(self.X)
        for line in UC:
            self.ax.plot(line[0],line[1],c='red',linewidth=1)
        
        if(self.scaleBar): super().addPlotScalebar()                            # Add a scale bar to the plot
        if(self.plotCaption):                                                   # Caption the image with Vbias and Iset
            plotCaption  = self.potentialType                                   # Caption image
            rx = abs(self.extent[1] - self.extent[0])
            ry = abs(self.extent[3] - self.extent[2])
            pos = [0.025*rx + self.extent[0],self.extent[3] - 0.1*ry]           # Put the caption in the top left
            super().addPlotCaption(plotCaption, pos)
                
        self.ax.set_xlim((left,right)); self.ax.set_ylim((bottom,top))          # Put back extent
    
    def updatePotential(self):
        func,params = self.potentialTypes[self.potentialType]
        a,X,V = func(params)
        
        x1,x2 = X
        self.extent = np.array([min(x1),max(x1),min(x2),max(x2)])
        
        cmap = self.cmaps[self.cmap][1]
        self.ax.imshow(V,extent=self.extent,cmap=cmap())
        self.ax.set_facecolor('black')
        
        self.a = a; self.X = X; self.V = V
        
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
            e[0].insert(0,self.potentialTypes[name][1][idx])
        
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
                
    def submitForm(self,name,skipParams=False):
        paramsChanged = False
        params = self.potentialTypes[name][1]
        if(not skipParams):
            params    = []
            oldparams = self.potentialTypes[name][1]
            for p,e in enumerate(self.forms[name]['entries']):
                try:
                    params.append(np.float32(e[0].get()))
                except:
                    self.updateHelpLabel("Error in form: All values must be numeric.")
                    return
                if(not np.float32(oldparams[p]) == params[p]):
                    paramsChanged = True
            
            self.potentialTypes[name][1] = params
        
        if(name == "SimParams"):
            if(self.sim.running["main"]):
                super().stop()
                cnt = 0
                while(self.sim.running["main"]):
                    print("waiting to stop")
                    cnt += 1
                    if(cnt > 5):
                        self.sim.running["main"] = False
                        self.sim.valid = False
                        print("An error occured - probably due to lack of memory, try reducing ks or N and run again...")
                        global_.main_running.clear()
                        break
                self.updateHelpLabel("Simulation stopped.")
                self.forms["SimParams"]['buttons'][0][0].configure(text="Run Simulation")
                self.forms["SimParams"]['buttons'][0][0].configure(fg_color='red')
            else:
                for simType,isRunning in self.sim.running.items():
                    if(isRunning):
                        waitStr = "Wait for " + simType + " to finish or stop it before running"
                        self.updateHelpLabel(waitStr)
                        return
                self.updateHelpLabel("Running!")
                self.forms["SimParams"]['buttons'][0][0].configure(fg_color=['#3B8ED0', '#1F6AA5'])
                self.forms["SimParams"]['buttons'][0][0].configure(text="STOP")
                self.run(params)
        elif(paramsChanged):
            self.forms["SimParams"]['buttons'][0][0].configure(fg_color='red')
            self.updateHelpLabel("Warning: Current simulation does not reflect the potential shown!")
        self.update(upd=[0])
        
    def run(self,args=[]):
        if(args==[]):
            params = [['ks','N','m_eff']]
            return params
        ks,N,m_eff = args
        self.updateHelpLabel("Running Sim...")
        self.sim = EPWE(int(ks),int(N),m_eff)
        
        func = lambda : self.sim.run(self.a,self.X,self.V,initiator=self)
        super().threadTask(func)
    
    def simFinished(self,success):
        self.forms["SimParams"]['buttons'][0][0].configure(text="Run Simulation")
        
        if(not success):
            self.updateHelpLabel("Cannot run sim while another process is running")
            return
        
        self.updateHelpLabel("Finished Running!")
        if(self.sim.valid):
            self.update()
        
    def progress(self,progressStr):
        self.updateHelpLabel(progressStr)
    ###########################################################################
    # Drawing STS Markers
    ###########################################################################
    def addx0Bind(self,name):
        if(self.bound): return
        self.bound = True
        if(name == self.ldosPanel.name):
            self.lcMarkSTSBind = self.canvas.get_tk_widget().bind('<Button-1>', self.addx0) # Default to LDOS Panel
        if(name == self.sweepPanel.name):
            self.lcMarkSTSBind = self.canvas.get_tk_widget().bind('<Button-1>', self.addx0_sweep) # Placing markers from sweep panel
            
        self.rcMarkSTSBind     = self.canvas.get_tk_widget().bind('<Button-3>', self.cancelx0)
        self.motionMarkSTSBind = self.canvas.get_tk_widget().bind('<Motion>',   self.placex0)
        
        self.updateHelpLabel("Place the LDOS marker on the image\n"
                             + "Left click to place\nRight click to cancel")
        
    def addx0Unbind(self):
        self.canvas.get_tk_widget().unbind('<Button-1>', self.lcMarkSTSBind)
        self.canvas.get_tk_widget().unbind('<Button-3>', self.rcMarkSTSBind)
        self.canvas.get_tk_widget().unbind('<Motion>',   self.motionMarkSTSBind)
        self.bound = False
        
        self.updateHelpLabel("")
        
    def cancelx0(self,event):
        self.addx0Unbind()
        self.currentstsPos = []
        self.update(upd=[0])
        
    def addx0(self,event):
        size = self.fig.get_size_inches()*self.fig.dpi                          # size in pixels
        x = event.x
        y = size[1] - event.y
        X = super()._getX(x)
        Y = super()._getY(y)
        
        stsPos = [X,Y]
        self.ldosPanel.setx0(stsPos)
        self.currentstsPos = []
        
        self.addx0Unbind()
        self.update(upd=[0])
        
    def addx0_sweep(self,event):
        size = self.fig.get_size_inches()*self.fig.dpi                          # size in pixels
        x = event.x
        y = size[1] - event.y
        X = super()._getX(x)
        Y = super()._getY(y)
        
        stsPos = [X,Y]
        self.sweepPanel.setx0(stsPos)
        self.currentstsPos = []
        
        self.addx0Unbind()
        self.update(upd=[0])
    
    def placex0(self,event):
        size = self.fig.get_size_inches()*self.fig.dpi                          # size in pixels
        x = event.x
        y = size[1] - event.y
        X = super()._getX(x)
        Y = super()._getY(y)
        
        self.currentstsPos = [X,Y]
        self.update(upd=[0])
       
    ###########################################################################
    # Misc
    ###########################################################################
    def save(self):
        if(not self.sim or (not self.sim.valid)):
            self.updateHelpLabel("Error, simulation not valid. Need to rerun before saving.")
            return
        
        default = 'sim.epwe'
        path = filedialog.asksaveasfilename(title="Save as",initialfile=default)
        if(not path.endswith('.epwe')): path += '.epwe'
        
        potType = self.potentialType
        params  = self.potentialTypes[potType][1]
        
        self.sim.save(path)
        
        pkldict = pickle.load(open(path,'rb'))
        pkldict["potparams"] = {}
        pkldict["potparams"]["type"]   = potType
        pkldict["potparams"]["params"] = params
        
        pkldict["simparams"] = self.potentialTypes["SimParams"][1]
        
        pickle.dump(pkldict,open(path,'wb'))
        
        self.updateHelpLabel("Save successful!")
    
    def load(self):
        path = self.browseFile()
        if(not path.endswith('.epwe')):
            self.updateHelpLabel("Error: Expected a .epwe file.")
            return
        
        self.sim.load(path)
        
        for panel in self.panels:
            panel.load()
        
        pkldict = pickle.load(open(path,'rb'))
        self.potentialType = pkldict["potparams"]["type"]
        self.potentialTypes[self.potentialType][1] = pkldict["potparams"]["params"]
        
        self.potentialTypes["SimParams"][1] = pkldict["simparams"]
        
        self.hideForm(self.potentialType)
        self.showForm(self.potentialType)
        
        self.hideForm("SimParams")
        self.showForm("SimParams")
        
        self.btn['Potential'].set(self.potentialType)
        self.forms["SimParams"]['buttons'][0][0].configure(fg_color=['#3B8ED0', '#1F6AA5'])
        
        self.update()
        
    def getUC(self,X):
        x1,x2 = X
        a1,a2 = self.a
        
        mid = int(len(x1)/2)
        u = []
        oo = np.array(x1[mid],x2[mid])
        oo = oo -a1/2 -a2/2
        u.append(oo + a1)
        u.append(u[0] + a2)
        u.append(u[1] - a1)
        u.append(u[2] - a2)
        u.append(u[0])
        
        UC = []
        for p in range(5):
            xx = u[p][0],u[p+1][0]
            yy = u[p][1],u[p+1][1]
            UC.append([xx,yy])
            if(p==3): break
        
        return  np.array(UC)
        
    def potselect(self,option):
        if(option == "Potential Type"):
            self.btn['Potential'].set(self.potentialType)
            return
        if(option != self.potentialType):
            self.hideForm(self.potentialType)
            self.potentialType = option
            self.showForm(option)
            self.forms["SimParams"]['buttons'][0][0].configure(fg_color='Red')
            self.update()
            
            if(self.sweepPanel.active): self.sweepPanel.updateParams()
    
    def openPanel(self,option):
        if(option == "Map Generator"):      self.mapsPanel.create()
        if(option == "Map Viewer"):         self.mapViewerPanel.create()
        if(option == "2D Bands"):           self.bs2DPanel.create()
        if(option == "3D Bands"):           self.bs3DPanel.create()
        if(option == "Rebuilt Potential"):  self.potentialPanel.create()
        if(option == "LDOS"):               self.ldosPanel.create()
        if(option == "Sweep"):              self.sweepPanel.create()
        if(option == "Fitting"):            self.fitPanel.create()
        
        self.btn['OpenPanel'].set("Open Panel")
        
    def overlay(self,option):
        if(not self.init):
            self.btn['Overlay'].set("Overlay")
            return
        if(option == "Caption"):    self.toggleCaption()
        if(option == "Scale Bar"):  self.toggleScaleBar()
        if(option == "Export PNG"): super().exportPNG()
        self.btn['Overlay'].set("Overlay")
    
    def export(self,option):
        if(option == "PNG"): super().exportPNG()
        if(option == "Pickle"): self.exportPickle()
        self.btn['Export'].set("Export")
    
    def exportPickle(self):
        pklDict = {"V": self.V,
                   "a": self.a,
                   "X": self.X,
                   "extent": self.extent}
        
        super().exportPickle(pklDict=pklDict,initialfile="potential")
        
    def toggleCaption(self):
        self.plotCaption = not self.plotCaption
        self.update(upd=[0])
    
    def toggleScaleBar(self):
        self.scaleBar = not self.scaleBar
        self.update(upd=[0])
        
    def quit(self):
        self.master.destroy()                                                   # Close the Tkinter GUI
        os._exit(00)                                                            # Easier to restart kernal than to figure out Tkinter memory leaks (think this is a problem with spyder, not tkinter)