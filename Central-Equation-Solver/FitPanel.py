# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 15:47:55 2022

@author: jced0001
"""

from Panel import Panel
import tkinter as tk
import customtkinter as ctk
import numpy as np
import nanonispy as nap
from scipy.signal import savgol_filter as savgol

from lmfit import Model, Parameters, fit_report
from lmfit.models import GaussianModel, ConstantModel

class FitPanel(Panel):
    curveIdx = -1
    componentIdx = -1
    forms = {}
    formActive = False
    curve = []
    ###########################################################################
    # Constructor
    ###########################################################################
    def __init__(self, master, width, height, dpi, mainPanel):
        super().__init__(master, width, height, dpi, mainPanel=mainPanel,length=8,btnSize=2)
        self.buttons()
        self.reset()
        
    ###########################################################################
    # Panel
    ###########################################################################
    def buttons(self):
        self.btn = {
            "Load":     ctk.CTkButton(self.master, text="Load .dat",  command=self.loadDat),
            "Add":      ctk.CTkComboBox(self.master,values=["Add Fit"], command=self.addFitCurve),
            "Edit":     ctk.CTkComboBox(self.master,values=["Edit Fit"],command=self.editForm),
            "Reset":    ctk.CTkButton(self.master, text="Reset",      command=self.reset),
            "Close":    ctk.CTkButton(self.master, text="Close",      command=self.destroy)
            }
        
        addValues=["Add Fit","Gaussian","Fermi-Dirac","Point-Spectrum","LDOS"]
        self.btn['Add'].configure(values=addValues,variable="Add Fit")
    
    def buttonHelp(self):
        helpStr = "Load in a curve to fit (.dat)"
        self.btn['Load'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Add a new component to the fit"
        self.btn['Add'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Edit one of the existing fit components"
        self.btn['Edit'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Start from scratch"
        self.btn['Reset'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
        helpStr = "Close this panel"
        self.btn['Close'].bind('<Enter>',lambda event, s=helpStr: self.updateHelpLabel(s))
        
    def special(self):
        # Fermi-Dirac Form
        params = []
        params.append(['A','Amin','Amax'])
        params.append(['x0','x0min','x0max'])
        params.append(['T'])
        params.append(['Tb','Tbmin','Tbmax'])
        self.buildForm(name="Fermi-Dirac", params=params)
        
        # Gaussian Form
        params = []
        params.append(['A','Amin','Amax'])
        params.append(['x0','x0min','x0max'])
        params.append(['sigma','sigmaMin','sigmaMax'])
        params.append(['c','cmin','cmax'])
        self.buildForm(name="Gaussian", params=params)
        
        # Point Spectrum Form
        params = []
        params.append(['A','Amin','Amax'])
        params.append(['c','cmin','cmax'])
        self.buildForm(name="Point-Spectrum", params=params)
        
        # LDOS Form
        # params = self.mainPanel.potentialTypes[self.mainPanel.potentialType][0](args=[])
        params = []
        params.append(['A','Amin','Amax'])
        params.append(['c','cmin','cmax'])
        params.append(['v0','v0min','v0max'])
        # params.append(['vCu','vCumin','vCumax'])
        params.append(['Tilt','Tmin','Tmax'])
        params.append(['expA','expC','expP'])
        self.buildForm(name="LDOS",params=params)
        
        l = [ctk.CTkLabel(self.master, text="GridPos"),6,self.pos + 2*2]
        l[0].grid(row=l[1],column=l[2],columnspan=1)
        l[0].configure(width=int(self.width/self.length),height=27)
        
        e = [ctk.CTkEntry(self.master),6,self.pos+2*3]
        e[0].grid(row=e[1],column=e[2],columnspan=1)
        e[0].configure(width=int(self.width/self.length),height=27)
        self.gridIdxE = e
        
    def removeSpecial(self):
        self.hideForm()
    
    def buildForm(self,name,params):
        self.forms[name] = {"labels"  : [],
                            "entries" : [],
                            "buttons" : []}
        idr = 0; row = 7
        for idr,p in enumerate(params):
            for idp,param in enumerate(p):
                self.forms[name]['labels'].append([ctk.CTkLabel(self.master, text=param),row+idr,self.pos + 2*idp])
                self.forms[name]['entries'].append([ctk.CTkEntry(self.master),row+idr,self.pos+2*idp+1])
        
        idr += 1
        self.forms[name]['buttons'] = []
        self.forms[name]['buttons'].append([ctk.CTkButton(self.master, text="submit", command=lambda n=name: self.submitForm(n)),row+idr,self.pos])
        self.forms[name]['buttons'].append([ctk.CTkButton(self.master, text="cancel", command=lambda n=name: self.cancelForm(n)),row+idr,self.pos])
        self.forms[name]['buttons'].append([ctk.CTkButton(self.master, text="remove", command=lambda n=name: self.removeForm(n)),row+idr,self.pos])
        
    ###########################################################################
    # Update and Plotting
    ###########################################################################
    def update(self):
        if(not self.mainPanel.init): return
        if(not self.active): return
        
        self.ax.cla()                                                           # Clear the axis
        self.plotSTS()                                                          # Plots the sts curve selected from the sts panel
        self.plotFit()                                                          # 
        self.ax.set_position([0.13, 0.1, 0.83, 0.83])                           # Leave room for axis labels and title
        
        self.canvas.figure = self.fig                                           # Assign the figure to the canvas
        self.canvas.draw()                                                      # Redraw the canvas with the updated figure
    
    def plotSTS(self):
        # self.validateCurveIdx()
        # if(self.curveIdx < 0):
        #     self.updateHelpLabel("Add at least one spectrum to the STSPanel to start fitting")
        #     return
        
        # Ex = self.mainPanel.ldosPanel.Ex.copy()
        # spectrum = self.mainPanel.ldosPanel.smoothedLDOS[self.curveIdx].copy()
        # self.curve = Ex,spectrum
        
        if(not len(self.curve)): return
        
        Ex,spectrum = self.curve
        self.ax.plot(Ex,spectrum,linewidth=1.5)
       
        self.ax.set_xlabel("Energy (eV)")
        self.ax.set_ylabel("LDOS (arb)")
        # self.ax.set_title("Curve Fitting")
    
    def plotFit(self):
        result = self.fit()
        if(not result): return
        print(fit_report(result))
        c = self.mainPanel.mplibColours[self.curveIdx + 1]
        fit = result.init_fit.reshape(self.curve[0].shape)
        self.ax.plot(self.curve[0], fit, '--', label='init fit', c=c, linewidth=0.75, alpha=0.25)
        fit = result.best_fit.reshape(self.curve[0].shape)
        self.ax.plot(self.curve[0], fit, '--', label='best fit', c=c, linewidth=1,    alpha=0.8)
        
    ###########################################################################
    # Curve Fitting...
    ###########################################################################
    def fit(self):
        model = 0
        pars = Parameters()
        if(not self.curve): return
        xx = self.curve[0]
        yy = self.curve[1]
        if("Fermi-Dirac" in self.fitDict):
            fermiEdge = []
            fermiDirac = self.fitDict["Fermi-Dirac"]
            for idx,fd in enumerate(fermiDirac):
                A  = fd[0]; Amin  = fd[1]; Amax  = fd[2]
                x0 = fd[3]; x0min = fd[4]; x0max = fd[5]
                Tb = fd[6]; Tbmin = fd[7]; Tbmax = fd[8]
                kT = 8.617e-5*fd[9]
                
                pars.add('FD' + str(idx) + '_A',  value=A,  min=Amin,  max=Amax)
                pars.add('FD' + str(idx) + '_x0', value=x0, min=x0min, max=x0max)
                pars.add('FD' + str(idx) + '_T',  value=-kT, vary=False)
                pars.add('FD' + str(idx) + '_Tb', value=Tb, min=Tbmin, max=Tbmax)
                pars.add('FD' + str(idx) + '_c',  value=A,  min=0,     max=2*A)
                
                fermiEdge.append(Model(self.fermiDiracFunc,prefix='FD' + str(idx) + '_'))
                if(not model): model = fermiEdge[-1]
                else:          model += fermiEdge[-1]
        
        if("Gaussian" in self.fitDict):
            gaussian = self.fitDict["Gaussian"]
            for idx,g in enumerate(gaussian):
                A  = g[0]; Amin  = g[1]; Amax  = g[2]
                x0 = g[3]; x0min = g[4]; x0max = g[5]
                sigma = g[6]; sigmaMin = g[7]; sigmaMax = g[8]
                c  = g[9]; cmin  = g[10]; cmax = g[11]
                
                pars.add('GAUSS' + str(idx) + '_center',    value=x0,    min=x0min,    max=x0max)
                pars.add('GAUSS' + str(idx) + '_amplitude', value=A,     min=Amin,     max=Amax)
                pars.add('GAUSS' + str(idx) + '_sigma',     value=sigma, min=sigmaMin, max=sigmaMax)
                pars.add('GAUSS' + str(idx) + '_c',         value=c,     min=cmin,     max=cmax)
                
                if(not model): model  = GaussianModel(prefix='GAUSS' + str(idx) + '_')
                else:          model += GaussianModel(prefix='GAUSS' + str(idx) + '_')
                
                model += ConstantModel(prefix='GAUSS' + str(idx) + '_')
        
        if("Point-Spectrum" in self.fitDict):
            pointSpec = self.fitDict["Point-Spectrum"]
            for idx,ps in enumerate(pointSpec):
                A = ps[0]; Amin = ps[1]; Amax = ps[2]
                c = ps[3]; cmin = ps[4]; cmax = ps[5]
                
                pars.add('PS' + str(idx) + '_A', value=A, min=Amin, max=Amax)
                pars.add('PS' + str(idx) + '_c', value=c, min=cmin, max=cmax)
                pars.add('PS' + str(idx) + '_datFile', value=idx, vary=False)
                
                if(not model): model  = Model(self.pointSpecFunc,prefix='PS' + str(idx) + '_')
                else:          model += Model(self.pointSpecFunc,prefix='PS' + str(idx) + '_')
        
        if("LDOS" in self.fitDict):
            LDOS = self.fitDict["LDOS"]
            for idx,ldos in enumerate(LDOS):
                A = ldos[0]; Amin = ldos[1]; Amax = ldos[2]
                c = ldos[3]; cmin = ldos[4]; cmax = ldos[5]
                
                v0 = ldos[6]; v0min = ldos[7];   v0max = ldos[8]
                T = ldos[9];  Tmin = ldos[10];   Tmax = ldos[11]
                
                expA = ldos[12] 
                expC = ldos[13]
                expP = ldos[14]
                # vCu = ldos[9]; vCumin = ldos[10]; vCumax = ldos[11]
                
                pars.add('PS' + str(idx) + '_A', value=A, min=Amin, max=Amax)
                pars.add('PS' + str(idx) + '_c', value=c, min=cmin, max=cmax)
                
                pars.add('PS' + str(idx) + '_v0', value=v0, min=v0min, max=v0max)
                # pars.add('PS' + str(idx) + '_vCu', value=vCu, min=vCumin, max=vCumax)
                pars.add('PS' + str(idx) + '_T', value=T, min=Tmin, max=Tmax)
                
                pars.add('PS' + str(idx) + '_expA', value=expA, vary=False)
                pars.add('PS' + str(idx) + '_expC', value=expC, vary=False)
                pars.add('PS' + str(idx) + '_expP', value=expP, vary=False)
                
                if(not model): model  = Model(self.LDOSFunc,prefix='PS' + str(idx) + '_')
                else:          model += Model(self.LDOSFunc,prefix='PS' + str(idx) + '_')
                
        if(model == 0): return 0
        
        model.eval(pars,x=xx)
        return model.fit(yy,pars,x=xx)
    
    ###########################################################################
    # Form Actions (Show, Submit, Cancel, Remove, etc.)
    ###########################################################################
    def showForm(self,name):
        if(self.formActive): self.hideForm()
        self.formActive = True
        
        for l in self.forms[name]['labels']:
            l[0].grid(row=l[1],column=l[2],columnspan=1)
            l[0].configure(width=int(self.width/self.length),height=27)
            
        for idx,e in enumerate(self.forms[name]['entries']):
            e[0].grid(row=e[1],column=e[2],columnspan=1)
            e[0].configure(width=int(self.width/self.length),height=27)
            if(self.componentIdx < 0): continue
            e[0].delete(0,tk.END)
            e[0].insert(0,self.fitDict[name][self.componentIdx][idx])
        
        for idx,b in enumerate(self.forms[name]['buttons']):
            # if(b[0].text == "remove" and self.componentIdx == -1): continue     # Only show the 'remove' button if we're editing a selection
            b[0].grid(row=b[1],column=b[2] + 2*idx,columnspan=2)
            b[0].configure(width=int(self.width/self.length),height=27)
            print("Showing button:",b)
            
    def hideForm(self,name=""):
        self.formActive = False
        names = [name]
        if(not names[0]): names = self.forms.keys()
        
        for name in names:
            for l in self.forms[name]['labels']:
                l[0].grid_forget()
                
            for e in self.forms[name]['entries']:
                e[0].grid_forget()
            
            for b in self.forms[name]['buttons']:
                b[0].grid_forget()
        
        self.btn['Add'].set("Add Fit")
        self.btn['Edit'].set("Edit Fit")
            
    def submitForm(self,name):
        filename = ""
        if(name == "Point-Spectrum"):
            filename = self.browseFile()
            if(not filename.endswith(".dat")):
                print("Excpecting .dat file")
                return
        
        params = []
        for e in self.forms[name]['entries']:
            try:
                params.append(np.float32(e[0].get()))
            except:
                self.updateHelpLabel("Error in form: All values must be numeric.")
                return
        
        if(filename): params.append(filename)
        
        if(self.componentIdx == -1):
            self.fitDict[name].append(params)
        else:
            self.fitDict[name][self.componentIdx] = params
        
        self.hideForm(name)
        self.update()
        self.updateEditButton()
        self.componentIdx = -1
        
    def cancelForm(self,name):
        self.hideForm(name=name)
        self.componentIdx = -1
        
    def removeForm(self,name):
        del self.fitDict[name][self.componentIdx]
        self.componentIdx = -1
        self.updateEditButton()
        self.hideForm(name=name)
        self.update()
        
    def editForm(self,name):
        if(name == "Edit Fit"): return
        name,index = name.split(" ")
        self.componentIdx = int(index)
        self.showForm(name=name)
        
    ###########################################################################
    # Custom Fitting Curves (not in lmfit)
    ###########################################################################
    def fermiDiracFunc(self,x, A, x0, T, Tb, c):
        """
        Fermi-Dirac function

        Parameters
        ----------
        x  : x-axis
        A  : Amplitude
        x0 : Onset
        T  : Temperature (K)

        Returns
        -------
        Fermi-Dirac Function

        """
        return A/(1+np.exp((x0-x)/(-8.617e-5*T*Tb))) + c
    
    def pointSpecFunc(self,x,A,c,datFile):
        datFile = self.fitDict['Point-Spectrum'][datFile][6]
        ps = self.getDIDV(datFile=datFile)
        PS = self.mainPanel.stsPanel.getReferenceForCurve(x,reference=ps)       # Need to do something like this to interpolate curves 
        return A*PS + c
    
    def LDOSFunc(self,x,A,c,v0,T,expA,expC,expP):
        func,params = self.mainPanel.potentialTypes[self.mainPanel.potentialType]
        params[4] = v0
        # params[6] = vCu
        params[6] = T
        
        a,X,V = func(params)
        self.mainPanel.sim.run(a,X,V)
        
        dE   = 0.06
        emin = np.min(self.curve[0])
        emax = np.max(self.curve[0])
        pts  = len(self.curve[0])
        
        x0   = self.mainPanel.ldosPanel.x0s[0]
        
        LDOS,Ex = self.mainPanel.sim.getLDOS(emin,emax,dE,int(pts),[x0])
        
        return A*LDOS + c + expA*np.exp(Ex*expP) + expC
        
    ###########################################################################
    # Misc Button Functions
    ###########################################################################
    def loadDat(self):
        filename = self.browseFile()
        if(filename.endswith(".dat")):
            self.curve = self.getDIDV(filename,normalise=True)
            self.update()
            return
        
        if(filename.endswith(".3ds")):
            header = {'Delay before measuring (s)':'0',                         # Not used, just stop nap from complaining
                      'Start time':'0',
                      'End time':'0',
                      'Comment':'none'}
            self.gridData = nap.read.Grid(filename, header_override=header)     # Read the file. see gridData.header.keys() and gridData.signals.keys() for structure
            sweep = self.gridData.signals['sweep_signal']
            gridIdx = int(self.gridIdxE[0].get())
            
            # x_idx = int((pos[0]/lxy[0]) * self.gridData.header['dim_px'][0])
            # y_idx = int((pos[1]/lxy[1]) * self.gridData.header['dim_px'][1])
            # indexes.append(np.array([x_idx,y_idx]))
            
            # filenum = self.gridData.header['dim_px'][0]*y_idx + x_idx
            
            y_idx = int(np.floor(gridIdx/self.gridData.header['dim_px'][0]))
            x_idx = gridIdx - self.gridData.header['dim_px'][0]*y_idx
            c = self.gridData.signals["LI Demod 1 X (A)"][y_idx][x_idx]
            print(gridIdx)
            print("len",c.shape)
            self.curve = self.getDIDV(curve=[sweep,c],ychannel="LI Demod 1 X (A)")
            self.update()
            return
        
        print("Excpecting .dat or .3ds file")
        self.updateHelpLabel("Error: Excpecting .dat or 3ds file")
    
    def getDIDV(self,datFile="",curve=[],xchannel="Bias calc (V)",ychannel="LI Demod 1 X (A)",normalise=False):
        V = 0; didv = 0
        if(datFile):
            dat = nap.read.Spec(datFile)
            V = dat.signals[xchannel]
            I = dat.signals[ychannel]
        elif(len(curve)):
            V = curve[0]
            I = curve[1]
        else:
            return V,didv
        
        dV = V[1] - V[0]
        
        didv = 0*I
        self.sg_pts  = 5
        self.sg_poly = 1
        if('Current' in ychannel):
            didv = savgol(I,self.sg_pts,self.sg_poly,deriv=1,delta=dV)
        
        if('Demod' in ychannel):
            didv = savgol(I,self.sg_pts,self.sg_poly,deriv=0)
        
        if(normalise): didv /= np.max(didv)
        
        return V,didv
    
    def getReferenceForCurve(self,x,reference):
        """
        This function is useful when the reference spectra is not exactly the 
        same range/number of points as the data. To return a valid reference, 
        the domain of the data must be within the domain of the refernce.
        Simple linear interpolation is used when the number of data points is
        greater than the number of points in the reference spectrum in the 
        overlapping region
        """
        try:
            return np.interp(x, reference[0], reference[1])
        except Exception as e:
            print(e)
            return 0
    
    def validateCurveIdx(self):
        numCurves = len(self.mainPanel.ldosPanel.smoothedLDOS)
        if(not numCurves): self.curveIdx = -1
        if(self.curveIdx >= numCurves): self.curveIdx = numCurves - 1
        if(self.curveIdx < 0 and numCurves): self.curveIdx = 0
        
    def addFitCurve(self,name):
        if(name == "Add Fit"): return
        if(not name in self.fitDict):
            self.fitDict[name] = []
        self.showForm(name=name)
        
    def updateEditButton(self):
        editValues = ["Edit Fit"]
        for key in self.fitDict.keys():
            for idx,p in enumerate(self.fitDict[key]):
                editValues.append(key + " " + str(idx))
        self.btn['Edit'].configure(values=editValues,variable="Edit Fit")
        
    def reset(self):
        self.fitDict = {}
        self.curve = []
        self.componentIdx = -1
        self.updateEditButton()
        self.update()