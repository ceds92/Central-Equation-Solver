# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 09:28:48 2022

@author: jced0001
"""

import customtkinter as ctk
from tkinter import filedialog
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib import colors
from matplotlib_scalebar.scalebar import ScaleBar
import global_
import threading
import pickle

class Panel():
    pos    = 0                                                                  # This panel's column position. init to zero. gets added to the RHS of the window when created
    active = False; cmap  = 0; shift = False                                    # Use alternate values - bound to the left shift key
    zunit = 1e-12; xunit = 1e-9                                                 # convert xy units to nm and z units to pm
    helpText = ""
    imprint  = False
    def __init__(self, master, width, height, dpi, mainPanel=None,length=4,btnSize=1,plotType='', scaleFactor=1):
        self.master = master
        self.width  = width
        self.height = height
        self.length = length
        self.btnSize = btnSize
        self.dpi    = dpi
        self.mainPanel = mainPanel
        
        if(not mainPanel): self.scaleFactor = scaleFactor
        else:              self.scaleFactor = self.mainPanel.scaleFactor
        
        # Set up canvas
        self.canvas = FigureCanvasTkAgg(master=master)
        self.canvas.get_tk_widget().configure(width=width*self.scaleFactor, height=height*self.scaleFactor)
        
        # Figure
        self.fig = plt.figure(figsize=(width/dpi,height/dpi),dpi=dpi)
        if(not plotType):     self.ax = self.fig.add_subplot(111)
        if(plotType == '3d'): self.ax = self.fig.add_subplot(111,projection='3d')
        
        # Misc init
        self._initCmaps()
    
    def initGlobs(self,name):
        exec("global_." + name + "_task = []")
        exec("global_." + name + "_running = threading.Event()")                # event to stop threads
        exec("self.task = global_." + name + "_task")
        exec("self.running = global_." + name + "_running")
        
    def addPlotCaption(self, plotCaption, pos, props=None, fontsize=16, ax=None):
        if(not props):
            props = dict(boxstyle='round',facecolor='white',alpha=0.5)          # Box properties. alpha = transparency
        if(not ax): ax = self.ax
        ax.text(pos[0],pos[1],plotCaption,bbox=props,fontsize=fontsize)
    
    def addPlotScalebar(self, fontsize=20, ax=None):
        # Add a scalebar to the sxm image
        font_properties = {"size": fontsize}
        scalebar = ScaleBar(dx=1, units='nm', length_fraction=.5,
        font_properties=font_properties, frameon=False, location=4, color='w')
        
        if(not ax): ax = self.ax
        ax.add_artist(scalebar)
        
    def addButtons(self):
        numCols = np.ceil(self.length/self.btnSize)
        row = 5; col = 0                                                        # Row=5 because rows 0-4 are taken up by the canvas
        self.master.rowconfigure(index=row,weight=1,minsize=40)
        for btn in self.btn.values():                                           # Now loop through each button and put it in the correct position. 4 per row
            btn.grid(row=row,column=self.btnSize*col+self.pos,columnspan=self.btnSize)
            btn.configure(width=int(self.width/5),height=27)
            col += 1
            if(col == numCols):
                col  = 0
                row += 1
                self.master.rowconfigure(index=row,weight=1,minsize=40)
                
    def removeButtons(self):                                                    # Remove all the buttons after destroying the canvas
        for btn in self.btn.values():
            btn.grid_forget()

    def special(self):
        pass
    def removeSpecial(self):
        pass
    
    def _helpLabel(self):
        self.helpLabel = ctk.CTkLabel(self.master,text="",justify=ctk.LEFT, wraplength=self.width - 5)
        for r in range(4):
            self.master.rowconfigure(index=r+14,weight=1,minsize=20)
        self.helpLabel.grid(row=14,column=self.pos,columnspan=self.length,rowspan=4,sticky='nsew')
    
    def removeHelpLabel(self):
        self.helpLabel.grid_forget()
    
    def updateHelpLabel(self,helpText):
        self.helpText = helpText
        self.helpLabel.configure(text=helpText)
    
    def buttonHelp(self):                                                       # Override this
        pass
    
    def create(self):                                                           # Displays the panel to the right of the previous one
        if(self.mainPanel):
            if(not self.mainPanel.init): return
        if(self.active): return                                                 # Do nothing if the panel is already active
        
        if(not self.mainPanel): self.pos = 0
        if(self.mainPanel): self.pos = self.mainPanel.getPos()
        
        for col in range(self.length):
            self.master.columnconfigure(index=col + self.pos,weight=1)
        
        self.canvas.get_tk_widget().grid(row=0,column=self.pos, columnspan=self.length,rowspan=4) # Put this panel after the last one (left to right)
        self.addButtons()                                                       # Display the buttons
        self.special()
        self._helpLabel()
        self.buttonHelp()
        self.active = True
        self.update()
        if(self.mainPanel):
            self.mainPanel.update()
            self.mainPanel.adjustWindowSize()
        
    def destroy(self):                                                          # Hide this canvas, it's panel is not active anymore
        self.canvas.get_tk_widget().grid_forget()
        self.removeButtons()                                                    # Also hide the buttons
        self.removeSpecial()
        self.removeHelpLabel()
        self.active = False
        self.mainPanel.update()
        self.mainPanel.reorderPanels(self.pos,self.length)
        self.mainPanel.adjustWindowSize()
        
    def _imprint(self):
        self.imprint = not self.imprint
        self.btn['Imprint'].configure(bg=['SystemButtonFace','red'][self.imprint])
        self.mainPanel.update()
    
    def _initCmaps(self):
        self.cmaps = {0 : ['viridis',lambda : "viridis"],                       # Colour map for all images. More at https://matplotlib.org/stable/gallery/color/colormap_reference.html
                      1 : ['plasma', lambda : "plasma"],
                      2 : ['inferno',lambda : "inferno"],
                      3 : ['magma',  lambda : "magma"],
                      4 : ['cividis',lambda : "cividis"],
                      5 : ['flame',  lambda : self.customCmap(cmap='flame')],
                      6 : ['Greys',  lambda : "Greys_r"],
                      7 : ['Blues',  lambda : "Blues_r"],
                      8 : ['Purples',lambda : "Purples_r"]
                      }
        
    def _cmap(self):
        self.cmap += 1
        if(self.cmap == len(self.cmaps)):
            self.cmap = 0
        self.btn['cmap'].configure(text=self.cmaps[self.cmap][0])
        self.update()
        
    def customCmap(self,cmap):
        nodes=[0, 0.2, 0.5, 0.8, 1]
        col = [[0, 0, 0, 255], [7, 0, 220, 255], [236, 0, 134, 255], [246, 246, 0, 255], [255, 255, 255, 255]]
        col = np.asarray(col) / 256
        return colors.LinearSegmentedColormap.from_list(cmap, list(zip(nodes, col)))
        
    def _getX(self, x):                                                         # Takes in x pixel position and returns the corresponding value on the x axis
        size = self.fig.get_size_inches()*self.fig.dpi # size in pixels
        pos = np.array(self.ax.get_position())
        lim = np.array(self.ax.get_xlim())
        dx = (lim[1] - lim[0])/((pos[1][0] - pos[0][0])*size[0])
        x = x - pos[0][0]*size[0]
        x = x*dx + lim[0]
        if(x < lim[0]): x = lim[0]
        if(x > lim[1]): x = lim[1]
        return x
    
    def _getY(self, y):                                                         # Takes in y pixel position and returns the corresponding value on the y axis
        size = self.fig.get_size_inches()*self.fig.dpi # size in pixels
        pos = np.array(self.ax.get_position())
        lim = np.array(self.ax.get_ylim())
        dy = (lim[1] - lim[0])/((pos[1][1] - pos[0][1])*size[1])
        y = y - pos[0][1]*size[1]
        y = y*dy + lim[0]
        if(y < lim[0]): y = lim[0]
        if(y > lim[1]): y = lim[1]
        return y
    
    def rotate(self,origin, point, angle):
        """
        Taken from:
        https://stackoverflow.com/questions/34372480/rotate-point-about-another-point-in-degrees-python
        Rotate a point counterclockwise by a given angle around a given origin.
    
        The angle should be given in radians.
        """
        ox, oy = origin
        px, py = point
    
        qx = ox + np.cos(angle) * (px - ox) - np.sin(angle) * (py - oy)
        qy = oy + np.sin(angle) * (px - ox) + np.cos(angle) * (py - oy)
        return np.array(qx, qy)
    
    def browseFile(self):
        path = filedialog.askopenfilename(title='Select File')                  # shows dialog box and return the path\filename
        return path
    
    def browseFolder(self):
        path = filedialog.askdirectory(title='Select Folder')                   # shows dialog box and return the path
        return path
    
    def saveFile(self,initialfile="",title="Save as"):
        path = filedialog.asksaveasfilename(title=title,initialfile=initialfile)
        return path
        
    def exportPNG(self,dpi=0):
        if(not self.init): return
        if(not dpi): dpi = self.dpi
        initialfile = "sim"
        
        path = filedialog.asksaveasfilename(title="Save as",initialfile=initialfile + '.png')
        if path == "":
            return None
        if(not path.endswith('.png') and not path.endswith('.svg')):
            path += '.png'
        self.fig.savefig(path,dpi=dpi,bbox_inches='tight',pad_inches=0)
    
    def exportPickle(self,pklDict,initialfile="pickle"):
        if(not initialfile.endswith('.pk')): initialfile += '.pk'
        path = self.saveFile(initialfile=initialfile)
        if path == "": return
        
        if(not path.endswith('.pk')): path += '.pk'
        pickle.dump(pklDict,open(path,'wb'))
        
    def threadTask(self,func):
        if self.running.is_set():
            self.updateHelpLabel("Simulation already running!")
            print("Something went wrong... sim already running")
            return
        
        self.running.set()
        t = threading.Thread(target=func)
        self.task = t
        t.start()
        
    def stop(self):
        if(self.running.is_set()):
            self.running.clear()
            self.task.join()