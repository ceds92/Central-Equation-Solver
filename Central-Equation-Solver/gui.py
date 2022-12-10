# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 09:24:00 2022

@author: jced0001
"""

import customtkinter as ctk
from MainPanel import MainPanel as mp
import ctypes

ctk.set_appearance_mode("Dark")                                                 # Modes: system (default), light, dark
ctk.set_default_color_theme("blue")                                             # Themes: blue (default), dark-blue, green

class App(ctk.CTk):
    WIDTH = 512
    HEIGHT = 512
    def __init__(self):
        super().__init__()
        self.title("EPWE")
        self.protocol("WM_DELETE_WINDOW", self.on_closing)
        
        dpi = self.winfo_fpixels('1i')
        try:
            scaleFactor = ctypes.windll.shcore.GetScaleFactorForDevice(0)/100   # Account for windows scale factor in display settings
        except:
            scaleFactor = 1                                                     # Might not work on mac - haven't tested
        self.mainPanel = mp(self, width=self.WIDTH, height=self.HEIGHT, dpi=dpi, scaleFactor=scaleFactor)
        self.geometry("%dx%d" % (self.WIDTH, 850))
    def on_closing(self, event=0):
        self.mainPanel.quit()
        
app = App()
app.mainloop()

