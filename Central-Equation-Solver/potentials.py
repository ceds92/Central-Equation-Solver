# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 16:59:47 2022

@author: jced0001
"""
import numpy as np

def muffin_primitive(args=[]):
    if(not len(args)):                                                          # This is used to define the input parameters
        params = []
        params.append(['L','r','rs'])                                           # Works best for GUI if done in groups of 3
        params.append(['v0'])
        return params
    
    L,r,rs,v0 = args                                                            # Unpack the arguments. (They will be in the order they are defined above)
    rs = int(rs)                                                                # Convert to int since inputs may be of type float
    
    a1 = np.array([L,0])                                                        # Define the primitive lattice vectors
    a2 = np.array([0,L])
    
    # Define the potential map...
    
    L1 = abs(a1[0])/2 + abs(a2[0])/2
    L2 = abs(a1[1])/2 + abs(a2[1])/2
    
    x1 = np.linspace(-L1,L1,rs)
    x2 = np.linspace(-L2,L2,rs)
    X1,X2 = np.meshgrid(x1,x2)
    R = np.sqrt(X1**2 + X2**2)
    
    V = np.array([R < r],dtype=np.float32)[0]
    
    return np.array([a1,a2]),np.array([x1,x2]),np.flipud(V)*v0                  # Return in this order: primitive lattice vectors, real space defined by [x1,x2], 2D potential map
    
def hexagonal_primitive(args=[]):
    if(not len(args)):
        params = []
        params.append(['L','r','rs'])
        params.append(['v0'])
        return params
    
    pi = np.pi
    
    L,r,rs,v0 = args
    rs = int(rs)
    a1 = np.array([L + L/2,-np.cos(pi/6)*L])
    a2 = np.array([L + L/2,+np.cos(pi/6)*L])
    
    L1 = abs(a1[0])/2 + abs(a2[0])/2
    L2 = abs(a1[1])/2 + abs(a2[1])/2
    
    x1 = np.linspace(-L1,L1,rs)
    x2 = np.linspace(-L2,L2,rs)
    X1,X2 = np.meshgrid(x1,x2)
    
    p = []
    p.append([-L/2,0])
    p.append([L/2,0])
    
    V = np.zeros_like(X1)
    for pp in p:
        R = np.sqrt((X1-pp[0])**2 + (X2+pp[1])**2)
        V += np.array([R < r],dtype=np.float32)[0]
    
    return np.array([a1,a2]),np.array([x1,x2]),np.flipud(V)*v0

def kagome_primitive(args=[]):
    if(not len(args)):
        params = []
        params.append(['L','r','rs'])
        params.append(['v0'])
        return params
    
    L,r,rs,v0 = args
    rs = int(rs)
    T = np.sqrt(3)*L/2

    a1 = np.array([2*L,0])
    a2 = np.array([1*L,2*T])
    
    L1 = abs(a1[0])/2 + abs(a2[0])/2
    L2 = abs(a1[1])/2 + abs(a2[1])/2
    
    x1 = np.linspace(-L1,L1,rs)
    x2 = np.linspace(-L2,L2,rs)
    X1,X2 = np.meshgrid(x1,x2)
    
    p = []
    p.append([-L/2,-T/2])
    p.append([L/2,-T/2])
    p.append([0,T - T/2])
    
    V = np.zeros_like(X1)
    for pp in p:
        R = np.sqrt((X1-pp[0])**2 + (X2+pp[1])**2)
        V += np.array([R < r],dtype=np.float32)[0]
    
    return np.array([a1,a2]),np.array([x1,x2]),V*v0
