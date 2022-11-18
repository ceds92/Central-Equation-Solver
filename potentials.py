# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 16:59:47 2022

@author: jced0001
"""
import numpy as np

def muffin_top(L,r,rs):
    a1 = np.array([L,0])
    a2 = np.array([0,L])
    
    L1 = abs(a1[0])/2 + abs(a2[0])/2
    L2 = abs(a1[1])/2 + abs(a2[1])/2
    
    x1 = np.linspace(-L1,L1,rs)
    x2 = np.linspace(-L2,L2,rs)
    X1,X2 = np.meshgrid(x1,x2)
    R = np.sqrt(X1**2 + X2**2)
    
    V = np.array([R < r],dtype=np.float32)[0]
    
    return [a1,a2],[x1,x2],V
    
def hexagonal_primitive(L,r,rs):
    pi = np.pi

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
    
    return [a1,a2],[x1,x2],V

def hexagonal_rec(L,r,rs):
    pi = np.pi

    a1 = np.array([2*L + 2*np.arcsin(pi/6)*L,0])
    a2 = np.array([0,2*np.arccos(pi/6)*L])
    
    L1 = abs(a1[0])/2 + abs(a2[0])/2
    L2 = abs(a1[1])/2 + abs(a2[1])/2
    
    x1 = np.linspace(-L1,L1,rs)
    x2 = np.linspace(-L2,L2,rs)
    X1,X2 = np.meshgrid(x1,x2)
    
    p = []
    p.append([-L/2,np.arccos(pi/6)*L/2])
    p.append([+L/2,np.arccos(pi/6)*L/2])
    p.append([-L/2-np.arcsin(pi/6)*L,-np.arccos(pi/6)*L/2])
    p.append([+L/2+np.arcsin(pi/6)*L,-np.arccos(pi/6)*L/2])
    
    V = np.zeros_like(X1)
    for pp in p:
        R = np.sqrt((X1-pp[0])**2 + (X2+pp[1])**2)
        V += np.array([R < r],dtype=np.float32)[0]
    
    return [a1,a2],[x1,x2],V

def kagome_primitive(L,r,rs):
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
    
    return [a1,a2],[x1,x2],V
    
def kagome_rec(L,r,rs):
    pi = np.pi

    a1 = np.array([L + 2*np.arcsin(pi/6)*L,0])
    a2 = np.array([0,4*np.arccos(pi/6)*L])
    
    L1 = abs(a1[0])/2 + abs(a2[0])/2
    L2 = abs(a1[1])/2 + abs(a2[1])/2
    
    x1 = np.linspace(-L1,L1,rs)
    x2 = np.linspace(-L2,L2,rs)
    X1,X2 = np.meshgrid(x1,x2)
    
    p = []
    p.append([0,0])
    
    p.append([-L*np.arcsin(pi/6),+L*np.arccos(pi/6)])
    p.append([+L*np.arcsin(pi/6),+L*np.arccos(pi/6)])
    p.append([-L*np.arcsin(pi/6),-L*np.arccos(pi/6)])
    p.append([+L*np.arcsin(pi/6),-L*np.arccos(pi/6)])
    
    p.append([-2*L*np.arcsin(pi/6),+2*L*np.arccos(pi/6)])
    p.append([+2*L*np.arcsin(pi/6),+2*L*np.arccos(pi/6)])
    p.append([-2*L*np.arcsin(pi/6),-2*L*np.arccos(pi/6)])
    p.append([+2*L*np.arcsin(pi/6),-2*L*np.arccos(pi/6)])
    
    
    V = np.zeros_like(X1)
    for pp in p:
        R = np.sqrt((X1-pp[0])**2 + (X2+pp[1])**2)
        V += np.array([R < r],dtype=np.float32)[0]
    
    return [a1,a2],[x1,x2],V
