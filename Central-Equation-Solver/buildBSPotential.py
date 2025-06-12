# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 15:12:06 2023

@author: jced0001
"""
import pickle
import numpy as np
import matplotlib.pyplot as plt
import nanonispy as nap

def _HATTriangle(p,side,theta,x1,x2,tiltV=-1,tiltA=0):
    T = np.sqrt(3)*side/2
    u =[]
    u.append(np.array([-side/2,-T/3]))
    u.append(np.array([side/2,-T/3]))
    u.append(np.array([0,+2*T/3]))
    
    theta = theta*np.pi/180
    for n in range(3):
        u[n] = rotate(np.array([0,0]),u[n],theta)
        u[n] += p
    X1,X2 = np.meshgrid(x1,x2)
    
    xx = np.array([u[0][0],u[1][0]])
    yy = np.array([u[0][1],u[1][1]])
    m,c = np.polyfit(xx,yy,1)
    line1 = m*x1 + c
    HAT = X2 > line1
    
    xx = np.array([u[1][0],u[2][0]])
    yy = np.array([u[1][1],u[2][1]])
    m,c = np.polyfit(yy,xx,1)
    line2 = m*x2 + c
    HAT = HAT & (X1 < line2.reshape(len(line2),1))
    
    xx = np.array([u[2][0],u[0][0]])
    yy = np.array([u[2][1],u[0][1]])
    m,c = np.polyfit(yy,xx,1)
    line3 = m*x2 + c
    HAT = HAT & (X1 > line3.reshape(len(line3),1))
    
    return np.flipud(np.array(HAT))

def _HATCircle(p,side,theta,x1,x2,tiltV=[1,1,1,1],t=0.1):
    # t = 0.1*side
    r = side/6
    
    u =[]
    u.append(np.array([0,0]))
    u.append(np.array([0,+2*r]))
    u.append(rotate([0,0],[0,+2*r],+2*np.pi/3))
    u.append(rotate([0,0],[0,+2*r],-2*np.pi/3))
    
    theta = theta*np.pi/180
    for n in range(4):
        u[n] = rotate(np.array([0,0]),u[n],theta)
        u[n] = u[n] + p
   
    X1,X2 = np.meshgrid(x1,x2)
    
    # hole = r - r/2
    
    t = t*r
    
    R = np.sqrt((X1 + u[0][0])**2 + (X2 + u[0][1])**2)
    # HAT = np.array((R > (r - t)) & (R < (r + t)),dtype=np.float32)*tiltV[0]
    HAT = np.array((R > (r - t)) & (R < (r)),dtype=np.float32)*tiltV[0]
    
    R = np.sqrt((X1 + u[1][0])**2 + (X2 + u[1][1])**2)
    # HAT += ((R > (r - t)) & (R < (r + t)) & (HAT == 0))*tiltV[1]
    HAT += ((R > (r - t)) & (R < (r)) & (HAT == 0))*tiltV[1]
    
    R = np.sqrt((X1 + u[2][0])**2 + (X2 + u[2][1])**2)
    # HAT += ((R > (r - t)) & (R < (r + t)) & (HAT == 0))*tiltV[2]
    HAT += ((R > (r - t)) & (R < (r)) & (HAT == 0))*tiltV[2]
    
    R = np.sqrt((X1 + u[3][0])**2 + (X2 + u[3][1])**2)
    # HAT += ((R > (r - t)) & (R < (r + t)) & (HAT == 0))*tiltV[3]
    HAT += ((R > (r - t)) & (R < (r)) & (HAT == 0))*tiltV[3]
    
    return np.array(HAT)
    
def rotate(origin, point, angle):
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
    return np.array([qx, qy])

def build(rs,A,B,vHAT,vCu,tilt,vs):
    rs = int(rs)
    sxmname  = 'C:/Users/jced0001/Development/Data/Highlights/Ag111/Helium/ncAFM/Bigstar_run6/220616_23-14-05_SPM_sbwatch__HAT-Cu-ncAFM-z-dep-6_-2.1785714285714288e-10_nm_0001.sxm'
    filename = 'C:/Users/jced0001/Development/Data/Highlights/Ag111/Helium/ncAFM/Bigstar_run6/220616_23-14-05_SPM_sbwatch__HAT-Cu-ncAFM-z-dep-6_-2.g80'
    pklDict  = pickle.load(open(filename,'rb'))
    
    sxm = nap.read.Scan(sxmname)
    im = sxm.signals['OC_M1_Freq._Shift']['forward']
    
    # rs = 501
    x1 = np.linspace(0,6e-9,rs)
    x2 = np.linspace(-6e-9,0,rs)
    X1,X2 = np.meshgrid(x1,x2)
    
    # Set up lattice vectors
    scan_angle = 37.5
    a1 = [B*2.8,(-7.0+scan_angle)*np.pi/180] # length,angle
    a2 = [B*2.8,(50.5+scan_angle)*np.pi/180]
    
    yFactor = 2.65/2.8
    a1 = np.array([a1[0]*np.cos(a1[1]),a1[0]*np.sin(a1[1])])*1e-9
    a2 = np.array([a2[0]*np.cos(a2[1]),a2[0]*np.sin(a2[1])])*1e-9
    
    V  = np.zeros((rs,rs))
    
    mainPanel = pklDict['MainPanel']
    
    molFiles = mainPanel['molFiles']    
    molPos = mainPanel['molPos']
    molRot = mainPanel['molRot']
    
    t1 = tilt
    
    t2 = 1 + t1
    t1 = 1 - t1
    t3 = 1
    
    tilt = []
    tilt.append([t3,t1,t2,t2])
    tilt.append([t3,t1,t2,t2])
    tilt.append([t2,t2,t2,t2])
    tilt.append([t3,t2,t2,t1])
    tilt.append([t3,t2,t2,t1])
    tilt.append([t2,t2,t2,t2])
    tilt.append([t3,t2,t1,t2])
    tilt.append([t3,t2,t2,t1])
    
    # plt.figure()
    
    n = 2
    skipHAT = [0,7,10,11,12,13]
    # skipHAT = []
    skipCu  = [0,10,11,15,16]
    
    # for a in np.linspace(-2,2,5):
    #     for b in np.linspace(-2,2,5):
            
    for a in range(1):
        for b in range(1):
            cnt = 0
            for n,pos in enumerate(molPos):
                p = np.array([pos[0],pos[1]/yFactor]) + a*a1 + b*a2
                if("HAT" in molFiles[n]):
                    if(n in skipHAT): continue
                    x = -p[0]
                    y = p[1]
                    V += vHAT*_HATCircle(p=[x,y], side=0.94e-9, theta=molRot[n] + 50, x1=x1, x2=x2,t=1,tiltV=tilt[cnt])
                    cnt += 1
                
                if("Cu" in molFiles[n]):
                    if(n - 14 in skipCu): continue
                    x = -p[0]
                    y = p[1]
                    R = np.sqrt((X1 + x)**2 + (X2 + y)**2)
                    V += vCu*(R < 0.05e-9)
    V += vs
    x1 = np.linspace(-A*3e-9,A*3e-9,rs)
    x2 = np.linspace(-A*3e-9,A*3e-9,rs)
    return np.array([A*a1*1e9,A*a2*1e9]),np.array([x1*1e9,x2*1e9]),V
    # extent = [0,6e-9,0,6e-9]
    # plt.imshow(V,extent=extent)
    # extent2 = [0,6e-9,0,6e-9/yFactor]
    # plt.imshow(im,extent=extent2,alpha=0.75)
    
    # cnt = 0
    # for n,pos in enumerate(molPos):
    #     p = np.array([pos[0],pos[1]/yFactor])
    #     if("HAT" in molFiles[n]):
    #         if(n in skipHAT): continue
    #         x = p[0]
    #         y = p[1]
    #         plt.annotate(str(cnt),[x,y])
    #         cnt += 1