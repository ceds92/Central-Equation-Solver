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

def HATCu_bigstar_CO(args=[]):
    import nanonispy as nap
    from PIL import Image
    import matplotlib.pyplot as plt
    
    if(not len(args)):
        params = []
        params.append(['rs'])                                                   # rs = real-space sampling
        params.append(['vs','v0'])                                              # vs = SS onset, vHAT = Potential for HAT, vCu = Potential for Cu
        return params
    
    rs,vs,v0 = args
    
    sxmname = 'C:/Users/jced0001/Development/Data/Highlights/Ag111/Helium/CO Functionalised/220610_14-52-03_SPM_sbwatch_0001.sxm'
    pngname = 'C:/Users/jced0001/Development/Data/Highlights/Ag111/Helium/CO Functionalised/BS_CO_14-52-03.png'
    
    sxm   = nap.read.Scan(sxmname)
    image = np.sum(np.asarray(Image.open(pngname)),2)
    
    lxy = sxm.header['scan_range']                                                  # Real length of scan (m)
    dxy = lxy/image.shape                                                           # Wrong way but doesn't matter for square image
    # dxy[1] *= 1.07
    dxy *= 1e9
    
    pts = []
    pts.append(np.array([175,115])) # [x,y]/[col,row]
    pts.append(np.array([352,115]))
    pts.append(np.array([263,259]))
    pts.append(np.array([86,259]))
    pts.append(np.array([175,115])) # pt 1 again
    pts = np.array(pts)
    
    direction = [1,-1,-1,1]
    k1 = np.arange(image.shape[0])
    k2 = np.arange(image.shape[1])
    K1,K2 = np.meshgrid(k1,k2)
    UC = np.zeros_like(image)
    for n in range(len(pts)-1):
        m  = pts[n+1][1] - pts[n][1]
        m /= pts[n+1][0] - pts[n][0]
        c  = pts[n][1] - m*pts[n][0]
        if(direction[n] > 0):
            UC += K2 > m*K1 + c
        else:
            UC += K2 < m*K1 + c
        plt.plot(pts[n][0],pts[n][1],'.',markersize=12)
    
    # image[UC < 4] = np.max(image[UC < 4])
    # image[image < 500] = np.max(image)
    # image = -image
    # image -= np.min(image)
    # image = image[pts[0][1] - 61:pts[2][1] + 61,pts[3][0]:pts[1][0]]
    
    # image[UC < 4] = np.min(image[UC < 4])
    # image -= np.min(image)
    # image = image[pts[0][1] - 61:pts[2][1] + 61,pts[3][0]:pts[1][0]]
    
    pores = image < 500
    image = -image
    image[UC < 4] = np.min(image[UC < 4])
    image -= np.min(image)
    image[pores] = np.min(image)
    image = image[pts[0][1] - 61:pts[2][1] + 61,pts[3][0]:pts[1][0]]

    V = v0*image/np.max(image)
    
    a1 = (pts[1] - pts[0])*dxy
    a2 = (pts[0] - pts[3])*dxy
    a2[1] = -a2[1]
    
    x1 = np.arange(pts[1][0] - pts[3][0])*dxy[0]
    x2 = np.arange(pts[2][1] - pts[0][1] + 122)*dxy[1]
    x1 = x1 - x1[-1]/2
    x2 = x2 - x2[-1]/2
    
    return np.array([a1,a2]),np.array([x1,x2]),V

def HATCu_ncAFM(args=[]):
    import buildBSPotential as bs
    if(not len(args)):
        params = []
        params.append(['rs','A','B'])                                               # L = Nearest neighbour distance, r = length of HAT side, rs = real-space sampling
        params.append(['vHAT','vCu'])                                           # t = thickness, vHAT = Potential for HAT, vCu = Potential for Cu
        params.append(['tilt','vs'])                                            # Tilt strength angle of HATs (0 to 1), vs = SS onset
        return params
    
    # rs,tilt,vHAT,vCu,vs = args
    
    return bs.build(*args)

def HATCu_bigstar(args=[]):
    if(not len(args)):
        params = []
        params.append(['L','r','rs'])                                           # L = Nearest neighbour distance, r = length of HAT side, rs = real-space sampling
        params.append(['t','vHAT','vCu'])                                       # t = thickness, vHAT = Potential for HAT, vCu = Potential for Cu
        params.append(['tilt','vs'])                                            # Tilt strength angle of HATs (0 to 1), vs = SS onset
        return params
    
    pi = np.pi

    # L,r,rs,v0,v1,v2,vCu,t1 = args
    L,r,rs,t,vHAT,vCu,t1,vs = args
    
    # t2 = (1+t1)/2
    
    t2 = 1 + t1
    t1 = 1 - t1
    t3 = 1
    # t3 = t2
    
    rs = int(rs)
    a1 = 2*np.array([L + L/2,-np.cos(pi/6)*L])
    a2 = 2*np.array([L + L/2,+np.cos(pi/6)*L])
    print("BS:")
    print(a1)
    print(a2)
    T = np.sqrt(3)*L/2
    # a1 = np.array([4*T,-2*L])
    # a2 = np.array([4*T,2*L])
    
    L1 = abs(a1[0])/2 + abs(a2[0])/2
    L2 = abs(a1[1])/2 + abs(a2[1])/2
    
    x1 = np.linspace(-L1,L1,rs)
    x2 = np.linspace(-L2,L2,rs)
    X1,X2 = np.meshgrid(x1,x2)
    
    p = []
    # tilt = np.array([t2,1,1,t1])
    tilt = np.array([t3,t2,t2,t1])
    p.append([+L/2,-np.cos(pi/6)*L,vHAT,+30,tilt.copy()])                         # Lower left
    # tilt = np.array([t2,t1,1,1])
    tilt = np.array([t3,t1,t2,t2])
    p.append([-L/2,-np.cos(pi/6)*L,vHAT,-30,tilt.copy()])                         # Lower right
    
    # tilt = np.array([t2,1,t1,1])
    tilt = np.array([t3,t2,t1,t2])
    p.append([+L/2,+np.cos(pi/6)*L,vHAT,+30,tilt.copy()])                         # Upper left
    # tilt = np.array([t2,1,t1,1])
    tilt = np.array([t3,t2,t1,t2])
    p.append([-L/2,+np.cos(pi/6)*L,vHAT,-30,tilt.copy()])                         # Upper right
    
    # tilt = np.array([t3,t3,t3,t3])
    tilt = np.array([t2,t2,t2,t2])
    p.append([-L/1,0,vHAT,+30,tilt.copy()])                                       # Inner left
    p.append([+L/1,0,vHAT,-30,tilt.copy()])                                       # Inner right
    
    # tilt = np.array([t2,t1,1,1])
    tilt = np.array([t3,t1,t2,t2])
    p.append([+2*L,0,vHAT,+30,tilt.copy()])                                       # Far left
    # tilt = np.array([t2,1,1,t1])
    tilt = np.array([t3,t2,t2,t1])
    p.append([-2*L,0,vHAT,-30,tilt.copy()])                                       # Far right
    p = np.array(p)
    
    V = np.zeros_like(X1) + vs
    for ppp in p:
        pp = ppp[0:2]
        v  = ppp[2]
        theta = ppp[3]
        tilt = ppp[4]
        V += v*_HATCircle(pp,r,theta,x1,x2,tiltV=tilt,t=t)
    
    cus = []
    cus.append(np.array([0,+np.cos(pi/6)*L]))
    cus.append(rotate([0,0],[0,+np.cos(pi/6)*L],pi/3))
    cus.append(rotate([0,0],[0,+np.cos(pi/6)*L],2*pi/3))
    cus.append(rotate([0,0],[0,+np.cos(pi/6)*L],3*pi/3))
    cus.append(rotate([0,0],[0,+np.cos(pi/6)*L],4*pi/3))
    cus.append(rotate([0,0],[0,+np.cos(pi/6)*L],5*pi/3))
    
    cus.append(np.array([-2*L + L/2,0]))
    cus.append(np.array([+2*L - L/2,0]))
    
    cus.append(rotate([0,0],[+2*L - L/2,0],-pi/3))
    cus.append(rotate([0,0],[+2*L - L/2,0],pi/3))
    
    hatcu = rotate([0,0],[0,+np.cos(pi/6)*L],4*pi/3) - np.array([L/2,0])
    cus.append(np.array([+2*L,0]) + hatcu)
    
    hatcu = rotate([0,0],[0,+np.cos(pi/6)*L],5*pi/3) - np.array([L/2,0])
    cus.append(np.array([+2*L,0]) + hatcu)
    for cu in cus:
        R  = np.sqrt((X1 + cu[0])**2 + (X2 + cu[1])**2)
        V += vCu*(R < r/10)
        
    return np.array([a1,a2]),np.array([x1,x2]),np.flipud(V)

def HATCu_bigstar_data(args=[]):
    if(not len(args)):
        params = []
        params.append(['L','r','rs'])                                           # L = Nearest neighbour distance, r = length of HAT side, rs = real-space sampling
        params.append(['t','vHAT','vCu'])                                       # t = thickness, vHAT = Potential for HAT, vCu = Potential for Cu
        params.append(['tilt','vs'])                                            # Tilt strength angle of HATs (0 to 1), vs = SS onset
        return params
    
    pi = np.pi

    L,r,rs,t,vHAT,vCu,t1,vs = args
    
    t2 = 1 + t1
    t1 = 1 - t1
    t3 = 1
    
    rs = int(rs)
    a1 = 2*np.array([L + L/2,-np.cos(pi/6)*L])
    a2 = 2*np.array([L + L/2,+np.cos(pi/6)*L])
    
    L1 = abs(a1[0])/2 + abs(a2[0])/2
    L2 = abs(a1[1])/2 + abs(a2[1])/2
    
    x1 = np.linspace(-L1,L1,rs)
    x2 = np.linspace(-L2,L2,rs)
    X1,X2 = np.meshgrid(x1,x2)
    
    p = []
    ox,oy = [0.1,-0.05]
    tilt = np.array([t3,t2,t2,t1])
    p.append([+L/2 - ox,-np.cos(pi/6)*L + oy,vHAT,+30,tilt.copy()])                         # Lower left
    
    ox,oy = [0.1,0]
    tilt = np.array([t3,t1,t2,t2])
    p.append([-L/2 - ox,-np.cos(pi/6)*L + oy,vHAT,-30,tilt.copy()])                         # Lower right
    
    ox,oy = [-0.1,-0.075]
    tilt = np.array([t3,t2,t1,t2])
    p.append([+L/2 - ox,+np.cos(pi/6)*L + oy,vHAT,+30,tilt.copy()])                         # Upper left
    
    ox,oy = [-0.1,0]
    tilt = np.array([t3,t2,t1,t2])
    p.append([-L/2 - ox,+np.cos(pi/6)*L - oy,vHAT,-30,tilt.copy()])               # Upper right
    
    ox,oy = [0,0.1]
    tilt = np.array([t2,t2,t2,t2])
    p.append([-L/1 + ox,0 + oy,vHAT,+30,tilt.copy()])                             # Inner left
    p.append([+L/1 - ox,0 - oy,vHAT,-30,tilt.copy()])                             # Inner right
    
    ox,oy = [0,0]
    tilt = np.array([t3,t1,t2,t2])
    p.append([+2*L - ox,0 + oy,vHAT,+30,tilt.copy()])                                       # Far right
    
    ox,oy = [-0.025,-0.08]
    tilt = np.array([t3,t2,t2,t1])
    p.append([-2*L - ox,0 + oy,vHAT,-30,tilt.copy()])                                       # Far left
    
    p = np.array(p)
    
    V = np.zeros_like(X1) + vs
    for ppp in p:
        pp = ppp[0:2]
        v  = ppp[2]
        theta = ppp[3]
        tilt = ppp[4]
        V += v*_HATCircle(pp,r,theta,x1,x2,tiltV=tilt,t=t)
    
    cus = []
    ox,oy = [0.1,-0.03]
    cus.append(np.array([0 - ox,+np.cos(pi/6)*L - oy]))                                   # Bottom hex
    
    oxy = np.array([0.075,0])
    cus.append(rotate([0,0],[0,+np.cos(pi/6)*L],pi/3) - oxy)
    
    oxy = np.array([-0.025,0.025])
    cus.append(rotate([0,0],[0,+np.cos(pi/6)*L],2*pi/3) - oxy)
    
    oxy = np.array([-0.1,-0.03])
    cus.append(rotate([0,0],[0,+np.cos(pi/6)*L],3*pi/3) - oxy)                        # Top
    
    oxy = np.array([-0.05,-0.075])
    cus.append(rotate([0,0],[0,+np.cos(pi/6)*L],4*pi/3) - oxy)
    
    oxy = np.array([0.09,-0.06])
    cus.append(rotate([0,0],[0,+np.cos(pi/6)*L],5*pi/3) - oxy)
    
    cus.append(np.array([-2*L + L/2,0]))
    cus.append(np.array([+2*L - L/2,0]))
    
    cus.append(rotate([0,0],[+2*L - L/2,0],-pi/3))
    cus.append(rotate([0,0],[+2*L - L/2,0],pi/3))
    
    hatcu = rotate([0,0],[0,+np.cos(pi/6)*L],4*pi/3) - np.array([L/2,0])
    cus.append(np.array([+2*L,0]) + hatcu)
    
    hatcu = rotate([0,0],[0,+np.cos(pi/6)*L],5*pi/3) - np.array([L/2,0])
    cus.append(np.array([+2*L,0]) + hatcu)
    for cu in cus:
        R  = np.sqrt((X1 + cu[0])**2 + (X2 + cu[1])**2)
        V += vCu*(R < r/10)
        
    return np.array([a1,a2]),np.array([x1,x2]),np.flipud(V)

def HATCu_smallstar(args=[]):
    if(not len(args)):
        params = []
        params.append(['L','r','rs'])                                           # L = Nearest neighbour distance, r = length of HAT side, rs = real-space sampling
        params.append(['t','vHAT','vCu'])                                       # t = thickness, vHAT = Potential for HAT, vCu = Potential for Cu
        params.append(['tilt','vs'])                                            # Tilt strength angle of HATs (0 to 1), vs = SS onset
        return params
    
    pi = np.pi

    L,r,rs,t,vHAT,vCu,t1,vs = args
    
    t2 = 1 + t1
    t1 = 1 - t1
    # t3 = 1
    t3 = t2
    
    rs = int(rs)
    # a1 = 2*np.array([L + L/2,-np.cos(pi/6)*L])
    # a2 = 2*np.array([L + L/2,+np.cos(pi/6)*L])
    
    T = np.sqrt(3)*L/2
    a1 = np.array([3*L,-2*T])
    a2 = np.array([3*L,2*T])
    print("SS:")
    print(a1)
    print(a2)
    
    L1 = abs(a1[0])/2 + abs(a2[0])/2
    L2 = abs(a1[1])/2 + abs(a2[1])/2
    
    x1 = np.linspace(-L1,L1,rs)
    x2 = np.linspace(-L2,L2,rs)
    X1,X2 = np.meshgrid(x1,x2)
    
    p = []
    tilt = np.array([1,t2,t1,t2])
    upper = [0,T+np.sqrt(3)*L/6]
    p.append([*upper,vHAT,60,tilt.copy()])                                      # Upper
    
    pos = rotate([0,0],upper,pi/3)
    tilt = np.array([1,t2,t1,t2])
    p.append([*pos,vHAT,0,tilt.copy()])                                         # Upper right
    
    pos = rotate([0,0],upper,2*pi/3)
    tilt = np.array([1,t2,t2,t1])
    p.append([*pos,vHAT,60,tilt.copy()])                                        # Bottom right
    
    pos = rotate([0,0],upper,3*pi/3)
    tilt = np.array([1,t2,t2,t1])
    p.append([*pos,vHAT,0,tilt.copy()])                                         # Bottom
    
    pos = rotate([0,0],upper,4*pi/3)
    tilt = np.array([1,t1,t2,t2])
    p.append([*pos,vHAT,60,tilt.copy()])                                        # Bottom left
    
    pos = rotate([0,0],upper,5*pi/3)
    tilt = np.array([1,t1,t2,t2])
    p.append([*pos,vHAT,0,tilt.copy()])                                         # Upper left
    p = np.array(p)
    
    V = np.zeros_like(X1) + vs
    for ppp in p:
        pp = ppp[0:2]
        v  = ppp[2]
        theta = ppp[3]
        tilt = ppp[4]
        V += v*_HATCircle(pp,r,theta,x1,x2,tiltV=tilt,t=t)
    
    cus = []
    cus.append(np.array([L,0]))
    cus.append(rotate([0,0],cus[0],pi/3))
    cus.append(rotate([0,0],cus[0],2*pi/3))
    cus.append(rotate([0,0],cus[0],3*pi/3))
    cus.append(rotate([0,0],cus[0],4*pi/3))
    cus.append(rotate([0,0],cus[0],5*pi/3))
    
    cus.append(np.array((a1)/2))
    cus.append(np.array((a2)/2))
    cus.append(np.array((a1+a2)/2))
    
    for cu in cus:
        R  = np.sqrt((X1 + cu[0])**2 + (X2 + cu[1])**2)
        V += vCu*(R < r/10)
    
    print("maxV:",np.max(V))
    return np.array([a1,a2]),np.array([x1,x2]),np.flipud(V)

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
    HAT = np.array((R > (r - t)) & (R < (r)),dtype=np.float32)*tiltV[0]
    # HAT = np.array((R > (r - t)) & (R < (r + t)),dtype=np.float32)*tiltV[0]
    
    R = np.sqrt((X1 + u[1][0])**2 + (X2 + u[1][1])**2)
    HAT += ((R > (r - t)) & (R < (r)) & (HAT == 0))*tiltV[1]
    # HAT += ((R > (r - t)) & (R < (r + t)) & (HAT == 0))*tiltV[1]
    
    R = np.sqrt((X1 + u[2][0])**2 + (X2 + u[2][1])**2)
    HAT += ((R > (r - t)) & (R < (r)) & (HAT == 0))*tiltV[2]
    # HAT += ((R > (r - t)) & (R < (r + t)) & (HAT == 0))*tiltV[2]
    
    R = np.sqrt((X1 + u[3][0])**2 + (X2 + u[3][1])**2)
    HAT += ((R > (r - t)) & (R < (r)) & (HAT == 0))*tiltV[3]
    # HAT += ((R > (r - t)) & (R < (r + t)) & (HAT == 0))*tiltV[3]
    
    # R = np.sqrt((X1 + u[0][0])**2 + (X2 + u[0][1])**2)
    # HAT = np.array((R > (r - hole)) & (R < (r + t/2)),dtype=np.float32)*tiltV[0]
    
    # R = np.sqrt((X1 + u[1][0])**2 + (X2 + u[1][1])**2)
    # HAT += ((R > (r - hole)) & (R < (r + t/2)) & (HAT == 0))*tiltV[1]
    
    # R = np.sqrt((X1 + u[2][0])**2 + (X2 + u[2][1])**2)
    # HAT += ((R > (r - hole)) & (R < (r + t/2)) & (HAT == 0))*tiltV[2]
    
    # R = np.sqrt((X1 + u[3][0])**2 + (X2 + u[3][1])**2)
    # HAT += ((R > (r - hole)) & (R < (r + t/2)) & (HAT == 0))*tiltV[3]
    
    return np.flipud(np.array(HAT))
    
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