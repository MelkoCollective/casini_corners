import numpy as np
from scipy import linalg
import time

start_time = time.time()

perx = pery = False # PBC or not along the x and y direction 
massterm = 0
for i in range(2,4):
    L= 8*i
    Lx = Ly = L
    Ns = Lx * Ly
    def site(x, y):
        return x + y*Lx
    K = np.zeros((Ns, Ns))
    for x in range(0, Lx):
        for y in range(0, Ly):
            K[site(x,y),site(x,y)]= 4.0 + (massterm ** (2))
            xp = (x+1)%Lx
            yp = (y+1)%Ly
            if (xp > x) or perx:
                K[site(xp, y), site(x, y)] = -1.0 
                K[site(x, y), site(xp, y)] = -1.0
            if (yp > y) or pery:
                K[site(x, yp), site(x, y)] = -1.0 
                K[site(x, y), site(x, yp)] = -1.0
    
    Eval,Evec = np.linalg.eigh(K) #, b=None, left=False, right=True, overwrite_a=False, overwrite_b=False, check_finite=True)
    
    P = 1./2. * np.matrix(Evec) * np.matrix(np.diag(np.sqrt(Eval))) * np.matrix(Evec.T)
    X = 1./2. * np.matrix(Evec) * np.matrix(np.diag(1. / np.sqrt(Eval))) * np.matrix(Evec.T)
    LxA , LyA = Lx , Ly/2
    sitesA = []
    for y in range(0, LyA):
        for x in range(0, LxA):
            sitesA.append(site(x, y))
    NsA=len(sitesA)
    Pred = np.zeros((NsA, NsA))
    Xred = np.zeros((NsA, NsA))
    # ............ Following step differs from Mathematica ...........        
    for a in range(0, NsA):
        for b in range(0, NsA):
            Pred[a,b] = P[sitesA[a] , sitesA[b]]
            Xred[a,b] = X[sitesA[a] , sitesA[b]]
    #....end of the step .......
     
    Csquared = (Xred.T).dot(Pred)
    Ev = np.sqrt(np.linalg.eigvals(Csquared))
        
    for ind in range(1, 3):
        renyi = ind /1
        Sn = 0.
        if renyi == 1:
            for j in range(0, NsA):
                if Ev[j] > 0.5:
                    Sn += (Ev[j]+1./2) * np.log(abs(Ev[j]+1./2.))-(Ev[j]-1./2.)*np.log(abs(Ev[j]-1./2)) 
            
        else:
            for j in range(0, NsA):
                if Ev[j] > 0.5:
                    Sn += -1/(1-renyi) * np.log((Ev[j]+1./2.) ** renyi -(Ev[j]-1./2.) ** renyi)
        print "%s %s %s %s" % (LxA, LyA, renyi, Sn)
        print("--- %s seconds ---" % (time.time() - start_time))

    
                
                       
            
            
        
                 
    
        
            
