import numpy as np
from scipy import linalg
import time
import dbm
import key_gen as kg 
import math

start_time = time.time()

perx = pery = False # PBC or not along the x and y direction 
massterm = 0
for s in range(3, 32):
    for Lx in range(int(math.ceil(s/2.)), s):
        Ly = s-Lx
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
   
        filename="Database_circle/"+"%02d_%02d"%(Lx,Ly)
        db = dbm.open(filename,'w')
        for key in db.keys():
            if db[key] == str(-99):
                clust = np.array(kg.decompress(key))
                print clust
                sitesA = []
                for y in range(0, Ly):
                    for x in range(0, Lx):
                        if clust[y, x] == 1:
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
                Sn = 0. 
                for j in range(0, NsA):
                    if Ev[j] > 0.5:
                        Sn += (Ev[j]+1./2) * np.log(abs(Ev[j]+1./2.))-(Ev[j]-1./2.)*np.log(abs(Ev[j]-1./2)) 
                db[key] = str("%.15f"%Sn)
        db.close()
        print filename , "done"
 
