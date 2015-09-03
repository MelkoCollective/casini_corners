

#from scipy import optimize, zeros, log, arange
import scipy as sp
import numpy as np

import bresenham

def circle_lattice(L): #borrowed from calculate.py
    '''
    Generates an array of points making up a rasterized disc with integer
    radius L.
    
    :param L: the radius of the disk
    '''
    coord_tuples = bresenham.raster_disk(0,0,L)
    unique_tuples = set(coord_tuples)

    coords = [sp.array(i) for i in unique_tuples]
    return coords

def circle_pointsin(rad, Lx, Ly): #This is Anushya Chandran's circle
    """
    rad: Radius
    Lx: total number of rows
    Ly: total number of columns
    
    Returns a Lx X Ly array with 1s where the circle is.
    """
    
    #First the canonical circle with center coordinates (Lx/2, Ly/2)
    (cols, rows) = np.meshgrid( np.arange(Ly),np.arange(Lx))
    
    distances = np.sqrt( (cols - Ly/2)**2 + (rows - Lx/2)**2 )
    
    mask = distances<rad
    #print(mask)
    mask2 = np.asarray(mask,dtype=int) #mask.astype(int)
    #print(mask2)
    
    #tmpmask = np.roll(mask, 0 - Ly/2, axis=1)
    #return np.roll(mask, 0 - Lx/2, axis=0)    
    return mask2


def cluster_cuts(Cx,Cy,Lx,Ly,latt):
 
   edge_counter = 0  #counts how many clusters have a unique edge
   for y in range (0,Ly-Cy):
      for x in range (0,Lx-Cx):
         sub = latt[y:y+Cy,x:x+Cx]
         #print sub #print lattice 
         if np.sum(sub) != 0 and np.sum(sub) != (Cx*Cy):
            edge_counter += 1
            filename=str(y)+"_"+str(x)+".dat"
            np.savetxt(filename, sub, delimiter=' ', fmt='%i') #save to plain text

   print edge_counter

def main():

   Cx=2;   #the linear dimensions of a cluster
   Cy=2;
   R = 3  #this is the radius of the circle

   Lx = 2*R+2*Cx+1  # this specifies an Lx x Ly lattice to embed the circle in
   Ly = 2*R+2*Cy+1  

   c = circle_lattice( R ) #This creates the Bresenham circle
   lattice = sp.zeros( (Ly,Lx), dtype = 'int' )
   for i in c:
       lattice[ i[0]+R+Cy,i[1]+R+Cx ] = 1 # + count  #This assigns '1' to the region A, offset by R
   np.savetxt('test1.out', lattice, delimiter=' ', fmt='%i') #save to plain text

   #cluster_cuts(Cx,Cy,Lx,Ly,lattice)  #generate all cuts

   #below, compare to the in-R circle

   d = circle_pointsin(R+1,Lx,Ly)
   #lattice2 = sp.zeros( (Ly,Lx))
   #for i in d:
   #    lattice2[ i[0]+R+Cy,i[1]+R+Cx ] = 1 # + count  #This assigns '1' to the region A, offset by R
   np.savetxt('test2.out', d, delimiter=' ', fmt='%i') #save to plain text


if __name__ == '__main__':    
    main()
