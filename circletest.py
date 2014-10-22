

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


def cluster_cuts(Lx,Ly,L,latt):
 
   edge_counter = 0  #counts how many clusters have a unique edge
   for y in range (0,L-Ly):
      for x in range (0,L-Lx):
         sub = latt[y:y+Ly,x:x+Lx]
         #print sub #print lattice 
         if np.sum(sub) != 0 and np.sum(sub) != (Lx*Ly):
            edge_counter += 1
            filename=str(y)+"_"+str(x)+".dat"
            np.savetxt(filename, sub, delimiter=' ', fmt='%i') #save to plain text

   print edge_counter

def main():

   Lx=5;   #the linear dimensions of a cluster
   Ly=5;
   R = 10  #this is the radius of the circle

   L = 2*R+2*Lx+1  # this specifies an LxL lattice to embed the circle in

   c = circle_lattice( R )
   print c[0]

   lattice = sp.zeros( (L,L), dtype = 'int' )

   #count = 0
   #print count
   for i in c:
       lattice[ i[0]+R+Lx,i[1]+R+Ly ] = 1 # + count  #This assigns '1' to the region A, offset by R
       #count += 1
       #print count

   cluster_cuts(Lx,Ly,L,lattice)

   np.savetxt('test.out', lattice, delimiter=' ', fmt='%i') #save to plain text


if __name__ == '__main__':    
    main()
