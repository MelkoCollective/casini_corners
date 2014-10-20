

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
    

def main():

   R = 8  #this is the radius of the circle
   L = 20  # this specifies an LxL lattice

   c = circle_lattice( R )
   print c[0]

   lattice = sp.zeros( (L,L), dtype = 'int8' )

   for i in c:
	   lattice[ i[0]+R,i[1]+R  ] = 1  #This assigns '1' to the region A, offset by R

   print lattice #optional

   np.savetxt('test.out', lattice, delimiter=' ', fmt='%i') #save to plain text


if __name__ == '__main__':    
    main()
