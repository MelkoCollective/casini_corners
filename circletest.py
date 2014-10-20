

#from scipy import optimize, zeros, log, arange
import scipy as sp
import numpy as np

import bresenham

def square_lattice(L):
    '''
    Generates an array of points making up a square LxL lattice.
    
    :param L: the dimension of the square lattice.
    '''
    x = sp.linspace(0,L,L+1)
    coord_arrays = sp.meshgrid(x,x)
    polygon = (sp.dstack(coord_arrays))
    polygon = sp.reshape(polygon,((L+1)**2,2))    
    return polygon


def circle_lattice(L):
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

   R = 2  #this is the radius of the circle
   L = 8  # this specifies an LxL lattice

   c = circle_lattice( R )
   print c[0]

   lattice = sp.zeros( (L,L), dtype = 'int8' )

   for i in c:
	   lattice[ i[0]+R,i[1]+R  ] = 1  #This assigns '1' to the region A, offset by R

   print lattice #optional

   np.savetxt('test.out', lattice, delimiter=' ', fmt='%i') #save to plain text


if __name__ == '__main__':    
    main()
