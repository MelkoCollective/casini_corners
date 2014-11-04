
#generates and empty database for m_n clusters
#from scipy import optimize, zeros, log, arange
import scipy as sp
import numpy as np
import dbm

import bresenham
import key_gen as fc

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

def square_lattice(L): # generates lattice coordinates inside a square
    coords = []
    for i in range (-L/2 , L/2):
       for j in range (-L/2 , L/2):
           coords.append([i,j])
    return coords


def cluster_cuts(Cx,Cy,Lx,Ly,latt):

   filename=str(Cx)+"_"+str(Cy)
   database = dbm.open(filename,'c')

   edge_counter = 0  #counts how many clusters have a unique edge
   for y in range (0,Ly-Cy):
      for x in range (0,Lx-Cx):
         sub = latt[y:y+Cy,x:x+Cx]
         #print sub #print lattice 
         if np.sum(sub) != 0 and np.sum(sub) != (Cx*Cy):
            edge_counter += 1
            filename=str(y)+"_"+str(x)+".dat"
            #np.savetxt(filename, sub, delimiter=' ', fmt='%i') #save to plain text
            cname = fc.compress(sub)
            #print cname
            #look in the database file Cx_Cy.db
            database[cname] = str(-99) #storing the entropy value under key cname

   database.close()
   print edge_counter


def main():

   Cx=7;   #the linear dimensions of a cluster
   Cy=1;
   R = 4  #this is the radius of the circle

   Lx = 2*R+2*Cx+1  # this specifies an Lx x Ly lattice to embed the circle in
   Ly = 2*R+2*Cy+1  

   c = square_lattice( R )

   lattice = sp.zeros( (Ly,Lx), dtype = 'int' )

   for i in c:
       lattice[ i[0]+R+Cy,i[1]+R+Cx ] = 1 # + count  #This assigns '1' to the region A, offset by R
   print lattice
   cluster_cuts(Cx,Cy,Lx,Ly,lattice)  #generate all cuts

   np.savetxt('test.out', lattice, delimiter=' ', fmt='%i') #save to plain text

   filename=str(Cx)+"_"+str(Cy)
   db = dbm.open(filename,'c')

   print db.keys()
   for key in db.keys():
      print(db[key])

   db.close()

if __name__ == '__main__':    
    main()
