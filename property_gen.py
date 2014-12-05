
#from scipy import optimize, zeros, log, arange
import scipy as sp
import numpy as np
import dbm
import math
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
    for i in range (-L , L):
       for j in range (-L , L):
           coords.append([i,j])
    return coords


def property(Cx,Cy,Lx,Ly,latt):

   filename="Database/"+"%02d_%02d"%(Cx,Cy)
   database = dbm.open(filename,'c')

   edge_counter = 0  #counts how many clusters have a unique edge
   pmn = 0   #calculates property for mxn cluster
   for y in range (0,Ly-Cy):
      for x in range (0,Lx-Cx):
         sub = latt[y:y+Cy,x:x+Cx]
         #print sub #print lattice 
         if np.sum(sub) != 0 and np.sum(sub) != (Cx*Cy):
            edge_counter += 1
            #fname=str(y)+"_"+str(x)+".dat"
            #np.savetxt(fname, sub, delimiter=' ', fmt='%i') #save to plain text
            cname = fc.compress(sub)
            flag = database.has_key(cname) #check database
            if not flag: 
                database[cname] = str(-99) #storing the entropy value under key cname
                print 'WARNING: DATABASE INCOMPLETE'
            elif database[cname] == str(-99): print 'WARNING: -99 USED IN NLCE' 
            else: pmn += float(database[cname]) #generating property from Database#
   database.close()
   #print edge_counter
   #print pmn
   return pmn
