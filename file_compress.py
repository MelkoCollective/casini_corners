#This is a simple test of a compression algorithm for the .dat files produced by circletest.py
import numpy as np
import dbm

import key_gen as kg

def main():

    #create the entangling cluster
    filename = 'temp.dat'  #this file was created by circletest.py
    fluster = np.loadtxt(filename, dtype=int)  #read file - could try dtype=bool
    print fluster
    print fluster.shape

    #compress the binary cluster data into a "name" or key
    cname = kg.compress(fluster)
    print len(cname), cname

    #attempt to recover the cluster from the compressed key; check properties
    clust = kg.decompress(cname)
    print clust.size
    print clust.shape
    print clust
    #save to file - can compare to temp2.dat
    np.savetxt("temp2.dat", clust, delimiter=' ', fmt='%i') #save to plain text

    #This is how the database will work
    database = dbm.open('m_n.clust','c')
    database[cname] = str(0.693) #storing the entropy value under key cname
    database.keys()         #get the file keys
    print database[cname]


if __name__ == '__main__':    
    main()
