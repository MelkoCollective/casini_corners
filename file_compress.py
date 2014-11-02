#This is a simple test of a compression algorithm for the .dat files produced by circletest.py
import numpy as np
import zlib
import binascii
import pickle
import dbm

def compress(fluster):
#takes a numpy array that is the cluster
#returns a hexadecimal key that is the compression of the *string* representing the cluster

    cluster = pickle.dumps(fluster) 
    #print len(cluster), binascii.hexlify(cluster)  #this is a very long key

    compressed = zlib.compress(cluster)
    #print 'Compressed   :', len(compressed), binascii.hexlify(compressed)

    return binascii.hexlify(compressed)  #convert to a hex number from binary


def decompress(cname):
#takes the hexadecimal name representing a cluster
#returns the cluster as a numpy array

    compressed = binascii.unhexlify(cname)  #convert back from hex to binary

    decompressed = zlib.decompress(compressed)
    #print 'Decompressed :', len(decompressed), decompressed

    decomp_array = pickle.loads(decompressed)

    return decomp_array


def main():

    #create the entangling cluster
    filename = 'temp.dat'  #this file was created by circletest.py
    fluster = np.loadtxt(filename, dtype=int)  #read file - could try dtype=bool
    print fluster
    print fluster.shape

    #compress the binary cluster data into a "name" or key
    cname = compress(fluster)
    print len(cname), cname

    #attempt to recover the cluster from the compressed key; check properties
    clust = decompress(cname)
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
