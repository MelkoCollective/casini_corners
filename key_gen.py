#This is set of functions to generate a key
import numpy as np
import zlib
import binascii
import pickle

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

