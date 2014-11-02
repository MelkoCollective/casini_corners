#This is a simple test of a compression algorithm for the .dat files produced by circletest.py
import numpy as np
import zlib
import binascii

def compress(fluster):
#takes a numpy array that is the cluster
#returns a hexadecimal key that is the compression of the *string* representing the cluster

    cluster = np.array_str(fluster)  #convert to a string
    #print cluster

    compressed = zlib.compress(cluster)
    print 'Compressed   :', len(compressed), binascii.hexlify(compressed)

    #e_file= str(binascii.hexlify(compressed)) + ".ent"  #this is the file name for the entropy values
    #print e_file

    return binascii.hexlify(compressed)  #convert to a hex number from binary


def decompress(cname):
#takes the hexadecimal name representing a cluster
#returns the cluster as a numpy array

    compressed = binascii.unhexlify(cname)  #convert back from hex to binary

    decompressed = zlib.decompress(compressed)
    #print 'Decompressed :', len(decompressed), decompressed

    #decomp_array = np.array(decompressed)  #convert back from a string to a numpy array

    return decompressed
    #return decomp_array


def main():

    filename = 'temp.dat'  #this file was created by circletest.py
    fluster = np.loadtxt(filename, dtype=int)  #read file - could try dtype=bool
    print fluster
    print fluster.shape

    cname = compress(fluster)
    print cname

    clust = decompress(cname)
    #print clust.shape
    print clust

    #np.savetxt("temp3.dat", clust, delimiter=' ', fmt='%i') #save to plain text

#    #the file should not exist here: test the exception handling
#    try:
#       fh = open(e_file, "r")
#    except IOError:
#       print "Error: FILE DOES NOT EXIST YET"
#    
#    
#    #here we create the file
#    fh = open(e_file, "w")
#    fh.write("your entropy value goes here")
#    fh.close()

if __name__ == '__main__':    
    main()
