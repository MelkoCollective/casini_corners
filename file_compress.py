#This is a simple test of a compression algorithm for the .dat files produced by circletest.py
import numpy as np
import zlib
import binascii

filename = 'temp.dat'  #this file was created by circletest.py

fluster = np.loadtxt(filename, dtype=int)  #read file - could try dtype=bool
cluster = np.array_str(fluster)

print cluster

compressed = zlib.compress(cluster)
print 'Compressed   :', len(compressed), binascii.hexlify(compressed)

e_file= str(binascii.hexlify(compressed)) + ".ent"  #this is the file name for the entropy values
print e_file

#the file should not exist here: test the exception handling
try:
   fh = open(e_file, "r")
except IOError:
   print "Error: FILE DOES NOT EXIST YET"


#here we create the file
fh = open(e_file, "w")
fh.write("your entropy value goes here")
fh.close()

decompressed = zlib.decompress(compressed)
print 'Decompressed :', len(decompressed), decompressed

