import numpy
import dbm
#import mxn_getdata
import mxn_weight
from frange import *
from mxn_order import *

#############################
# User settings

order_min = 2
order_max = 16 
order_step = 1
order = Arithmetic()
#############################
#for r in range(2,3):
#r = 1 #radius 4

R = 2  #Raduis 2, diameter 4

filename = "Results_%02d"%(R,)
print "Writing file ",filename
f = open(filename, 'a')
total = 0
w = {} # weights
missing = [] # The list of missing data files

clusters = []

for I in frange(order_min,order_max+0.01,order_step):
    for m,n in order.clusters(I):
        filename= "Database/"+"%02d_%02d"%(m,n)
        try:
            db = dbm.open(filename)
            db.close()
        except dbm.error:
            print "Can't open %s" %(filename)
            missing.append(filename)

        if len(missing) == 0:
            w = mxn_weight.weight(m,n,R,w) # performs cluster weight calculations

        #Embedding factor (1 for squares, 2 for rectangles):
            Lc = 1
            if m != n: Lc = 2   #see line 24, mxn_weight.py
            # cannot use total += w['%02d%02d'%(m,n)] or else W0202 somehow gets changed every iteration
            total = total + Lc*w['%02d%02d'%(m,n)]

    # Save result to file
    #f.write("%d %.15f"%(I,total)+'\n')

# Show all required data files
#print "The following data files are required:"
#for r in required:
#    print "  ",r 

# If any missing data
    if len(missing) > 0:
        print "The following data files were not found:"
        for m in missing:
            print "  ",m
    f.close()
