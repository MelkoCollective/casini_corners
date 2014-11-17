import dbm
import math
import key_gen as kg

'''filename=str(8)+"_"+str(2)
try:
    db = dbm.open(filename)
except dbm.error:
    print "Can't open %s" % (filename)

'''

for s in range(3, 6):
    for Lx in range(int(math.ceil(s/2.)), s):
        Ly = s-Lx
        filename="Database/"+"%02d_%02d"%(Lx,Ly)
        print filename
        db = dbm.open(filename)
        for key in db.keys():
            print kg.decompress(key), db[key]
        db.close()      
