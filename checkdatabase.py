import dbm
import math

for s in range(3, 9):
    for Lx in range(int(math.ceil(s/2.)), s):
        Ly = s-Lx
        filename=str(Lx)+"_"+str(Ly)
        print filename
        db = dbm.open(filename,'c')
        for key in db.keys():
            print db[key]
        db.close()
         
