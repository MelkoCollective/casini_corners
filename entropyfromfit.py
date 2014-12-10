import numpy as np
from numpy import matrix,linalg,mean
from math import log,sqrt
import sys

lmax = 4
lmin = 3

results_template = "Results_%02d"

def linfit(x_list,y_list): #includes linear fitting as a special case
    Y = matrix([[y_] for y_ in y_list])
    X = matrix([[x_,1] for x_ in x_list])

    (Q,R) = linalg.qr(X)

    (m,b) = [float(item) for item in linalg.solve(R,Q.transpose()*Y)]

    rms_err = sqrt(mean([float(item)**2 for item in Y-X*matrix([[m],[b]])]))

    return (m,b,rms_err)

entropy = []
for i in range (1,9):
    R=2*i
    fname = results_template%(R,)
    f = open(fname, 'r')
    invorder = []
    EE = []
    for lines in f.readlines()[(lmin-3):(lmax-3)]:
        line = lines.split()
        invorder.append(1/float(line[0]))
        EE.append(float(line[1]))
    f.close()
    m,b,err = linfit(invorder,EE)
    entropy.append(b)

f = open("entropy",'w')
for i in range (1,9):
    R = 2*i
    f.write("%02d %.20f\n"%(R,entropy[i-1]))
f.close()

