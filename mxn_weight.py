import numpy
import scipy as sp
import property_gen as pg
def weight(m,n,r,w):
    R = 4*r  #this is the radius of the circle

    Lx = 2*R+2*m+1  # this specifies an Lx x Ly lattice to embed the circle in
    Ly = 2*R+2*n+1

    #c = pg.square_lattice( R )
    c = [[0,0]]

    lattice = sp.zeros( (Ly,Lx), dtype = 'int' )

    for i in c:
        lattice[ i[0]+R+n,i[1]+R+m ] = 1 # + count  #This assigns '1' to the region A, offset by R
    
    
    w_mxn_name = '%02d%02d'%(m,n)


    # First term in weight of mxn is property of mxn
    w[w_mxn_name] = pg.property(m,n,Lx,Ly,lattice)
    #if (m != n): w[w_mxn_name] *= 2
    print w[w_mxn_name],

    wformula = "W%02d%02d=P%02d%02d"%(m,n,m,n)

    for y in range (1,n+1):
        for x in range (y,m+1):
            if (y < n or x < m) and x > 1:
                if x > n: coeff = (m-x+1)*(n-y+1)  # drop last term otherwise get negative coeff
                else: coeff = (m-x+1)*(n-y+1)+(m-y+1)*(n-x+1)
                
                if x==y: coeff = coeff/2   # different coefficents for squares

                wformula += "%+d*W%02d%02d"%(-coeff,x,y)

                w[w_mxn_name] -= coeff * w['%02d%02d'%(x,y)]
                
    print w[w_mxn_name]

    print wformula

    return w
