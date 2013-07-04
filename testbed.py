'''
Created on Jun 18, 2013

@author: ghwatson
'''
# from __future__ import division

# from scipy import optimize
# from scipy import optimize
# import scipy as sp
# from calculate import calculate_entropy, generate_square_lattice, calculate_correlations
# from maple import MapleLink

if __name__ == '__main__':

# WRITE MATRIX TO FILE ----------------------------------------

#     f1 = open('x.txt','w')
#     f2 = open('p.txt','w')
#      
#     n = 1
#     saved_correlations = {}
#     maple_link = MapleLink("/Library/Frameworks/Maple.framework/Versions/12/bin/maple -tu")
#  
# #     sizes = sp.array([1])
# #     entropies = sp.zeros(sizes.shape)
#     polygon = generate_square_lattice(7)
#  
#     X,P,saved_correlations = calculate_correlations(polygon,maple_link,20,saved_correlations,True)
#      
#     for i in xrange(X.shape[0]):
#         for j in xrange(X.shape[1]):
#             pass
#             f1.write(repr(X[i,j]) + ' ')
#             f2.write(repr(P[i,j]) + ' ')
#         f1.write('\n')
#         f2.write('\n')
#      
#      
#     f1.close()
#     f2.close()
#         

# This code tests simple_quadosc precision----------------------------------------------------------------
    
#     from tool.box import simple_quadosc
#     from sympy import mpmath as mpm
#     mpm.mp.dps = 30
#     
#     omega = 0.5
#     
#     def func(x):
#         return mpm.cos(omega*x)
#     
#     interval = [0,2*mpm.pi]
#     precision = 40
#     osc = 'cos'
#     
#     result = simple_quadosc(func, interval, omega, osc, precision)
# #     print result
#     print '%.40f' % result
#     pass


# This code tests parsing in strings -------
 
# 
#     import re
#     py_str = "blah"
#     def replacer(matchobj):
#         newstr = "mpf('{0}')".format(matchobj.group(0))
#         return newstr
#     py_str = re.sub(r"[-+]?[0-9]*\.?[0-9]+",replacer,py_str)
#             
# ------------------------------------------------------------------------


    # THis integral is solvable in maple...but not in mpmath.  Try different algorithm?
    # Maple gives 6.34499346714920200531149676632 at 30 precision.
    # Scipy gives 6.3449934671492025 with its regular precision (float).
    # Gauss-leg   6.34493993891195322762533970189 which is very inaccurate.
    
#     from mpmath import *
#     import scipy
#     from scipy import integrate
#     mp.dps = 30
#     mystr = "mpf('2')/(mpf('6')-mpf('2')*cos(y))**(mpf('1')/mpf('2'))*ellipk(mpf('2')/(mpf('3')-cos(y)))"
#     # This string is 2/sqrt(6-2*cos(y))*EllipticK(sqrt(2)/sqrt(3-cos(y)))  where the elliptic integral has k as its argument.
#     def func(y):
#         out = mpf(2)/(mpf(6)-mpf(2)*cos(y))**(mpf(1)/mpf(2))*ellipk(mpf('2')/(mpf('3')-cos(y)))
#         return out
#     
#     def funcfloat(y):
#         out = eval(mystr)
#         return float(out)
#     
#     myint = quad(func,[0,pi],verbose=True)
#     print myint
#     
#     myint = quad(funcfloat,[0,pi])
#     print myint
#     
#     mp.dps = 160
#     
#     myint = quad(func,[0,pi])
#     print myint
#     
#     mp.dps = 30
#     myint = quad(func,[0,pi],method='gauss-legendre')
#     print myint
#     
#     myint = quad(func,[0,pi],maxdegree=20)
#     print myint
# 
#     
#     
#     
#     spint = integrate.quad(funcfloat,0,float(pi))
#     print spint



# test of main + my tinkering----------------------------------------------------

#     import sympy
#     import scipy as sp
#     from calculate import Calculate
#     from maple import MapleLink
#       
#       
#     precision = 20
#     sympy.mpmath.mp.dps = precision
#     # The Renyi index.
#     n = 1
#     # Storage for correlations to pass to larger lattice sizes for optimization.
#     saved_correlations = {}
#                        
#     # Initiate MapleLink.
#     maple_link = MapleLink("/Library/Frameworks/Maple.framework/Versions/12/bin/maple -tu")
#       
#     # Get the entropy for many different sizes of a polygon.
#     sizes = sp.linspace(10,100,10)
#     sizes = sp.array([5])
#     entropies = sp.zeros(sizes.shape)
#           
#     for count,L in enumerate(sizes):
#         # The array of points in the polygon defining region V.
#         polygon = Calculate.square_lattice(L)
#               
#         print "Working on lattice size L={0}...".format(L)
#               
#         # Calculate the entropy
#         X,P,saved_correlations = Calculate.correlations(polygon,maple_link,precision,saved_correlations,True)
#          
#         # Check symmetry of XP
# #         XP = X*P
# #         delta = (XP - XP.transpose())
# #         delta.simplify()
# #         print delta
#           
#         # Write matrix to text. -------------------
#           
#         # WRITE MATRIX TO FILE ----------------------------------------
#   
#         f1 = open('x.txt','w')
#         f2 = open('p.txt','w')
#             
#         for i in xrange(X.shape[0]):
#             for j in xrange(X.shape[1]):
#                 pass
#                 f1.write(repr(X[i,j]) + ' ')
#                 f2.write(repr(P[i,j]) + ' ')
#             f1.write('\r\n')
#             f2.write('\r\n')
#             
#             
#         f1.close()
#         f2.close()
#            
#   
# # ------------------------------------------------------------------------
#           
#           
#         entropies[count] = Calculate.entropy(X,P,n,precision,True)
#             
#     print entropies
#          
        
# Test maple querying ------------------------------------------------------

#     from maple import MapleLink
#     import pexpect
#     maple_dir = "/Library/Frameworks/Maple.framework/Versions/12/bin/maple -tu"
#       
#     maple_link = MapleLink("/Library/Frameworks/Maple.framework/Versions/12/bin/maple -tu")
#       
# #     X = "int(1/sqrt(4-2*cos(x)-2*cos(y)),x=0..Pi) assuming y>=0;"
# #     X2 = maple_link.raw_query(X)
# #     X3 = "evalf(int({0},y=0..Pi),30);".format(X2)
# #     X4 = maple_link.raw_query(X3)
# #     X5 = "evalf(int(2/(6-2*cos(y))^(1/2)*EllipticK(2^(1/2)/(3-cos(y))^(1/2)),y=0..Pi),30);"
# #     X6 = maple_link.raw_query(X5)
# #     X7 = "evalf(x-3+x,30);"
# #     X8 = maple_link.raw_query(X7)
#     X9 = "evalf(int(cos(y)+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0+0,y=0..Pi/2),30);"
#   
# #     X10 = maple_link.raw_query(X9)
# #     X11 = maple_link.raw_query("y=0..Pi/2),30);")
#  
#     X12 = "evalf(int(2*cos(y)/(6-2*cos(y))^(1/2)*EllipticK(2^(1/2)/(-cos(y)+3)^(1/2))-4*cos(y)/(6-2*cos(y))^(1/2)*(-1/2+1/2*cos(y))*EllipticK(2^(1/2)/(-cos(y)+3)^(1/2))-4*cos(y)/(6-2*cos(y))^(1/2)*(3/2-1/2*cos(y))*EllipticE(2^(1/2)/(-cos(y)+3)^(1/2)),y=0..Pi),30)"
#          
#       
#     child = pexpect.spawn(maple_dir,maxread=100000)
#       
#     
#     child.expect('#-->')
#       
#     child.sendline(X12)
#     child.expect('#-->')
#       
#     if child.buffer != '':
#         child.expect('#-->')
#         out = child.before #this required 2 primings of the expect.
#           
#         while child.buffer != '':
#             child.expect('#-->')
#             out+=child.before
#           
#     else:
#         out = child.before
#   
#       
#     out = out[out.find(';')+1:].strip()
#     out = ''.join(out.split('\r\n'))
#     out = out.replace('\\','')
#       
#     b = float(out)
#     pass

# Test for precision leaks in the integrals.------------------------------------------------------------------------

# Compare to maple for some precisions.

# L=2 case.

#     from maple import MapleLink
#     import pexpect
#     maple_dir = "/Library/Frameworks/Maple.framework/Versions/12/bin/maple -tu"
#     import sympy.mpmath
#     from sympy.mpmath import cos,extraprec,pi,mpf,ellipk,ellipe
#     from maple import maple_EllipticK,maple_EllipticE
#        
#     maple_link = MapleLink("/Library/Frameworks/Maple.framework/Versions/12/bin/maple -tu")
#     sympy.mpmath.mp.dps = 30
# 
#     [i,j] = [1,1]
#         
# 
#         
# 
# 
#     # Symbolically solve inner integral, using Maple:
#     phi_str = "cos({0}*x)/sqrt(2*(1-cos(x))+2*(1-cos(y)))".format(int(i))
#     phi_integ_str = "int({0},x=0..Pi) assuming y >= 0;".format(phi_str)
#     pi_str = "cos({0}*x)*sqrt(2*(1-cos(x))+2*(1-cos(y)))".format(int(i))
#     pi_integ_str = "int({0},x=0..Pi) assuming y >= 0;".format(pi_str)
#     inner_phi_str = maple_link.query(phi_integ_str)
#     inner_pi_str = maple_link.query(pi_integ_str)
# #             inner_phi_str = maple_link.raw_query(phi_integ_str)
# #             inner_pi_str = maple_link.raw_query(pi_integ_str)
# 
#     # Create function using maple output. #TODO: switch out the eval for a parser when possible. eval is dangerous.
#     def phi_inner_integral(y):
#         out = eval(inner_phi_str)
#         return out*cos(int(j)*y)
#     def pi_inner_integral(y):
#         out = eval(inner_pi_str)
#         return out*cos(int(j)*y)
#      
#     # Perform the outer integrals.
#     phi_integ = sympy.mpmath.quad(extraprec(1000)(phi_inner_integral),[0,pi])
#     pi_integ = sympy.mpmath.quad(extraprec(1000)(pi_inner_integral),[0,pi])
#     
#     print phi_integ
#     print pi_integ

# ----Read in matrix and do eigenvalue stuff.----------------------------------


#     import sympy
#     from sympy import Matrix
#     from sympy.mpmath import mpf
#     import scipy as sp
#     from calculate import calculate_entropy, generate_square_lattice, calculate_correlations
#     from maple import MapleLink
#     
#     n = 1
#     precision = 30
#     mysize = 6
#     sympy.mpmath.mp.dps = precision
#     
#     X = sympy.zeros(mysize**2,mysize**2)
#     P = sympy.zeros(mysize**2,mysize**2)
# 
# #     f1 = open('x.txt','r')
# #     f2 = open('p.txt','r')
#     func = lambda x: x != '\r\n'
#     
# #     f1 = open ( 'x.txt' , 'r')
# #     for line in f1:
# #         test = line.split(' ')
# #         test = filter(func,test)
# # #         map(float,test)
# #          
# # #         map(float,test)
# #         print test
# #         test = [float(i) for i in test]
# #         pass
# # #         float(test)
# 
# 
#     f1 = open ( 'x.txt' , 'r')
#     l1 = [ map(float,filter(func,line.split(' '))) for line in f1 ]
#     X = Matrix(l1)
#     f2 = open ( 'p.txt' , 'r')
#     l2 = [ map(mpf,filter(func,line.split(' '))) for line in f2 ]
#     P = Matrix(l2)
#        
# #     for i in xrange(X.shape[0]):
# #         for j in xrange(X.shape[1]):
# #             pass
# #             X[i,j] = f1.read()
# #             P[i,j] = f2.read()
# #         f1.write('\r\n')
# #         f2.write('\r\n')
#        
#        
#     f1.close()
#     f2.close()
#       
#     
#     # ------------------------------------------------------------------------
#      
#     entropy = calculate_entropy(X,P,n,precision,True)
#            
#     print entropy

    # Initiate MapleLink.
#     from maple import MapleLink
#     from calculate import Calculate
#     precision = 20
#     saved_correlations = None
#     
#     maple_link = MapleLink("/Library/Frameworks/Maple.framework/Versions/12/bin/maple -tu")
#       
#     polygon = Calculate.square_lattice(5)
# 
#     X,P = Calculate.correlations(polygon,maple_link,precision,saved_correlations,True)
#     
#     pass



# precision datacollecting----------------------------------------------------------------------------
# 
#     from scipy import linalg
#     import scipy as sp
#     from math import log
#     from calculate import calculate_entropy, generate_square_lattice, calculate_correlations
#     from maple import MapleLink
#     
#     from itertools import chain
#     from sympy.mpmath import sqrt
#     import sympy
# 
#     
#     # Hardcoded presets.
#     n = 1
#     polygon = generate_square_lattice(2)
#     maple_link = MapleLink("/Library/Frameworks/Maple.framework/Versions/12/bin/maple -tu")
#     
#     # Run the code repeatedly over range of precisions.
# #     precisions = [15,17,20,23,25,28,30,32,35,37,40]
#     precisions = xrange(15,50)
#     
#     # File to save to.
#     f = open('precision.txt','w')
#     
#     eigvals = []
#     for precision in precisions:
# #         sympy.mpmath.mp.dps = precision
#         # Calculate the entropy
#         X,P,precomps = calculate_correlations(polygon, maple_link, precision,{},verbose=True)
#         
#         # This snippet is just from calculate_entropy -----------------
#         # Get the eigenvalues of sqrt(XP)
#         XP = X*P
#         v = XP.eigenvals()
#         sqrt_eigs = []
#         for eig, mult in v.iteritems():
#             sqrt_eigs.append([sqrt(sympy.re(sympy.N(eig,precision)))] * mult)
#         sqrt_eigs = list(chain.from_iterable((sqrt_eigs)))
#         
#         print sqrt_eigs
#         
#         #Scipy eigenvalues.
#         #TODO: see if this works, given mpmath structure.
# #         XPnum = sympy.matrix2numpy(XP)
# #         speigs = linalg.eigvals(XPnum)
# #         sqrtspeigs = sp.sqrt(speigs)
# #         sqrt_eigs = sqrtspeigs
#         
#         # Check that the eigenvalues are well-defined.
#         for eig in sqrt_eigs:
#             if eig <= 0.5:
#                 raise ValueError("At least one of the eigenvalues of sqrt(XP) is below 0.5!")
# 
#         # Calculate entropy.
#         S_n = 0
#         if n == 1:
#             for vk in sqrt_eigs:
#                 S_n += ((vk + 0.5)*log(vk + 0.5) - (vk - 0.5)*log(vk - 0.5))
#         else:
#             for vk in sqrt_eigs:
#                 S_n += log((vk + 0.5)**n - (vk - 0.5)**n)
#             S_n *= 1./(n-1)
#             
#         
#         # -------------------------------------------------------------
#         
#         # save the accumulated eigenvalues, and entropies to file.
#         
#         print 'writing data for precision {0}'.format(float(precision))
#         floateigs = [round(float(x),5) for x in sqrt_eigs]
#         f.write(repr(precision) + ' |||| ' + repr(round(float(S_n),5)) + '  ||||  ' + repr(floateigs))
#         f.write('\n')
#         
#     f.close()
#     
#     print 'done!'

# Get linear chain values ----------------------------------------------------------------------

    from maple import MapleLink
    from calculate import Calculate
    import scipy as sp
    precision = 20
    n = 1
    saved_correlations = None
    
    maple_link = MapleLink("/opt/maple13/bin/maple -tu")
      
    polygons = []
    polygons.append(sp.array([[0,0]]))
    polygons.append(sp.array([[0,0],[1,0]]))
    polygons.append(sp.array([[0,0],[1,0],[2,0]]))
    polygons.append(sp.array([[0,0],[1,0],[2,0],[3,0]]))

    ents = []
    for polygon in polygons:
        # Get this polygon's entropy.
        X,P = Calculate.correlations(polygon,maple_link,precision,saved_correlations,True)
        ents.append(Calculate.entropy(X, P, n, precision))
    
    print ents
