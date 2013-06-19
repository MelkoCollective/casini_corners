'''
Created on May 27, 2013

@author: ghwatson
'''

# from scipy import optimize
from scipy import optimize
import scipy as sp
from calculate import calculate_entropy, generate_square_lattice, calculate_correlations
from maple import MapleLink

if __name__ == '__main__':
    
    #TODO: Remove these hardcodings, and make them args. (n, maple_link, resolution)
    
    # The Renyi index.
    n = 1
    # Storage for correlations to pass to larger lattice sizes for optimization.
    saved_correlations = {}
        
    #TODO: Replace maplelink with sympy once they have elliptic integrals fully worked in.
    #      Pull request here: 
    #            https://github.com/sympy/sympy/pull/2165
    #      Discussion on adding elliptic integrals here: 
    #            http://colabti.org/irclogger/irclogger_log/sympy?date=2013-06-09#l151
    #      My original query here:
    #            http://stackoverflow.com/questions/16991997/decomposing-integral-into-elliptic-integrals/17013218?noredirect=1#comment24619933_17013218
    
    # Initiate MapleLink.
    maple_link = MapleLink("/Library/Frameworks/Maple.framework/Versions/12/bin/maple -tu")

    # Get the entropy for many different sizes of a polygon.
#     sizes = sp.linspace(10,100,10)
    sizes = sp.array([2,4,7,10,15,20,30,40])
    sizes = sp.array([2,3,4,5,6])
    entropies = sp.zeros(sizes.shape)
    
    for count,L in enumerate(sizes):
        # The array of points in the polygon defining region V.
        polygon = generate_square_lattice(L)
        #DEBUGGING: hardwiring in desired L
#         polygon = generate_square_lattice(2)
        
        print "Working on lattice size L={0}...".format(L)
        
        # Calculate the entropy
        X,P,saved_correlations = calculate_correlations(polygon,maple_link,35,saved_correlations,True)
        entropies[count] = calculate_entropy(X,P,n,True)

    #TODO: BELOW NOT YET TESTED IN ANY WAY.
    
    # Take all results and perform a fitting to find s_n:
    print "Performing fit..."
    
    def func_to_fit(L,c0,c1,c2,c3,s_n):
        # Note that the coefficient's names are not the same as Casini's notation.
        return c0 + c1*L + c2*(1./L) + c3*(1./L**2) - 4*s_n*sp.log(L)
    
    popt, pcov = optimize.curve_fit(func_to_fit,sizes,entropies)

    s_n = popt[4]
    print("Corner coefficient is %d",s_n)
