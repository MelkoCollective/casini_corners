'''
Created on May 27, 2013

@author: ghwatson
'''
# TODO: Replace maplelink with sympy once they have elliptic integrals fully worked in.
#      Pull request here: 
#            https://github.com/sympy/sympy/pull/2165
#      Discussion on adding elliptic integrals here: 
#            http://colabti.org/irclogger/irclogger_log/sympy?date=2013-06-09#l151
#      My original query here:
#            http://stackoverflow.com/questions/16991997/decomposing-integral-into-elliptic-integrals/17013218?noredirect=1#comment24619933_17013218

import sys

from scipy import optimize, zeros, linspace, log

from calculate import Calculate
from maple import MapleLink

if __name__ == '__main__':
    
    precision = int(sys.argv[1])
    n = int(sys.argv[2])
    maple_dir = sys.argv[3]
    
    # Storage to avoid repeated computation of ij correlations.
    saved_correlations = {}
        
    # Initiate MapleLink
    maple_link = MapleLink(maple_dir)

    # Get the entropy for many different sizes of square lattice.
    sizes = linspace(10,100,10)
    entropies = zeros(sizes.shape)
    for count,L in enumerate(sizes):
        print "Working on lattice size L={0}...".format(L)
        polygon = Calculate.square_lattice(L)
        X,P,saved_correlations = Calculate.correlations(polygon,maple_link,precision,saved_correlations,True)
        entropies[count] = Calculate.entropy(X,P,n,precision, True, True)

    #TODO: BELOW NOT YET TESTED IN ANY WAY.
    # Fitting.
    print "Performing fit..."
    
    def func_to_fit(L,c0,c1,c2,c3,s_n):
        # Note that the coefficient's names are not the same as Casini's notation.
        return c0 + c1*L + c2*(1./L) + c3*(1./L**2) - 4*s_n*log(L)
    
    popt, pcov = optimize.curve_fit(func_to_fit,sizes,entropies)

    s_n = popt[4]
    print("Corner coefficient is %d",s_n)
