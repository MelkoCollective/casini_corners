#!/usr/bin/python
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

from scipy import optimize, zeros, linspace, log, array, arange

from calculate import Calculate
from maple import MapleLink
import pickle

def main():
    precision = int(sys.argv[1])
    n = int(sys.argv[2])
    maple_dir = sys.argv[3]
    fit_to = int(sys.argv[4])
    load_pickle = bool(int(sys.argv[5]))
    save_pickle = bool(int(sys.argv[6]))
    
    # Create an instance of the class used for computation.
    calc = Calculate()
    
    # Storage to avoid repeated computation of ij correlations. If desired,
    # load in pickled correlators from a previous run.
    saved_correlations = {}
    if load_pickle is True:
        pkl_file = open('corrs.pkl', 'rb')
        saved_correlations = pickle.load(pkl_file)
        pkl_file.close()
        
    # Initiate MapleLink
    maple_link = MapleLink(maple_dir,precision)
    
    # Print presets
    print "START-------------------------------------"
    print "------------------------------------------"
    print "PARAMETERS:"
    print "precision: {0}".format(precision)
    print "Renyi index: {0}".format(n)
    print "fit_to: {0}".format(fit_to)
    print "------------------------------------------"

    # Get the entropy for many different sizes of square lattice.
    sizes = arange(1,fit_to)
    entropies = zeros(sizes.shape)
    for count,L in enumerate(sizes):
        print "---------------Working on lattice size L={0}...".format(L)
        polygon = calc.square_lattice(L)
        X,P,saved_correlations = calc.correlations(polygon,maple_link,precision,saved_correlations,True)
        entropies[count] = calc.entropy(X,P,n,precision, True, True)
        
        # Save the 'saved' correlations to file for later use.
        if save_pickle is True:
            output = open('corrs.pkl','wb')
            pickle.dump(saved_correlations,output,-1)
            output.close()
        
    #TODO: BELOW NOT YET TESTED IN ANY WAY.
    # Fitting.
    print "Performing fit..."
    
    def func_to_fit(L,c0,c1,c2,c3,s_n):
        # Note that the coefficient's names are not the same as Casini's notation.
        return c0 + c1*(L) + c2*(1./(L)) + c3*(1./(L)**2) - s_n*log(L)
    
    p0 = [1,1,1,1,0.05]
    popt, pcov = optimize.curve_fit(func_to_fit,sizes,entropies,p0)

    s_n = popt[4]
    print "Corner coefficient is: " 
    print s_n
    print "covariances are: "
    print pcov


if __name__ == '__main__':    
    main()
