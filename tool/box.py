'''
Created on Jun 20, 2013

Generic tools.

@author: ghwatson
'''

import sympy
from sympy import mpmath as mpm

def simple_quadosc(func, interval, omega, osc, precision):
    '''
    Splits up an integral weighted by cos or sin using the frequency. Note: 
    This function does not multiply the integrand by the trig function.
    
    :param func:        the integrand.
    :param interval:    integral bounds.
    :param omega:       period of trig function.
    :param osc:         'cos' or 'sin' for the weighting trig function.
    '''
    
    with mpm.workdps(precision):
        
        width = mpm.pi/mpm.mpf(omega)
        
        # Initialize leftmost zero. #TODO: Is there a more efficient way of going about this? I feel there must be.
        if osc == 'cos':
            offset = mpm.pi/(mpm.mpf(2)*mpm.mpf(omega))
        elif osc == 'sin':
            offset = mpm.mpf(0)
        else:
            raise ValueError("osc must be 'cos' or 'sin'!")
 
        if interval[0] >= offset:
            steps = mpm.ceil((interval[0]-offset)/width)
            zero = offset + steps*width
        elif interval[0] < offset:
            steps = (mpm.floor((offset - interval[0])/width))
            zero = offset + steps*width
                
        if (False == (interval[0] <= offset <= interval[1])) and (steps == 0):
            integral = mpm.quad(func,interval)
        else:
            integral = 0
            integral += mpm.quad(func,[interval[0],zero])
            while zero + width < interval[1]:
                integral += mpm.quad(func,[zero,zero+width])
                zero += width
            integral += mpm.quad(func,[zero,interval[1]])
    
    return integral