'''
Created on May 29, 2013

A module for querying Maple.

@author: ghwatson
'''

import pexpect
from sympy.mpmath import ellipe, ellipf, ellipk, asin

def maple_EllipticE(z,k=None):
    '''
    A wrapper to compute the mpmath complete/incomplete elliptic integral of
    the second kind via the maple definition.
    
    :param z: as defined in maple docs.  If computing complete version, then
    this becomes k.
    :param k: as defined in maple docs.
    '''
    # We deal with the complete and incomplete cases.
    if k is None:
        k = z     # If we only receive 1 argument, then we treat it as k, not z.
        m = k*k
        out = ellipe(m)
    else:
        m = k*k
        phi = asin(z)
        out = ellipe(phi,m)
    return out


def maple_EllipticF(z,k):
    '''
    A wrapper to compute the mpmath incomplete elliptic integral of the first
    kind using the maple definition.
    
    :param z: as defined in maple docs.
    :param k: as defined in maple docs.
    '''
    m = k*k
    phi = asin(z)
    out = ellipf(phi,m)
    return out

def maple_EllipticK(k):
    '''
    A wrapper to compute the mpmath complete elliptic integral of the first
    kind using the maple definition.
    
    :param k: as defined in maple docs.
    '''
    m = k*k
    out = ellipk(m)
    return out

## -----------------------

class MapleLink:
    '''
    A class for Python to communicate with Maple.
    
    Disclaimer: In no way is this currently robustly/generically designed as it
    was created quickly for solving a very particular integral.  Review/revise
    before using elsewhere.
    '''
    #TODO: code this class's structure differently. This dict stuff is ugly/hardcodey.
    
    child = None        
    last_query = None
    replace_dict = {"EllipticF":"maple_EllipticF","EllipticE":"maple_EllipticE","EllipticK":"maple_EllipticK","^":"**"}  
    
    def __init__(self,maple_dir = None):
        self.connect_to(maple_dir)
    
    def connect_to(self, maple_dir):
        '''
        Spawn instance of maple.
        '''
        self.child = pexpect.spawn(maple_dir)
        self.child.expect('#--')

    def disconnect(self):
        '''
        Kill the instance.
        '''
        #TODO: test if this is the right pid code.
        self.child.kill(0)
    
    def query(self,X):
        '''
        Gives a Python-friendly copy of the raw query.
        '''
        maple_str = self.raw_query(X)
        py_str = self._parse_to_python(maple_str)
        return py_str
    
    def raw_query(self,X):
        '''
        Send the string X to maple, and get the output string.
        Note: not entirely raw in that it performs some minor formatting
        (ex: line break characters).
        '''        
        if self.child is None:
            raise ValueError("There is no connection.")
        
        self.child.sendline(X)
        self.child.expect('#--')
        out = self.child.before
        out = out[out.find(';')+1:].strip()
        out = ''.join(out.split('\r\n'))
        out = out.replace('\\','')
        
        self.last_query = out
        return out
      
    def _parse_to_python(self,maple_str):
        '''
        Prepares a maple string for evaluation by Python.          
        NOTE: many cases are  not supported here.  This parser is very
              primitive in that it just searches for strings and replaces
              them.
        '''
          
        # Perform replacements with the dictionary.
        py_str = maple_str
        for maple_entry,py_entry in self.replace_dict.iteritems():
            py_str = py_str.replace(maple_entry,py_entry)
         
        # Add more parsing functionality here if needs grow.
         
        return py_str



    
