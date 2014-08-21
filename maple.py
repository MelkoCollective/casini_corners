'''
Created on May 29, 2013

A module for querying Maple.

@author: ghwatson
'''

import pexpect
from sympy.mpmath import ellipe, ellipf, ellipk, asin
from math_parse import math_parser

import types
def copy_func(f, name=None):
    return types.FunctionType(f.func_code, f.func_globals, name or f.func_name,
        f.func_defaults, f.func_closure)

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
    
    This class also contains a parser (from math_parse) that can be used to
    evaluate the strings that Maple outputs. For example, if it outputs the
    string '2+4*x' after performing a computation, then the user can use the 
    eval_stack function to evaluate the expression.
    
    Disclaimer: In no way is this currently robustly/generically designed as it
    was created quickly for solving a very particular integral.  Review/revise
    before using elsewhere.
    '''
    
    child = None        
    last_query = None
    last_good_query = None
    replace_dict = {"EllipticF":"maple_EllipticF","EllipticE":"maple_EllipticE","EllipticK":"maple_EllipticK"}  
    
    def __init__(self,maple_dir = None,precision=15,var_names=['x','y']):
        self.connect_to(maple_dir)
        # Instantiate a parser, with the maple functions added to the grammar.
        extra_binary = {'maple_EllipticF':maple_EllipticF}
        extra_unary = {'maple_EllipticK':maple_EllipticK}
        extra_dyn = {'maple_EllipticE':maple_EllipticE}
        
        self.parser = math_parser(var_names, extra_fn=extra_unary, extra_bin_fn=extra_binary, extra_dynamic_fn=extra_dyn,precision=precision)
        
    
    def connect_to(self, maple_dir):
        '''
        Spawn instance of maple.
        '''
        self.child = pexpect.spawn(maple_dir,timeout=3000,maxread=100000)
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
        self.last_good_query = py_str
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
        self.child.expect('#-->')
        
        # If the buffer is not empty, then that means the input took up more
        # than 1 line. Additional work is required.
        if self.child.buffer != '':
            self.child.expect('#-->')
            out = self.child.before #this required 2 primings of the expect.
            out = ''.join(out.split('\r\n')) 
            out = ''.join(out.split('\r'))
            out = out.replace('\\','')
            out = out.strip() # Assuming no important whitespace on ends?
            while self.child.buffer != '':
                self.child.expect('#-->')
                out+=self.child.before
                out = ''.join(out.split('\r\n'))
                out = ''.join(out.split('\r'))
                out = out.replace('\\','')
                out = out.strip()
#                 print out
        else:
            out = self.child.before
                
        out = out[out.find(';')+1:].strip()
        out = ''.join(out.split('\r\n'))
        out = out.replace('\\','')
        
        self.last_query = out
        return out
    
    def parse(self,string):
        '''
        Parse a string. The internal parser will then create a parsed stack
        which can be evaluated given values for the variables.
        '''
        return self.parser.parse(string)
    
    def get_parsed_stack(self):
        '''
        Get the parsed stack from the internal parser.
        '''
        return self.parser.get_stack()
    
    def eval_stack(self,stack,var_dict):
        '''
        Evaluate a string generated by Maple using the internal parser.
        '''
        return self.parser.value_of_stack(stack,var_dict)
      
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
         
        return py_str
