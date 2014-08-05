'''
Created on July 2, 2014

Original code by paul McGuire (see old header below). Wrapped in a class and 
altered to allow variables defined at evaluation time, and the ability to parse
functions with more than one argument.

Strings can be evaluated to arbitrary precision, with the use of sympy.mpmath.

Note: The implementation only handles unary and binary functions, but can be
easily expanded.

BUGS:
- There is a bug where if you define a function starting with 'E' or 'PI',
  then the parser will mess up.
- a string like '-(-X)' where X is some expression(s) does not currently work.

@author: ghwatson
'''

# fourFn.py
#
# Demonstration of the pyparsing module, implementing a simple 4-function expression parser,
# with support for scientific notation, and symbols for e and pi.
# Extended to add exponentiation and simple built-in functions.
# Extended test cases, simplified _pushFirst method.
#
# Copyright 2003-2006 by Paul McGuire
#
from pyparsing import Literal,CaselessLiteral,Word,Combine,Group,Optional,\
    ZeroOrMore,Forward,nums,alphas,delimitedList,MatchFirst
import math
import operator
from functools import partial
from itertools import chain
from copy import deepcopy
import sympy.mpmath as mpm

# Parses mathematical expressions. The user is capable of expanding the grammar
# if desired. Currently only handles unary and binary functions.
class math_parser:
    
    def __init__( self, var_names=None, extra_opn=None, extra_fn=None, extra_bin_fn=None,extra_dynamic_fn=None,precision=None ):
        if precision:
            mpm.mp.dps = precision
        self.mpfing = True # Flag to switch between casting numbers with mpf or float
        
        # For parsing
        self.exprStack = [] # The parsed results, ordered for recursive
                            # evaluation.
        self.bnf = None     # The grammatical structure (Backus-Naur form)
        self.results = None # The parsed results
                
        # Correspondences between strings and mathematical syntax. The user
        # can tweak these, but this provides some basic syntax.
        if var_names:
            self.var_names = var_names
        else: # Default labels
            self.var_names = ['x','y','z']
        self.epsilon = 1e-12
        self.opn = { "+" : operator.add,
                "-" : operator.sub,
                "*" : operator.mul,
                "/" : operator.truediv,
                "^" : operator.pow }
        self.fn  = { "sin" : mpm.sin,
                "cos" : mpm.cos,
                "tan" : mpm.tan,
                "abs" : abs,
                "trunc" : lambda a: int(a),
                "round" : round,
                "sgn" : lambda a: abs(a)>self.epsilon and cmp(a,0) or 0}
        # binary functions
        self.bin_fn = {}
        # functions that can take on a variable number of arguments.
        self.dynamic_fn = {}
    
        # The user can append additional elements to the syntax dictionaries
        if extra_opn:
            self.opn = dict(chain(self.opn.iteritems(), extra_opn.iteritems()))
        if extra_fn:
            self.fn = dict(chain(self.fn.iteritems(), extra_fn.iteritems()))
        if extra_bin_fn:
            self.bin_fn = dict(chain(self.bin_fn.iteritems(), extra_bin_fn.iteritems()))
        if extra_dynamic_fn:
            self.dynamic_fn = dict(chain(self.dynamic_fn.iteritems(), extra_dynamic_fn.iteritems()))
    
    
    # INTERNAL FUNCTIONS (user functions further below)
    
    
    def _pushFirst(self, strg, loc, toks ):
        self.exprStack.append( toks[0] )
    def _pushUMinus(self, strg, loc, toks ):
        if toks and toks[0]=='-': 
            self.exprStack.append( 'unary -' )
            #~ exprStack.append( '-1' )
            #~ exprStack.append( '*' )
    
    def _BNF(self):
        """
        expop   :: '^'
        multop  :: '*' | '/'
        addop   :: '+' | '-'
        integer :: ['+' | '-'] '0'..'9'+
        atom    :: PI | E | real | fn '(' expr ')' | '(' expr ')'
        factor  :: atom [ expop factor ]*
        term    :: factor [ multop factor ]*
        expr    :: term [ addop term ]*
        """
        if not self.bnf:
            point = Literal( "." )
            e     = CaselessLiteral( "E" )
            fnumber = Combine( Word( "+-"+nums, nums ) + 
                               Optional( point + Optional( Word( nums ) ) ) +
                               Optional( e + Word( "+-"+nums, nums ) ) )
            ident = Word(alphas, alphas+nums+"_$")
         
            plus  = Literal( "+" )
            minus = Literal( "-" )
            mult  = Literal( "*" )
            div   = Literal( "/" )
            lpar  = Literal( "(" ).suppress()
            rpar  = Literal( ")" ).suppress()
#             comma = Literal( "," ).suppress()
            comma = Literal( "," )
            addop  = plus | minus
            multop = mult | div
            expop = Literal( "^" )
            pi    = CaselessLiteral( "PI" )
            var_list = [Literal(i) for i in self.var_names]
            
            expr = Forward()
            arg_func = Forward()
            or_vars = MatchFirst(var_list)
#             atom = (Optional("-") + ( pi | e | fnumber | ident + lpar + delimitedList(Group(expr)) + rpar | or_vars ).setParseAction( self._pushFirst ) | ( lpar + delimitedList(Group(expr)).suppress() + rpar ) ).setParseAction(self._pushUMinus) 
            atom = (Optional("-") + ( pi | e | fnumber | ident + lpar + arg_func + rpar | or_vars ).setParseAction( self._pushFirst ) | \
                     ( lpar + arg_func.suppress() + rpar ) ).setParseAction(self._pushUMinus) 
            
#             expr + ZeroOrMore( "," + expr )
            # by defining exponentiation as "atom [ ^ factor ]..." instead of "atom [ ^ atom ]...", we get right-to-left exponents, instead of left-to-righ
            # that is, 2^3^2 = 2^(3^2), not (2^3)^2.
            factor = Forward()
            factor << atom + ZeroOrMore( ( expop + factor ).setParseAction( self._pushFirst ) )
            
            term = factor + ZeroOrMore( ( multop + factor ).setParseAction( self._pushFirst ) )
            expr << term + ZeroOrMore( ( addop + term ).setParseAction( self._pushFirst ) )
            arg_func << expr + ZeroOrMore( (comma + expr).setParseAction( self._pushFirst))
            self.bnf = expr
        return self.bnf
    
    # Evaluate the string s. var_list is a dictionary containing the variable values
    # corresponding to the string segments in s. Ex:
    # var_list = {'x':x,'y':y}  'x+y' -> x+y
    # var_list = {'x1':x,'x2':y,'x3':z} 'x1+x2+x3' -> x+y+z
    def _evaluateStack(self, s, var_list={}):
        # Evaluate.
        op = s.pop()
        if op == ',':
            op = s.pop()
        if op == 'unary -':
            return -self._evaluateStack( s,var_list )
        if op in "+-*/^":
            op2 = self._evaluateStack( s,var_list )
            op1 = self._evaluateStack( s,var_list )
            return self.opn[op]( op1, op2 )
        elif op == "PI":
            return math.pi # 3.1415926535
        elif op == "E":
            return math.e  # 2.718281828
        elif op in self.fn:
            return self.fn[op]( self._evaluateStack( s,var_list ) )
        elif op in self.bin_fn:
            op2 = self._evaluateStack( s,var_list )
            if s is not []:
                op1 = self._evaluateStack( s,var_list )
                return self.bin_fn[op]( op1, op2 )
            else:
                return self.bin_fn[op](op2)
        elif op in self.dynamic_fn:
            ops = []
            # This if clause will parse through a comma separated list of arguments
            # until there is no more.
            if s[-1] is ',':
                s.pop()
                while True:
                    ops.append(self._evaluateStack(s, var_list))
                    if s[-1] != ',':
                        break
                    else:
                        s.pop()
            ops.append(self._evaluateStack(s, var_list))
            return self.dynamic_fn[op](*reversed(ops))
        elif op in var_list:
            return var_list[op]
        elif op[0].isalpha():
            return 0
        elif not self.mpfing:
            return float(op)
        else:
            return mpm.mpf(op)
    
    
    # USER FUNCTIONS
        
    def get_stack(self):
        return deepcopy(self.exprStack)
        
    def parse(self,string):
        self.exprStack = []
        self.results = self._BNF().parseString( string )
        
    # Just a wrapper for _evaluateStack (see above) which stops 
    # the input's copy of the stack from being emptied upon evaluation.
    def value_of_stack(self,stack,var_list={}):
        local_stack = deepcopy(stack)
        return self._evaluateStack(local_stack, var_list)
        
if __name__ == "__main__":
    
    # binary functions
    bin_fn = {'sum_this':operator.add}
    
    # variable values
    x = 0.4
    y = 0.5
    var_list = {'x':x,'y':y}
    
    # Instantiate parser, with the added binary function
    parser = math_parser(extra_bin_fn = bin_fn)

    def test( s, expVal ):
        parser.parse(s)
        results = parser.results
        exprStack = parser.exprStack
        val = parser.value_of_stack(exprStack,var_list)
 
        if val == expVal:
            print s, "=", val, results, "=>", exprStack
        else:
            print s+"!!!", val, "!=", expVal, results, "=>", exprStack
            
        
    test("sum_this(x,x+y)",x+x+y)
    test( "9", 9 )
    test( "-9", -9 )
    test( "--9", 9 )
    test( "-E", -math.e )
    test( "9 + 3 + 6", 9 + 3 + 6 )
    test( "9 + 3 / 11", 9 + 3.0 / 11 )
    test( "(9 + 3)", (9 + 3) )
    test( "(9+3) / 11", (9+3.0) / 11 )
    test( "9 - 12 - 6", 9 - 12 - 6 )
    test( "9 - (12 - 6)", 9 - (12 - 6) )
    test( "2*3.14159", 2*3.14159 )
    test( "3.1415926535*3.1415926535 / 10", 3.1415926535*3.1415926535 / 10 )
    test( "PI * PI / 10", math.pi * math.pi / 10 )
    test( "PI*PI/10", math.pi*math.pi/10 )
    test( "PI^2", math.pi**2 )
    test( "round(PI^2)", round(math.pi**2) )
    test( "6.02E23 * 8.048", 6.02E23 * 8.048 )
    test( "e / 3", math.e / 3 )
    test( "sin(PI/2)", math.sin(math.pi/2) )
    test( "trunc(E)", int(math.e) )
    test( "trunc(-E)", int(-math.e) )
    test( "round(E)", round(math.e) )
    test( "round(-E)", round(-math.e) )
    test( "E^PI", math.e**math.pi )
    test( "2^3^2", 2**3**2 )
    test( "2^3+2", 2**3+2 )
    test( "2^9", 2**9 )
    test( "sgn(-2)", -1 )
    test( "sgn(0)", 0 )
    test( "sgn(0.1)", 1 )


"""
Test output:
>pythonw -u fourFn.py
sum_this(x,x+y) = 1.3 ['sum_this', ['x'], ['x', '+', 'y']] => ['x', 'x', 'y', '+', 'sum_this']
9 = 9.0 ['9'] => ['9']
9 + 3 + 6 = 18.0 ['9', '+', '3', '+', '6'] => ['9', '3', '+', '6', '+']
9 + 3 / 11 = 9.27272727273 ['9', '+', '3', '/', '11'] => ['9', '3', '11', '/', '+']
(9 + 3) = 12.0 [] => ['9', '3', '+']
(9+3) / 11 = 1.09090909091 ['/', '11'] => ['9', '3', '+', '11', '/']
9 - 12 - 6 = -9.0 ['9', '-', '12', '-', '6'] => ['9', '12', '-', '6', '-']
9 - (12 - 6) = 3.0 ['9', '-'] => ['9', '12', '6', '-', '-']
2*3.14159 = 6.28318 ['2', '*', '3.14159'] => ['2', '3.14159', '*']
3.1415926535*3.1415926535 / 10 = 0.986960440053 ['3.1415926535', '*', '3.1415926535', '/', '10'] => ['3.1415926535', '3.1415926535', '*', '10', '/']
PI * PI / 10 = 0.986960440109 ['PI', '*', 'PI', '/', '10'] => ['PI', 'PI', '*', '10', '/']
PI*PI/10 = 0.986960440109 ['PI', '*', 'PI', '/', '10'] => ['PI', 'PI', '*', '10', '/']
PI^2 = 9.86960440109 ['PI', '^', '2'] => ['PI', '2', '^']
6.02E23 * 8.048 = 4.844896e+024 ['6.02E23', '*', '8.048'] => ['6.02E23', '8.048', '*']
e / 3 = 0.90609394282 ['E', '/', '3'] => ['E', '3', '/']
sin(PI/2) = 1.0 ['sin', 'PI', '/', '2'] => ['PI', '2', '/', 'sin']
trunc(E) = 2 ['trunc', 'E'] => ['E', 'trunc']
E^PI = 23.1406926328 ['E', '^', 'PI'] => ['E', 'PI', '^']
2^3^2 = 512.0 ['2', '^', '3', '^', '2'] => ['2', '3', '2', '^', '^']
2^3+2 = 10.0 ['2', '^', '3', '+', '2'] => ['2', '3', '^', '2', '+']
2^9 = 512.0 ['2', '^', '9'] => ['2', '9', '^']
sgn(-2) = -1 ['sgn', '-2'] => ['-2', 'sgn']
sgn(0) = 0 ['sgn', '0'] => ['0', 'sgn']
sgn(0.1) = 1 ['sgn', '0.1'] => ['0.1', 'sgn']
>Exit code: 0
"""