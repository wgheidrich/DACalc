
import sys
import io
from contextlib import redirect_stdout

from pydoc import pager

import ply.lex as lex
import ply.yacc as yacc

import dacalc.concretenumber as cn
from dacalc.concretenumber import ConcreteNumber as CN
from dacalc.unitparser import parse as uparse




#######
# data structures
#######

univar_func = {
    "sqrt": cn.sqrt,
    "sin": cn.sin,
    "cos": cn.cos,
    "tan": cn.tan,
    "asin": cn.asin,
    "acos": cn.acos,
    "atan": cn.atan,
    "abs": cn.fabs,
    "ln": cn.ln,
    "log2": cn.log2,
    "log10": cn.log10,
}

bivar_func = {
    "atan2": cn.atan2,
    "log": cn.log,
    "pow": cn.pow
}

# start with empty variables
variables = {}

# user defined unit names (only used for printing help)
user_units = []

# solution count for the last analysis command 
soln_ctr = 0


#######
# lexer
#######


# tokens
tokens = (
    'NEWLINE',
    'HELP',
    'ASSIGN',
    'USE',
    'DEF',
    'LOAD',
    'ANALYZE',
    'LPAR',
    'RPAR',
    'LBRAC',
    'RBRAC',
    'COMMA',
    'ADD',
    'SUB',
    'MUL',
    'DIV',
    'EXP',
    'EXPS',
    'ROOT',
    'RT',
    'INT',
    'NUMBER',
    'UNIVARFUNC',
    'BIVARFUNC',
    'IDENTIFIER',
    'UNITS',
    'STRING'
)


# single character tokens
t_HELP =   r'\?'
t_ASSIGN = r'='
t_LPAR =   r'\('
t_RPAR =   r'\)'
t_LBRAC =  r'\{'
t_RBRAC =  r'\}'
t_COMMA =  r','
t_ADD =    r'\+'
t_SUB =    r'-'
t_MUL =    r'\*'
t_DIV =    r'/'
t_RT =     r'\u221a'


# variable/constant names and builtin functions
def t_IDENTIFIER(t):
    r'[a-zA-Z_][a-zA-Z0-9_]*'
    if t.value in univar_func:
        t.type = "UNIVARFUNC"
        t.value = univar_func[t.value]
    elif t.value in bivar_func:
        t.type = "BIVARFUNC"
        t.value = bivar_func[t.value]
    elif t.value == "use":
        t.type = 'USE'
    elif t.value == "def":
        t.type = "DEF"
    elif t.value == "import":
        t.type = "LOAD"
    elif t.value == "analyze":
        t.type = "ANALYZE"
    return t

# signed integer exponent (i.e. leading '^')
def t_EXP(t):
    r'\^[+-]?\d+'
    t.value = int(t.value[1:])
    return t

# signed integer exponent in unicode superscript
def t_EXPS(t):
    r'(\u207a|\u207b)?((\u2070|\xb9|\xb2|\xb3|\u2074|\u2075|\u2076|\u2077|\u2078|\u2079)+)'
    t.value = cn.sup_int(t.value)
    return t

# signed integer root (i.e. leading '%')
def t_ROOT(t):
    r'\%-?\d+'
    t.value = int(t.value[1:])
    return t

def t_NUMBER(t):
    r'(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?'
    t.value = float(t.value)
    if t.value.is_integer():
        t.type = "INT"
        t.value = int(t.value)
    return t

# a unit string, bounded by brackets
def t_UNITS(t):
    r'\[[^\]]*\]'
    t.value = t.value[1:-1]
    return t

def t_STRING(t):
    r'"[^"]*"'
    t.value = t.value[1:-1]
    return t

def t_comment(t):
    r'\#[^\n]*'
    pass

# Define a rule so we can track line numbers
def t_NEWLINE(t):
    r'[;\n]+'
    t.lexer.lineno += len(t.value)
    return t
    
# A string containing ignored characters (spaces and tabs)
t_ignore  = ' \t'

# Error handling rule
def t_error(t):
    print("Illegal character:\n  \"... "+t.value[0:10]+"...\"\n       ^",
          file = sys.stderr)
    t.lexer.skip(1)

    

# Build the lexer
da_lexer = lex.lex(optimize=1,lextab='dacalc.dalextab')


def dbg_lexer(s):
    da_lexer.input(s)
    for tok in iter(lex.token, None):
        print(repr(tok.type), repr(tok.value))




##########
# parser
##########


precedence = (
    ('left','ADD','SUB'),
    ('left','MUL','DIV'),
    ('right','UMINUS'),
    ('left','EXP','ROOT'),
    ('right','UNITS')
)

def p_lineseq(t):
    '''
    lines : lines NEWLINE line
          | line
    '''
    t[0] = t[1] if type(t[1]) == str else ""
    if len(t) == 4 and type(t[3]) == str:
        t[0] + t[3]


def p_oneline(t):
    '''
    line : typestatement
         | usestatement
         | unitdef
         | helpline
         | loadstatement
         | analyzestatement
         | empty
    '''
    t[0] = t[1]

def p_empty(t):
    'empty :'
    pass

def p_helpline(t):
    '''
    helpline : HELP
             | HELP IDENTIFIER
    '''
    if len(t) == 2:
        t[0] =  '''
  DA calculator -- a calculator for dimensional analysis
  
  Dimensional analysis means calculating with concrete values, i.e. values
  that have units attached, and thus have physical meanings. For details, try
        
        ? <topic>
  
  where <topic> is one of (any unique prefix will do):

        basics:       basic calculations
        commands:     all commands
        operators:    built-in arithmetic operators
        func:         built-in functions (short list)
        functions:    built-in functions (details)
        constants:    built-in constants
        units:        currently defined units
        variables:    currently defined variables
        '''
    elif t[2] == 'basics'[:len(t[2])]:
        t[0] = '''
  At the core of dacalc is a scientific calculator that keeps tracks of
  units and dimensions, and performs automatic unit conversion. Here is
  a simple example for computing a force and a torque in different unit
  systems:

        Welcome to the dimensional analysis calculator!
        Try "?" for help...
        DA > m = 5[kg]
        m = 5 [kg]
        DA > a = 10[m/s^2]
        a = 10 [m s⁻²]
        DA > F = m*a
        F = 50 [N]
        DA > l = 10[in]
        l = 254 [mm]
        DA > T = F*l
        T = 12.7 [J]
        DA > [lbf in] T
        112.404 [lbf in]
        DA > 
        
   Refer to the 'commands' help page for more details for more options
   on custom units and result display.

   Of course dacalc also supports all the usual math library functions
   (see help page for 'functions'). 

        DA > r = 5[cm]
        r = 50 [mm]
        DA > alpha = 60[deg]
        alpha = 1.0472 [rad]
        DA > l = r*sin(alpha)
        l = 43.3013 [mm]
        DA > 
        DA > beta = atan2(3[mm],5[cm])
        beta = 59.9282 [mrad]
        DA > [deg] beta
        3.43363 [deg]        
        '''
    elif t[2] == 'commands'[:len(t[2])]:
        t[0] = '''
  Apart from arithmetic expressions, variable assignments, and unit conversions
  (see help page 'basics') several special commands are supported. They are:

        use :      select output units and other choices for output formatting
        
        def :      define a new unit

        import :   load a dacalc script
        
        analyze :  perform dimensional analysis on a set of variables

  The specifics of the individual commands are as follows:

  * use:
        The use command controls the output of values. It has three different
        variants: 
        
        use <system> :  choose a specific unit system for outputing all
                        quantities. Supported systems are SI and US, as
                        well as SI_base and US_base. The latter two only
                        utilize base units, so that for example a force
                        would be displayed as [kg m s^-2] instead of [N]

        use [unit] :    choose a default unit for quantities of a specific
                        dimension. Examples:

                        DA > use [lb]
                        Using [lb]
                        DA > 5[kg]
                        11.0231 [lb]

                        DA > use [N m]
                        Using [N m]
                        DA > 2[kg] * _g * 1[m]
                        19.6133 [N m]

        use <digits> :  where digits is a non-negative integer. This selects
                        the number of displayed digits in numerical values.

  * def:
        The def command introduces new units.

        def unit expr : where "unit" is the name of the new unit, and "expr"
                        is an artithmetic expression that calculates the value
                        of the new unit expressed in known qunatities. For
                        example, to define the light year as a unit:

                        DA > def ly 365.25[d] * _c
                        New unit ly = 9.46073e+15 [m]
                        DA > .2[ly]
                        1.89215e+15 [m]
  
  * import:
        Imports a dacalc script. All variables, units and settings made in
        the script will become available. Syntax:

        import file :   where file is  a file name enclosed in single (') or
                        double (") quotes

  * analyze:
        Perform dimensional analysis of multiple variables to try to match
        a target unit.

        analyze [unit] varlist :
                        here, "unit" is a target unit, which could be [] for
                        dimensionless quantities

                        varlist is a comma-separated list of variables and
                        constants, enclosed in brackets {}.

                        analyze tries to match the target units by combining
                        different powers of the provided variables. The
                        possible solutions are stored in variables named
                        "solutionN", where N is a number.

                        A simple example for finding the equation for
                        potential energy is given below; please check the
                        example notebooks for more detailed examples.

                        DA > m = 3[kg]
                        m = 3 [kg]
                        DA > h =  1[m]
                        h = 1 [m]
                        DA > analyze [J] {m,h,_g}
                        solution0 = m * h * _g = 29.42 [J]

        '''        
    elif t[2] == 'operators'[:len(t[2])]:
        t[0] = '''
  The following operators are supported (a and b are concrete values):
        
        a + b:        addition
        a - b:        subtraction
        a * b:        multiplication
        a / b:        division
        - a:          negation
        a ^ i:        i-th power of a (i is an int, see below)
        a % i:        i-th root of a
  
  In the last two operators, <i> refers to an integer constant, i.e. this
  cannot be a calculated value but has to be specified directly as a^-2,
  b%2 etc. Alternatively, it is  possible to directly enter unicode
  superscript and root characters, e.g.

        a\u207b\xb2           is identical to a ^ -2, and
        \xb3\u221aa           is identical to a % 3

  
  The order of precedence of these operations matches standard math
  conventions, i.e. the input

        a + -b / c * d^2

  is equivalent to

        a + (((-b) / c) * (d^2))

        '''
    elif t[2] == 'func'[:len(t[2])]:
        t[0] = "  The following functions are supported:\n"
        for n,f in univar_func.items():
            t[0] += '\t' + n + "(a)" + '\n'
        for n,f in bivar_func.items():
            t[0] += '\t' + n + "(a,b)" + '\n'
        t[0] += '\n'
    elif t[2] == 'functions'[:len(t[2])]:
        t[0] = "  The following functions are supported:\n"
        for n,f in univar_func.items():
            t[0] += "\n  " + n + "(.):\n" + f.__doc__ + '\n' + '-'*75 + '\n'
        for n,f in bivar_func.items():
            t[0] += "\n  " + n + "(.,.):\n" + f.__doc__ + '\n' + '-'*75 + '\n'
        t[0] += '\n'
    elif t[2] == 'constants'[:len(t[2])]:
        t[0] = '''
  Constants are predefined values that cannot be changed. To distinguish
  them from variables, they all start with an underscore '_', i.e. the
  constant 'pi' can be accessed as '_pi'.

  List of constants:
  
 '''
        ctr = 0
        for c in cn.const:
            t[0] += "   "+"{:<5} : {:<15}".format(c[0],c[2])
            ctr+= 1
            if ctr%3 == 0:
                t[0] += '\n ' # indentation fudging
        t[0] += '\n'
    elif t[2] == 'units'[:len(t[2])]:
        t[0] = "  SI units:\n"
        ctr = 0
        for u in cn.si_u:
            t[0] += "    {:<16}".format(u[2]+" ["+u[0]+"],")
            ctr+= 1
            if ctr%4 == 0:
                t[0] += '\n'
        t[0] += "\n\n  SI prefixes (can be used with any SI unit):\n"
        ctr = 0
        for p in CN.SI_PREFIXES.items():
            t[0] += "    %s: %1.0e" % p
            ctr+= 1
            if ctr%5 == 0:
                t[0] += '\n'
        t[0] += "\n  US units (not that some of these differ from UK units):\n"
        ctr = 0
        for u in cn.us_u:
            t[0] += "    {:<20}".format(u[2]+" ["+u[0]+"],")
            ctr+= 1
            if ctr%3 == 0:
                t[0] += '\n'
        t[0] += "\n\n  Other units (often outdated):\n"
        ctr = 0
        for u in cn.misc_u:
            t[0] += "    {:<20}".format(u[2] + " [" + u[0] + "],")
            ctr+= 1
            if ctr%3 == 0:
                t[0] += '\n'
        t[0] += "\n\n  User defined units:\n"
        ctr = 0
        for u in user_units:
            t[0] += "    {:<10}".format(u)
            ctr+= 1
            if ctr%5 == 0:
                t[0] += '\n'
        t[0] += '\n'
    elif t[2] == "vars"[:len(t[2])]:
        t[0] = '''
  Variables are user-specified quantities that can be defined and changed.
  The names of variables must start with a letter (unlike constants, which
  start with an underscore '_').
  
  Currently defined variables and their values:

'''
        for n,v in variables.items():
            t[0] += '\t' + ' ' + n + ' :\t' + str(v) + '\n'
        t[0] += '\n'
    else:
        print("Unknown help topic",file = sys.stderr)
        t[0] = None


def p_use_units(t):
    'usestatement : USE UNITS'
    newu = uparse(t[2])
    CN.use(newu,t[2])
    print("Using ["+t[2]+"]")

    
def p_use_system(t):
    'usestatement : USE IDENTIFIER'
    CN.use_system(t[2])
    print("Using",t[2],"system")
   
    
def p_use_precision(t):
    'usestatement : USE INT'
    CN.set_precision(t[2])
    print("Using precision of",t[2],"digits")

    
def p_def(t):
    'unitdef : DEF IDENTIFIER expr'
    if t[2] in user_units:
        # just overwriting a previous user-defined unit
        cn.add_unit_to_dict(t[2],t[3])
        print("Redefined unit",t[2],'=',t[3])
    elif t[2] in CN.units.keys():
        # if the unit exists but is not user-defined then this
        # is an attempt to redefine a built-in unit
        print("Cannot redefine built-in unit",file = sys.stderr)
    else:
        # new unit
        cn.add_unit_to_dict(t[2],t[3])
        user_units.append(t[2])
        print("New unit",t[2],'=',t[3])

        
def p_load(t):
    'loadstatement : LOAD STRING'
    f = open(t[2])
    script = f.read()
    f.close()
    new_lexer = da_lexer.clone()
    with io.StringIO() as buf, redirect_stdout(buf):
        # only error messages will be reported
        da_parser.parse(script,lexer=new_lexer)


def var_lookup(name):
    'look up the value of a variable or constant'
    result = 0
    if name[0] == '_':
        try:
            result = CN.const[name[1:]]
        except LookupError:
            print("Unknown constant '%s'" % name,file = sys.stderr)
    else:
        try:
            result = variables[name]
        except LookupError:
            print("Unknown variable '%s'" % name,file = sys.stderr)
    return result


def p_analyze(t):
    'analyzestatement : ANALYZE UNITS LBRAC varlist RBRAC'
    global soln_ctr

    target = uparse(t[2])
    val_list = [var_lookup(var) for var in t[4]]
    solns = CN.analyze(val_list,target)
    
    # delete old solution variables
    for i in range(0,soln_ctr):
        variables.pop("solution"+str(i),None)
    soln_ctr = 0
    for s in solns:
        sol_str = ""
        sol = CN()
        for name,val in zip(t[4],s):
            if val[1] != 0:
                sol *= val[0]**val[1]
                sol_str += ' * '+name
                if val[1]!= 1:
                    sol_str+= '^'+str(val[1])
        var_str = "solution"+str(soln_ctr)
        sol_str = CN.make_super(sol_str)
        print(var_str+' ='+sol_str[2:]+' = '+sol.__str__(t[2]))
        variables[var_str] = sol
        soln_ctr+= 1
    if len(solns) == 0:
        print("No solution found for analysis problem",file = sys.stderr)
        

def p_varlist(t):
    '''
    varlist : varlist COMMA IDENTIFIER
            | IDENTIFIER
    '''
    if len(t) == 2:
        t[0] = [t[1]]
    else:
        t[0] = t[1]
        t[0].append(t[3])
    
    
def p_tstatement_typecast(t):
    'typestatement : UNITS statement'
    t[0] = t[2]
    print(t[2].__str__(t[1]))

    
def p_tstatement_plain(t):
    'typestatement : statement'
    t[0] = t[1]
    print(t[1])

    
def p_statement_assign(t):
    'statement : IDENTIFIER ASSIGN expr'
    t[0] = t[3]
    if t[1][0] == '_':
        print("Cannot define new constants...",file = sys.stderr)
    else:
        variables[t[1]] = t[3] # store new / updated variable
        print(t[1],end=' = ')

        
def p_statement_expr(t):
    'statement : expr'
    t[0] = t[1]

    
def p_expr_binop(t):
    '''
    expr : expr ADD expr
    expr : expr SUB expr
    expr : expr MUL expr
    expr : expr DIV expr
    '''
    if t[2] == '+'  : t[0] = t[1] + t[3]
    elif t[2] == '-': t[0] = t[1] - t[3]
    elif t[2] == '*': t[0] = t[1] * t[3]
    elif t[2] == '/': t[0] = t[1] / t[3]

def p_expr_unitmul(t):
    'expr : expr UNITS'
    t[0] = t[1] * uparse(t[2])

def p_expr_exp(t):
    '''
    expr : expr EXP %prec EXP
         | expr EXPS %prec EXP
    '''
    t[0] = t[1]**t[2]

def p_expr_root(t):
    'expr : expr ROOT %prec ROOT'
    t[0] = cn.root(t[1],t[2])

# alternative way to do roots by prefix superscript and root symbol
def p_expr_root_sup(t):
    '''
    expr : EXPS RT expr %prec ROOT
    '''
    if len(t) == 4:
        t[0] = cn.root(t[3],t[1])
    else:
        t[0] = cn.root(t[2],2)
    
def p_expr_uminus(t):
    'expr : SUB expr %prec UMINUS'
    t[0] = -t[2]

def p_expr_univar(t):
    'expr : UNIVARFUNC LPAR expr RPAR'
    t[0] = t[1](t[3])

def p_expr_bivar(t):
    'expr : BIVARFUNC LPAR expr COMMA expr RPAR'
    t[0] = t[1](t[3],t[5])
    
def p_expr_group(t):
    'expr : LPAR expr RPAR'
    t[0] = t[2]

def p_expr_const(t):
    '''
    expr : INT
    expr : NUMBER
    '''
    t[0] = t[1]*CN.u("")


def p_expr_var(t):
    'expr : IDENTIFIER'
    t[0] = var_lookup(t[1])

    
def p_error(t):
    print("Syntax error at '%s'" % t.value,file = sys.stderr)


da_parser = yacc.yacc(optimize=1,tabmodule='dacalc.daparsetab')

def parse(s):
    return da_parser.parse(s,lexer=da_lexer)

def main():
    import readline

    print('Welcome to the dimensional analysis calculator!')
    print('Try "?" for help...')
    while True:
        try:
            s = input('DA > ')
        except EOFError:
            break
        try:
            help_msg = parse(s)
            if type(help_msg) == str and help_msg != "":
                pager(help_msg)
        except TypeError as te:
            print(te,file = sys.stderr)
        except Exception as exc:
            print(exc,file = sys.stderr)

    print("\nGoodbye...")


if __name__ == '__main__':
    main()

