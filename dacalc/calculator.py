
import sys
import io
from contextlib import redirect_stdout

from pydoc import pager

import ply.lex as lex
import ply.yacc as yacc

import dacalc.concretenumber as cn
from dacalc.concretenumber import ConcreteNumber as CN
from dacalc.unitparser import parse as uparse
import dacalc.conversions as conv
import dacalc.expression as ex

import matplotlib.pyplot as plt
import numpy


#######
# data structures
#######

# built-in functions
functions = {
    "sqrt": cn.sqrt,
    "sin": cn.sin,
    "cos": cn.cos,
    "tan": cn.tan,
    "asin": cn.asin,
    "acos": cn.acos,
    "atan": cn.atan,
    "atan2": cn.atan2,
    "abs": cn.fabs,
    "pow": cn.pow,
    "exp": cn.exp,
    "log": cn.log,
    "log2": cn.log2,
    "log10": cn.log10,
    "Re": cn.Re,
    "Im": cn.Im,
    "phase": cn.phase,
    "conj": cn.conj,
    "rect": cn.rect,
    "SI_to_Gauss": conv.SI_to_Gauss,
    "SI_to_ESU": conv.SI_to_ESU,
    "SI_to_EMU": conv.SI_to_EMU,
    "Gauss_to_SI": conv.Gauss_to_SI,
    "ESU_to_SI": conv.ESU_to_SI,
    "EMU_to_SI": conv.EMU_to_SI,
}

# user defined functions
# each entry is a tuple of a parameter list and an Expression
user_functions = {}

# start with just the constants
variables = {'_'+n: v for n, v in CN.const.items()}

# user defined unit names (only used for printing help)
user_units = []

# solution count for the last analysis command
soln_ctr = 0

# whether we are in complex mode
complex_mode = False


#######
# lexer
#######


# tokens
tokens = (
    'NEWLINE',
    'HELP',
    'ASSIGN',
    'USE',
    'OUTPUT',
    'DEF',
    'LOAD',
    'IMAGE',
    'DIM',
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
    'INT',
    'NUMBER',
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
t_DIV =    r'/'
t_EXP =    r'\^'


# this can be either MUL or EXP in Python syntax depending on the number of *s
def t_MUL(t):
    r'\*\*?'
    if len(t.value) > 1:
        t.type = "EXP"
        t.value = '^'
    return t


# root is either % or unicode 221a
def t_ROOT(t):
    r'[%\u221a]'
    t.value = '\u221a'  # used for output
    return t


# variable/constant names and builtin functions
def t_IDENTIFIER(t):
    r'[a-zA-Z_][a-zA-Z0-9_]*'
    if t.value == "use":
        t.type = 'USE'
    elif t.value == "output":
        t.type = 'OUTPUT'
    elif t.value == "def":
        t.type = "DEF"
    elif t.value == "import":
        t.type = "LOAD"
    elif t.value == "image":
        t.type = "IMAGE"
    elif t.value == "dim":
        t.type = "DIM"
    elif t.value == "analyze":
        t.type = "ANALYZE"
    return t


# signed integer exponent in unicode superscript
def t_EXPS(t):
    (r'(\u207a|\u207b)?'
     r'((\u2070|\xb9|\xb2|\xb3|\u2074|\u2075|\u2076|\u2077|\u2078|\u2079)+)')
    t.value = cn.sup_int(t.value)
    return t


def t_NUMBER(t):
    r'((\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?)j?'
    if t.value[-1] == 'j':
        # this is a complex number
        if complex_mode:
            t.value = complex(t.value)
        else:
            raise(TypeError("Please enable complex number mode."))
    else:
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
t_ignore = ' \t'


# Error handling rule
def t_error(t):
    print("Illegal character:\n  \"... " + t.value[0:10] + "...\"\n       ^",
          file=sys.stderr)
    t.lexer.skip(1)


# Build the lexer
da_lexer = lex.lex(optimize=1, lextab='dacalc.dalextab')


def dbg_lexer(s):
    da_lexer.input(s)
    for tok in iter(lex.token, None):
        print(repr(tok.type), repr(tok.value))


##########
# parser
##########


precedence = (
    ('left', 'ADD', 'SUB'),
    ('left', 'MUL', 'DIV'),
    ('right', 'UMINUS'),
    ('left', 'EXP', 'EXPS', 'ROOT'),
    ('right', 'UROOT'),
    ('right', 'UNITS')
)


def p_lineseq(t):
    '''
    lines : lines NEWLINE line
          | line
    '''
    # concatenate all the results for the lines into a list
    if len(t) == 4:
        t[0] = t[1]
        t[0].append(t[3])
    else:
        t[0] = [t[1]]


def p_oneline(t):
    '''
    line : typestatement
         | usestatement
         | outputstatement
         | unitdef
         | helpline
         | loadstatement
         | imgstatement
         | stringstatement
         | dimstatement
         | analyzestatement
         | newfunc
         | empty
    '''
    if isinstance(t[1],CN):
        variables['_'] = t[1]
    if t[1] is not None:
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

        use :      select output units

        output :   formatting options for output

        def :      define a new unit

        dim :      show the dimensions of an expression

        import :   load a dacalc script

        analyze :  perform dimensional analysis on a set of variables

  The specifics of the individual commands are as follows:

  * use:
        The use command controls the output of values. It has different
        variants:

        use <system> :  choose a specific unit system for outputing all
                        quantities. Supported systems are SI and US, as
                        well as SI_base and US_base. The latter two only
                        utilize base units, so that for example a force
                        would be displayed as [kg m s^-2] instead of [N]

                        DA calc also supports obsolete systmes based on
                        centimetre-gram-second including just base
                        units (CGS_base), Gauss units (CGS or CGS_Gauss),
                        CGS plus electrostatic units (CGS_ESU) and CGS plus
                        electromagentic units (CGS_EMU). Selection of any of
                        these also implies half dimensions (see "use halfdim"
                        below).

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

        use complex :   switch from real numbers to complex numbers. This
                        cannot be undone, except by restarting the calculator.

        use halfdim :   enable half dimensions for length and mass as needed
                        by CGS unit systems. Once enabled, it cannot be
                        disabled, except by restarting the calculator.

  * output:
        Various options for output formatting

        output plain|unicode|html : text format for the output.
                        "plain" uses plain text, with exponents shown as a^-2.
                        "unicode" uses unicode superscript to format exponents
                        "html" uses HTML superscripts
                        while "plain" and "unicode" produces output that can
                        be copy/pasted and used as input, this is not the
                        case for "html"

        output <num> :  where num is a non-negative integer. This selects
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

  * dim:
        Display the dimensions of a quantity

        dim expr:       where expr is any expression

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
        a ^ e:        e-th power of a (see below for valid values of e)
        i \u221a a:        i-th root of a (i is an integer)
        \u221a a:          square root of a

  The power and root operators also have alternative syntax variants:

        a\u207b\xb2           (unicode superscript) is identical to a^-2
        a**-2         is also identical to a^-2 (Python syntax for powers)
        3%a           is identical to 3\u221aa

  Note that the unicode superscript syntax is really only available for
  integer constants since there is no unicode symbol for superscript
  decimal points. Also note that all power and root operations are only
  valid if they produce integer dimensions. This means that floating point
  exponents are only available for dimensionless quantities. So for example,

        2^.4          is valid, as is
        3\u221a1[L]        since the cube root of a volume is a length, but
        2[km]^.4      is not valid since length^.4 is physically ill-defined.

  The order of precedence of these operations matches standard math
  conventions, i.e. the input

        a + -b / c * d^2

  is equivalent to

        a + (((-b) / c) * (d^2))

  The unit operator [.] has higher priority than any other operator.

  Note that the power operator has a higher priority than the unary minus,
  i.e. -2^2 is equivalent to -(2^2) instead of (-2)^2. This choice was made
  in order to be compatible with Python expressions.

  Also note the interaction between the power and root operators: \u221a-2^2
  gets evaluated as \u221a(-(2^2)), which is the complex number (0+2j).
        '''
    elif t[2] == 'func'[:len(t[2])]:
        t[0] = "  The following functions are supported:\n"
        for n, f in functions.items():
            t[0] += '\t' + n + "()" + '\n'
        t[0] += '\n'
    elif t[2] == 'functions'[:len(t[2])]:
        t[0] = "  The following functions are supported:\n"
        for n, f in functions.items():
            t[0] += "\n  " + n + "():\n" + f.__doc__ + '\n' + '-' * 75 + '\n'
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
            t[0] += "  {:<5} : {:<15}".format(c[0], c[2])
            ctr += 1
            if ctr % 3 == 0:
                t[0] += '\n  '
        t[0] += '''

  Also, a single '_' refers to the value of the last expression, e.g.

    DA > 2[kg] * _g
    19.6133 [N]

    DA > [lbf] _
    4.40925 [lbf]

  Statements that may produce multiple solutions (e.g. the 'analyze' statement)
  define constants named '_0', '_1', and so forth for their results. These
  'constants' are defined only until the next such statement overwrites them.

  '''
    elif t[2] == 'units'[:len(t[2])]:
        t[0] = "  SI units:\n"
        ctr = 0
        for u in cn.si_u:
            t[0] += "    {:<16}".format(u[2] + " [" + u[0] + "],")
            ctr += 1
            if ctr % 4 == 0:
                t[0] += '\n'
        t[0] += "\n\n  SI prefixes (can be used with any SI unit):\n"
        ctr = 0
        for p in CN.SI_PREFIXES.items():
            t[0] += "    %s: %1.0e" % p
            ctr += 1
            if ctr % 5 == 0:
                t[0] += '\n'
        t[0] += "\n  US units (not that some of these differ from UK units):\n"
        ctr = 0
        for u in cn.us_u:
            t[0] += "    {:<20}".format(u[2] + " [" + u[0] + "],")
            ctr += 1
            if ctr % 3 == 0:
                t[0] += '\n'
        t[0] += "\n\n  Other units (often outdated):\n"
        ctr = 0
        for u in cn.misc_u:
            t[0] += "    {:<20}".format(u[2] + " [" + u[0] + "],")
            ctr += 1
            if ctr % 3 == 0:
                t[0] += '\n'
        t[0] += "\n\n  User defined units:\n"
        ctr = 0
        for u in user_units:
            t[0] += "    {:<10}".format(u)
            ctr += 1
            if ctr % 5 == 0:
                t[0] += '\n'
        t[0] += '\n'
    elif t[2] == "vars"[:len(t[2])]:
        t[0] = '''
  Variables are user-specified quantities that can be defined and changed.
  The names of variables must start with a letter (unlike constants, which
  start with an underscore '_').

  Currently defined variables and their values:

'''
        for n, v in variables.items():
            t[0] += '\t' + ' ' + n + ' :\t' + str(v) + '\n'
        t[0] += '\n'
    else:
        print("Unknown help topic", file=sys.stderr)
        t[0] = None


def p_use_units(t):
    'usestatement : USE UNITS'
    newu = uparse(t[2])
    CN.use(newu, t[2])
    print("Using [" + t[2] + "]")


def p_use_system(t):
    'usestatement : USE IDENTIFIER'
    if t[2] == "complex":
        import cmath
        cn.cnmath = cmath
        global complex_mode
        complex_mode = True
        print("Switching to complex number mode")
    elif t[2] == "halfdim":
        CN.half_dimensions = True
        print("Enabling CGS-style half dimensions for mass and length")
    else:
        CN.use_system(t[2])
        print("Using", t[2], "system")


def p_use_precision(t):
    'usestatement : USE INT'
    CN.set_precision(t[2])
    print("Using precision of", t[2], "digits")


def p_output_format(t):
    '''
    outputstatement : OUTPUT IDENTIFIER
                    | OUTPUT INT
    '''
    if type(t[2]) == int:
        CN.set_precision(t[2])
        print("Using precision of", t[2], "digits")
    elif type(t[2]) == str:
        CN.output_mode = t[2]
        print("New output mode:", t[2])
    else:
        print("Syntax error at", t[2], file=sys.stderr)


def p_def(t):
    'unitdef : DEF IDENTIFIER expr'
    # evaluate expressions with variables and check that all variables
    # are defined
    val = t[3].eval(variables) if isinstance(t[3], ex.Expression) else t[3]
    if isinstance(val, ex.Expression):
        raise(NameError("def statement has undefined variable(s):\n\t"
                        + val.get_undefined().keys()))

    if t[2] in user_units:
        # just overwriting a previous user-defined unit
        cn.add_unit_to_dict(t[2], val)
        print("Redefined unit", t[2], '=', val)
    elif t[2] in CN.units.keys():
        # if the unit exists but is not user-defined then this
        # is an attempt to redefine a built-in unit
        print("Cannot redefine built-in unit", file=sys.stderr)
    else:
        # new unit
        cn.add_unit_to_dict(t[2], val)
        user_units.append(t[2])
        print("New unit", t[2], '=', val)


def p_load(t):
    'loadstatement : LOAD STRING'
    with open(t[2]) as f:
        script = f.read()
    new_lexer = da_lexer.clone()
    with io.StringIO() as buf, redirect_stdout(buf):
        # only error messages will be reported
        da_parser.parse(script, lexer=new_lexer)


def p_image(t):
    'imgstatement : IMAGE STRING'
    t[0] = plt.imread(t[2])


def p_string(t):
    'stringstatement : STRING'
    # substitute string literals before output
    print(t[1].encode('raw_unicode_escape').decode('unicode_escape'))


def p_dim(t):
    '''
    dimstatement : DIM expr
                 | DIM UNITS
    '''
    if type(t[2]) == str:
        val = uparse(t[2])
    else:
        val = t[2].eval(variables) if isinstance(t[2], ex.Expression) else t[2]
        if isinstance(val, ex.Expression):
            raise(NameError("dim statement has undefined variable(s):\n\t"
                            + val.get_undefined().keys()))
    print(val.dimensionstr())


def p_analyze(t):
    'analyzestatement : ANALYZE UNITS LBRAC exprlist RBRAC'
    global soln_ctr

    # target dimension for the anlysis
    target = uparse(t[2])

    # check that each expression in the list is actually a variable
    # name, and the variable is fully defined (has no free variables
    # itself)
    var_list = []
    val_list = []
    for var in t[4]:
        if not isinstance(var, ex.VarExpr):
            raise(TypeError("Expected variable name, got " + str(var)))
        var_list.append(var.name)
        val = var.eval(variables)
        if isinstance(val, ex.Expression):
            raise(TypeError("Can only analyze fully defined expressions"))
        val_list.append(val)

    # perform the analysis
    solns = CN.analyze(val_list, target)

    # delete old solution variables
    for i in range(0, soln_ctr):
        variables.pop("solution" + str(i), None)
    soln_ctr = 0
    # define new variables for all solutions we found.
    for s in solns:
        sol_str = ""
        sol = CN()
        for name, val in zip(var_list, s):
            if val[1] != 0:
                sol *= val[0]**val[1]
                sol_str += ' * ' + name
                if val[1] != 1:
                    sol_str += '^' + str(val[1])
        var_str = "_" + str(soln_ctr)
        sol_str = CN.make_super(sol_str)
        print(var_str + ' =' + sol_str[2:] + ' = ' + sol.__str__(t[2]))
        variables[var_str] = sol
        soln_ctr += 1
    if len(solns) == 0:
        print("No solution found for analysis problem", file=sys.stderr)


def p_newfunc(t):
    'newfunc : IDENTIFIER LPAR exprlist RPAR ASSIGN expr'
    if t[1] in functions.keys():
        raise TypeError("Cannot redefine built-in function "+t[1])
    if t[6] is None:
        raise(TypeError("Assignment needs right hand side"))
    
    # check that each expression in the parameter list is actually a
    # variable name
    param_list = []
    for var in t[3]:
        if isinstance(var, ex.VarExpr):
            param_list.append(var.name)
        else:
            raise(TypeError("Expected variable name, got "+var))

    # prevent substitution of parameters here to implement proper
    # argument scoping
    scope_vars = dict(variables)
    for param in param_list:
        scope_vars.pop(param,True)
    t[0] = t[6].eval(scope_vars) \
        if isinstance(t[6], ex.Expression) else t[6]

    # we don't allow free variables
    for v in t[0].get_undefined():
        if v not in param_list:
            raise NameError("Undefined variable: " + v)

    user_functions[t[1]] = (param_list, t[0])
    print("New function " + t[1] + '(' + ', '.join(a for a in param_list)
          + ') = ' + str(t[0]))


def p_tstatement_typecast(t):
    'typestatement : UNITS statement'
    t[0] = t[2]
    if isinstance(t[0], ex.Expression):
        print("Can only unit-convert fully defined quantities",
              file=sys.stderr)
    else:
        print(t[2].__str__(t[1]))


def p_tstatement_plain(t):
    'typestatement : statement'
    t[0] = t[1]
    print(t[1])


def p_statement_assign(t):
    'statement : IDENTIFIER ASSIGN expr'
    if t[3] is None:
        raise(TypeError("Assignment needs right hand side"))

    t[0] = t[3].eval(variables) \
        if isinstance(t[3], ex.Expression) else t[3]

    # we don't allow free variables
    if isinstance(t[0], ex.Expression):
        for v in t[0].get_undefined():
            raise NameError("Undefined variable: " + v)

    if t[1][0] == '_':
        print("Cannot define new constants...", file=sys.stderr)
    else:
        variables[t[1]] = t[0]  # store new / updated variable
        print(t[1], end=' = ')


def p_statement_expr(t):
    'statement : expr'
    t[0] = t[1].eval(variables) \
        if isinstance(t[1], ex.Expression) else t[1]
    # check for free variables
    if isinstance(t[0], ex.Expression):
        for v in t[0].get_undefined():
            raise NameError("Undefined variable: " + v)
    
    


def p_expr_binop(t):
    '''
    expr : expr ADD expr
         | expr SUB expr
         | expr MUL expr
         | expr DIV expr
         | expr EXP expr
         | expr ROOT expr
    '''
    if isinstance(t[1], ex.Expression) or isinstance(t[3], ex.Expression):
        # dealing with partially defined expression -- create expression tree
        t[0] = ex.BinOpExpr(t[2], t[1], t[3])
    else:
        # otherwise, evaluate directly
        if t[2] == '+':
            t[0] = t[1] + t[3]
        elif t[2] == '-':
            t[0] = t[1] - t[3]
        elif t[2] == '*':
            t[0] = t[1] * t[3]
        elif t[2] == '/':
            t[0] = t[1] / t[3]
        elif t[2] == '^':
            t[0] = t[1] ** t[3]
        elif t[2] == '\u221a':
            t[0] = cn.root(t[3], t[1])


# exponent as superscipt
def p_expr_exp(t):
    '''
    expr : expr EXPS %prec EXP
    '''
    if isinstance(t[1], ex.Expression):
        t[0] = ex.BinOpExpr('^', t[1], t[2])
    else:
        t[0] = t[1]**t[2]


# unary root == square root
def p_expr_uroot(t):
    'expr : ROOT expr %prec UROOT'
    if isinstance(t[2], ex.Expression):
        t[0] = ex.UniOpExpr(t[1], t[2])
    else:
        t[0] = cn.root(t[2], 2)


# unary minus
def p_expr_uminus(t):
    'expr : SUB expr %prec UMINUS'
    if isinstance(t[2], ex.Expression):
        t[0] = ex.UniOpExpr(t[1], t[2])
    else:
        t[0] = -t[2]


def p_expr_unitmul(t):
    'expr : expr UNITS'
    arg2 = uparse(t[2])
    if isinstance(t[1], ex.Expression):
        t[0] = ex.BinOpExpr('*', t[1], arg2)
    else:
        t[0] = t[1] * arg2


def p_expr_func(t):
    'expr : IDENTIFIER LPAR exprlist RPAR'
    # first, evaluate all parameters
    param_vals = [p.eval(variables) if isinstance(p, ex.Expression) else p
                  for p in t[3]]

    if t[1] in functions:
        # builtin function
        t[0] = functions[t[1]](*param_vals)
    elif t[1] in user_functions:
        # user defined function
        par_names, body = user_functions[t[1]]
        if len(par_names) != len(param_vals):
            raise(TypeError("Wrong number of arguments for function " + t[1] +
                            "({:0} instead of {:1})".format(len(param_vals),
                                                            len(par_names))))
        params = {n: v for n, v in zip(par_names, param_vals)}
        t[0] = body.eval(params)
    else:
        raise(NameError("Unknown function "+t[1]))


def p_expr_list(t):
    '''
    exprlist : exprlist COMMA expr
             | expr
    '''
    if len(t) == 2:
        t[0] = [t[1]]
    else:
        t[0] = t[1]
        t[0].append(t[3])


def p_expr_group(t):
    'expr : LPAR expr RPAR'
    t[0] = t[2]


def p_expr_const(t):
    '''
    expr : INT
    expr : NUMBER
    '''
    t[0] = CN(t[1])


def p_expr_var(t):
    'expr : IDENTIFIER'
    t[0] = ex.VarExpr(t[1])


def p_error(t):
    if t is None:
        raise(SyntaxError("Syntax error: premature end of line"))
    else:
        raise(SyntaxError("Syntax error at '%s'" %  t.value))


da_parser = yacc.yacc(optimize=1, tabmodule='dacalc.daparsetab')


def parse(s):
    return da_parser.parse(s, lexer=da_lexer)


def main():
    from dacalc.pygment_lexer import DALexer
    from prompt_toolkit import PromptSession
    from prompt_toolkit.shortcuts import prompt
    from prompt_toolkit.lexers import PygmentsLexer
    
    if len(sys.argv) > 1:
        # do some more sophisticated argument parsing here in the future
        if len(sys.argv) == 2:
            f = open(sys.argv[1], 'r')
            script = f.read()
            f.close()
            parse(script)
        else:
            print("Expecting at most 1 commandline argument", file=sys.stderr)
            exit(1)
    else:
        # interactive mode
        print('Welcome to the dimensional analysis calculator!')
        print('Try "?" for help...')
        session = PromptSession(lexer=PygmentsLexer(DALexer))
        while True:
            try:
                s = session.prompt('DA > ')
            except EOFError:
                break
            try:
                results = parse(s)
                for item in results:
                    if isinstance(item, str) and item != "":
                        pager(item)
                    elif isinstance(item, numpy.ndarray):
                        plt.imshow(item)
                        plt.show()
#                        plt.close()
            except TypeError as te:
                print(te, file=sys.stderr)
            except Exception as exc:
                print(exc, file=sys.stderr)

        print("\nGoodbye...")


if __name__ == '__main__':
    main()
