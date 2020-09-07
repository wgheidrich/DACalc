
import sys
import ply.lex as lex
import ply.yacc as yacc

import dacalc.concretenumber as cn
from dacalc.concretenumber import ConcreteNumber as CN


#######
# lexer
#######

# tokens
tokens = (
    'DIV',
    'ONE',
    'SYMBOL',
    'EXP'
)


# single character rules
t_DIV    = r'/'

# symbol name is a letter string
t_SYMBOL = r'[a-zA-Z_]+'

# constant 1
def t_ONE(t):
    r'1'
    t.value = "" # units of dimensionless entity
    return t
    
# signed integer exponent (i.e. leading '^')
def t_EXP(t):
    r'\^-?\d+'
    t.value = int(t.value[1:])
    return t

# signed integer exponent in unicode superscript
def t_EXPUNI(t):
    r'(\u207a|\u207b)?((\u2070|\xb9|\xb2|\xb3|\u2074|\u2075|\u2076|\u2077|\u2078|\u2079)+)'
    t.value = cn.sup_int(t.value)
    t.type = "EXP"
    return t


# Ignored characters
t_ignore = " \t"

def t_newline(t):
    r'\n+'
    t.lexer.lineno += t.value.count("\n")
    
def t_error(t):
    print('Illegal character in unit string:\n  \"... '+t.value+"\"\n       ^",
          file = sys.stderr)
    t.lexer.skip(1)
    
    
# Build the lexer (prepare for optimizing)
ulexer = lex.lex(optimize=1,lextab='dacalc.unitlextab')

# debug function for lexer
def lex_debug(s):
    ulexer.input(s)
    for tok in iter(lex.token, None):
        print(repr(tok.type), repr(tok.value),file = sys.stderr)



########
# parser
########


def p_unit_string_empty(t):
    "unit_string :"
    t[0] = CN.units[""]
    
def p_unit_simple(t):
    "unit_string : unit_seq"
    t[0] = t[1]

def p_unit_divided(t):
    "unit_string : unit_seq DIV unit_seq"
    t[0] = t[1] / t[3]

def p_unitseq(t):
    "unit_seq : unit_seq unit"
    t[0] = t[1]*t[2]

def p_unitseq_single(t):
    "unit_seq : unit"
    t[0] = t[1]

def p_unit(t):
    """
    unit : ONE
    unit : SYMBOL
    unit : SYMBOL EXP
    """
    try:
        t[0] = CN.units[t[1]]
    except KeyError:
        print("Unknown unit: ", t[1],file = sys.stderr)
        t[0] = CN.units[""]
    if len(t) == 3:
        t[0] = t[0]**t[2]

def p_error(p):
    print("Syntax error in unit string",file = sys.stderr)
    print("  ["+current_unit_string+"]",file = sys.stderr)
    print(' '*(p.lexpos+3)+'^',file = sys.stderr)

# bulid paerse in optimize mode
uparser = yacc.yacc(optimize=1,tabmodule='dacalc.unitparsetab')

# copy of the currently processed unit sring, for error messages
current_unit_string = ""

def parse(s):
    global current_unit_string
    current_unit_string = s
    return uparser.parse(s,lexer=ulexer)


if __name__ == '__main__':
    import readline

    while True:
        try:
            s = input('calc > ')
        except EOFError:
            break
        uparser.parse(s,lexer=u_lexer)

        print("\nGoodbye...")
