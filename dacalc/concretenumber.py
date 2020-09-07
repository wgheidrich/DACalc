
import math
import sys


##########
# Helper Stuff for superscript syntax
##########


# unicode superscript tables for integers
superscript_chars = {
    u'0': u'\u2070',
    u'1': u'\xb9',
    u'2': u'\xb2',
    u'3': u'\xb3',
    u'4': u'\u2074',
    u'5': u'\u2075',
    u'6': u'\u2076',
    u'7': u'\u2077',
    u'8': u'\u2078',
    u'9': u'\u2079',
    u'+': u'\u207a',
    u'-': u'\u207b',
    u'^': ''          # fake character to simplify sting conversion
}

superscript_ints = {
    u'\u2070' : 0,
    u'\xb9'   : 1,
    u'\xb2'   : 2,
    u'\xb3'   : 3,
    u'\u2074' : 4,
    u'\u2075' : 5,
    u'\u2076' : 6,
    u'\u2077' : 7,
    u'\u2078' : 8,
    u'\u2079' : 9,
}


def sup_int(s):
    if s[0] == superscript_chars['-']:
        return -sup_int(s[1:])
    elif s[0] == superscript_chars['+']:
        return sup_int(s[1:])
    else:
        val = 0
        for char in s:
            val = 10 * val + superscript_ints[char]
        return val



##########
# ConcreteNumber class
##########

class ConcreteNumber:
    # dimesnions of a dimensionless variable
    # the entries are in order of:
    #   - mass [kg]
    #   - length [m]
    #   - time [s]
    #   - current [A]
    #   - temperature [K]
    #   - amount of substance [mol]
    #   - luminance [cd]
    #   - angle or solid angle [rad or sr] as pseudo-units 
    DIMENSION_LESS = (0,0,0,0,0,0,0,0)
    
    # dimensions of each SI base unit (and for rad and sr)
    KG =   (1,0,0,0,0,0,0,0)
    M =    (0,1,0,0,0,0,0,0)
    S =    (0,0,1,0,0,0,0,0)
    A =    (0,0,0,1,0,0,0,0)
    K =    (0,0,0,0,1,0,0,0)
    MOL =  (0,0,0,0,0,1,0,0)
    CD =   (0,0,0,0,0,0,1,0)
    RAD =  (0,0,0,0,0,0,0,1)
    SR =   (0,0,0,0,0,0,0,2)

    # all the recognized SI prefixes used for parsing
    SI_PREFIXES = {
        "y": 1e-24, # yocto
        "z": 1e-21, # zepto
        "a": 1e-18, # atto
        "f": 1e-15, # femto
        "p": 1e-12, # pico
        "n": 1e-9,  # nano
        "u": 1e-6,  # micro
        "m": 1e-3,  # milli
        "c": 1e-2,  # centi
        "d": 1e-1,  # deci
        "da": 1e1,  # deca
        "h": 1e2,   # hecto
        "k": 1e3,   # kilo
        "M": 1e6,   # mega
        "G": 1e9,   # giga
        "T": 1e12,  # tera
        "P": 1e15,  # peta
        "E": 1e18,  # exa
        "Z": 1e21,  # zetta
        "Y": 1e24,  # yotta
    }

    # SI prefixes for output (only powers of 3 as well as the empty prefix),
    SI_OUT_PREF = [
        (1e-24, "y"),
        (1e-21, "z"),
        (1e-18, "a"),
        (1e-15, "f"),
        (1e-12, "p"),
        ( 1e-9, "n"),
        ( 1e-6, "u"),
        ( 1e-3, "m"),
        (    1, "" ),
        (  1e3, "k"),
        (  1e6, "M"),
        (  1e9, "G"),
        ( 1e12, "T"),
        ( 1e15, "P"),
        ( 1e18, "E"),
        ( 1e21, "Z"),
        ( 1e24, "Y")
    ]

    # dictionary of constants
    const = {}

    # dictionary of recognized units
    units = {}
    
    # number formatting string
    format_str = "{:.6g}"

    # whether to use superscript for exponents
    use_super = True
    
    # whether to use SI prefixes for output of SI units
    # (single units only, not composite ones like [N m])
    use_si_prefix = True
    
    # preferred units for output
    # each key is a dimension tuple
    # each entry is a pair (ConcreteNumber unit, string uname bool prefix)
    # uname is the string representation of the unit
    # if prefix == True then SI prefixes are allowed to modify the unit
    preferred_units = {}

    # base units
    # the base units for the 7 dimensions. Used to piece together unit
    # stings for quantitiews that don't mach in "preferred units"
    base_units = {}
    
    @classmethod
    def use(cls,unit,uname,prefix=False):
        '''
        specify a specific unit to be used for output of variables
        ConcreteNumbers of the respective dimension
        
        Parameters:
            unit   : a ConcreteNumber representign the unit
            uname  : unit string to be displayed when using this unit
            prefix : boolean specifying whether to use SI prefixes

        Returns:
            None
        '''
        cls.preferred_units[unit.units] = (unit,uname,prefix)
    
    @classmethod
    def use_system(cls,sys):
        '''
        Select a new defaule unit system for showing outputs
        
        Parameters:
            sys : name of the unit system (e.g. "si", or "us")

        Returns:
            None
        '''
        cls.preferred_units.clear()
        if sys == "SI_base":
            # only the base units
            # meaning that preferred unit remain empty
            cls.base_units = unit_systems["si"]
        elif sys == "SI":
            # set base units 
            cls.base_units = unit_systems["si"]
            # add all SI-class units as preferences
            for u in si_u:
                cls.use(u[1],u[0],True)
        elif sys == "US_base":
            # only the base units
            cls.base_units = unit_systems["us"]
        elif sys == "US":
            # set base units
            cls.base_units = unit_systems["us"]
            # add US units from a separate list of deault units
            for u in default_us_u:
                cls.use(u[1],u[0])
        else:
            raise TypeError("Unknown unit system",sys)
        

    @classmethod
    def set_precision(cls,prec):
        '''
        specify the precision (significant digits) in an output
        
        Parameters:
            prec : integer
        
        Returns:
            None
        '''
        cls.format_str = "{:."+str(prec)+"g}"

    @classmethod
    def make_super(cls,s):
        '''
        convert well formatted unit string with exponents in '^' notation to
        superscript if we are in unicode mode
        
        Parameters:
            s : unit string
        
        Returns:
            prettified unit string
        '''
        if cls.use_super:
            return "".join(superscript_chars.get(c,c) for c in s)
        else:
            return s

    @classmethod
    def analyze(cls,var,target):
        '''
        perform a dimensional analysis of a set of variables -
        ie. try to matche the target dimensions by combining different
        powers of the provided quantities 
        
        Parameters:
            - var:      a list of ConcreteNumbers
        
            - traget:   ConcreteNumber representing the target dimensions

        Returns:
            a list of solutions
            each solution is a list of terms, and each term is a tuple
            of one of the input quantities and an exponent (in the same
            order provided in the vars argument)
        '''
        def tuple_rank(t):
            maxval = 0
            td = 0
            lex = 0
            for i in t:
                ia = math.fabs(i)
                if ia> maxval:
                    maxval = ia
                td += ia
                lex = lex/10.0 + ia + (0 if i>= 0 else 0.000001)
            return (maxval,td,lex)
        
        def mk_indices(num_var):
            idx = [tuple()]
            for i in range (0,num_var):
                idx = [j+(k,) for k in range(-3,4) for j in idx]
            list.sort(idx,key=tuple_rank)
            return idx
        
        def is_multiple(exp1,exp2):
            mult = 0
            for e1,e2 in zip(exp1,exp2):
                if e2!= 0:
                    mult = e1 // e2
            for e1,e2 in zip(exp1,exp2):
                if e1 != mult * e2:
                    return False
            return True
            
        idx_list = mk_indices(len(var))
        sol_exps = []
        for idx in idx_list:
            soln = ConcreteNumber()
            for i,v in zip(idx,var):
                soln *= v**i
            if soln.units == target.units:
                redundant = False
                for s in sol_exps:
                    redundant |= is_multiple(idx,s)
                if not redundant:
                    sol_exps.append(idx)
        # convert exponent list to the actual solution format
        return [[(v,e) for v,e in zip(var,exp)] for exp in sol_exps]

    
    @classmethod
    def u(cls,unit_str):
        '''
        constructor from a unit string

        Parameters:
            unit_str:   a unit string to be parsed by unitparser
        
        Returns:
            a new ConcreteNumber object representing the specified units
        '''
        return dacalc.unitparser.parse(unit_str)

    
    def __init__(self,val=1,unit_tuple=(0,0,0,0,0,0,0,0)):
        '''
        create new ConcreteNumber object
        
        Parameters:
            val :        numerical value
            unit_tuple : an 8-tuple specifying the exponents for the 8
                         base dimensions (kg,m,s,A,K,mol,cd,rad)
        '''
        self.value = val
        self.units = unit_tuple
        

    def unitstr(self):
        '''
        make a string describing the dimensions in the currently chosen
        unit system
        
        Parameters:
            None
        
        Returns:
            unit sting representing the value using only base units
        '''
        s = ""
        # process all but the radians
        for u,exp in zip(ConcreteNumber.base_units,self.units[:-1]):
            if exp!= 0:
                s += u[0]
                if exp!= 1:
                    s += "^" + str(exp)
                s += ' '
        # deal with radians and steradians
        rad = self.units[-1]
        if rad != 0:
            if rad == 2:
                s += "sr"
            elif rad == -2:
                s += "sr^-1"
            else:
                s+= "rad" + (("^"+str(rad)) if rad != 1 else "")
        return ConcreteNumber.make_super(s.strip())
        
    
    def check_dimensionless(self):
        '''
        check if the value is dimension-less, raise TypeError if not
        
        Parameters:
            None
        
        Returns:
            None
        '''
        if self.units[:-1] != ConcreteNumber.DIMENSION_LESS[:-1]:
            raise TypeError("Expecting dimensionless value, got ["
                            + self.unitstr() + "]")

        
    def check_dim_match(self,other,warnrad=True):
        '''
        check if the units match with the units of another ConcreteNumber
        raise TypeError if not
        
        Mismatches in radians and steradians do not lead to an exception,
        but possibly a warning message.
        
        Parameters:
            other    : the other ConcreteNumber
            warnrad  : whether to pritn a waring message if radians mismatch
        
        Returns:
            None
        '''
        if self.units[:-1] != other.units[:-1]:
            raise TypeError("Dimension mismatch between [" + self.unitstr()
                                + "] and [" + other.unitstr() + "]")
        if warnrad and self.units[-1] != other.units[-1]:
            print("WARNING: radians don't match in [" + self.unitstr()
                            + "] and [" + other.unitstr() + "]",
                  file = sys.stderr)
            

    def convert_to_base(self):
        '''
        convert value from the internal representation (SI base) to
        whatever other base system is selected via the base_units
        class variable
        
        Parameters:
            None
        
        Returns:
            floating point value corresponding to self.value in base units
        '''
        new_u = ConcreteNumber(1) # start dimensionless
        for i in range(0,7):
            new_u *= ConcreteNumber.base_units[i][1] ** self.units[i]
        return self.value / new_u.value
    
        
    def __str__(self,display_units = None):
        '''
        string reprsentation of the ConcreteNumber
        
        If display_units is provided, it represents a unit conversion
        request. If this conversion is not possible due to a dimension
        mismatch, then two unit strings are produced, e.g.:

            ConcreteNumber.u("m/s").__str__("m")

        will produce the output
        
            1 [Hz] [m]
        
        i.e. the residual in dimenstion (in this case s^-1) will be
        represented in the best fitting units, followed by the desired
        units.

        If display_units is None, and  dimensions of the object match
        those of a preferred unit, we convert the value to that unit.
        Otherwise we represent the value in the chosen base_units.
        
        Parameters:
            display_units : unit string for representing the value
                            (default: None)
        
        Returns:
            string
        '''
        # first we check if specific display units have been provided
        if display_units != None:
            du = ConcreteNumber.u(display_units)
            return str(self/du)+' ['+ConcreteNumber.make_super(display_units)+']'
        
        # dimensionless variables are just the numerical values
        if self.units == ConcreteNumber.DIMENSION_LESS:
            return ConcreteNumber.format_str.format(self.value)

        # for everything else we first check if one of the preferred
        # units matches
        try:
            pref_unit = ConcreteNumber.preferred_units[self.units]
            converted = self/pref_unit[0]
            if ConcreteNumber.use_si_prefix and pref_unit[2]:
                # calculate the SI prefix
                (mult,pref) = ConcreteNumber.SI_OUT_PREF[-1]
                for (m,p) in ConcreteNumber.SI_OUT_PREF:
                    if converted.value< 900*m:
                        (mult,pref) = (m,p)
                        break
                return (ConcreteNumber.format_str.format(converted.value/mult)
                        + " [" + pref+pref_unit[1]+ "]")
            else:
                return (ConcreteNumber.format_str.format(converted.value)
                        + " [" + pref_unit[1] + "]")
        except:
            # if there are no preferred units, we just use the
            # standard unit string conversion
            return (ConcreteNumber.format_str.format(self.convert_to_base())
                    + " [" + self.unitstr() + "]")

    
    ### arithmetic ops
    
    def __add__(self,other):
        '''
        addition operator for two ConcreteNumbers
        
        The dimensions between both ConcreteNumbers must match
        
        Parameters:
            other : a second ConcreteNumber with matching dimensions
        
        Returns:
            sum as ConcreteNumber
        '''
        if not isinstance(other, ConcreteNumber):
            return NotImplemented
        self.check_dim_match(other)
        return ConcreteNumber(self.value+other.value, self.units)

            
    def __sub__(self,other):
        '''
        subtraction operator for two ConcreteNumbers
        
        The dimensions between both ConcreteNumbers must match
        
        Parameters:
            other : a second ConcreteNumber with matching dimensions
        
        Returns:
            difference as ConcreteNumber
        '''
        if not isinstance(other, ConcreteNumber):
            return NotImplemented
        self.check_dim_match(other)
        return ConcreteNumber(self.value-other.value, self.units)


    def __neg__(self):
        '''
        unary minus operator
        
        Parameters:
            None
        
        Returns:
            negative of self
        '''
        return ConcreteNumber(-self.value, self.units)
    

    def __mul__(self,other):
        '''
        multiplication operator
        
        Parameters:
            other : a second ConcreteNumber, an int or a float
        
        Returns:
            product as a ConcreteNumber
        '''
        if isinstance(other, ConcreteNumber):
            return ConcreteNumber(self.value*other.value,
                   tuple(u1+u2 for (u1,u2) in zip(self.units,other.units)))
        elif type(other) == float or type(other) == int:
            return ConcreteNumber( self.value*other, self.units )
        else:
            return NotImplemented

            
    def __rmul__(self,other):
        '''
        right-sided multiplication operator
        
        Parameters:
            other : an int or a float
        
        Returns:
            product as a ConcreteNumber
        '''
        # replicating is a bit faster than adding another function call
        # by calling __mul__
        if isinstance(other, ConcreteNumber):
            return ConcreteNumber(self.value*other.value,
                   tuple(u1+u2 for (u1,u2) in zip(self.units,other.units)))
        elif type(other) == float or type(other) == int:
            return ConcreteNumber( self.value*other, self.units )
        else:
            return NotImplemented

    
    def __truediv__(self,other):
        '''
        division operator
        
        Parameters:
            other : a second ConcreteNumber, an int or a float
        
        Returns:
            division as a ConcreteNumber
        '''
        if isinstance(other, ConcreteNumber):
            return ConcreteNumber(self.value/other.value,
                                  tuple(u1-u2 for (u1,u2) in zip(self.units,other.units)))
        elif type(other) == float or type(other) == int:
            return ConcreteNumber( self.value/other, self.units )
        else:
            return NotImplemented
       
    
    def __rtruediv__(self,other):
        '''
        right-sided division operator
        
        Parameters:
            other : an int or a float
        
        Returns:
            product as a ConcreteNumber
        '''
        if type(other) == float or type(other) == int:
            return ConcreteNumber(val=other)/self
        else:
            return NotImplemented

    
    def __pow__(self,exp):
        '''
        power operator
        
        Parameters:
            exp : integer exponent
        
        Returns:
            self ^ exp
        '''
        if type(exp) == int:
            return ConcreteNumber(self.value**exp,
                                  tuple(u*exp for u in self.units))
        else:
            return NotImplemented
    



    
    
### math functions using ConcreteNumber

def root(val,exp):
    '''
    The <exp>-th root
    
    Parameters:
        val:    a concrete value, where the units must appear in powers of <exp>
        exp:    integer exponent

    Returns:
        the <exp>-th root of <val>
    '''
    param: val 
    if type(exp) == int:
        # check if unit exponents are divisible by the root exponent
        for u in val.units[:-1]:
            if (u%exp) != 0:
                raise TypeError("Cannot take the " + str(exp) + "-root of ["
                                + str(val.unitstr()) + "]")
        return ConcreteNumber(math.pow(val.value,1/exp),
                              tuple(int(u/exp) for u in val.units))
    else:
        return NotImplemented

        
def sqrt(val):
    '''
    Square root
    
    Parameters:
        val:    a concrete value where all SI base units appear in even powers
    
    Returns:
        the square root of <val>
    '''
    return root(val,2)
 

def sin(angle):
    '''
    Sine
    
    Parameters:
        angle:  the angle (in radians) with units of [] or [rad]
    
    Returns:
        the sine of <angle> in units of []
    '''
    angle.check_dimensionless()
    val = math.sin(angle.value)
    if angle.units == ConcreteNumber.SR: # steradians become radians
        return ConcreteNumber(val,ConcreteNumber.RAD)
    else: # everything else becomes dimensionless
        return ConcreteNumber(val)

    
def cos(angle):
    '''
    Cosine
    
    Parameters:
        angle:  the angle (in radians) with units of [] or [rad]
    
    Returns:
        the cosine of <angle> in units of []
    '''
    angle.check_dimensionless()
    val = math.cos(angle.value)
    if angle.units == ConcreteNumber.SR: # steradians become radians
        return ConcreteNumber(val,ConcreteNumber.RAD)
    else: # everything else becomes dimensionless
        return ConcreteNumber(val)


def tan(angle):
    '''
    Tangens
    
    Parameters:
        angle:  the angle (in radians) with units of [] or [rad]
    
    Returns:
        the tangens of <angle> in units of []
    '''
    angle.check_dimensionless()
    val = math.tan(angle.value)
    if angle.units == ConcreteNumber.SR: # steradians become radians
        return ConcreteNumber(val,ConcreteNumber.RAD)
    else: # everything else becomes dimensionless
        return ConcreteNumber(val)


def asin(val):
    '''
    Arc Sine
    
    Parameters:
        val:    a scalar value (units of [])
    
    Returns:
        the arc sine of <val> in units of [rad]
    '''
    val.check_dimensionless()
    angle = math.asin(val.value)
    if val.units == ConcreteNumber.RAD: # radians become steradians
        return ConcreteNumber(angle,ConcreteNumber.SR)
    else: # everything else becomes radians
        return ConcreteNumber(angle,ConcreteNumber.RAD)

    
def acos(val):
    '''
    Arc Cosine
    
    Parameters:
        val:  a scalar value (units of [])
    
    Returns:
        the arc cosine of <val> in units of [rad]
    '''
    val.check_dimensionless()
    angle = math.acos(val.value)
    if val.units == ConcreteNumber.RAD: # radians become steradians
        return ConcreteNumber(angle,ConcreteNumber.SR)
    else: # everything else becomes radians
        return ConcreteNumber(angle,ConcreteNumber.RAD)

    
def atan(val):
    '''
    Arc Tangens
    
    Parameters:
        val:    a scalar value (units of [])
    
    Returns:
        the arc tangens of <val> in units of [rad]
    '''
    val.check_dimensionless()
    angle = math.atan(val.value)
    if val.units == ConcreteNumber.RAD: # radians become steradians
        return ConcreteNumber(angle,ConcreteNumber.SR)
    else: # everything else becomes radians
        return ConcreteNumber(angle,ConcreteNumber.RAD)

    
def atan2(val1,val2):
    '''
    Arc Tangens with two arguments
    
    Parameters:
        val1:   a concrete number
        val2:   a concrete number with the same units as <val1>
    
    Returns:
        the arc tangens of <val1>/<val2> in units of [rad]
    '''
    val1.check_dim_match(val2,warnrad=False)
    angle = math.atan2(val1.value,val2.value)
    if val1.units == ConcreteNumber.RAD: # radians become steradians
        return ConcreteNumber(angle,ConcreteNumber.SR)
    else: # everything else becomes radians
        return ConcreteNumber(angle,ConcreteNumber.RAD)

    
def fabs(val):
    '''
    Absolute value
    
    Parameters:
        val:    a concrete number
    
    Returns:
        the absolute value of <val>, with the same units as <val>
    '''
    return ConcreteNumber(val=math.fabs(val.value),unit_tuple=val.units)


def log(val,base):
    '''
    General logarithm
    
    Parameters:
        val:    a scalar value (units of [])
        base:   a float, int, or scalar value (units of [])
    
    Returns:
        the logarithm base <base> of <val>
    '''
    val.check_dimensionless()
    if type(base) == float or type(base) == int:
        base.check_dimensionless()
        return ConcreteNumber(val=math.log(val.value,base))
    else:
        return ConcreteNumber(val=math.log(val.value,base.value))


def ln(val):
    '''
    Logarithm base 2

    Parameters:
        val:    a scalar value (units of [])

    Returns:
        natural logarithm of <val>
    '''
    val.check_dimensionless()
    return ConcreteNumber(val=math.log(val.value))
    

def log2(val):
    '''
    Logarithm base 2

    Parameters:
        val:    a scalar value (units of [])

    Returns:
        logarithm base 2 of <val>
    '''
    val.check_dimensionless()
    return ConcreteNumber(val=math.log2(val.value))


def log10(val):
    '''
    Logarithm base 10

    Parameters:
        val:    a scalar value (units of [])

    Returns:
        logarithm base 10 of <val>
    '''
    val.check_dimensionless()
    return ConcreteNumber(val=math.log10(val.value))


def pow(val,exp):
    '''
    Power

    Parameters:
        val:    a scalar value (units of [])
        exp:    a scalar value (units of [])

    Returns:
        <val> raised to the power of <exp>
    '''
    if type(exp) == int:
        return val**exp
    elif type(exp) == float:
        val.check_dimensionless()
        return ConcreteNumber(val=math.pow(val.value,exp))
    elif isinstance(exp, ConcreteNumber):
        val.check_dimensionless()
        exp.check_dimensionless()
        return ConcreteNumber(val=math.pow(val.value,exp.value))        
    else:
        return NotImplemented


# add a single unit to the unit dictionary (check for name clashes)
def add_unit_to_dict(name,newunit):
    if name in ConcreteNumber.units.keys():
        print("Warning: unit " + name
              + "already exists (overwriting by new definition).",
              file = sys.stderr)
    ConcreteNumber.units[name] = newunit

# add a new unit both to a unit group list and the global unit
# dictionary; define all SI prefixes if needed
def new_unit(unitgroup,unit,si_prefixes=False):
    unitgroup.append(unit)
    add_unit_to_dict(unit[0],unit[1])
    if si_prefixes:
        for prefix, multiplier in ConcreteNumber.SI_PREFIXES.items():
            add_unit_to_dict(prefix+unit[0], multiplier*unit[1])
        


# the folllwing are lists of unit groups, which consist of tuples
# (name string, ConcreteNumber representation, description string)

# add units of a dimensionless variable
add_unit_to_dict("",ConcreteNumber(1,ConcreteNumber.DIMENSION_LESS))

# with the class defined we can do a circular import of the
# unitparser, so that definign the units becomes a bit less verbose
import dacalc.unitparser

# list of all )base+derived) SI units (without prefix)
si_u = []
# we define each new unit and directly add it, so that we can use each
# new unit right away
new_unit(si_u,("g",   ConcreteNumber(1e-3,ConcreteNumber.KG), "gram"), True)
new_unit(si_u,("m",   ConcreteNumber(1,ConcreteNumber.M),     "metre"), True)
new_unit(si_u,("s",   ConcreteNumber(1,ConcreteNumber.S),     "second"), True)
new_unit(si_u,("A",   ConcreteNumber(1,ConcreteNumber.A),     "Ampere"), True)
new_unit(si_u,("K",   ConcreteNumber(1,ConcreteNumber.K),     "Kelvin"), True)
new_unit(si_u,("mol", ConcreteNumber(1,ConcreteNumber.MOL),   "Mole"), True)
new_unit(si_u,("cd",  ConcreteNumber(1,ConcreteNumber.CD),    "Candela"), True)
new_unit(si_u,("rad", ConcreteNumber(1,ConcreteNumber.RAD),   "radian"))

# SI derived units
new_unit(si_u,("sr",  ConcreteNumber.u("rad^2"),      "steradian"))
new_unit(si_u,("N",   ConcreteNumber.u("kg m / s^2"), "Newton"), True)
new_unit(si_u,("Pa",  ConcreteNumber.u("N / m^2"),    "Pascal"), True)
new_unit(si_u,("J",   ConcreteNumber.u("N m"),        "Joule"), True)
new_unit(si_u,("W",   ConcreteNumber.u("J / s"),      "Watt"), True)
new_unit(si_u,("C",   ConcreteNumber.u("s A"),        "Coulomb"), True)
new_unit(si_u,("V",   ConcreteNumber.u("W / A"),      "Volt"), True)
new_unit(si_u,("F",   ConcreteNumber.u("C / V"),      "Farad"), True)
new_unit(si_u,("Ohm", ConcreteNumber.u("V / A"),      "Ohm"), True)
new_unit(si_u,("S",   ConcreteNumber.u("Ohm^-1"),     "Siemens"), True)
new_unit(si_u,("Wb",  ConcreteNumber.u("V s"),        "Weber"), True)
new_unit(si_u,("T",   ConcreteNumber.u("Wb / m^2"),   "Tesla"), True)
new_unit(si_u,("H",   ConcreteNumber.u("Wb / A"),     "Henry"), True)
new_unit(si_u,("lm",  ConcreteNumber.u("cd sr"),      "Lumen"), True)
new_unit(si_u,("lx",  ConcreteNumber.u("cd / m^2"),   "Lux"), True)
new_unit(si_u,("Bq",  ConcreteNumber.u("1 / s"),      "Bequerel"), True)
#  after Bequerel, so that Hertz is default unit for 1/s
new_unit(si_u,("Hz",  ConcreteNumber.u("1 / s"),      "Hertz"), True) 
new_unit(si_u,("Gy",  ConcreteNumber.u("J / kg"),     "Gray"), True)
new_unit(si_u,("Sv",  ConcreteNumber.u("Gy"),         "Sievert"), True)
new_unit(si_u,("kat", ConcreteNumber.u("mol / s"),    "Katal"), True)


#new_unit(si_u,("sr",  units["rad"]*units["rad"],              "steradian"))
#new_unit(si_u,("N",   units["kg"]*units["m"]*units["s"]**-2,  "Newton"), True)
#new_unit(si_u,("Pa",  units["N"]/units["m"]**2,               "Pascal"), True)
#new_unit(si_u,("J",   units["N"]*units["m"],                  "Joule"), True)
#new_unit(si_u,("W",   units["J"]/units["s"],                  "Watt"), True)
#new_unit(si_u,("C",   units["s"]*units["A"],                  "Coulomb"), True)
#new_unit(si_u,("V",   units["W"]/units["A"],                  "Volt"), True)
#new_unit(si_u,("F",   units["C"]/units["V"],                  "Farad"), True)
#new_unit(si_u,("Ohm", units["V"]/units["A"],                  "Ohm"), True)
#new_unit(si_u,("S",   units["Ohm"]**-1,                       "Siemens"), True)
#new_unit(si_u,("Wb",  units["V"]*units["s"],                  "Weber"), True)
#new_unit(si_u,("T",   units["Wb"]*units["m"]**-2,             "Tesla"), True)
#new_unit(si_u,("H",   units["Wb"]/units["A"],                 "Henry"), True)
#new_unit(si_u,("lm",  units["cd"]*units["sr"],                "Lumen"), True)
#new_unit(si_u,("lx",  units["cd"]/units["m"]**2,              "Lux"), True)
#new_unit(si_u,("Bq",  units["s"]**-1,                         "Bequerel"), True)
#new_unit(si_u,("Hz",  units["s"]**-1,                         "Hertz"), True) #  after Bequerel, so that Hertz is default unit for 1/s
#new_unit(si_u,("Gy",  units["J"]/units["kg"],                 "Gray"), True)
#new_unit(si_u,("Sv",  units["Gy"],                            "Sievert"), True)
#new_unit(si_u,("kat", units["mol"]/units["s"],                "Katal"), True)



# miscellaneaous units (including some outdated/dperecated ones)
misc_u = []
new_unit(misc_u,("min",    60*ConcreteNumber.u("s"),             "minute"))
new_unit(misc_u,("h",      60*ConcreteNumber.u("min"),           "hour"))
new_unit(misc_u,("d",      24*ConcreteNumber.u("h"),             "day"))
new_unit(misc_u,("deg",    math.pi/180.0*ConcreteNumber.u("rad"),"degree"))
new_unit(misc_u,("arcmin", ConcreteNumber.u("deg")/60.0,         "arc minute"))
new_unit(misc_u,("arcsec", ConcreteNumber.u("arcmin")/60.0,      "arc second"))
new_unit(misc_u,("NM",     1852*ConcreteNumber.u("m"),         "nautical mile"))
new_unit(misc_u,("knot",   ConcreteNumber.u("NM / h"),           "knot"))
new_unit(misc_u,("ha",     ConcreteNumber.u("hm^2"),             "hectare"))
new_unit(misc_u,("L",      ConcreteNumber.u("dm^3"),             "litre"), True)
new_unit(misc_u,("t",      1e3*ConcreteNumber.u("kg"),           "tonne"))
new_unit(misc_u,("Da" ,    1.66053904020e-27*ConcreteNumber.u("kg"), "Dalton"))
new_unit(misc_u,("dyn",    1e-5*ConcreteNumber.u("N"),           "Dyne"))
new_unit(misc_u,("rpm" ,   2*math.pi*ConcreteNumber.u("rad/min"),"rpm"))
new_unit(misc_u,("eV",     1.602176634e-19*ConcreteNumber.u("J"),"electron volt"),True)
new_unit(misc_u,("erg",    1e-7*ConcreteNumber.u("J"),           "erg"))
new_unit(misc_u,("cal",    4.184*ConcreteNumber.u("J"),          "calorie"))
new_unit(misc_u,("kcal",   4184*ConcreteNumber.u("J"),           "kilocalorie"))
new_unit(misc_u,("torr",   133.3224*ConcreteNumber.u("Pa"),      "Torr"))
new_unit(misc_u,("atm",    101325*ConcreteNumber.u("Pa"),        "atmosphere"))
new_unit(misc_u,("bar",    1e5*ConcreteNumber.u("Pa"),           "Bar"),True)
new_unit(misc_u,("nt",     ConcreteNumber.u("cd / m^2"),         "nit"))
new_unit(misc_u,("ct",     .2*ConcreteNumber.u("g"),             "carat"))
new_unit(misc_u,("degC",   ConcreteNumber.u("K"),             "deg. Celsius"))
new_unit(misc_u,("degF",   (5.0/9.0)*ConcreteNumber.u("K"),   "deg. Farenheit"))

#new_unit(misc_u,("min",    60*units["s"],               "minute"))
#new_unit(misc_u,("h",      60*units["min"],             "hour"))
#new_unit(misc_u,("d",      24*units["h"],               "day"))
#new_unit(misc_u,("deg",    math.pi/180.0*units["rad"],  "degree"))
#new_unit(misc_u,("arcmin", units["deg"]/60.0,           "arc minute"))
#new_unit(misc_u,("arcsec", units["arcmin"]/60.0,        "arc second"))
#new_unit(misc_u,("NM",     1852*units["m"],             "nautical mile"))
#new_unit(misc_u,("knot",   units["NM"]/units["h"],      "knot"))
#new_unit(misc_u,("ha",     units["hm"]**2,              "hectare"))
#new_unit(misc_u,("L",      units["dm"]**3,              "litre"), True)
#new_unit(misc_u,("t",      1e3*units["kg"],             "tonne"))
#new_unit(misc_u,("Da" ,    1.66053904020e-27*units["kg"], "Dalton"))
#new_unit(misc_u,("dyn",    1e-5*units["N"],             "Dyne"))
#new_unit(misc_u,("rpm" ,   2*math.pi*units["rad"]/units["min"], "rpm"))
#new_unit(misc_u,("eV",     1.602176634e-19*units["J"],  "electron volt"),True)
#new_unit(misc_u,("erg",    1e-7*units["J"],             "erg"))
#new_unit(misc_u,("cal",    4.184*units["J"],            "calorie"))
#new_unit(misc_u,("kcal",   4184*units["J"],             "kilocalorie"))
#new_unit(misc_u,("torr",   133.3224*units["Pa"],        "Torr"))
#new_unit(misc_u,("atm",    101325*units["Pa"],          "atmosphere"))
#new_unit(misc_u,("bar",    1e5*units["Pa"],             "Bar"),True)
#new_unit(misc_u,("nt",     units["cd"]/units["m"]**2,   "nit"))
#new_unit(misc_u,("ct",     .2*units["g"],               "carat"))
#new_unit(misc_u,("degC",   units["K"],                  "deg. Celsius"))
#new_unit(misc_u,("degF",   (5.0/9.0)*units["K"],        "deg. Farenheit"))


# list of US/imperial units (including some older ones)
# Note: different for the UK/imperial ones especially on the volumes!
us_u = []
# US length
new_unit(us_u,("in",   25.4*ConcreteNumber.u("mm"),            "inch"))
new_unit(us_u,("ft",   12*ConcreteNumber.u("in"),              "foot"))
new_unit(us_u,("yd",   36*ConcreteNumber.u("in"),              "yard"))
new_unit(us_u,("mi",   63360*ConcreteNumber.u("in"),           "mile"))
new_unit(us_u,("mil",  ConcreteNumber.u("in")/1000,            "mil"))
new_unit(us_u,("thou", ConcreteNumber.u("mil"),                "thou"))
new_unit(us_u,("uin",  ConcreteNumber.u("mil")/1000,           "microinch"))
# US area
new_unit(us_u,('acre', (43560*(1200.0/3937)**2)*ConcreteNumber.u("m^2"),"acre"))
# US speed
new_unit(us_u,("mph",  ConcreteNumber.u("mi / h"),             "mph"))
# US weight
new_unit(us_u,("lb",   0.45359237*ConcreteNumber.u("kg"),      "pound"))
new_unit(us_u,("ton",  2000*ConcreteNumber.u("lb"),            "ton"))
new_unit(us_u,("oz",   ConcreteNumber.u("lb")/16,              "ounze"))
new_unit(us_u,("gr",   ConcreteNumber.u("lb")/7000,            "grain"))
new_unit(us_u,("ozt",  480*ConcreteNumber.u("gr"),             "troy_ounce"))
new_unit(us_u,("lbt",  12*ConcreteNumber.u("ozt"),             "troy_pound"))
# US force
new_unit(us_u,("lbf",  9.80665*ConcreteNumber.u("lb m / s^2"), "pound_force"))
# US volume (fluid -- didn't bother with dry volumes)
new_unit(us_u,("gal",  3.785411784*ConcreteNumber.u("L"),      "gallon"))
new_unit(us_u,("qt",   ConcreteNumber.u("gal")/4,              "quart"))
new_unit(us_u,("pt",   ConcreteNumber.u("qt")/2,               "pint"))
new_unit(us_u,("cp",   ConcreteNumber.u("pt")/2,               "cup"))
new_unit(us_u,("floz", ConcreteNumber.u("cp")/8,               "fluid_ounce"))
new_unit(us_u,("Tbsp", ConcreteNumber.u("floz")/2,             "tablespoon"))
new_unit(us_u,("tsp",  ConcreteNumber.u("Tbsp")/3,             "teaspoon"))
# US miscellaneous or outdated
new_unit(us_u,("slug", ConcreteNumber.u("lbf s^2 / ft"),       "slug"))
new_unit(us_u,("pdl",  ConcreteNumber.u("lb ft / s^2"),        "poundal"))
new_unit(us_u,("psi",  ConcreteNumber.u("lbf / in^2"),         "PSI"))
new_unit(us_u,("fc",   ConcreteNumber.u("lm / ft^2"),          "foot_candle"))
new_unit(us_u,("lbmol",453.59237*ConcreteNumber.u("mol"),      "pound_mole"))
new_unit(us_u,("degR", (5.0/9.0)*ConcreteNumber.u("K"),        "deg. Rankine"))


# because the US system has multiple units for each dimension, we need
# to pick some specific ones as default units
default_us_u =[
    ("lb",      ConcreteNumber.u("lb")),
    ("ft",      ConcreteNumber.u("ft")), # areas and volumes are ft^2 and ft^3
    ("s",       ConcreteNumber.u("s")),
    ("Hz",      ConcreteNumber.u("Hz")),
    ("A",       ConcreteNumber.u("A")),
    ("degR",    ConcreteNumber.u("degR")), # degree Rankine is the default
    ("lbmol",   ConcreteNumber.u("lbmol")),
    ("cd",      ConcreteNumber.u("cd")),
    ("lbf",     ConcreteNumber.u("lbf")),
    ("lbf ft",  ConcreteNumber.u("ft lbf")),
    ("mph",     ConcreteNumber.u("mph")),
    ("psi",     ConcreteNumber.u("psi")),
    ("fc",      ConcreteNumber.u("fc"))
]
    



# base units for different unit systems
# (list of the seven base units for each of: mass, length, time,
# current, temperature, substance amount, luminous intensity)
# the order of these entries matters, and should align with the dimensions

unit_systems = {
    "si" : [
        ("kg",   ConcreteNumber.u("kg"),     "kilogram"),
        ("m",    ConcreteNumber.u("m"),      "metre"),
        ("s",    ConcreteNumber.u("s"),      "second"),
        ("A",    ConcreteNumber.u("A"),      "Ampere"),
        ("K",    ConcreteNumber.u("K"),      "Kelvin"),
        ("mol",  ConcreteNumber.u("mol"),    "mole"),
        ("cd",   ConcreteNumber.u("cd"),     "candela")
    ],

    "us" : [
        ("lb",    ConcreteNumber.u("lb"),    "pound"),
        ("ft",    ConcreteNumber.u("ft"),    "foot"),
        ("s",     ConcreteNumber.u("s"),     "second"),
        ("A",     ConcreteNumber.u("A"),     "Ampere"),
        ("degR",  ConcreteNumber.u("degR"),  "deg. Rankine"),
        ("lbmol", ConcreteNumber.u("lbmol"), "pound-mole"),
        ("cd",    ConcreteNumber.u("cd"),    "candela")
    ]
}


# initialize the unit system to SI
ConcreteNumber.use_system("SI")



# some mathematical and physical constants
# (this is the version that contains help text)
const = [
    ("pi", math.pi*ConcreteNumber.u(""),                     "\u03C0"),
    ("e",  math.e*ConcreteNumber.u(""),                      "Euler c."),
    ("zeroC",273.15*ConcreteNumber.u("K"),                   "0 [degC] abs."),
    ("zeroF",255.372*ConcreteNumber.u("K"),                  "0 [degF] abs."),
    ("c",  299792458*ConcreteNumber.u("m / s"),              "speed of light"),
    ("G",  6.67408e-11*ConcreteNumber.u("m^3 / kg s^2"),     "gravit. c." ),
    ("g",  9.80665*ConcreteNumber.u("m / s^2"),              "std. gravity"),
    ("h",  6.62607004e-34*ConcreteNumber.u("J s"),           "Planck c."),
    ("ec", 1.602176634e-19*ConcreteNumber.u("C"),            "elem. charge"),
    ("alpha",0.0072973525693*ConcreteNumber.u(""),           "fine struct.c."),
    ("m_e",9.1093837015e-31*ConcreteNumber.u("kg"),          "electron mass"),
    ("m_p",1.67262192369e-27*ConcreteNumber.u("kg"),         "proton mass"),
    ("N_A",6.02214076e23/ConcreteNumber.u("mol"),            "Avogadro c."),
    ("F",  96485.3321233100184*ConcreteNumber.u("C / mol"),  "Faraday c."),
    ("R",  8.31446261815324*ConcreteNumber.u("J / K mol"),   "univ. gas c."),
    ("k",  1.380649e-23*ConcreteNumber.u("J / K"),           "Boltzmann c."),
    ("mu_0",1.25663706212e-6*ConcreteNumber.u("H / m"),      "vac. permeab."),
    ("eps_0",8.8541878128e-12*ConcreteNumber.u("F / m"),     "vac. permit."),
    ("PI", 2.0678338141336296e-15*ConcreteNumber.u("Wb"),    "mag. flux quant")
]

# add pre-defined constants into class dictionary for easier lookup
# and more compact python API syntax
for (c,v,_) in const:
    ConcreteNumber.const[c] = v


