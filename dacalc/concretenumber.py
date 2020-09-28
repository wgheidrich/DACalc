
import math
import sys
import re


# configuration for flake8
# ignore extra white space errors since we use tabular formatting a lot
# noqa: E201,E222


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
    u'\u2070': 0,
    u'\xb9':   1,
    u'\xb2':   2,
    u'\xb3':   3,
    u'\u2074': 4,
    u'\u2075': 5,
    u'\u2076': 6,
    u'\u2077': 7,
    u'\u2078': 8,
    u'\u2079': 9,
    u'\u207a': "+",
    u'\u207b': "-",
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

    # dimensions of a dimensionless variable
    # the entries are in order of:
    #   - mass [kg]
    #   - length [m]
    #   - time [s]
    #   - current [A]
    #   - temperature [K]
    #   - amount of substance [mol]
    #   - luminance [cd]
    #   - angle or solid angle [rad or sr] as pseudo-units
    #
    # NOTE: we use fixed point representation with the LSB representing 0.5
    #       (this is necessary to represent CGS fractional dimensions)
    DIMENSION_LESS = (0, 0, 0, 0, 0, 0, 0, 0)

    # dimensions of each SI base unit (and for rad and sr)
    KG =   (2, 0, 0, 0, 0, 0, 0, 0)
    M =    (0, 2, 0, 0, 0, 0, 0, 0)
    S =    (0, 0, 2, 0, 0, 0, 0, 0)
    A =    (0, 0, 0, 2, 0, 0, 0, 0)
    K =    (0, 0, 0, 0, 2, 0, 0, 0)
    MOL =  (0, 0, 0, 0, 0, 2, 0, 0)
    CD =   (0, 0, 0, 0, 0, 0, 2, 0)
    RAD =  (0, 0, 0, 0, 0, 0, 0, 2)
    SR =   (0, 0, 0, 0, 0, 0, 0, 4)

    # all the recognized SI prefixes used for parsing
    SI_PREFIXES = {
        "y": 1e-24,      # yocto
        "z": 1e-21,      # zepto
        "a": 1e-18,      # atto
        "f": 1e-15,      # femto
        "p": 1e-12,      # pico
        "n": 1e-9,       # nano
        "u": 1e-6,       # micro as ascii alternative
        "\u03BC": 1e-6,  # micro as greek mu
        "m": 1e-3,       # milli
        "c": 1e-2,       # centi
        "d": 1e-1,       # deci
        "da": 1e1,       # deca
        "h": 1e2,        # hecto
        "k": 1e3,        # kilo
        "M": 1e6,        # mega
        "G": 1e9,        # giga
        "T": 1e12,       # tera
        "P": 1e15,       # peta
        "E": 1e18,       # exa
        "Z": 1e21,       # zetta
        "Y": 1e24,       # yotta
    }

    # SI prefixes for output (only powers of 3 as well as the empty prefix),
    SI_OUT_PREF = [
        (1e-24, "y"),
        (1e-21, "z"),
        (1e-18, "a"),
        (1e-15, "f"),
        (1e-12, "p"),
        ( 1e-9, "n"),
        ( 1e-6, "\u03BC"),
        ( 1e-3, "m"),
        (    1, ""),
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

    # text formatting mode ("plain", "unicode", or "html")
    output_mode = "unicode"

    # enable CGS fractional units
    fractional_units = True

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
    def use(cls, unit, uname, prefix=False):
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
        cls.preferred_units[unit.units] = (unit, uname, prefix)

    @classmethod
    def use_system(cls, sys):
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
                cls.use(u[1], u[0], True)
        elif sys == "US_base":
            # only the base units
            cls.base_units = unit_systems["us"]
        elif sys == "US":
            # set base units
            cls.base_units = unit_systems["us"]
            # add US units from a separate list of deault units
            for u in default_us_u:
                cls.use(u[1], u[0])
        elif sys == "CGS_base":
            # only base units
            cls.base_units = unit_systems["cgs"]
        elif sys == "CGS" or sys == "CGS_Gauss":
            cls.base_units = unit_systems["cgs"]
            for u in default_cgs_gauss_u:
                cls.use(u[1], u[0])
        elif sys == "CGS_ESU":
            cls.base_units = unit_systems["cgs"]
            for u in default_cgs_esu_u:
                cls.use(u[1], u[0])
        elif sys == "CGS_EMU":
            cls.base_units = unit_systems["cgs"]
            for u in default_cgs_emu_u:
                cls.use(u[1], u[0])
        else:
            raise TypeError("Unknown unit system", sys)

    @classmethod
    def set_precision(cls, prec):
        '''
        specify the precision (significant digits) in an output

        Parameters:
            prec : integer

        Returns:
            None
        '''
        cls.format_str = "{:." + str(prec) + "g}"

    @classmethod
    def make_super(cls, s):
        '''
        convert well formatted unit string with exponents into target format
        (plain text, unicode, or html)

        Parameters:
            s : unit string

        Returns:
            prettified unit string
        '''
        def plainify(match):
            s = match.group(0)
            return "^" + "".join(str(superscript_ints.get(c, c)) for c in s)

        # convert to plain first
        # remove unicode
        s = re.sub(r'(\u207a|\u207b)?((\u2070|\xb9|\xb2|\xb3|\u2074|\u2075|\u2076|\u2077|\u2078|\u2079)+)', plainify, s)
        # remove HTML tags
        s = re.sub(r'<[^>]*>', '', s)

        if cls.output_mode == "html":
            return re.sub(r"\^([+-]?\d+(.5)?)", r"<sup>\1</sup>", s)
        elif (cls.output_mode == "unicode" and
              not ConcreteNumber.fractional_units):
            return "".join(superscript_chars.get(c, c) for c in s)
        else:
            return s

    @classmethod
    def analyze(cls, var, target):
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
                if ia > maxval:
                    maxval = ia
                td += ia
                lex = lex / 10.0 + ia + (0 if i >= 0 else 0.000001)
            return (maxval, td, lex)

        def mk_indices(num_var):
            idx = [tuple()]
            for i in range(0, num_var):
                idx = [j + (k,) for k in range(-3, 4) for j in idx]
            list.sort(idx, key=tuple_rank)
            return idx

        def is_multiple(exp1, exp2):
            mult = 0
            for e1, e2 in zip(exp1, exp2):
                if e2 != 0:
                    mult = e1 // e2
            for e1, e2 in zip(exp1, exp2):
                if e1 != mult * e2:
                    return False
            return True

        idx_list = mk_indices(len(var))
        sol_exps = []
        for idx in idx_list:
            soln = ConcreteNumber()
            for i, v in zip(idx, var):
                soln *= v**i
            if soln.units == target.units:
                redundant = False
                for s in sol_exps:
                    redundant |= is_multiple(idx, s)
                if not redundant:
                    sol_exps.append(idx)
        # convert exponent list to the actual solution format
        return [[(v, e) for v, e in zip(var, exp)] for exp in sol_exps]

    def __init__(self, val=1, units=(0, 0, 0, 0, 0, 0, 0, 0)):
        '''
        create new ConcreteNumber object

        Parameters:
            val :       numerical value
            units :     - EITHER an 8-tuple specifying the exponents for the
                          8 SI base units  (kg,m,s,A,K,mol,cd,rad),
                        - OR a unit string to be parsed
        '''
        if type(units) == str:
            u = dacalc.unitparser.parse(units)
            self.value = val * u.value
            self.units = u.units
        else:
            self.value = val
            self.units = units

    def dimensionstr(self):
        '''
        the dimensions shown in a unit-independent fashion

        Returns:
            string
        '''
        strs = ("mass", "length", "time", "current", "temperature",
                "substance", "luminance")
        dims = ""
        for d, e in zip(strs, self.units):
            if e != 0:
                ei = e / 2
                dims += ' ' + d + '^' + str(int(ei) if ei.is_integer() else ei)
        return ConcreteNumber.make_super(dims[1:])

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
        for u, fexp in zip(ConcreteNumber.base_units, self.units[:-1]):
            h = fexp / 2
            exp = int(h) if h.is_integer() else h
            if exp != 0:
                s += u[0]
                if exp != 1:
                    s += "^" + str(exp)
                s += ' '
        # deal with radians and steradians
        rad = int(self.units[-1] / 2)
        if rad != 0:
            if rad == 2:
                s += "sr"
            elif rad == -2:
                s += "sr^-1"
            else:
                s += "rad" + (("^" + str(rad)) if rad != 1 else "")
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

    def check_dim_match(self, other, warnrad=True):
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
                  file=sys.stderr)

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
        new_u = ConcreteNumber(1)  # start dimensionless
        for i in range(0, 7):
            ef = self.units[i] / 2
            if ef.is_integer():
                e = int(ef)
            else:
                if ConcreteNumber.fractional_units:
                    e = ef
                else:
                    print("Fractional units detected --",
                          "please enable fractional unit mode",
                          file=sys.stderr)
                    raise(TypeError)
            new_u *= ConcreteNumber.base_units[i][1] ** e
        return self.value / new_u.value

    def __str__(self, display_units=None):
        '''
        string reprsentation of the ConcreteNumber

        If display_units is provided, it represents a unit conversion
        request. If this conversion is not possible due to a dimension
        mismatch, then two unit strings are produced, e.g.:

            ConcreteNumber(1,"m/s").__str__("m")

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
        if display_units is not None:
            du = ConcreteNumber(1, display_units)
            return str(self / du) + ' [' + \
                ConcreteNumber.make_super(display_units) + ']'

        # dimensionless variables are just the numerical values
        if self.units == ConcreteNumber.DIMENSION_LESS:
            return ConcreteNumber.format_str.format(self.value)

        # for everything else we first check if one of the preferred
        # units matches
        try:
            pref_unit = ConcreteNumber.preferred_units[self.units]
            converted = self / pref_unit[0]
            if ConcreteNumber.use_si_prefix and pref_unit[2]:
                # calculate the SI prefix
                (mult, pref) = ConcreteNumber.SI_OUT_PREF[-1]
                for (m, p) in ConcreteNumber.SI_OUT_PREF:
                    if converted.value < 900 * m:
                        (mult, pref) = (m, p)
                        break
                return(ConcreteNumber.format_str.format(converted.value / mult)
                       + " [" + pref + pref_unit[1] + "]")
            else:
                return (ConcreteNumber.format_str.format(converted.value)
                        + " [" + pref_unit[1] + "]")
        except KeyError:
            # if there are no preferred units, we just use the
            # standard unit string conversion
            return (ConcreteNumber.format_str.format(self.convert_to_base())
                    + " [" + self.unitstr() + "]")

    def __repr__(self):
        '''
        string representation
        '''
        return self.__str__()

    #
    # arithmetic ops
    #

    def __add__(self, other):
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
        return ConcreteNumber(self.value + other.value, self.units)

    def __sub__(self, other):
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
        return ConcreteNumber(self.value - other.value, self.units)

    def __neg__(self):
        '''
        unary minus operator

        Parameters:
            None

        Returns:
            negative of self
        '''
        return ConcreteNumber(-self.value, self.units)

    def __mul__(self, other):
        '''
        multiplication operator

        Parameters:
            other : a second ConcreteNumber, an int or a float

        Returns:
            product as a ConcreteNumber
        '''
        if isinstance(other, ConcreteNumber):
            return ConcreteNumber(self.value * other.value,
                                  tuple(u1 + u2 for (u1, u2) in
                                        zip(self.units, other.units)))
        elif type(other) == float or type(other) == int:
            return ConcreteNumber(self.value * other, self.units)
        else:
            return NotImplemented

    def __rmul__(self, other):
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
            return ConcreteNumber(self.value * other.value,
                                  tuple(u1 + u2 for (u1, u2) in
                                        zip(self.units, other.units)))
        elif type(other) == float or type(other) == int:
            return ConcreteNumber(self.value * other, self.units)
        else:
            return NotImplemented

    def __truediv__(self, other):
        '''
        division operator

        Parameters:
            other : a second ConcreteNumber, an int or a float

        Returns:
            division as a ConcreteNumber
        '''
        if isinstance(other, ConcreteNumber):
            return ConcreteNumber(self.value / other.value,
                                  tuple(u1 - u2 for (u1, u2) in
                                        zip(self.units, other.units)))
        elif type(other) == float or type(other) == int:
            return ConcreteNumber(self.value / other, self.units)
        else:
            return NotImplemented

    def __rtruediv__(self, other):
        '''
        right-sided division operator

        Parameters:
            other : an int or a float

        Returns:
            product as a ConcreteNumber
        '''
        if type(other) == float or type(other) == int:
            return ConcreteNumber(val=other) / self
        else:
            return NotImplemented

    def __pow__(self, exp):
        '''
        power operator

        Parameters:
            exp : integer exponent

        Returns:
            self ^ exp
        '''
        err_msg = "Illegal exponent " + str(exp) + " for dimension " + \
            '[' + ConcreteNumber.dimensionstr(self) + ']'

        if type(exp) == int or (type(exp) == float and exp.is_integer()):
            # integer exponents
            return ConcreteNumber(self.value**int(exp),
                                  tuple(u * exp for u in self.units))
        elif type(exp) == float:
            # check that we have exactly half dimensions
            if not (ConcreteNumber.fractional_units and (exp * 2).is_integer()):
                raise(TypeError(err_msg))
            dim_t = tuple(int(u * exp) if (u * exp).is_integer()
                          else u * exp for u in self.units)
            # even with fractional units enabled, only half dimensions
            # are allowed, and only for mass and length
            if type(dim_t[0]) != int or type(dim_t[1]) != int:
                raise TypeError(err_msg)
            for u in dim_t[2:]:
                if type(u) != int or u % 2 != 0:
                    raise TypeError(err_msg)

            return ConcreteNumber(self.value**exp, dim_t)
        else:
            return NotImplemented


#
# math functions using ConcreteNumber
#

def root(val, exp):
    '''
    The <exp>-th root

    Parameters:
        val:    a concrete value, where the units must be powers of <exp>
        exp:    integer exponent

    Returns:
        the <exp>-th root of <val>
    '''
    err_msg = "Cannot take the " + str(exp) + "-root of [" \
        + str(val.unitstr()) + "]"
    if type(exp) == int:
        # tentative dimensions of result, subject to sanity tests
        dim_t = tuple(int(u / exp) for u in val.units)

        if ConcreteNumber.fractional_units:
            # even with fractional units enabled, only half dimensions
            # are allowed, and only for mass and length
            if type(dim_t[0]) != int or type(dim_t[1]) != int:
                raise TypeError(err_msg)
            for u in dim_t[2:]:
                if u % 2 != 0:
                    raise TypeError(err_msg)
        else:
            # no fractional units -> all result must be even in fixed point
            for u in dim_t:
                if type(u) != int or u % 2 != 0:
                    raise TypeError(err_msg)
        return ConcreteNumber(math.pow(val.value, 1.0 / exp), dim_t)
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
    return root(val, 2)


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
    if angle.units == ConcreteNumber.SR:  # steradians become radians
        return ConcreteNumber(val, ConcreteNumber.RAD)
    else:  # everything else becomes dimensionless
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
    if angle.units == ConcreteNumber.SR:  # steradians become radians
        return ConcreteNumber(val, ConcreteNumber.RAD)
    else:  # everything else becomes dimensionless
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
    if angle.units == ConcreteNumber.SR:  # steradians become radians
        return ConcreteNumber(val, ConcreteNumber.RAD)
    else:  # everything else becomes dimensionless
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
    if val.units == ConcreteNumber.RAD:  # radians become steradians
        return ConcreteNumber(angle, ConcreteNumber.SR)
    else:  # everything else becomes radians
        return ConcreteNumber(angle, ConcreteNumber.RAD)


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
    if val.units == ConcreteNumber.RAD:  # radians become steradians
        return ConcreteNumber(angle, ConcreteNumber.SR)
    else:  # everything else becomes radians
        return ConcreteNumber(angle, ConcreteNumber.RAD)


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
    if val.units == ConcreteNumber.RAD:  # radians become steradians
        return ConcreteNumber(angle, ConcreteNumber.SR)
    else:  # everything else becomes radians
        return ConcreteNumber(angle, ConcreteNumber.RAD)


def atan2(val1, val2):
    '''
    Arc Tangens with two arguments

    Parameters:
        val1:   a concrete number
        val2:   a concrete number with the same units as <val1>

    Returns:
        the arc tangens of <val1>/<val2> in units of [rad]
    '''
    val1.check_dim_match(val2, warnrad=False)
    angle = math.atan2(val1.value, val2.value)
    if val1.units == ConcreteNumber.RAD:  # radians become steradians
        return ConcreteNumber(angle, ConcreteNumber.SR)
    else:  # everything else becomes radians
        return ConcreteNumber(angle, ConcreteNumber.RAD)


def fabs(val):
    '''
    Absolute value

    Parameters:
        val:    a concrete number

    Returns:
        the absolute value of <val>, with the same units as <val>
    '''
    return ConcreteNumber(val=math.fabs(val.value), unit_tuple=val.units)


def log(val, base):
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
        return ConcreteNumber(val=math.log(val.value, base))
    else:
        return ConcreteNumber(val=math.log(val.value, base.value))


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


def pow(val, exp):
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
        if ConcreteNumber.fractional_units and (exp * 2).is_integer():
            return val**exp
        else:
            val.check_dimensionless()
            return ConcreteNumber(val=math.pow(val.value, exp))
    elif isinstance(exp, ConcreteNumber):
        val.check_dimensionless()
        exp.check_dimensionless()
        return ConcreteNumber(val=math.pow(val.value, exp.value))
    else:
        return NotImplemented


# add a single unit to the unit dictionary (check for name clashes)
def add_unit_to_dict(name, newunit):
    if name in ConcreteNumber.units.keys():
        print("Warning: unit", name,
              " already exists (overwriting by new definition).",
              file=sys.stderr)
    ConcreteNumber.units[name] = newunit


# add a new unit both to a unit group list and the global unit
# dictionary; define all SI prefixes if needed
def new_unit(unitgroup, unit, si_prefixes=False):
    unitgroup.append(unit)
    add_unit_to_dict(unit[0], unit[1])
    if si_prefixes:
        for prefix, multiplier in ConcreteNumber.SI_PREFIXES.items():
            add_unit_to_dict(prefix + unit[0], multiplier * unit[1])


# the folllwing are lists of unit groups, which consist of tuples
# (name string, ConcreteNumber representation, description string)

# add units of a dimensionless variable
add_unit_to_dict("", ConcreteNumber(1, ConcreteNumber.DIMENSION_LESS))

# with the class defined we can do a circular import of the
# unitparser, so that defining the units becomes a bit less verbose
import dacalc.unitparser    # noqa: E402

#
# list of all )base+derived) SI units (without prefix)
#
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
new_unit(si_u,("sr",  ConcreteNumber(1,"rad^2"),      "steradian"))
new_unit(si_u,("N",   ConcreteNumber(1,"kg m / s^2"), "Newton"), True)
new_unit(si_u,("Pa",  ConcreteNumber(1,"N / m^2"),    "Pascal"), True)
new_unit(si_u,("J",   ConcreteNumber(1,"N m"),        "Joule"), True)
new_unit(si_u,("W",   ConcreteNumber(1,"J / s"),      "Watt"), True)
new_unit(si_u,("C",   ConcreteNumber(1,"s A"),        "Coulomb"), True)
new_unit(si_u,("V",   ConcreteNumber(1,"W / A"),      "Volt"), True)
new_unit(si_u,("F",   ConcreteNumber(1,"C / V"),      "Farad"), True)
new_unit(si_u,("Ohm", ConcreteNumber(1,"V / A"),      "Ohm"), True)
new_unit(si_u,("S",   ConcreteNumber(1,"Ohm^-1"),     "Siemens"), True)
new_unit(si_u,("Wb",  ConcreteNumber(1,"V s"),        "Weber"), True)
new_unit(si_u,("T",   ConcreteNumber(1,"Wb / m^2"),   "Tesla"), True)
new_unit(si_u,("H",   ConcreteNumber(1,"Wb / A"),     "Henry"), True)
new_unit(si_u,("lm",  ConcreteNumber(1,"cd sr"),      "Lumen"), True)
new_unit(si_u,("lx",  ConcreteNumber(1,"cd / m^2"),   "Lux"), True)
new_unit(si_u,("Bq",  ConcreteNumber(1,"1 / s"),      "Bequerel"), True)
#  after Bequerel, so that Hertz is default unit for 1/s
new_unit(si_u,("Hz",  ConcreteNumber(1,"1 / s"),      "Hertz"), True)
new_unit(si_u,("Gy",  ConcreteNumber(1,"J / kg"),     "Gray"), True)
new_unit(si_u,("Sv",  ConcreteNumber(1,"Gy"),         "Sievert"), True)
new_unit(si_u,("kat", ConcreteNumber(1,"mol / s"),    "Katal"), True)

#
# miscellaneaous units (including some outdated/dperecated ones)
#
misc_u = []
new_unit(misc_u,("min",    ConcreteNumber(60,"s"),             "minute"))
new_unit(misc_u,("h",      ConcreteNumber(60,"min"),           "hour"))
new_unit(misc_u,("d",      ConcreteNumber(24,"h"),             "day"))
new_unit(misc_u,("deg",    ConcreteNumber(math.pi/180.0,"rad"), "degree"))
new_unit(misc_u,("arcmin", ConcreteNumber(1 / 60.0,"deg"),     "arc minute"))
new_unit(misc_u,("arcsec", ConcreteNumber(1 / 60.0,"arcmin"),  "arc second"))
new_unit(misc_u,("NM",     ConcreteNumber(1852,"m"),           "nautical mile"))
new_unit(misc_u,("knot",   ConcreteNumber(1,"NM / h"),         "knot"))
new_unit(misc_u,("ha",     ConcreteNumber(1,"hm^2"),           "hectare"))
new_unit(misc_u,("L",      ConcreteNumber(1,"dm^3"),           "litre"), True)
new_unit(misc_u,("t",      ConcreteNumber(1e3,"kg"),           "tonne"))
new_unit(misc_u,("Da" ,    ConcreteNumber(1.66053904020e-27,"kg"), "Dalton"))
new_unit(misc_u,("rpm" ,   ConcreteNumber(2*math.pi,"rad/min"),"rpm"))
new_unit(misc_u,("eV",     ConcreteNumber(1.602176634e-19,"J"),
                 "electron volt"),True)
new_unit(misc_u,("cal",    ConcreteNumber(4.184,"J"),          "calorie"))
new_unit(misc_u,("kcal",   ConcreteNumber(4184,"J"),           "kilocalorie"))
new_unit(misc_u,("torr",   ConcreteNumber(133.3224,"Pa"),      "Torr"))
new_unit(misc_u,("atm",    ConcreteNumber(101325,"Pa"),        "atmosphere"))
new_unit(misc_u,("bar",    ConcreteNumber(1e5,"Pa"),           "Bar"),True)
new_unit(misc_u,("nt",     ConcreteNumber(1,"cd / m^2"),       "nit"))
new_unit(misc_u,("ct",     ConcreteNumber(.2,"g"),             "carat"))
new_unit(misc_u,("degC",   ConcreteNumber(1,"K"),              "deg. Celsius"))
new_unit(misc_u,("degF",   ConcreteNumber(5.0 / 9.0,"K"),     "deg. Farenheit"))

#
# list of US/imperial units (including some older ones)
# Note: different for the UK/imperial ones especially on the volumes!
#
us_u = []
# US length
new_unit(us_u,("in",   ConcreteNumber(25.4,"mm"),            "inch"))
new_unit(us_u,("ft",   ConcreteNumber(12,"in"),              "foot"))
new_unit(us_u,("yd",   ConcreteNumber(26,"in"),              "yard"))
new_unit(us_u,("mi",   ConcreteNumber(63360,"in"),           "mile"))
new_unit(us_u,("mil",  ConcreteNumber(1e-3,"in"),            "mil"))
new_unit(us_u,("thou", ConcreteNumber(1,"mil"),              "thou"))
new_unit(us_u,("uin",  ConcreteNumber(1e-3,"mil"),           "microinch"))
# US area
new_unit(us_u,('acre', ConcreteNumber(43560*(1200.0/3937)**2,"m^2"),"acre"))
# US speed
new_unit(us_u,("mph",  ConcreteNumber(1,"mi / h"),           "mph"))
# US weight
new_unit(us_u,("lb",   ConcreteNumber(0.45359237,"kg"),      "pound"))
new_unit(us_u,("ton",  ConcreteNumber(2000,"lb"),            "ton"))
new_unit(us_u,("oz",   ConcreteNumber(1 / 16,"lb"),          "ounze"))
new_unit(us_u,("gr",   ConcreteNumber(1 / 7000,"lb"),        "grain"))
new_unit(us_u,("ozt",  ConcreteNumber(480,"gr"),             "troy_ounce"))
new_unit(us_u,("lbt",  ConcreteNumber(12,"ozt"),             "troy_pound"))
# US force
new_unit(us_u,("lbf",  ConcreteNumber(9.80665,"lb m / s^2"), "pound_force"))
# US volume (fluid -- didn't bother with dry volumes)
new_unit(us_u,("gal",  ConcreteNumber(3.785411784,"L"),       "gallon"))
new_unit(us_u,("qt",   ConcreteNumber(.25,"gal"),             "quart"))
new_unit(us_u,("pt",   ConcreteNumber(.5,"qt"),               "pint"))
new_unit(us_u,("cp",   ConcreteNumber(.5,"pt"),               "cup"))
new_unit(us_u,("floz", ConcreteNumber(.125,"cp"),             "fluid_ounce"))
new_unit(us_u,("Tbsp", ConcreteNumber(.5,"floz"),             "tablespoon"))
new_unit(us_u,("tsp",  ConcreteNumber(1 / 3.0,"Tbsp"),        "teaspoon"))
# US miscellaneous or outdated
new_unit(us_u,("slug", ConcreteNumber(1,"lbf s^2 / ft"),      "slug"))
new_unit(us_u,("pdl",  ConcreteNumber(1,"lb ft / s^2"),       "poundal"))
new_unit(us_u,("psi",  ConcreteNumber(1,"lbf / in^2"),        "PSI"))
new_unit(us_u,("fc",   ConcreteNumber(1,"lm / ft^2"),         "foot_candle"))
new_unit(us_u,("lbmol",ConcreteNumber(453.59237,"mol"),       "pound_mole"))
new_unit(us_u,("degR", ConcreteNumber(5.0 / 9.0,"K"),         "deg. Rankine"))

# because the US system has multiple units for each dimension, we need
# to pick some specific ones as default units
default_us_u = [
    (u, ConcreteNumber(1,u)) for u in ["lb", "ft", "s", "Hz", "A", "degR",
                                       "lbmol", "cd", "lbf", "ft lbf", "mph",
                                       "psi", "fc"]
]

#
# legacy CGS units including Gauss, EMU and ESU units
# (excluding SI units from above)
#
cgs_u = []
new_unit(cgs_u,("dyn",    ConcreteNumber(1e-5,"N"),          "Dyne"),True)
new_unit(cgs_u,("erg",    ConcreteNumber(1e-7,"J"),          "erg"),True)
new_unit(cgs_u,("statW",  ConcreteNumber(1,"erg/s"),         "stat Watt"),True)
new_unit(cgs_u,("Ba",     ConcreteNumber(1,"dyn/cm^2"),      "barye"),True)
new_unit(cgs_u,("Gal",    ConcreteNumber(1,"cm/s^2"),        "Galileo"),True)
new_unit(cgs_u,("P",      ConcreteNumber(1,"g/cm s"),        "poise"),True)
new_unit(cgs_u,("St",     ConcreteNumber(1,"cm^2/s"),        "Stokes"),True)

new_unit(cgs_u,("Fr",     ConcreteNumber(math.sqrt(1e-9),
                                         (1,3,-2,0,0,0,0,0)),"Franklin"),
         True)
new_unit(cgs_u,("statC",  ConcreteNumber(1,"Fr"),            "stat Coulomb"),
         True)
new_unit(cgs_u,("abC",    ConcreteNumber(math.sqrt(1e-5),
                                         (1,1,0,0,0,0,0,0)),"ab Coulomb"),
         True)
new_unit(cgs_u,("statA",  ConcreteNumber(1,"Fr/s"),          "stat Ampere"),
         True)
new_unit(cgs_u,("abA",    ConcreteNumber(1,"abC/s"),         "ab Ampere"),True)
new_unit(cgs_u,("Bi",     ConcreteNumber(1,"abA"),           "Biot"),True)
new_unit(cgs_u,("statV",  ConcreteNumber(1,"erg/Fr"),        "stat Volt"),True)
new_unit(cgs_u,("abV",    ConcreteNumber(1,"abA cm/s"),      "ab Volt"),True)
new_unit(cgs_u,("statOhm",ConcreteNumber(1,"statV/statA"),   "stat Ohm"),True)
new_unit(cgs_u,("abOhm",  ConcreteNumber(1,"abV/abA"),       "ab Ohm"),True)
new_unit(cgs_u,("statF",  ConcreteNumber(1,"Fr/statV"),      "stat Faraday"),
         True)
new_unit(cgs_u,("abF",    ConcreteNumber(1,"abC/abV"),       "ab Faraday"),
         True)
new_unit(cgs_u,("statWb", ConcreteNumber(1,"statV s"),       "stat Weber"),True)
new_unit(cgs_u,("statT",  ConcreteNumber(1,"statWb/cm^2"),   "stat Tesla"),True)
new_unit(cgs_u,("Mx",     ConcreteNumber(1,"statV cm"),      "Maxwell"),True)
new_unit(cgs_u,("G",      ConcreteNumber(1,"Mx/cm^2"),       "Gauss"),True)
new_unit(cgs_u,("Oe",     ConcreteNumber(1,"G"),             "Oersted"),True)
new_unit(cgs_u,("pole",   ConcreteNumber(1,"dyn/Oe"),        "unit pole"),True)
new_unit(cgs_u,("statH",  ConcreteNumber(1,"erg/statA^2"),   "stat Henry"),True)
new_unit(cgs_u, ("abH",   ConcreteNumber(1,"cm"),            "ab Henry"),True)
new_unit(cgs_u,("D",      ConcreteNumber(1e-18,"statC cm"),  "Debye"),True)

new_unit(cgs_u,("ph",    ConcreteNumber(1,"lm/cm^2"),        "phot"))
new_unit(cgs_u,("sb",    ConcreteNumber(1,"cd/cm^2"),        "stilb"),True)
new_unit(cgs_u,("Lb",    ConcreteNumber(1/math.pi,"cd/cm^2"),"Lambert"),True)

# default units for CGS-Gauss
default_cgs_gauss_u = [
    (u, ConcreteNumber(1,u)) for u in ["g", "cm", "s", "Hz", "Gal", "dyn",
                                       "erg", "erg/s", "Ba", "P", "St",
                                       "ph", "Lb",
                                       "Fr", "Fr/s", "statV", "statV/cm",
                                       "Fr/cm^2", "G","erg/G"]
]

# default units for CGS-ESU
default_cgs_esu_u = [
    (u, ConcreteNumber(1,u)) for u in ["g", "cm", "s", "Hz", "Gal", "dyn",
                                       "erg", "erg/s", "Ba", "P", "St",
                                       "ph", "Lb",
                                       "statC", "statA", "statV", "statV/cm",
                                       "statC/cm^2", "statT", "statC cm^2",
                                       "statA/cm", "statWb"]
]

# default units for CGS-EMU
default_cgs_emu_u = [
    (u, ConcreteNumber(1,u)) for u in ["g", "cm", "s", "Hz", "Gal", "dyn",
                                       "erg", "erg/s", "Ba", "P", "St",
                                       "ph", "Lb",
                                       "abC", "Bi", "abV", "abV/cm",
                                       "abC/cm^2", "abC cm", "G", "Bi cm^2",
                                       "Oe", "Mx", "abOhm", "abOhm cm", "abF"]
]

# base units for different unit systems
# (list of the seven base units for each of: mass, length, time,
# current, temperature, substance amount, luminous intensity)
# the order of these entries matters, and should align with the dimensions

unit_systems = {
    "si" : [
        ("kg",   ConcreteNumber(1,"kg"),     "kilogram"),
        ("m",    ConcreteNumber(1,"m"),      "metre"),
        ("s",    ConcreteNumber(1,"s"),      "second"),
        ("A",    ConcreteNumber(1,"A"),      "Ampere"),
        ("K",    ConcreteNumber(1,"K"),      "Kelvin"),
        ("mol",  ConcreteNumber(1,"mol"),    "mole"),
        ("cd",   ConcreteNumber(1,"cd"),     "candela")
    ],

    "us" : [
        ("lb",    ConcreteNumber(1,"lb"),    "pound"),
        ("ft",    ConcreteNumber(1,"ft"),    "foot"),
        ("s",     ConcreteNumber(1,"s"),     "second"),
        ("A",     ConcreteNumber(1,"A"),     "Ampere"),
        ("degR",  ConcreteNumber(1,"degR"),  "deg. Rankine"),
        ("lbmol", ConcreteNumber(1,"lbmol"), "pound-mole"),
        ("cd",    ConcreteNumber(1,"cd"),    "candela")
    ],

    # although the CGS system has alternative units for current, these
    # are dimensionally different. Also the Ampere is still used in many
    # practical CGS settings
    "cgs" : [
        ("g",    ConcreteNumber(1,"g"),      "gram"),
        ("cm",   ConcreteNumber(1,"cm"),     "centimetre"),
        ("s",    ConcreteNumber(1,"s"),      "second"),
        ("A",    ConcreteNumber(1,"A"),      "Ampere"),
        ("K",    ConcreteNumber(1,"K"),      "Kelvin"),
        ("mol",  ConcreteNumber(1,"mol"),    "mole"),
        ("cd",   ConcreteNumber(1,"cd"),     "candela")
    ]
}


# initialize the unit system to SI
ConcreteNumber.use_system("SI")

#
# some mathematical and physical constants
# (this is the version that contains help text)
#
const = [
    ("pi",   ConcreteNumber(math.pi,                ""),      "\u03C0"),
    ("e",    ConcreteNumber(math.e,                 ""),      "Euler c."),
    ("zeroC",ConcreteNumber(273.15,                 "K"),     "0 [degC] abs."),
    ("zeroF",ConcreteNumber(255.372,                "K"),     "0 [degF] abs."),
    ("c",    ConcreteNumber(299792458,              "m / s"), "speed of light"),
    ("G",    ConcreteNumber(6.67408e-11,     "m^3 / kg s^2"), "gravit. c."),
    ("g",    ConcreteNumber(9.80665,                "m/s^2"), "std. gravity"),
    ("h",    ConcreteNumber(6.62607004e-34,         "J s") ,  "Planck c."),
    ("ec",   ConcreteNumber(1.602176634e-19,        "C"),     "elem. charge"),
    ("alpha",ConcreteNumber(0.0072973525693,        ""),      "fine struct.c."),
    ("m_e",  ConcreteNumber(9.1093837015e-31,       "kg"),    "electron mass"),
    ("m_p",  ConcreteNumber(1.67262192369e-27,      "kg"),    "proton mass"),
    ("N_A",  ConcreteNumber(6.02214076e23,          "1/mol"), "Avogadro c."),
    ("F",    ConcreteNumber(96485.3321233100184,    "C/mol"), "Faraday c."),
    ("R",    ConcreteNumber(8.31446261815324,     "J/K mol"), "univ. gas c."),
    ("k",    ConcreteNumber(1.380649e-23,           "J / K"), "Boltzmann c."),
    ("mu_0", ConcreteNumber(1.25663706212e-6,       "H / m"), "vac. permeab."),
    ("eps_0",ConcreteNumber(8.8541878128e-12,       "F / m"), "vac. permit."),
    ("PI",   ConcreteNumber(2.0678338141336296e-15, "Wb"),    "mag. flux quant")
]

# add pre-defined constants into class dictionary for easier lookup
# and more compact python API syntax
for (c,v,_) in const:
    ConcreteNumber.const[c] = v
