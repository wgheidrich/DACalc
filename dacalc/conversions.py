from dacalc.concretenumber import ConcreteNumber as CN
import dacalc.concretenumber as cn
from dacalc.unitparser import parse as uparse
import math as m
import sys


c = CN.const["c"] * CN(1e-10, "s/cm")
sq_eps4pi =  cn.sqrt(4 * m.pi * CN.const["eps_0"])
sq_mu4pi =   cn.sqrt(4 * m.pi * CN.const["mu_0"])
sq_4pi_eps = cn.sqrt(4 * m.pi / CN.const["eps_0"])
sq_4pi_mu =  cn.sqrt(4 * m.pi / CN.const["mu_0"])

# SI units for electromagnetism quantities
Q_SI =     CN(1, "C")
phi_e_SI = CN(1, "V m")
I_SI =     CN(1, "A")
V_SI =     CN(1, "V")
E_SI =     CN(1, "V/m")
D_SI =     CN(1, "C/m^2")
p_SI =     CN(1, "C m")
m_SI =     CN(1, "A m^2")
B_SI =     CN(1, "T")
H_SI =     CN(1, "A/m")
phi_m_SI = CN(1, "Wb")
R_SI =     CN(1, "Ohm")
rho_SI =   CN(1, "Ohm m")
C_SI =     CN(1, "F")
L_SI =     CN(1, "H")


CGSconversionFactors = {
    "Gauss": {
        Q_SI.units:           c * CN(1e-1,        "Fr/C"),
        phi_e_SI.units:       c * CN(m.pi * 4e-1, "Fr/V m"),
        I_SI.units:           c * CN(1e-1,        "Fr/A s"),
        V_SI.units:       1 / c * CN(1e8,         "statV/V"),
        E_SI.units:       1 / c * CN(1e6,         "statV m/V cm"),
        D_SI.units:           c * CN(1e-5,        "Fr m^2/C cm^2"),
        p_SI.units:           c * CN(1e19,        "D/C m"),
        m_SI.units:               CN(1e3,         "erg/G A m^2"),
        B_SI.units:               CN(1e4,         "G/T"),
        H_SI.units:               CN(m.pi * 4e-3, "Oe m/A"),
        phi_m_SI.units:           CN(1e8,         "Mx/Wb"),
        R_SI.units:   1 / c / c * CN(1e9,         "s/Ohm cm"),
        rho_SI.units: 1 / c / c * CN(1e11,        "s/Ohm m"),
        C_SI.units:       c * c * CN(1e-9,        "cm/F"),
        L_SI.units:   1 / c / c * CN(1e9,         "s^2/cm H")
    },

    "ESU": {
        Q_SI.units:           c * CN(1e-1,        "statC/C"),
        phi_e_SI.units:       c * CN(m.pi * 4e-1, "statC/V m"),
        I_SI.units:           c * CN(1e-1,        "statA/A"),
        V_SI.units:       1 / c * CN(1e8,         "statV/V"),
        E_SI.units:       1 / c * CN(1e6,         "statV m/V cm"),
        D_SI.units:           c * CN(1e-5,        "statC m^2/C cm^2"),
        p_SI.units:           c * CN(10,          "statC cm /C m"),
        m_SI.units:           c * CN(1e3,         "statC cm^2/A m^2"),
        B_SI.units:       1 / c * CN(1e4,         "statT/T"),
        H_SI.units:           c * CN(m.pi * 4e-3, "statA m/A cm"),
        phi_m_SI.units:   1 / c * CN(1e8,         "statWb/Wb"),
        R_SI.units:   1 / c / c * CN(1e9,         "s/Ohm cm"),
        rho_SI.units: 1 / c / c * CN(1e11,        "s/Ohm m"),
        C_SI.units:       c * c * CN(1e-9,        "cm/F"),
        L_SI.units:   1 / c / c * CN(1e9,         "s^2/cm H")
    },

    "EMU": {
        Q_SI.units:     CN(1e-1,        "abC/C"),
        phi_e_SI.units: CN(1e-1,        "abC/V m"),
        I_SI.units:     CN(1e-1,        "Bi/A"),
        V_SI.units:     CN(1e8,         "abV/V"),
        E_SI.units:     CN(1e6,         "abV m/V cm"),
        D_SI.units:     CN(1e-5,        "abC m^2/C cm^2"),
        p_SI.units:     CN(10,          "abC cm/C m"),
        m_SI.units:     CN(1e3,         "Bi cm^2/A m^2"),
        B_SI.units:     CN(1e4,         "G/T"),
        H_SI.units:     CN(m.pi * 4e-3, "Oe m/A"),
        phi_m_SI.units: CN(1e8,         "Mx/Wb"),
        R_SI.units:     CN(1e9,         "abOhm/Ohm"),
        rho_SI.units:   CN(1e11,        "abOhm cm/Ohm m"),
        C_SI.units:     CN(1e-9,        "abF/F"),
        L_SI.units:     CN(1e9,         "abH/H")
    }
}


def fromSI(cn, target_sys):
    '''
    convert from SI system to one of the CGS subsystems

    If the units are not related to electromagnetic quanitties, the input and
    output are the same; otherwise a conversion will be attmepted, which
    changes the internal representation (dimension) of the quantity.

    Parameters:
        cn :         a concrete number in SI units

        target_sys : the system to convert to (one of "Gauss", "ESU", or "EMU")

    Returns:
        concrete number that is valid in the chosen target CGS subsystem
    '''
    if cn.units[3] == 0:
        # no Amperes in the unit -- there fore we are not dealing
        # with an electromagnetic qantity, and the input is already
        # compatible with CGS
        return cn
    else:
        try:
            return CGSconversionFactors[target_sys][cn.units] * cn
        except LookupError:
            print("Don't know how to convert", cn, "from SI to", target_sys,
                  file=sys.stderr)
        return cn


def toSI(cn, target, source_sys):
    '''
    convert from one of the CGS subsystems to SI

    If the units are not related to electromagnetic quanitties, the input and
    output are the same; otherwise a conversion will be attmepted, which
    changes the internal representation (dimension) of the quantity to use
    Amperes instead of half-dimensions.

    Since the CGS subsytems are dimnesionally ambiguous, an SI target dimension
    to convert to must be provided.

    Parameters:
        cn :         a concrete number in the source CGS subsystem

        target :     concrete number of the same dimesnions we want to
                     convert to. Either in the form of a concrete number
                     object, or in the form of a unit string

        source_sys : the CGS subsystem in which the source is expressed
                     (one of  "Gauss", "ESU", or "EMU")
    '''
    if type(target) == str:
        target = uparse(target)
    if cn.units == target.units:
        return cn
    else:
        try:
            result = cn / CGSconversionFactors[source_sys][target.units]
            if result.units == target.units:
                return result
        except LookupError:
            pass
        print("Don't know how to convert", cn, " from ", source_sys, "to SI",
              file=sys.stderr)
        return cn


def SI_to_Gauss(cn):
    '''
    Convert SI quantities to equivalent quantities in CGS Gaussian units

    This is a conversion of the internal representation of the quantity to
    different dimensions.

    If the dimensions of the provided concrete number is valid in CGS
    Gaussian units, (i.e if the value does not have dimensions of Ampere)
    the input will be returned unaltered. Only if the input dimensions contain
    Ampere, a conversion is attempted (typically this will result in an entity
    with half-dimensions for mass and lengh.)

    Parameters:
        cn : a concrete number in SI units

    Returns:
        concrete number that is valid in CGS Gauss units.
    '''
    return fromSI(cn, "Gauss")


def Gauss_to_SI(cn, target):
    '''
    Convert CGS Gaussian quantities to equivalent quantities in SI units

    This is a conversion of the internal representation of the quantity to
    different dimensions.

    Because there is a large amount of dimensional ambiguity in Gaussian
    units (e.g. between electric an magnetic field, time and resistivity
    etc.), the target dimensions must be provided in the form of a concrete
    number (the numerical value of this quantity will be ignored).

    If input and target dimensions already match, the representation
    in Gaussian and SI units is identical, and the input is returned
    unaltered. Otherwise, a conversion is attempted.

    Parameters:
        cn : a concrete number in CGS Gaussian units

    Returns:
        concrete number that is valid in SI units.
    '''
    return toSI(cn, target, "Gauss")


def SI_to_ESU(cn):
    '''
    Convert SI quantities to equivalent quantities in CGS electrostatic units

    This is a conversion of the internal representation of the quantity to
    different dimensions.

    If the dimensions of the provided concrete number is valid in CGS ESU
    units, (i.e if the value does not have dimensions of Ampere) the input
    will be returned unaltered. Only if the input dimensions contain Ampere,
    a conversion is attempted (typically this will result in an entity
    with half-dimensions for mass and lengh.)

    Parameters:
        cn : a concrete number in SI units

    Returns:
        concrete number that is valid in CGS ESU units.
    '''
    return fromSI(cn, "ESU")


def ESU_to_SI(cn, target):
    '''
    Convert CGS electrostatic quantities to equivalent quantities in SI units

    This is a conversion of the internal representation of the quantity to
    different dimensions.

    Because there is a large amount of dimensional ambiguity in ESU
    units, the target dimensions must be provided in the form of a concrete
    number (the numerical value of this quantity will be ignored).

    If input and target dimensions already match, the representation
    in ESU and SI units is identical, and the input is returned
    unaltered. Otherwise, a conversion is attempted.

    Parameters:
        cn : a concrete number in CGS ESU units

    Returns:
        concrete number that is valid in SI units.
    '''
    return toSI(cn, target, "ESU")


def SI_to_EMU(cn):
    '''
    Convert SI quantities to equivalent quantities in CGS electromagnetic units

    This is a conversion of the internal representation of the quantity to
    different dimensions.

    If the dimensions of the provided concrete number is valid in CGS EMU
    units, (i.e if the value does not have dimensions of Ampere) the input
    will be returned unaltered. Only if the input dimensions contain Ampere,
    a conversion is attempted (typically this will result in an entity
    with half-dimensions for mass and lengh.)

    Parameters:
        cn : a concrete number in SI units

    Returns:
        concrete number that is valid in CGS EMU units.
    '''
    return fromSI(cn, "EMU")


def EMU_to_SI(cn, target):
    '''
    Convert CGS electromagnetic quantities to equivalent quantities in SI units

    This is a conversion of the internal representation of the quantity to
    different dimensions.

    Because there is a large amount of dimensional ambiguity in EMU
    units, the target dimensions must be provided in the form of a concrete
    number (the numerical value of this quantity will be ignored).

    If input and target dimensions already match, the representation
    in EMU and SI units is identical, and the input is returned
    unaltered. Otherwise, a conversion is attempted.

    Parameters:
        cn : a concrete number in CGS EMU units

    Returns:
        concrete number that is valid in SI units.
    '''
    return toSI(cn, target, "EMU")
