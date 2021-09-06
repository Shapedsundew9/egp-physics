"""Manages endpoint type interactions.

Endpoint types are identified by a signed 16-bit value or
a fully qualified name.

Erasmus GP types have values < 0.
An invalid ep type has the value -32768
NoneType == 0
All other types values are > 0
"""


from json import load
from os.path import dirname, join
from enum import IntEnum
from logging import NullHandler, getLogger


_logger = getLogger(__name__)
_logger.addHandler(NullHandler())


# Load type data
with open(join(dirname(__file__), "data/ep_types.json"), "r") as file_ptr:
    ep_type_lookup = load(file_ptr)
    ep_type_lookup['v2n'] = {int(k): v for k, v in ep_type_lookup['v2n'].items()}
    ep_type_lookup['n2v'] = {k: int(v) for k, v in ep_type_lookup['n2v'].items()}
    ep_type_lookup['instanciation'] = {int(k): v for k, v in ep_type_lookup['instanciation'].items()}


INVALID_EP_TYPE_NAME = 'egp_invalid_type'
INVALID_EP_TYPE_VALUE = -32768
UNKNOWN_EP_TYPE_NAME = 'egp_unknown_type'
UNKNOWN_EP_TYPE_VALUE = -32767
ep_type_lookup['n2v'][INVALID_EP_TYPE_NAME] = INVALID_EP_TYPE_VALUE
ep_type_lookup['v2n'][INVALID_EP_TYPE_VALUE] = INVALID_EP_TYPE_NAME
ep_type_lookup['instanciation'][INVALID_EP_TYPE_VALUE] = [None] * 5
ep_type_lookup['n2v'][UNKNOWN_EP_TYPE_NAME] = UNKNOWN_EP_TYPE_VALUE
ep_type_lookup['v2n'][UNKNOWN_EP_TYPE_VALUE] = UNKNOWN_EP_TYPE_NAME
ep_type_lookup['instanciation'][UNKNOWN_EP_TYPE_VALUE] = [None] * 5


class inst(IntEnum):
    """EP type 'instanciation' value list index."""

    PACKAGE = 0 # (str) package name
    VERSION = 1 # (str) package version number
    MODULE = 2  # (str) module name
    NAME = 3    # (str) object name
    PARAM = 4   # (bool or None)


class vtype(IntEnum):
    """Validation type to use in validate()."""

    EP_TYPE_INT = 0
    EP_TYPE_STR = 1
    INSTANCE_STR = 2
    OBJECT = 3
    TYPE_OBJECT = 4


def import_str(ep_type_int):
    """Return the import string for ep_type_int.

    Args
    ----
    ep_type_int (int): A valid ep_type value.

    Returns
    -------
    (str): The import e.g. 'from numpy import float32 as numpy_float32'
    """
    i = ep_type_lookup['instanciation'][ep_type_int]
    if i[inst.MODULE] is None:
        return 'None'
    package = '' if i[inst.PACKAGE] is None else i[inst.PACKAGE] + '.'
    return f'from {package}{i[inst.MODULE]} import {i[inst.NAME]} as {i[inst.MODULE]}_{i[inst.NAME]}'


# TODO: Move this into validate or a function for asint & asstr too. Get rid of circular import.
# Import all types.
# If a type does not exist on this system remove it (all instance will be reteated as INVALID)
for ep_type_int, data in tuple(filter(lambda x: x[1][inst.MODULE] is not None, ep_type_lookup['instanciation'].items())):
    try:
        exec(import_str(ep_type_int))
    except ModuleNotFoundError:
        _logger.warning(f"Module '{data[inst.MODULE]}' was not found. '{data[inst.NAME]}' will be treated as an INVALID type.")
        del ep_type_lookup['n2v'][ep_type_lookup['v2n'][ep_type_int]]
        del ep_type_lookup['instanciation'][ep_type_int]
        del ep_type_lookup['v2n'][ep_type_int]
    else:
        _logger.info(import_str(ep_type_int))


# Must be defined after the imports
EP_TYPE_NAMES = set(ep_type_lookup['n2v'].keys())
EP_TYPE_VALUES = set(ep_type_lookup['v2n'].keys())


def validate(obj, vt=vtype.EP_TYPE_INT):
    """Validate an object as an EP type.

    The type of validation is specified by vt and may take one of 5 values:
        vtype.EP_TYPE_INT: object is an int and represents an EP type.
        vtype.EP_TYPE_STR: object is a str and represents an EP type.
        vtype.INSTANCE_STR: object is a valid python code str instanciating an object.
        vtype.OBJECT: object is an object
        vtype.TYPE_OBJECT: object is a type object of the type to be validated.

    Args
    ----
    object (object): See description above.
    vt (IntEnum): See description above.

    Returns
    -------
    (bool) True if the type is defined else false.
    """
    if vt == vtype.TYPE_OBJECT:
        return fully_qualified_name(obj()) in EP_TYPE_NAMES
    if vt == vtype.OBJECT:
        return fully_qualified_name(obj) in EP_TYPE_NAMES
    if vt == vtype.INSTANCE_STR:
        # TODO: try on nameerror do the inport
        return fully_qualified_name(eval(obj)) in EP_TYPE_NAMES
    if vt == vtype.EP_TYPE_STR:
        return obj != INVALID_EP_TYPE_NAME and obj in EP_TYPE_NAMES
    return obj != INVALID_EP_TYPE_VALUE and obj in EP_TYPE_VALUES


def asint(obj, vt=vtype.EP_TYPE_STR):
    """Return the EP type value for an object.

    The interpretation of the object is specified by vt and may take one of 5 values:
        vtype.EP_TYPE_INT: object is an int and represents an EP type.
        vtype.EP_TYPE_STR: object is a str and represents an EP type.
        vtype.INSTANCE_STR: object is a valid EP type name str.
        vtype.OBJECT: object is an object
        vtype.TYPE_OBJECT: object is a type object of the type to be validated.

    Args
    ----
    object (object): See description above.
    vt (IntEnum): See description above.

    Returns
    -------
    (int) The EP type of the object (may be egp.invalid_type)
    """
    if vt == vtype.TYPE_OBJECT:
        return ep_type_lookup['n2v'].get(fully_qualified_name(obj()), INVALID_EP_TYPE_VALUE)
    if vt == vtype.OBJECT:
        return ep_type_lookup['n2v'].get(fully_qualified_name(obj), INVALID_EP_TYPE_VALUE)
    if vt == vtype.INSTANCE_STR:
        return ep_type_lookup['n2v'].get(fully_qualified_name(eval(obj)), INVALID_EP_TYPE_VALUE)
    if vt == vtype.EP_TYPE_STR:
        return ep_type_lookup['n2v'].get(obj, INVALID_EP_TYPE_VALUE)


def asstr(obj, vt=vtype.EP_TYPE_INT):
    """Return the EP type string for an object.

    The interpretation of the object is specified by vt and may take one of 5 values:
        vtype.EP_TYPE_INT: object is an int and represents an EP type.
        vtype.EP_TYPE_STR: object is a str and represents an EP type.
        vtype.INSTANCE_STR: object is a valid EP type name str.
        vtype.OBJECT: object is an object
        vtype.TYPE_OBJECT: object is a type object of the type to be validated.

    Args
    ----
    object (object): See description above.
    vt (IntEnum): See description above.

    Returns
    -------
    (str) The EP type of the object (may be egp.invalid_type)
    """
    if vt == vtype.TYPE_OBJECT:
        ep_type_name = fully_qualified_name(obj())
        return ep_type_name if ep_type_name in EP_TYPE_NAMES else INVALID_EP_TYPE_NAME
    if vt == vtype.OBJECT:
        ep_type_name = fully_qualified_name(obj)
        return ep_type_name if ep_type_name in EP_TYPE_NAMES else INVALID_EP_TYPE_NAME
    if vt == vtype.INSTANCE_STR:
        ep_type_name = fully_qualified_name(eval(obj))
        return ep_type_name if ep_type_name in EP_TYPE_NAMES else INVALID_EP_TYPE_NAME
    if vt == vtype.EP_TYPE_INT:
        return ep_type_lookup['v2n'].get(obj, INVALID_EP_TYPE_NAME)


def fully_qualified_name(obj):
    """Return the fully qualified type name for obj.

    Args
    ----
    obj (object): Any object

    Returns
    -------
    (str): Fully qualified type name.
    """
    return obj.__class__.__module__ + '_' + obj.__class__.__qualname__


def compatible(a, b):
    """If EP type a is compatible with gc type b return True else False.

    a and b must be of the same type.

    TODO: Define what compatible means. For now it means 'exactly the same type'.

    Args
    ----
    a (int or str): A valid ep_type value or name.
    b (int or str): A valid ep_type value or name.

    Returns
    -------
    (bool): True if a and b are compatible.
    """
    return a == b


def type_str(ep_type_int):
    """Return the type string for ep_type_int.

    Args
    ----
    ep_type_int (int): A valid ep_type value.

    Returns
    -------
    (str): The type string e.g. 'int' or 'str'
    """
    i = ep_type_lookup['instanciation'][ep_type_int]
    return i[inst.NAME] if i[inst.MODULE] is None else f'{i[inst.MODULE]}_{i[inst.NAME]}'


def instance_str(ep_type_int, param_str=''):
    """Return the instanciation string for ep_type_int.

    Args
    ----
    ep_type_int (int): A valid ep_type value.
    param_str (str): A string to be used as instanciation parameters.

    Returns
    -------
    (str): The instanciation e.g. numpy_float32(<param_str>)
    """
    inst_str = type_str(ep_type_int)
    if ep_type_lookup['instanciation'][ep_type_int][inst.PARAM] is not None:
        inst_str += f'({param_str})'
    return inst_str


def interface_definition(xputs, vt=vtype.TYPE_OBJECT):
    """Create an interface definition from xputs.

    Used to define the inputs or outputs of a GC from an iterable
    of types, objects or EP definitions.

    Args
    ----
    xputs (iterable(object)): Object is of the type defined by vt.
    vt (vtype): The interpretation of the object is specified by vt and may take one of 5 values:
        vtype.EP_TYPE_INT: object is an int and represents an EP type.
        vtype.EP_TYPE_STR: object is a str and represents an EP type.
        vtype.INSTANCE_STR: object is a valid EP type name str.
        vtype.OBJECT: object is an object
        vtype.TYPE_OBJECT: object is a type object of the type to be validated.

    Returns
    -------
    tuple(list(ep_type_int), list(ep_type_int), list(int)): A list of the xputs as EP type in value
        format, a list of the EP types in xputs in value format in ascending order and a list of
        indexes into it in the order of xputs.
    """
    xput_eps = tuple((asint(x, vt) for x in xputs))
    xput_types = sorted(set(xput_eps)),
    return xput_eps, xput_types, [xput_types.index(x) for x in xputs]
