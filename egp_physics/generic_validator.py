"""Generic Genetic Code Validation Schema.

The generic_validator() provides validation functions for _GC derived classes.
Validation is in the context of the individual GC i.e. does not consider if
dependenceis, like gca or ancestora exist etc.

Normalization functions are provided for Fields that can be defined within the
context of the individual e.g. signature, input_types etc.

The SCHEMA is a 'catch all' schema and will be modified based on the specific
_GC class being validated.
"""

from datetime import datetime
from hashlib import sha256
from json import load
from os.path import dirname, join
from pprint import pformat
from random import choice, getrandbits

from .utils.base_validator import BaseValidator

with open(join(dirname(__file__), "formats/generic_gc_format.json"), "r") as file_ptr:
    SCHEMA = load(file_ptr)


_SIGN = (1, -1)


def random_reference():
    """Quickest way to get a random (enough) reference."""
    return getrandbits(63) * choice(_SIGN)


def define_signature(gc):
    """Define the signature of a genetic code.

    The signature for a codon GC is slightly different to a regular GC.

    Args
    ----
    gc(dict): Must at least be an cGC.

    Returns
    -------
    (str): Lowercase hex SHA256 string.
    """
    # NOTE: This needs to be very specific and stand the test of time!
    string = pformat(gc['graph'], indent=0, sort_dicts=True, width=65535, compact=True) + str(gc['gca']) + str(gc['gcb'])

    # If it is a codon glue on the mandatory definition
    if "generation" in gc and gc["generation"] == 0:
        if "meta_data" in gc and "function" in gc["meta_data"]:
            string += gc["meta_data"]["function"]["python3"]["0"]["inline"]
            if 'code' in gc["meta_data"]["function"]["python3"]["0"]:
                string += gc["meta_data"]["function"]["python3"]["0"]["code"]
    return sha256(string.encode()).hexdigest()


class generic_validator(BaseValidator):
    """Generic validator of _GC classes."""

    # TODO: Make errors ValidationError types for full disclosure
    # https://docs.python-cerberus.org/en/stable/customize.html#validator-error

    def _check_with_valid_created(self, field, value):

        if isinstance(value, datetime):
            date_time_obj = value
        else:
            try:
                date_time_obj = datetime.strptime(value, "%Y-%m-%dT%H:%M:%S.%fZ")
            except ValueError:
                self._error(
                    field, "Created date-time is not valid. Unknown error parsing.")
                return

        if date_time_obj > datetime.utcnow():
            self._error(field, "Created date-time cannot be in the future.")

    def _check_with_valid_type_input_index(self, field, value):
        # TODO: Check type indices make sense
        pass

    def _check_with_valid_type_output_index(self, field, value):
        # TODO: Check type indices make sense
        pass

    def _check_with_valid_inline(self, field, value):
        # TODO: Check right number of return parameters and arguments
        pass

    def _check_with_valid_callable(self, field, value):
        # TODO: Check right number of return parameters and arguments. Check arguments all have default=None.
        pass

    def _normalize_default_setter_set_ref(self, document):
        """Fast way to get a unique (enough) reference."""
        return random_reference()

    def _normalize_default_setter_set_signature(self, document):
        """Define the signature of a genetic code."""
        return define_signature(self.document)

    def _normalize_default_setter_set_input_types(self, document):
        # Gather all the input endpoint types. Reduce in a set then order the list.
        inputs = []
        for row in document["graph"].values():
            inputs.extend((ep[2] for ep in filter(lambda x: x[0] == 'I', row)))
        return sorted(set(inputs))

    def _normalize_default_setter_set_output_types(self, document):
        # Gather all the output endpoint types. Reduce in a set then order the list.
        return sorted(set((ep[2] for ep in document["graph"].get("O", tuple()))))

    def _normalize_default_setter_set_input_indices(self, document):
        # Get the type list then find all the inputs in order & look them up.
        type_list = self._normalize_default_setter_set_input_types(document)
        inputs = []
        for row in document["graph"].values():
            inputs.extend((ep for ep in filter(lambda x: x[0] == 'I', row)))
        return bytes((type_list.index(ep[2]) for ep in sorted(inputs, key=lambda x: x[1])))

    def _normalize_default_setter_set_output_indices(self, document):
        # Get the type list then find all the inputs in order & look them up.
        type_list = self._normalize_default_setter_set_output_types(document)
        bytea = (type_list.index(ep[2]) for ep in sorted(document["graph"].get("O", tuple()), key=lambda x: x[1]))
        return bytes(bytea)

    def _normalize_default_setter_set_num_inputs(self, document):
        return len(document["graph"].get("I", tuple()))

    def _normalize_default_setter_set_num_outputs(self, document):
        return len(document["graph"].get("O", tuple()))

    def _normalize_default_setter_set_opt_num_codons(self, document):
        return 1 if document['gca'] is None and document['gcb'] is None else 0

    def _normalize_default_setter_set_created(self, document):
        return datetime.utcnow()


"""
Some benchmarking on SHA256 generation
======================================
Python 3.8.5

>>> def a():
...     start = time()
...     for _ in range(10000000): int(sha256("".join(string.split()).encode()).hexdigest(), 16)
...     print(time() - start)
...
>>> a()
8.618626356124878
>>> def b():
...     start = time()
...     for _ in range(10000000): int.from_bytes(sha256("".join(string.split()).encode()).digest(), 'big')
...     print(time() - start)
...
>>> b()
7.211490631103516
>>> def c():
...     start = time()
...     for _ in range(10000000): sha256("".join(string.split()).encode()).hexdigest()
...     print(time() - start)
...
>>> c()
6.463267803192139
>>> def d():
...     start = time()
...     for _ in range(10000000): sha256("".join(string.split()).encode()).digest()
...     print(time() - start)
...
>>> d()
6.043259143829346
>>> def e():
...     start = time()
...     for _ in range(10000000): {sha256("".join(string.split()).encode()).digest(): "Test"}
...     print(time() - start)
...
>>> e()
6.640311002731323
>>> def f():
...     start = time()
...     for _ in range(10000000): {int.from_bytes(sha256("".join(string.split()).encode()).digest(), 'big'): "Test"}
...     print(time() - start)
...
>>> f()
7.6320412158966064
>>> def g():
...     start = time()
...     for _ in range(10000000): {sha256("".join(string.split()).encode()).hexdigest(): "Test"}
...     print(time() - start)
...
>>> g()
7.144319295883179
>>> def h1():
...     start = time()
...     for _ in range(10000000): getrandbits(256)
...     print(time() - start)
...
>>> h1()
1.0232288837432861
>>> def h2():
...     start = time()
...     for _ in range(10000000): getrandbits(128)
...     print(time() - start)
...
>>> h2()
0.8551476001739502
>>> def h3():
...     start = time()
...     for _ in range(10000000): getrandbits(64)
...     print(time() - start)
...
>>> h3()
0.764052152633667
>>> def i():
...     start = time()
...     for _ in range(10000000): getrandbits(256).to_bytes(32, 'big')
...     print(time() - start)
...
>>> i()
2.038336753845215
"""


"""
Some Benchmarking on hashing SHA256
===================================
Python 3.8.5

>>> a =tuple( (getrandbits(256).to_bytes(32, 'big') for _ in range(10000000)))
>>> b =tuple( (int(getrandbits(63)) for _ in range(10000000)))
>>> start = time(); c=set(a); print(time() - start)
1.8097834587097168
>>> start = time(); d=set(b); print(time() - start)
1.0908379554748535
"""
