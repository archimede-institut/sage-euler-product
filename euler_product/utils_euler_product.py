"""
Utils_euler_product utilities  for Euler Product

utils_euler_product.py defines  functions
Main Engines

AUTHORS:

- Olivier Ramaré (2023-01-008) : initial version
- Dominique Benielli(2023-02_15) :
    Aix Marseille Université,
    Integration as SageMath package.
    Cellule de developpement Institut Archimède

WARNING:

    Needs Sage version at least 9.0
    CAREFUL, this is Python 3 code!

EXAMPLES::

    sage: from euler_product.utils_euler_product import  LatticeInvariantClasses

"""
# *****************************************************************************
#       Copyright (C) 2023 Olivier Ramare
#       < olivier . ramare @univ-amu.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from __future__ import print_function, absolute_import

from builtins import sum as add
from builtins import max
import numpy as np

from sage.functions.other import floor
from sage.arith.misc import euler_phi
from sage.arith.misc import gcd
from sage.arith.misc import divisors
from sage.arith.functions import lcm
from sage.rings.real_mpfi import RealIntervalField
from sage.rings.real_mpfi import RealIntervalFieldElement
from sage.rings.real_mpfr import RealNumber
from sage.rings.real_mpfr import RealField
from sage.rings.complex_interval_field import ComplexIntervalField
from sage.rings.complex_mpfr import ComplexField
from sage.modular.dirichlet import DirichletGroup
from sage.arith.misc import moebius
from sage.misc.misc_c import prod
from sage.functions.log import log
from sage.rings.integer_ring import ZZ
from sage.functions.transcendental import hurwitz_zeta
from sage.sets.primes import Primes
from sage.modules.free_module_element import vector


def nb_common_digits(a, b):
    r"""
    Returns -1 if floor(a) != floor(b).

    INPUT:

    - ``a`` -- float
        first float to compare.

    - ``b`` -- float
        second float to compare.

    OUTPUT:

    int
        Returns -1 if floor(a) != floor(b), or the number of common digits.

    EXAMPLES::

        sage: from euler_product.utils_euler_product import nb_common_digits
        sage: import numpy as np
        sage: nb_common_digits(1.33333, 1.334444)
        2
        sage: nb_common_digits(1.33333, 2.334444)
        -1
        sage: nb_common_digits(1.33333, np.inf)
        -1
        sage: nb_common_digits(np.inf, np.nan)
        -1
    """
    #  Returns -1 if floor(a) != floor(b)
    #  This is tailored for positive inputs
    if a is np.nan or b is np.nan or a is np.inf or b is np.inf:
        return (-1)
    if (isinstance(a, RealIntervalFieldElement) and isinstance(b, RealIntervalFieldElement)) or \
       (isinstance(a, RealNumber) and isinstance(b , RealNumber)):
        if a.is_NaN() or b.is_NaN() or a.is_infinity() or b.is_infinity():
            return (-1)
    if floor(a) != floor(b):
        return (-1)
    a = a - floor(a)
    b = b - floor(b)
    nb = 0
    while floor(10 * a) == floor(10 * b):
        a = 10 * a - floor(10 * a)
        b = 10 * b - floor(10 * b)
        nb += 1
    return(nb)


def laTeX_for_number(w, how_many, nb_block_sto_cut):
    r"""
    Return a character string representing the real number ``w`` made of its
    integer part followed by every decimal up to the``how_many`` -th decimals,
    where every block of 5 decimal is separated by \'\\\\,\',
    and every succession of ``how_many`` blocks is separated by \'\\n\'.
    The string has a \`&\` after the decimal point and ends with
    the string \`\\\\cdots\`.

    INPUT:

    - ``w`` -- float
        ``w`` is a real number with a (short) integer part and a floating point.

    - ``how_many`` -- int
        number of decimals, separated every 5 of them by \'\\,\'
        and every block of ``nb_block_sto_cut``, on a different line. '\\cdots' ends the string.

    - ``nb_block_sto_cut`` -- int
        See above.

    OUTPUT:

    str
        a character string int(w).separated_decimals where separated_decimals is LaTeX formatted
        version of the decimal expansion of ``w``, see the description of the function.

    EXAMPLES::

        sage: from euler_product.utils_euler_product import laTeX_for_number
        sage: laTeX_for_number(22.01234567812345, 100, 8)
        '22.&01234\\,56781\\,235\\cdots'

    """
    thelen = 5
    listchar = list(str(w))
    begpos = listchar.index('.') + 1
    listchar.insert(begpos, '&')

    if len(listchar) > how_many + begpos:
        listchar = listchar[0:how_many + begpos + 1]

    n = thelen + begpos + 1
    while n < len(listchar):
        if (n - begpos) / (thelen + 1) % nb_block_sto_cut == 0:
            listchar.insert(n, '\n')
        else:
            listchar.insert(n, '\\,')
        n += 1 + thelen

    listchar.append('\\cdots')
    return ''.join(listchar)


def sub_group_generated(n, q):
    """
    Return the frozenset of the multiplicative subgroup generated by the powers of ``n`` modulo ``q``.
    It is expected that ``n`` and ``q`` are coprime.

    INPUT:

    - ``n`` -- int
        an integer, expected to be coprime to ``q``.

    - ``q`` -- int
        a positive integer.

    OUTPUT:

    frozenset
        immutable set of the powers of ``n`` modulo ``q``

    EXAMPLES::

        sage: from euler_product.utils_euler_product import sub_group_generated
        sage: sub_group_generated(5, 3)
        frozenset({1, 2})
    """
    return frozenset(n**j % q for j in range(1, euler_phi(q) + 1))


class LatticeInvariantClasses():
    """
    This class takes a modulus ``q`` (i.e. a positive integer) and has two named accessors,
    ``the_SG_tuple`` and ``the_Class_tuple``.
    The SG tuple is the list of the multiplicative subgroups of :math:`(\mathbb{Z}/q\mathbb{Z})^*` that are generated
    by a single element.
    The Class tuple is the list of Lattice Invariant classes, namely the partition of :math:`(\mathbb{Z}/q\mathbb{Z})^*`
    made by the smallest non-empty intersections of elements of ``the_SG_tuple``.

    EXAMPLES::

        sage: from euler_product.utils_euler_product import LatticeInvariant
        sage: LatticeInvariant(30)
        ((frozenset({1}),
          frozenset({1, 11}),
          frozenset({1, 19}),
          frozenset({1, 29}),
          frozenset({1, 7, 13, 19}),
          frozenset({1, 17, 19, 23})),
         (frozenset({1}),
          frozenset({11}),
          frozenset({19}),
          frozenset({29}),
          frozenset({7, 13}),
          frozenset({17, 23})))

    """
    def __init__(self):
        pass

    def __call__(self, q):
        """summary for __call__ for usage of Class like function

        INPUT:

        - ``q`` -- int
          A positive integer.

        OUTPUT:

        pair
            pair of the tuples: the_SG_tuple, the_Class_tuple.

        EXAMPLES::

            sage: from euler_product.utils_euler_product import LatticeInvariant
            sage: LatticeInvariant(10)
            ((frozenset({1}), frozenset({1, 9}), frozenset({1, 3, 7, 9})),
             (frozenset({1}), frozenset({9}), frozenset({3, 7})))
            sage: lat = LatticeInvariant
            sage: lat.the_SG_tuple == LatticeInvariant(10)[0]
            True
            sage: lat.the_Class_tuple == LatticeInvariant(10)[1]
            True
            sage: lat.the_Class_tuple == LatticeInvariant(20)[1]
            False

        """
        if hasattr(self, 'q') and q == self.q:
            return (self.the_SG_tuple, self.the_Class_tuple)
        else:
            the_SG_list = []
            for n in filter(lambda w: gcd(w, q) == 1, range(1, q)):
                my_sub = sub_group_generated(n=n, q=q)
                #  print my_sub
                if my_sub not in the_SG_list:
                    the_SG_list.append(my_sub)
            the_SG_list.sort(key=len)
            #  the_SG_tuple = tuple(the_SG_list)

            #  Then get the classes:
            #  Create a copy:
            the_Class_list = the_SG_list.copy()

            for n in range(0, len(the_Class_list)):
                aux = the_Class_list[n]
                for m in range(n + 1, len(the_Class_list)):
                    if aux.issubset(the_SG_list[m]):
                        the_Class_list[m] = frozenset(the_Class_list[m].difference(aux))
            self.q = q
            #  Create immutable tuples from the mutable lists:
            self.the_SG_tuple = tuple(the_SG_list)
            self.the_Class_tuple = tuple(the_Class_list)
        return (self.the_SG_tuple, self.the_Class_tuple)


LatticeInvariant = LatticeInvariantClasses()

#  BaseTupledname = namedtuple('BaseComponentStructure', ['the_SG_tuple', 'the_Class_tuple', 'nb_class', 'the_exponent',
#                                                            'phi_q', 'character_group', 'invertibles', 'invariant_characters'])


class ComponentStructure():
    """
    This class takes a positive integer ``q`` and creates the following list of accessors:

    - ``phi_q``: the value of the Euler-phi function at ``q``.
    - ``the_exponent``: the exponent of the multiplicative group :math:`(\mathbb{Z}/q\mathbb{Z})^*`.
    - ``character_group``: the group of Dirichlet characters modulo ``q``, see this function for its description.
    - ``invertibles``: the tuple of the integers between 1 and ``q`` that are prime to ``q``.
    - ``the_SG_tuple`` and ``the_Class_tuple`` as in the class LatticeInvariantClass.
    - ``nb_class``: the number of Lattice Invariant classes.
    - ``invariant_characters``: given a subgroup in ``the_SG_tuple``, the tuple of the characters that leaves this subgroup invariant is created.
      ``invariant_characters`` is this list of tuples, arranged as in ``the_SG_tuple``.
    - ``getr_A_Kt``: a method used only for ``get_CA_Km`` and ``get_CA_Km_F_sur_H``.
        The coefficient C(A,K,m, F/H) are a sum on a variable t of s(F/H,m/t) times a function of t, say f(t).
        The lattice class A in given by its index ``ind_A`` in ``the_Class_tuple``, the subgroup K is given by its index ``ind_K`` in ``the_SG_tuple``.
        The function ``get_r_A_K_t`` answers a dictionary which to every ``(ind_A, ind_K, t)`` associates this f(t) (with the moebius factor).
        The list of ``t`` is of course limited and given as the input parameter of ``get_r_A_K_t``.
        This is the list of elements that form a divisor-closed subset of integers.
        This list is the same as the list of necessary values of ``m``.
    - ``get_CA_Km``: a method used for ``get_vs``.
        The coefficient C(A,K,m) are a sum on a variable t of a function of the value computed by ``getr_A_K_t``.
        The lattice class A in given by its index ``ind_A`` in ``the_Class_tuple``, the subgroup K is given by its index ``ind_K`` in ``the_SG_tuple``.
        The function ``get_CA_Km`` answers a dictionary which to every ``(ind_A, ind_K, m)`` associates this value.
    - ``get_CA_Km_F_sur_H``: a method used for ``get_euler_products``.
        The coefficient C(A,K,m, F/H) are a sum on a variable t of s(F/H, m/t) times a function of the value computed by ``getr_A_K_t``.
        The lattice class A in given by its index ``ind_A`` in ``the_Class_tuple``, the subgroup K is given by its index ``ind_K`` in ``the_SG_tuple``.
        The function ``get_CA_Km_F_sur_H`` answers a dictionary which to every ``(ind_A, ind_K, m)`` associates this value.
        When ``F == 1`` and ``H == 1-X``, the output of ``get_CA_Km_F_sur_H`` is the same as the one of ``get_CA_Km``.
    - ``get_L_values``: a method used only for ``get_gamma``.
    - ``get_gamma``: outputs the tuple defined in (22) of the :doc:`corresponding paper<../tutorial/LoeschianConstant-NS-04-MCOMP>`.
        For every cyclic subgroup :math:`G_0` in ``the_SG_tuple``,
        we compute :math:`\sum_{\chi\in G_0^\perp} \log L_P(t*s, \chi)`, where :math:`L_P(x,\chi)` is the L-series associated to :math:`\chi`,
        save that we remove the Euler factors for primes below ``P==big_p``.
        The output is the list of these values computed with ``prec`` correct binary digits.

    EXAMPLES::

        sage: from euler_product.utils_euler_product import ComponentStructure
        sage: structure = ComponentStructure(3)

    """
    def __init__(self, q):
        self._get_structure(q=q)
        self.q = q

    def _get_structure(self, q):
        """summary for _get_structure

        INPUT:

        - ``q`` -- int
            the modulus, i.e. a positive integer.

        OUTPUT:

        namedtuple
            see description of the class ``ComponentStructure``.
        """
        (the_SG_tuple, the_Class_tuple) = LatticeInvariant(q)
        #  size of our matrices:
        nb_class = len(the_SG_tuple)
        #  Getting the exponent:
        the_exponent = 1
        for i in range(0, nb_class):
            the_exponent = lcm(the_exponent, len(the_SG_tuple[i]))
        #  Now the character group:
        character_group = DirichletGroup(q)
        #  The range of invertible classes
        invertibles = tuple(n for n in filter(lambda w: gcd(w, q) == 1, range(1, q)))
        if q == 1:
            invertibles = (1,)
        #  For each cyclic subgroup G0, we need the list of characters for G0perp:
        invariant_characters = tuple(tuple(ind_e for ind_e in range(0, len(character_group))
                                           if the_SG_tuple[n].issubset(character_group[ind_e].kernel()))
                                     for n in range(0, nb_class))
        phi_q = euler_phi(q)
        self.the_SG_tuple = the_SG_tuple
        self.the_Class_tuple = the_Class_tuple
        self.nb_class = nb_class
        self.the_exponent = the_exponent
        self.phi_q = phi_q
        self.character_group = character_group
        self.invertibles = invertibles
        self.invariant_characters = invariant_characters
        self.q = q
        #  return (the_SG_tuple, the_Class_tuple, nb_class, the_exponent, phi_q,
        #        character_group, invertibles, invariant_characters)

    def getr_A_Kt(self, my_indices):
        """
        This method is used only for ``get_CA_Km`` and ``get_CA_Km_F_sur_H``. The coefficient C(A,K,m, F/H) are a sum on a variable t of s(F/H,m/t) times a function of t, say f(t).
        The lattice class A in given by its index ``ind_A`` in ``the_Class_tuple``, the subgroup K is given by its index ``ind_K`` in ``the_SG_tuple``.
        The function ``get_r_A_K_t`` answers a dictionary which to every ``(ind_A, ind_K, t)`` associates this f(t) (with the moebius factor).
        The list of ``t`` is of course limited and given as the input parameter of ``get_r_A_K_t``. This is the list of elements that form a divisor-closed subset of integers.
        This list is the same as the list of necessary values of ``m``.

        INPUT:

        - ``my_indices`` -- list
            list of indices (positive integers) ``t``. It should be divisor-closed (and include 1) and ordered increasingly.

        OUTPUT:

        dictionary
            output is a the dictionary (ind_A, ind_K, t) --> value, see above.

        EXAMPLES::

            sage: from euler_product.utils_euler_product import ComponentStructure
            sage: structure = ComponentStructure(3)
            sage: structure.getr_A_Kt([1, 2, 3, 4, 6])
            {(0, 0, 1): 1/2,
            (0, 0, 2): 0,
            (0, 0, 3): -1/2,
            (0, 0, 4): 0,
            (0, 0, 6): 0,
            (0, 1, 1): 0,
            (0, 1, 2): -1,
            (0, 1, 3): 0,
            (0, 1, 4): 0,
            (0, 1, 6): 1,
            (1, 0, 1): -1/2,
            (1, 0, 2): 0,
            (1, 0, 3): 1/2,
            (1, 0, 4): 0,
            (1, 0, 6): 0,
            (1, 1, 1): 1,
            (1, 1, 2): 0,
            (1, 1, 3): -1,
            (1, 1, 4): 0,
            (1, 1, 6): 0}


        """
        #  (theSGTuple, theClassTuple, nb_classes, theExponent,
        #  phi_q, characterGroup, invertibles, invariantCharacters) = structure
        r_A_Kt = {}
        for ind_a in range(0, self.nb_class):
            #  TRICK ! The subgroup generated by A is at the same index!
            sgr_A = self.the_SG_tuple[ind_a]
            for ind_k in range(0, self.nb_class):
                k_tuple = self.the_SG_tuple[ind_k]
                for t_i in my_indices:
                    r_A_Kt[ind_a, ind_k, t_i] = 0
                    for SG in self.the_SG_tuple:
                        if k_tuple.issubset(SG):
                            L_set = {x**t_i % self.q for x in SG}
                            if L_set == sgr_A:
                                #  In Python 3, len(..)/len(...) is a real number, say 2.0
                                #  and the value of moebius becomes ... 1!
                                r_A_Kt[ind_a, ind_k, t_i] += moebius(t_i) * moebius(int(len(SG) / len(k_tuple))) * len(k_tuple) / self.phi_q
        return r_A_Kt

    def get_CA_Km(self, my_indices):
        """
        ``get_CA_Km`` is a method used for ``get_vs``. The coefficient C(A,K,m) are a sum on a variable t of a function of the value computed by ``getr_A_K_t``.
        The lattice class A in given by its index ``ind_A`` in ``the_Class_tuple``, the subgroup K is given by its index ``ind_K`` in ``the_SG_tuple``.
        The function ``get_CA_Km`` answers a dictionary which to every ``(ind_A, ind_K, m)`` associates this value.

        INPUT:

        - ``my_indices`` -- [int]
            list of indices (positive integers) ``m``. It should be divisor-closed (and include 1) and ordered increasingly.

        OUTPUT:

        dictionary
            outputs the dictionary (ind_A, ind_K, m) --> value, see above.

        EXAMPLES::

            sage: from euler_product.utils_euler_product import ComponentStructure
            sage: from collections import OrderedDict
            sage: structure = ComponentStructure(3)
            sage: OrderedDict(structure.get_CA_Km([1, -4, 4, 2, -4, 1]))  # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
            OrderedDict([((0, 0, 1), 1/2), ((0, 0, -4), 1/2), ((0, 0, 4), 1/2), ((0, 0, 2), 1/2), ((0, 1, 1), 0), ((0, 1, -4), -1), ((0, 1, 4), -1), ((0, 1, 2), -1),
                ((1, 0, 1), -1/2), ((1, 0, -4), -1/2), ((1, 0, 4), -1/2), ((1, 0, 2), -1/2), ((1, 1, 1), 1), ((1, 1, -4), 1), ((1, 1, 4), 1), ((1, 1, 2), 1)])

        """
        # (theSGTuple, theClassTuple, nb_classes, theExponent,
        # phi_q, characterGroup, invertibles, invariantCharacters) = structure

        r_A_Kt = self.getr_A_Kt(my_indices=my_indices)
        C_A_Km = {}
        for ind_a in range(0, self.nb_class):
            for ind_k in range(0, self.nb_class):
                for m in my_indices:
                    C_A_Km[ind_a, ind_k, m] = 0
                    for t in divisors(m):
                        C_A_Km[ind_a, ind_k, m] += r_A_Kt[ind_a, ind_k, t]

        return C_A_Km

    def get_L_values(self, m, big_p, CIF, CF):
        """
        for every Dirichlet character :math:`\chi` modulo ``q``, we compute
        the L-series :math:`L_P(m, \chi)` associated to :math:\chi:,
        save that we remove the Euler factors for primes below ``P==big_p``.
        The output is the list of these values computed with ``prec`` correct binary digits.

        INPUT:


        - ``m`` -- [ComplexIntervalFieldElement]
            the point where the L-series are computed. The real part should be > 1 .

        - ``big_p`` -- int
            a positive integer. The Euler products are computed for primes above ``big_p``.

        - ``CIF`` -- Complex Interval Field
            [description]

        - ``CF`` -- Complex Field
            not used. Only ``CR.prec`` is used?

        OUTPUT:

        tuple
            the tuple of the values of :math:`L_P(m,\chi)`, where :math:`\chi` varies on the Dirichlet characters,
            values computed with ``prec`` correct binary digits.

        EXCEPTIONS:

            ValueError parameter ``m`` not in ``CIF``

        EXAMPLES::

            sage: from euler_product.utils_euler_product import ComponentStructure
            sage: structure = ComponentStructure(10)
            sage: CIF = ComplexIntervalField(200)
            sage: CF = ComplexField(200 + 1)
            sage: m = CIF(2)
            sage: structure.get_L_values(m, 200, CIF, CF)
            (1.0007481024252386196893654501571877025514323183079093676480?,
            1.0000226377974809104806639790897095274193344466859037418898? - 0.0000131408916900437454874106515694589606441168219958035059?*I,
            0.9999899240511933872962748479693199723956317768469030922497? + 4.14392795471732850815599881400588351007002717820829591633?e-63*I,
            1.0000226377974809104806639790897095274193344466859037418898? + 0.0000131408916900437454874106515694589606441168219958035059?*I)
            sage: m = CIF(2.1)
            sage: structure.get_L_values(m, 200, CIF, CF)
            (1.0004029274879933694024714910876346995176209724918492580239?,
            1.000013330852742794601876671961697811977029714891503324800? - 7.7538957108934769520297959484934618269499996602296768?e-6*I,
            0.9999947644552454506994437910117325481790746758589349726959? + 3.15595539279556818499806833653635488946195252544140357784?e-63*I,
            1.000013330852742794601876671961697811977029714891503324800? + 7.7538957108934769520297959484934618269499996602296768?e-6*I)

        """
        #  m belongs to CIF
        #  (theSGTuple, theClassTuple, nb_classes, theExponent,
        #  phi_q, characterGroup, invertibles, invariantCharacters) = structure
        if not m in CIF:
            raise ValueError("m parameter must belongs to CIF parameter")
        if m in ZZ and m > 1:
           preci = CF.prec()
           RIF = RealIntervalField(preci)
           RF = RealField(preci)
           CField = ComplexField(preci)
           CG = self.character_group.change_ring(CField)
           hurwitz_values = tuple(CIF(hurwitz_zeta(s=ZZ(m),
                                                x=CIF(a / self.q))._eval_self(RF)) / CIF(self.q)**m
                               for a in self.invertibles)  # type: ignore
        else:
            CG = self.character_group.change_ring(CF)
            hurwitz_values = tuple(CIF(hurwitz_zeta(s=m,
                                                x=CIF(a / self.q))) / CIF(self.q)**m
                               for a in self.invertibles)  # type: ignore
        aux0 = [[1 - CIF(e(p)) / CIF(p)**m
                for p in filter(lambda w: (w in Primes()), range(2, big_p))]
                for e in CG]
        aux1 = [prod(v) for v in aux0]
        aux2 = [sum([CIF(e(self.invertibles[ind_invert])) * hurwitz_values[ind_invert]  # type: ignore
                    for ind_invert in range(0, self.phi_q)]) for e in CG]

        res = tuple(aux1[ind_e] * aux2[ind_e] for ind_e in range(0, self.phi_q))
        return res

    def get_CA_Km_F_sur_H(self, my_indices, coeff_sf, coeff_sh):
        """
        ``get_CA_Km_F_sur_H``: a method used for ``get_euler_products`. The coefficient C(A,K,m, F/H) are a sum on a variable t of s(F/H, m/t) times a function of the value computed by ``getr_A_K_t``.
        The lattice class A in given by its index ``ind_A`` in ``the_Class_tuple``, the subgroup K is given by its index ``ind_K`` in ``the_SG_tuple``.
        The function ``get_CA_Km_F_sur_H`` answers a dictionary which to every ``(ind_A, ind_K, m)`` associates this value.
        When ``F == 1`` and ``H == 1-X``, the output of ``get_CA_Km_F_sur_H`` is the same as the one of ``get_CA_Km``.

        INPUT:

        - ``my_indices`` -- list[int]
            list of indices (positive integers) ``m``. It should be divisor-closed (and include 1) and ordered increasingly.

        - ``coeff_sf`` -- list[float]
            the list of the sum of the m-th power of the inverses of the roots of F.

        - ``coeff_sh`` -- [type]
            the list of the sum of the m-th power of the inverses of the roots of H.

        OUTPUT

        dictionary
            outputs the dictionary (ind_A, ind_K, m) --> value, see above.

        EXAMPLES:

            sage: from euler_product.utils_euler_product import ComponentStructure
            sage: structure = ComponentStructure(3)
            sage: structure.get_CA_Km_F_sur_H([1, 2, 3, 4, 5, 6], [1], [1, 0, -1])   # doctest: +NORMALIZE_WHITESPACE
            {(0, 0, 1): 0,
            (0, 0, 2): 1,
            (0, 0, 3): 0,
            (0, 0, 4): 1,
            (0, 0, 5): 0,
            (0, 0, 6): 0,
            (0, 1, 1): 0,
            (0, 1, 2): 0,
            (0, 1, 3): 0,
            (0, 1, 4): -2,
            (0, 1, 5): 0,
            (0, 1, 6): 0,
            (1, 0, 1): 0,
            (1, 0, 2): -1,
            (1, 0, 3): 0,
            (1, 0, 4): -1,
            (1, 0, 5): 0,
            (1, 0, 6): 0,
            (1, 1, 1): 0,
            (1, 1, 2): 2,
            (1, 1, 3): 0,
            (1, 1, 4): 2,
            (1, 1, 5): 0,
            (1, 1, 6): 0}
        """
        #  my_indices should be divisor-closed (and include 1) and ordered
        # (the_SG_tuple, the_Class_tuple, nb_classes, the_exponent,
        #  phi_q, character_group, invertibles, invariant_characters) = structure
        #  print(q, structure, my_indices, coeffsF, coeffsH)
        assert isinstance(my_indices, list), "my_indice must be a list"
        r_A_Kt = self.getr_A_Kt(my_indices)  # same indices as my_indices is divisor-closed
        max_index = my_indices[len(my_indices) - 1]
        s_f = get_vector_sf(coeff_sf, max_index + 1)
        s_h = get_vector_sf(coeff_sh, max_index + 1)

        CAKmF_sur_H = {}
        for ind_a in range(0, self.nb_class):
            for ind_k in range(0, self.nb_class):
                for m in my_indices:
                    CAKmF_sur_H[ind_a, ind_k, m] = 0
                    for t in divisors(m):
                        #  print("add",  rAKt[ind_A, ind_K, t]*(sH[m/t] - sF[m/t]))
                        CAKmF_sur_H[ind_a, ind_k, m] += r_A_Kt[ind_a, ind_k, t] * (s_h[m / t] - s_f[m / t])

        return CAKmF_sur_H

    def get_gamma(self, t, s, big_p, prec):
        """
        Outputs the tuple defined in (5.1) of the corresponding paper: for every cyclic subgroup :math:`G_0` in ``the_SG_tuple``,
        we compute :math:`\sum_{\chi\in G_0^\perp} \log L_P(t*s, \chi)`, where :math:`L_P(x,\chi)` is the L-series associated to :math:`\chi`,
        save that we remove the Euler factors for primes below ``P==big_p``.
        The output is the list of these values computed with ``prec`` correct binary digits.

        INPUT:

        - ``t`` -- int
            the L-series are computed at ``t*s``.

        - ``s`` -- float
            the L-series are computed at ``t*s``. The separation of ``t`` and ``s`` is only for readability of the code.

        - ``big_p`` -- int
            a positive integer. The Euler products are computing for primes larger than ``big_p``.

        - ``prec`` -- int
            number of correct binary digits in the output.

        OUTPUT:

        tuple
            the list of values of :math:`\sum_{\chi\in G_0^\perp} \log L_P(t*s, \chi)`, see the function description.


        EXAMPLES::

            sage: from euler_product.utils_euler_product import ComponentStructure
            sage: structure  = ComponentStructure(5)
            sage: structure.invariant_characters
            ((0, 1, 2, 3), (0, 2), (0,))
            sage: structure.get_gamma(1, 1.2, 20, 100)
            (0.412058674847838475387476473?, 0.3959326495526308567412224144?, 0.4113672762131896194520237806?)

        """
        # (theSGTuple, theClassTuple, nb_classes, theExponent,
        #             phi_q, characterGroup, invertibles, invariantCharacters) = structure
        CIF = ComplexIntervalField(prec)
        CF = ComplexField(prec + 1)
        if s * t * log(big_p) > (prec + 10) * log(2):
            one = RealIntervalField(prec + 10)(1 - 2**(- prec - 9), 1 + 2**(- prec - 9))
            l_values = (one, ) * len(self.character_group)
        else:
            m = CIF(t * s)
            l_values = self.get_L_values(m, big_p, CIF, CF)
        return vector([log(CIF(prod([l_values[ind_e] for ind_e in self.invariant_characters[ind_G0]])).real())for ind_G0 in range(0, self.nb_class)])


def get_vector_sf(coeffs_f, how_many):
    """
    A polynomial F is given by its list of coefficients, the first one being 1.
    The output is the list :math:`s_F(m)` for m less than ``how_many``, where :math:`s_F(m)`
    is the sum of the m-th power of the inverses of the roots of F.

    INPUT:

    - ``coeffs_f`` -- list[float]
        coefficients of the polynomial f, starting by 1.

    - ``how_many`` -- int
        number of computed coefficients.

    OUTPUT:

    list
        list des coefficient s_f(m) over ``m <= how_many``.

    EXAMPLES::

            sage: from euler_product.utils_euler_product import get_vector_sf
            sage: get_vector_sf([1, -1], 5)
            [1, 1, 1, 1, 1]
            sage: get_vector_sf([1, 1, 1], 10)
            [2, -1, -1, 2, -1, -1, 2, -1, -1, 2]
    """
    ann_i = coeffs_f + ((how_many - len(coeffs_f)) * [0])
    s_f = [0 for i in range(0, how_many)]
    s_f[0] = len(coeffs_f) - 1  # = degrees of F
    for k in range(1, how_many):
        s_f[k] = -k * ann_i[k] - add(ann_i[i] * s_f[k - i] for i in range(1, k))
    return s_f

def get_beta(F):
    """
    Outputs the maximum of 1 and of the inverse of the norm of the non-zero roots of the polynomial ``F``.

    INPUT:

    - ``F`` -- pol
        a polynomial with RealField coefficients.

    OUTPUT:

    float
       the maximum of 1 and of the inverse of the norm of the non-zero roots of ``F``.

    EXAMPLES::

        sage: from euler_product.utils_euler_product import get_beta
        sage: R0 = RealField(30)
        sage: R0X = R0['x']
        sage: (x,) = R0X._first_ngens(1)
        sage: F0 = R0X(1 - x^2)
        sage: get_beta(F0)
        1
    """
    my_roots = F.roots(multiplicities=False)
    # my_root must be never 0
    return max(1, max([1 / abs(c) for c in my_roots if c != 0], default=1))


def get_beta_rough(coeffs_f):
    """
    Outputs the maximum of 1 and of the sum of the norm of the coefficients of the polynomial ``F``,
    which is precisely given as the list ``coeffs_f``. This is intended to be an easy upper bound when the function
    ``get_beta`` takes too much time.

    INPUT:

    - ``coeffs_f`` -- float
        a list of floats, supposedly representing a polynomial ``F``.

    OUTPUT:

    float
        Outputs the maximum of 1 and of the sum of the norm of the elements of ``coeffs_f``.

    EXAMPLES::

        sage: from euler_product.utils_euler_product import get_beta_rough
        sage: get_beta_rough([1, 3, 4])
        8

    """
    return max(1, add(abs(c) for c in coeffs_f))
