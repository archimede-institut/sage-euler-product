"""
Utils for Euler Product

utils_euler_product.py defines  functions
Main Engines

AUTHORS:

- Olivier Ramré (2023-01-008) : initial version
- Dominique BENIELLI (2023-02_15) :
  AMU University <dominique.benielli@univ-amu.fr>,
  Integration as  SageMath package. 
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
import sys
from collections import namedtuple
from timeit import default_timer as timer

from sage.symbolic.function import GinacFunction, BuiltinFunction
from sage.functions.other import ceil, floor
from sage.arith.misc import prime_divisors
from sage.arith.misc import euler_phi
from sage.arith.misc import gcd
from sage.arith.misc import sigma
from sage.arith.misc import divisors
from sage.arith.misc import prime_factors
from sage.arith.functions import lcm
from sage.rings.real_mpfi import RealIntervalField
from sage.rings.complex_interval_field import ComplexIntervalField
from sage.rings.complex_mpfr import ComplexField
from sage.modular.dirichlet import DirichletGroup
from sage.arith.misc import moebius


def nb_common_digits(a, b):
    """summary for nb_common_digits
    Returns -1 if floor(a) != floor(b)
    
    INPUT:

    - ''a'' -- [float]
        [first float to compare]

    - ''b'' -- [float]
        [second float to compare]
        
    
    OUTPUT:
    
    [int]
        [Returns -1 if floor(a) != floor(b), or number of common digit]
    
    """
    #Returns -1 if floor(a) != floor(b)
    # This is tailored for positive inputs
    if floor(a) != floor(b):
        return(-1)
    a = a - floor(a)
    b = b - floor(b)
    nb = 0
    while floor(10*a)==floor(10*b):
        a = 10*a - floor(10*a)
        b = 10*b - floor(10*b)
        nb += 1
    return(nb)

def laTeX_for_number(w, how_many, nb_block_sto_cut):
    """summary for laTeX_for_number
    Return a character string int(w).digits where digits concerns the first
    
    INPUT:

    - ''w'' -- [float]
        [w is a real number with a (short) integer part and a floating point]

    - ''how_many'' -- [int]
        [number of decimal,decimals, separated every 5 of them by \'\,\' 
        et every block of ''nb_block_sto_cut'', on a different line. '\cdots' ends the string]

    - ''nb_block_sto_cut'' -- [type]
        [description]

    OUTPUT:

    [str]
        [a character string int(w).digits where digits concerns the first]
        
        
    EXAMPLES:
    
        sage: from euler_product.utils_euler_product import laTeX_for_number
        sage: laTeX_for_number(22.01234567812345, 100, 8) 
        '22.&01234\\,56781\\,235\\cdots'
    
    """
    thelen = 5
    listchar = list(str(w))
    begpos = listchar.index('.') + 1
    listchar.insert(begpos, '&')

    if len(listchar) > how_many + begpos:
        listchar = listchar[0:how_many + begpos +1]

    n = thelen + begpos + 1
    while n < len(listchar):
        if (n - begpos)/(thelen +1) % nb_block_sto_cut == 0:
            listchar.insert(n, '\n')
        else:
            listchar.insert(n, '\\,')
        n += 1 + thelen

    listchar.append('\\cdots')
    return ''.join(listchar)

def sub_group_generated(n, q):
    """summary for sub_group_generated
    Main Engines function
    
    INPUT:

    - ''n'' -- int
        [description]

    - ''q'' -- int
        [description]

    OUTPUT:

    frozenset
        immutable set of 
        
    EXAMPLES:
    
        sage: from euler_product.utils_euler_product import sub_group_generated
        sage: sub_group_generated(5, 3) 
        frozenset({4, 7})
    """
    return frozenset(n^j % q for j in range(1, euler_phi(q)+1))

class LatticeInvariantClasses():
    """
    Summary for LatticeInvariantClasses
    
    EXAMPLE:
    
        sage: from euler_product.utils_euler_product import LatticeInvariant
        sage: LatticeInvariant(30)
        ((frozenset({0, 2, 3, 4, 5, 6, 7, 9}),
          frozenset({0, 1, 2, 3, 4, 5, 6, 15}),
          frozenset({3, 8, 9, 10, 12, 13, 14, 15}),
          frozenset({5, 8, 9, 10, 11, 12, 14, 15}),
          frozenset({16, 18, 19, 20, 21, 22, 23, 25}),
          frozenset({16, 17, 18, 20, 21, 22, 23, 27}),
          frozenset({16, 17, 18, 19, 20, 21, 22, 31}),
          frozenset({21, 24, 25, 26, 27, 28, 30, 31})),
         (frozenset({0, 2, 3, 4, 5, 6, 7, 9}),
          frozenset({0, 1, 2, 3, 4, 5, 6, 15}),
          frozenset({3, 8, 9, 10, 12, 13, 14, 15}),
          frozenset({5, 8, 9, 10, 11, 12, 14, 15}),
          frozenset({16, 18, 19, 20, 21, 22, 23, 25}),
          frozenset({16, 17, 18, 20, 21, 22, 23, 27}),
          frozenset({16, 17, 18, 19, 20, 21, 22, 31}),
          frozenset({21, 24, 25, 26, 27, 28, 30, 31})))
        
    """
    def __init__(self):
        pass

    def __call__(self, q):
        """summary for __call__ for usage of Class like function

        INPUT:

        - ''q'' -- [type]
            [description]

        OUTPUT:

        [tuple]
            [tuple of tuple, the_SG_tuple, the_Class_tuple]
            
        EXAMPLES:
        
            sage: from euler_product.utils_euler_product import LatticeInvariant
            sage: LatticeInvariant(10)
            ((frozenset({0, 2, 3, 5}),
              frozenset({0, 1, 2, 7}),
              frozenset({3, 4, 5, 6}),
              frozenset({8, 10, 11, 13})),
             (frozenset({0, 2, 3, 5}),
              frozenset({0, 1, 2, 7}),
              frozenset({3, 4, 5, 6}),
              frozenset({8, 10, 11, 13})))
             sage: LatticeInvariant(10)
            ((frozenset({0, 2, 3, 5}),
              frozenset({0, 1, 2, 7}),
              frozenset({3, 4, 5, 6}),
              frozenset({8, 10, 11, 13})),
            (frozenset({0, 2, 3, 5}),
             frozenset({0, 1, 2, 7}),
             frozenset({3, 4, 5, 6}),
             frozenset({8, 10, 11, 13})))
             sage: LatticeInvariant(20)
            ((frozenset({0, 2, 3, 4, 5, 6, 7, 9}),
              frozenset({0, 1, 2, 4, 5, 6, 7, 11}),
              frozenset({0, 1, 2, 3, 4, 5, 6, 15}),
              frozenset({1, 8, 10, 11, 12, 13, 14, 15}),
              frozenset({3, 8, 9, 10, 12, 13, 14, 15}),
              frozenset({5, 8, 9, 10, 11, 12, 14, 15}),
              frozenset({16, 18, 19, 20, 21, 22, 23, 25}),
              frozenset({16, 17, 18, 20, 21, 22, 23, 27})),
             (frozenset({0, 2, 3, 4, 5, 6, 7, 9}),
              frozenset({0, 1, 2, 4, 5, 6, 7, 11}),
              frozenset({0, 1, 2, 3, 4, 5, 6, 15}),
              frozenset({1, 8, 10, 11, 12, 13, 14, 15}),
              frozenset({3, 8, 9, 10, 12, 13, 14, 15}),
              frozenset({5, 8, 9, 10, 11, 12, 14, 15}),
              frozenset({16, 18, 19, 20, 21, 22, 23, 25}),
              frozenset({16, 17, 18, 20, 21, 22, 23, 27})))
        
        """
        if  hasattr(self, 'q') and  q == self.q:
            return (self.the_SG_tuple, self.the_Class_tuple)
        else:
            self.q = q
            the_SG_list = []
            for n in filter(lambda w: gcd(w, q) == 1, range(1, q)):
                my_sub = sub_group_generated(n, q)
                # print my_sub
                if my_sub not in the_SG_list:
                    the_SG_list.append(my_sub)

            the_SG_list.sort(key=len)
            the_SG_tuple = tuple(the_SG_list)

            ## Then get the classes:
            ## Create a copy:
            the_Class_list = the_SG_list.copy()

            for n in range(0, len(the_Class_list)):
                aux = the_Class_list[n]
                for m in range(n + 1, len(the_Class_list)):
                    if aux.issubset(the_SG_list[m]):
                        the_Class_list[m] = frozenset(the_Class_list[m].difference(aux))

            ## Create immutable tuples from the mutable lists:
            self.the_SG_tuple = tuple(the_SG_list)
            self.the_Class_tuple = tuple(the_Class_list)
        return (self.the_SG_tuple, self.the_Class_tuple)

LatticeInvariant = LatticeInvariantClasses()

BaseComponentStructure = namedtuple('BaseComponentStructure', ['the_SG_tuple', 'the_Class_tuple', 'nb_class', 'the_exponent',
                                                            'phi_q', 'character_group', 'invertibles', 'invariant_characters'])


class ComponentStructure(BaseComponentStructure):
    """summary for ComponentStructure inherit from  ''BaseComponentStructure''
    
    EXAMPLE:

        sage: from from euler_product.utils_euler_product import ComponentStructure
        sage: my_structure = ComponentStructure(30)

    """
    __slots__ = ()  # " the_SG_tuple, the_Class_tuple, nb_class, the_exponent, phi_q, character_group, invertibles, invariant_characters):

    def __init__(self, q):
        super(ComponentStructure, self).__init__()
        self._get_structure(q)

    def _get_structure(self, q):
        """summary for _get_structure

        INPUT:

        - ''q'' -- int
            [description]

        OUTPUT:

        namedtuple
            [description]
        """
        (the_SG_tuple, the_Class_tuple) = LatticeInvariant(q)
        ## size of our matrices:
        nb_class = len(the_SG_tuple)
        ## Getting the exponent:
        the_exponent = 1
        for i in range(0, nb_class):
            the_exponent = lcm(the_exponent, len(the_SG_tuple[i]))
        ## Now the character group:
        character_group = DirichletGroup(q)
        ## The range of invertible classes
        invertibles = tuple(n for n in filter(lambda w: gcd(w, q)==1, range(1, q)))
        if q == 1:
            invertibles = (1)
        ## For each cyclic subgroup G0, we need the list of characters for G0perp:
        invariant_characters = tuple(tuple(ind_e for ind_e in range(0, len(character_group))
                                        if the_SG_tuple[n].issubset(character_group[ind_e].kernel()))
                                                                for n in range(0, nb_class))
        
        self.the_SG_tuple = the_SG_tuple
        self.the_Class_tuple = the_Class_tuple
        self.nb_class = nb_class
        self.the_exponent = the_exponent
        self.phi_q = euler_phi(q)
        self.characters_group = character_group
        self.invertibles = invertibles
        self.invariant_characters = invariant_characters 
        self.q = q
        
    def getr_A_Kt(self, q,  my_indices):
        """ summary for getr_A_Kt

        INPUT:

        - ''q'' -- int
            [description]



        - ''my_indices'' -- int
            [description]

        OUTPUT:

        [set]
            [description]
        """
        # (theSGTuple, theClassTuple, nbclasses, theExponent,
        # phiq, characterGroup, invertibles, invariantCharacters) = structure
        r_A_Kt = {}
        for ind_A in range(0, self.nb_class):
            # TRICK ! The subgroup generated by A is at the same index!
            sgrA = self.the_SG_tuple[ind_A]
            for ind_K in range(0, self.nb_classes):
                K = self.the_SG_tuple[ind_K]
                for t in my_indices:
                    r_A_Kt[ind_A, ind_K, t] = 0
                    for L in self.the_SG_tuple:
                        if K.issubset(L):
                            Lt = {x^t % q for x in L}
                            if Lt == sgrA:
                                ## In Python 3, len(..)/len(...) is a real number, say 2.0
                                ## and the value of moebius becomes ... 1!
                                r_AK_t[ind_A, ind_K, t] += moebius(t) * moebius(int(len(L)/len(K)))*len(K)/self.phi_q

        return r_A_Kt

    def get_CA_Km(self, q, my_indices):
        """summary for get_CA_Km

        INPUT:

        - ''q'' -- [int]
            [description]

        - ''my_indices'' -- [int]
            [description]

        OUTPUT:

        [set]
            [description]
        """
        # (theSGTuple, theClassTuple, nbclasses, theExponent,
        # phiq, characterGroup, invertibles, invariantCharacters) = structure

        r_A_Kt = self.getr_A_Kt(q, my_indices)
        C_A_Km = {}
        for ind_A in range(0, self.nb_classes):
            for ind_K in range(0, self.nb_classes):
                for m in my_indices:
                    C_A_Km[ind_A, ind_K, m] = 0
                    for t in divisors(m):
                        C_A_Km[ind_A, ind_K, m] += r_A_Kt[ind_A, ind_K, t]

        return C_A_Km
    
##################################################
#######  Computing Gamma  ########################
##################################################

    def get_L_values(self, q, m, big_p, CIF, CF):
        """summary for get_L_values
        Computing Gamma

        INPUT:

        - ''q'' -- [int]
            [description]

        - ''m'' -- [int]
            [description]

        - ''big_pCIF'' -- [type]
            [description]

        - ''CF'' -- [type]
            [description]

        OUTPUT:

        [type]
            [description]
        """
        # m belongs to CIF
        # (theSGTuple, theClassTuple, nbclasses, theExponent,
        # phiq, characterGroup, invertibles, invariantCharacters) = structure

        CG = characterGroup.change_ring(CF)
        hurwitzvalues = tuple(hurwitz_zeta(s = m, x = CIF(a/q))/CIF(q)^m for a in structure.invertibles)

        aux0 = [[1-CIF(e(p))/CIF(p)^m
                for p in filter(lambda w: (w in Primes()), range(2, bigP))]
                for e in CG]
        aux1 = [prod(v) for v in aux0]
        aux2 = [sum([CIF(e(invertibles[ind_invert])) * hurwitzvalues[ind_invert]
                    for ind_invert in range(0, structure.phi_q)]) for e in CG]

        res = tuple(aux1[ind_e]*aux2[ind_e] for ind_e in range(0, structure.phi_q))
        return res

    def get_CA_Km_F_sur_H(self, q, my_indices, coeff_sf, coeff_sh):
        """summary for get_CA_Km_F_sur_H

        INPUT:

        - ''q'' -- [type]
            [description]

        - ''structure'' -- [type]
        
        
            [description]

        - ''my_indices'' -- [type]
            [description]

        - ''coeff_sf'' -- [type]
            [description]

        - ''coeff_sh'' -- [type]
            [description]

        OUTPUT

        [type]
            [description]
            
        EXAMPLES:
        
            sage: get_CA_Km_F_sur_H[1, -4, 4, 2, -4, 1], 11)
            [0, 4, 2, 2, 2, 3, 1, 0, -8, -22, -53]
                
        """
        # myindices should be divisor-closed (and include 1) and ordered
        #(the_SG_tuple, the_Class_tuple, nb_classes, the_exponent,
        # phi_q, character_group, invertibles, invariant_characters) = structure
        #print(q, structure, myindices, coeffsF, coeffsH)
        r_A_Kt = self.getr_A_Kt(q, my_indices)  # same indices as my_indices is divisor-closed
        max_index = my_indices[len(my_indices) - 1] 
        s_f = get_vector_sf(coeff_sf, max_index+1)
        s_h = get_vector_sf(coeff_sh, max_index+1)
        
        CAKmFsurH = {}
        for ind_a in range(0, self.nb_classes):
            for ind_k in range(0, self.nb_classes):
                for m in my_indices:
                    CAKmFsurH[ind_a, ind_k, m] = 0
                    for t in divisors(m):
                        #print("add",  rAKt[ind_A, ind_K, t]*(sH[m/t] - sF[m/t]))
                        CAKmFsurH[ind_a, ind_k, m] += r_A_Kt[ind_a, ind_k, t] * (s_h[m/t] - s_f[m/t])

        return C_A_KmFsurH

# checkGetLvalues(30, 2, 200, 212)
def get_gamma(q, t, structure, s, big_p, prec):
    """summary for get_gamma

    INPUT:

    - ''q'' -- [type]
        [description]

    - ''t'' -- [type]
        [description]

    - ''structure'' -- [type]
        [description]

    - ''s'' -- [type]
        [description]

    - ''big_p'' -- [type]
        [description]

    - ''prec'' -- [type]
        [description]

    OUTPUT:

    [type]
        [description]
        
    EXAMPLE:
    
         sage: check_get_L_values(30, 2, 200, 212)
    """
    #(theSGTuple, theClassTuple, nbclasses, theExponent,
    #             phiq, characterGroup, invertibles, invariantCharacters) = structure
    CIF = ComplexIntervalField(prec)
    CF = ComplexField(prec + 1)
    ##
    if s*t*log(big_p) > (prec + 10) * log(2):
        one = RealIntervalField(prec + 10)(1 - 2^(-prec-9), 1 + 2^(-prec - 9))
        Lvalues = (one, ) * len(characterGroup)
    else:
        m = CIF(t * s)
        Lvalues = GetLvalues(q, m, structure, big_p, CIF, CF)
    #print(q, t, s, bigP, prec, Lvalues)
    return vector([log(CIF(prod([Lvalues[ind_e] for ind_e in invariantCharacters[ind_G0]])).real())for ind_G0 in range(0, structure.nb_class)])

# myCIF = ComplexIntervalField(200)
# [real(u) for u in GetGamma(30, myCIF(2), GetStructure(30), 200, myCIF)]
################################################
############  Witt Decomposition  ##############
################################################

def get_vector_sf(coeffs_f, how_many):
    """

    INPUT:

    - ''coeffs_f'' -- [type]
        coeeficient f

    - ''how_many'' -- [type]
        [description]

    OUTPUT

    list 
        list des coefficient f
    """

    ann_i = coeffs_f  + ((how_many - len(coeffs_f)) * [0])
    s_f = [0 for i in range(0, how_many)]
    s_f[0] = len(coeffs_f)-1  # = degres of F
    for k in range(1, how_many):
        s_f[k] = -k * ann_i[k] - add(ann_i[i] * s_f[k-i] for i in range(1, k))
    return s_f


def get_vector_bf(coeffs_f, how_many):
    """get_vector_bf give the vectoc of coefficient bf
    Not used in the main program
    
    INPUT:

    - ''coeffs_f'' -- [type]
        [description]

    - ''how_many'' -- int
        number of coeeficients bf

    OUTPUT

    [list]
        [description]
        
    EXAMPLE:

        sage: from euler_product.utils_euler_product import get_vector_bf
        sage: get_vector_bf([1, -4, 4, 2, -4, 1], 11)
        [0, 4, 2, 2, 2, 3, 1, 0, -8, -22, -53]
        
    """
    
    b_f = [0 for i in range(0, how_many)] ## bf[0] is not used
    s_f = get_vector_sf(coeffs_f, how_many)
    for k in range(1, how_many):
        b_f[k] = add(moebius(k/d)*s_f[d] for d in divisors(k))/k
    return b_f



# strut = GetStructure(30)
# GetLvalues(30 ,1 ,strut,2, 200, 212)
def check_get_L_values(q, m, big_p, prec):
    """AI is creating summary for checkGetLvalues

    INPUT:

    - ''q'' -- [type]
        [description]

    - ''m'' -- [type]
        [description]

    - ''big_p'' -- [type]
        [description]

    - ''prec'' -- [type]
        [description]

    OUTPUT:

    [type]
        [description]
        
    EXAMPLE:

        sage: from euler_product.utils_euler_product import check_get_L_values
        sage: check_get_L_values(30, 2, 200, 212)
    """
    structure = ComponentStructure(q)
    #(theSGTuple, theClassTuple, nbclasses, theExponent,
    # phiq, characterGroup, invertibles, invariantCharacters) = structure
    CIF = ComplexIntervalField(prec)   
    Lval = GetLvalues(m , structure , big_p, CIF)
    return tuple(CIF(u) for u in [Lval[index] / prod([1 - characterGroup[index](p)/p^s
                                        for p in filter(lambda w: (w in Primes()), range(2, big_p))])
                                for index in range(0, structure.phi_q)])

    
def get_beta(F):
    """AI is creating summary for GetBeta

    INPUT:

    - ''F'' -- [type]
        [description]

    OUTPUT:

    [type]
        [description]
    """
    my_roots = F.roots(multiplicities=False)
    if len(my_roots) == 0:
        return 1
    else:
        return max(1, max([1/abs(c) for c in my_roots]))

def get_BetaRough(coeffs_f):
    """summary for get_BetaRough

    INPUT:

    - ''coeffs_f'' -- [type]
        [description]

    OUTPUT:

    [type]
        [description]
    """
    return max(1, max([abs(c) for c in coeffs_f]))