"""
Utils for Euler Product

utils_euler_product.py defines  functions
Main Engines

AUTHORS:

- Olivier Ramaré (2023-01-008) : initial version
- Dominique Benielli(2023-02_15) :
  AMU University <dominique.benielli@univ-amu.fr>,
  Integration as  SageMath package. 
  Cellule de developpement Institut Archimède

WARNING:

    Needs Sage version at least 9.0
    CAREFUL, this is Python 3 code!
  
EXAMPLES::

    sage: from euler_product.utils_euler_product import  LatticeInvariantClasses
    # *****************************************************************************
#       Copyright (C) 2023 Olivier Ramare
#       < olivier . ramare @univ-amu.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************
"""
from __future__ import print_function, absolute_import

from builtins import sum as add
from builtins import max
from sage.functions.other import floor
from sage.arith.misc import euler_phi
from sage.arith.misc import gcd
from sage.arith.misc import divisors
from sage.arith.functions import lcm
from sage.rings.real_mpfi import RealIntervalField
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
    r"""summary for nb_common_digits
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
    r"""summary for laTeX_for_number
    Return a character string int(w).digits where digits concerns the first
    
    INPUT:

    - ''w'' -- [float]
        [w is a real number with a (short) integer part and a floating point]

    - ''how_many'' -- [int]
        [number of decimal,decimals, separated every 5 of them by \'\\,\' 
        et every block of ''nb_block_sto_cut'', on a different line. '\cdots' ends the string]

    - ''nb_block_sto_cut'' -- [int]
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

        - ''q'' -- [int]
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
            the_SG_list = []
            for n in filter(lambda w: gcd(w, q) == 1, range(1, q)):
                my_sub = sub_group_generated(n=n, q=q)
                # print my_sub
                if my_sub not in the_SG_list:
                    the_SG_list.append(my_sub)

            the_SG_list.sort(key=len)
            # the_SG_tuple = tuple(the_SG_list)

            ## Then get the classes:
            ## Create a copy:
            the_Class_list = the_SG_list.copy()

            for n in range(0, len(the_Class_list)):
                aux = the_Class_list[n]
                for m in range(n + 1, len(the_Class_list)):
                    if aux.issubset(the_SG_list[m]):
                        the_Class_list[m] = frozenset(the_Class_list[m].difference(aux))
            self.q = q
            ## Create immutable tuples from the mutable lists:
            self.the_SG_tuple = tuple(the_SG_list)
            self.the_Class_tuple = tuple(the_Class_list)
        return (self.the_SG_tuple, self.the_Class_tuple)

LatticeInvariant = LatticeInvariantClasses()
# 
# BaseTupledname = namedtuple('BaseComponentStructure', ['the_SG_tuple', 'the_Class_tuple', 'nb_class', 'the_exponent',
#                                                            'phi_q', 'character_group', 'invertibles', 'invariant_characters'])

class ComponentStructure:
    """creating summary for ComponentStructure
    
    
    EXAMPLES:
    
        sage: from euler_product.utils_euler_product import ComponentStructure
        sage: structure = ComponentStructure(3)
        
    """
    def __init__(self, q):
        self._get_structure(q=q)
        self.q = q
    
    def _get_structure(self, q):
        """summary for _get_structure

        INPUT:

        - ''q'' -- int
            [description]

        OUTPUT:

        namedtuple
            [description]
        """
        (the_SG_tuple, the_Class_tuple) =  LatticeInvariant(q)
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
        #return (the_SG_tuple, the_Class_tuple, nb_class, the_exponent, phi_q, 
        #        character_group, invertibles, invariant_characters)
    
    def getr_A_Kt(self, my_indices):
        """ summary for getr_A_Kt

        INPUT:


        - ''my_indices'' -- list
            [list of indices]

        OUTPUT:

        [set]
            [description]
            
        EXAMPLES:
        
            sage: from euler_product.utils_euler_product import ComponentStructure
            sage: structure = ComponentStructure(3)            
            sage: structure.getr_A_Kt([1, -4, 4, 2, -4, 1])
            {(0, 0, 1): 0, (0, 0, -4): 0, (0, 0, 4): 0, (0, 0, 2): 0}   
            
        """
        # (theSGTuple, theClassTuple, nb_classes, theExponent,
        # phi_q, characterGroup, invertibles, invariantCharacters) = structure
        r_A_Kt = {}
        for ind_a in range(0, self.nb_class):
            # TRICK ! The subgroup generated by A is at the same index!
            sgr_A = self.the_SG_tuple[ind_a]
            for ind_k in range(0, self.nb_class):
                k_tuple = self.the_SG_tuple[ind_k]
                for t_i in my_indices:
                    r_A_Kt[ind_a, ind_k, t_i] = 0
                    for SG in self.the_SG_tuple:
                        if k_tuple.issubset(SG):
                            L_set = {x**t_i % self.q for x in SG}
                            if L_set == sgr_A:
                                ## In Python 3, len(..)/len(...) is a real number, say 2.0
                                ## and the value of moebius becomes ... 1!
                                r_A_Kt[ind_a, ind_k, t_i] += moebius(t_i) * moebius(int(len(SG)/len(k_tuple)))*len(k_tuple)/self.phi_q
        return r_A_Kt

    def get_CA_Km(self, my_indices):
        """summary for get_CA_Km

        INPUT:

        - ''my_indices'' -- [int]
            [description]

        OUTPUT:

        [set]
            [description]
            
        EXAMPLES:
        
            sage: from euler_product.utils_euler_product import ComponentStructure
            sage: structure = ComponentStructure(3)
            sage: structure.get_CA_Km([1, -4, 4, 2, -4, 1])
            {(0, 0, 1): 0, (0, 0, -4): 0, (0, 0, 4): 0, (0, 0, 2): 0}
            
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
        """summary for get_L_values
        Computing Gamma

        INPUT:


        - ''m'' -- [ComplexIntervalFieldElement>]
            [description]

        - ''big_p'' -- [type]
            [description]
    
        - ''CIF'' -- [type]
            [description]
                        
        - ''CF'' -- [type]
            [description]

        OUTPUT:

        [type]
            [description]
            
        EXAMPLES:
        
            sage: from euler_product.utils_euler_product import ComponentStructure
            sage: structure = ComponentStructure(3)
            sage: myCIF = ComplexIntervalField(200)
            sage: structure.get_L_values(m, big_p, CIF, CF)
            
            
            
        """
        # m belongs to CIF
        # (theSGTuple, theClassTuple, nb_classes, theExponent,
        # phi_q, characterGroup, invertibles, invariantCharacters) = structure

        CG = self.character_group.change_ring(CF)
        if m in ZZ:
            m_new = ZZ(m)
        else:
            m_new = m
        hurwitz_values = tuple(CIF(hurwitz_zeta(s=m_new, 
                                                x=CIF(a / self.q))._complex_mpfi_(CF)) / CIF(self.q)**m 
                            for a in self.invertibles)  # type: ignore

        aux0 = [[1-CIF(e(p))/CIF(p)**m
                for p in filter(lambda w: (w in Primes()), range(2, big_p))]
                for e in CG]
        aux1 = [prod(v) for v in aux0]
        aux2 = [sum([CIF(e(self.invertibles[ind_invert])) * hurwitz_values[ind_invert] # type: ignore
                    for ind_invert in range(0, self.phi_q)]) for e in CG]

        res = tuple(aux1[ind_e]*aux2[ind_e] for ind_e in range(0, self.phi_q))
        return res

    def get_CA_Km_F_sur_H(self, my_indices, coeff_sf, coeff_sh):
        """summary for get_CA_Km_F_sur_H

        INPUT:

        - ''structure'' -- [type]
            [description]

        - ''my_indices'' -- [list]
            [list of indices]

        - ''coeff_sf'' -- [type]
            [description]

        - ''coeff_sh'' -- [type]
            [description]

        OUTPUT

        [type]
            [description]
            
        EXAMPLES:
        
            sage: from euler_product.utils_euler_product import ComponentStructure
            sage: structure = ComponentStructure(3)
            sage: structure.get_CA_Km_F_sur_H([1, -4, 4, 2, -4, 1], 11)
            [0, 4, 2, 2, 2, 3, 1, 0, -8, -22, -53]
                
        """
        # my_indices should be divisor-closed (and include 1) and ordered
        #(the_SG_tuple, the_Class_tuple, nb_classes, the_exponent,
        # phi_q, character_group, invertibles, invariant_characters) = structure
        # print(q, structure, my_indices, coeffsF, coeffsH)
        assert isinstance(my_indices, list), "my_indice must be a list"
        r_A_Kt = self.getr_A_Kt(my_indices)  # same indices as my_indices is divisor-closed
        max_index = my_indices[len(my_indices) - 1] 
        s_f = get_vector_sf(coeff_sf, max_index+1)
        s_h = get_vector_sf(coeff_sh, max_index+1)
        
        CAKmF_sur_H = {}
        for ind_a in range(0, self.nb_class):
            for ind_k in range(0, self.nb_class):
                for m in my_indices:
                    CAKmF_sur_H[ind_a, ind_k, m] = 0
                    for t in divisors(m):
                        #print("add",  rAKt[ind_A, ind_K, t]*(sH[m/t] - sF[m/t]))
                        CAKmF_sur_H[ind_a, ind_k, m] += r_A_Kt[ind_a, ind_k, t] * (s_h[m/t] - s_f[m/t])

        return CAKmF_sur_H

    def get_gamma(self, t, s, big_p, prec):
        """summary for get_gamma

        INPUT:

        - ''t'' -- [type]
            [description]

        - ''s'' -- [type]
            [description]

        - ''big_p'' -- [int]
            [description]

        - ''prec'' -- [int]
            [number of digits]

        OUTPUT:

        [type]
            [description]
            
        EXAMPLE:
        
            sage: from euler_product.utils_euler_product import ComponentStructure
            sage: structure  = ComponentStructure(3)
            sage: structure.get_gamma(30, 2, 200, 212)
            (0)
            sage: check_get_L_values(30, 2, 200, 212)
            
        """
        #(theSGTuple, theClassTuple, nb_classes, theExponent,
        #             phi_q, characterGroup, invertibles, invariantCharacters) = structure
        CIF = ComplexIntervalField(prec)
        CF = ComplexField(prec + 1)
        ##
        if s*t*log(big_p) > (prec + 10) * log(2):
            one = RealIntervalField(prec + 10)(1 - 2**(-prec-9), 1 + 2**(-prec - 9))
            l_values = (one, ) * len(self.character_group)
        else:
            m = CIF(t * s)
            l_values = self.get_L_values(m, big_p, CIF, CF)
        # print(q, t, s, bigP, prec, l_values)
        return vector([log(CIF(prod([l_values[ind_e] for ind_e in self.invariant_characters[ind_G0]])).real())for ind_G0 in range(0, self.nb_class)])

        
# myCIF = ComplexIntervalField(200)
# [real(u) for u in GetGamma(30, myCIF(2), GetStructure(30), 200, myCIF)]
################################################
############  Witt Decomposition  ##############
################################################

def get_vector_sf(coeffs_f, how_many):
    """

    INPUT:

    - ''coeffs_f'' -- [type]
        coefficient f

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
    """get_vector_bf give the vector of coefficient bf
    Not used in the main program
    
    INPUT:

    - ''coeffs_f'' -- [type]
        [description]

    - ''how_many'' -- int
        number of coefficients bf

    OUTPUT

    [list]
        [description]
        
    EXAMPLES:

        sage: from euler_product.utils_euler_product import get_vector_bf
        sage: get_vector_bf([1, -4, 4, 2, -4, 1], 11)
        [0, 4, 2, 2, 2, 3, 1, 0, -8, -22, -53]
        
    """
    
    b_f = [0 for i in range(0, how_many)] ## bf[0] is not used
    s_f = get_vector_sf(coeffs_f, how_many)
    for k in range(1, how_many):
        b_f[k] = add(moebius(k/d) * s_f[d] for d in divisors(k)) / k  # type: ignore
    return b_f



# strut = GetStructure(30)
# GetLvalues(30 ,1 ,strut,2, 200, 212)
def check_get_L_values(q, m, s, big_p, prec):
    """summary for check_get_L_values

    INPUT:

    - ''q'' -- [type]
        [description]

    - ''m'' -- [type]
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

        sage: from euler_product.utils_euler_product import check_get_L_values
        sage: check_get_L_values(30, 1, 2, 200, 212)
        ([+infinity .. +infinity],
         [.. NaN ..] + [.. NaN ..]*I,
         [.. NaN ..] + [.. NaN ..]*I,
         [.. NaN ..] + [.. NaN ..]*I,
         [.. NaN ..] + [.. NaN ..]*I,
         [.. NaN ..] + [.. NaN ..]*I,
         [.. NaN ..] + [.. NaN ..]*I,
         [.. NaN ..] + [.. NaN ..]*I)
         
    """
    structure = ComponentStructure(q)
    #(theSGTuple, theClassTuple, nb_classes, theExponent,
    # phi_q, characterGroup, invertibles, invariantCharacters) = structure
    CIF = ComplexIntervalField(prec)   
    CF = ComplexField(prec + 1)
    l_val = structure.get_L_values(m=m, big_p=big_p, CIF=CIF, CF=CF) # type: ignore
    return tuple(CIF(u) for u in [l_val[index] / prod([1 - structure.character_group[index](p)/p**s
                                        for p in filter(lambda w: (w in Primes()), 
                                                        range(2, big_p))])
                                for index in range(0, structure.phi_q)])


def get_beta(F):
    """summary for get_beta

    INPUT:

    - ''F'' -- [type]
        [description]

    OUTPUT:

    [type]
        [description]
    """
    my_roots = F.roots(multiplicities=False)
    # my_root must be never 0
    return max(1, max([1/abs(c) for c in my_roots if c != 0.0], default=1))

def get_beta_rough(coeffs_f):
    """summary for get_BetaRough

    INPUT:

    - ''coeffs_f'' -- [type]
        [description]

    OUTPUT:

    [type]
        [description]
    """
    return max(1, max(abs(c) for c in coeffs_f))
