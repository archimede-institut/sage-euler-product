"""
Lattice Invariant

lattice_invariant_euler_product.py defines  class.
Computing some structural invariants


AUTHORS:

- Olivier Ramaré (2023-01-008) : initial version
- Dominique Benielli (2023-02_15) :
  Aix Marseille Université ,
  Integration as  SageMath package.
  Cellule de developpement Institut Archimède

...WARNING:

  Needs Sage version at least 9.0
  CAREFUL, this is Python 3 code!

EXAMPLES::

    sage: from euler_product.lattice_invariant_euler_products import get_euler_products

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
from timeit import default_timer as timer

from sage.functions.other import ceil, floor
from sage.sets.primes import Primes
from sage.modules.free_module_element import vector
from sage.functions.log import exp
from sage.functions.log import log
from sage.misc.misc_c import prod
from sage.misc.flatten import flatten
from sage.arith.misc import prime_factors
from sage.arith.misc import prime_divisors
from sage.rings.integer import Integer
from sage.rings.real_mpfi import RealIntervalField
from sage.rings.real_mpfr import RealField
from euler_product.utils_euler_product import ComponentStructure
from euler_product.utils_euler_product import get_beta
from euler_product.utils_euler_product import laTeX_for_number
from euler_product.utils_euler_product import nb_common_digits


def get_vs(q, s, nb_decimals=100, big_p=100, verbose=2, with_laTeX=0, digits_offset=10):
    """summary for get_vs
    Computes an approximate value of the list (zeta(s; q, A))
    for A in the lattice-invariant classes.
    We compute directly what happens for primes < big_p.
    Main Engines GetVs(12, 2, 100, 110)

    ... WARNING ...
        Don't try me and avoid to select a prime number for big_p.

    INPUT:

    - ``q`` -- [type]
        [description]

    - ``s`` -- [type]
        [description]

    - ``nb_decimals`` -- [type]
        [description], by default 100

    - ``big_p`` -- int, optional
        [description], by default 100

    - ``verbose`` -- int, optional
        [description], by default 2

    - ``with_laTeX`` -- int, optional
        [description], by default 0

    - ``digits_offset`` : int, optional
        [description], by default 10

    OUTPUT:

    [type]
        [description]

    EXAMPLES::

        sage: from euler_product.lattice_invariant_euler_products import get_vs
        sage: get_vs(8, 3, 100)  # doctest: +NORMALIZE_WHITESPACE
        Computing the structural invariants ...  done.
        Computing big m ... Computing the finite product for p < 100 ...  done.
        done: we use big_m = 18 .
        Building indices ... done: there are 5 summands.
        -------------------
        For p + 8ZZ in frozenset({1})
        the product of 1 / (1 - p^{-3}) is between
        1.00022487189858709876213790850059727367160406698000012383809358224349804195751158838756078343541771445946357
        and
        1.00022487189858709876213790850059727367160406698000012383809358224349804195751158838756078343541771445946376
        (Obtained:  105  correct decimal digits)
        -------------------
        For p + 8ZZ in frozenset({3})
        the product of 1 / (1 - p^{-3}) is between
        1.03941995442465275247765320826895776436688562576128093245373042938586210850475647986883816951846918262944993
        and
        1.03941995442465275247765320826895776436688562576128093245373042938586210850475647986883816951846918262945008
        (Obtained:  103  correct decimal digits)
        -------------------
        For p + 8ZZ in frozenset({5})
        the product of 1 / (1 - p^{-3}) is between
        1.00859929667035271781603021669885964068116773568938335524379250743431508655959698706547733008351071600945961
        and
        1.00859929667035271781603021669885964068116773568938335524379250743431508655959698706547733008351071600945975
        (Obtained:  105  correct decimal digits)
        -------------------
        For p + 8ZZ in frozenset({7})
        the product of 1 / (1 - p^{-3}) is between
        1.00305724526111080016727462667660715366818826631992203857221002140263476460120525161534802433614351034528702
        and
        1.00305724526111080016727462667660715366818826631992203857221002140263476460120525161534802433614351034528716
        (Obtained:  105  correct decimal digits)
        ((frozenset({1}), frozenset({3}), frozenset({5}), frozenset({7})),
        ((1.00022487189858709876213790850059727367160406698000012383809358224349804195751158838756078343541771445946357,
        1.00022487189858709876213790850059727367160406698000012383809358224349804195751158838756078343541771445946376),
        (1.03941995442465275247765320826895776436688562576128093245373042938586210850475647986883816951846918262944993,
        1.03941995442465275247765320826895776436688562576128093245373042938586210850475647986883816951846918262945008),
        (1.00859929667035271781603021669885964068116773568938335524379250743431508655959698706547733008351071600945961,
        1.00859929667035271781603021669885964068116773568938335524379250743431508655959698706547733008351071600945975),
        (1.00305724526111080016727462667660715366818826631992203857221002140263476460120525161534802433614351034528702,
        1.00305724526111080016727462667660715366818826631992203857221002140263476460120525161534802433614351034528716)))

    TESTS::

        sage: from euler_product.lattice_invariant_euler_products import get_vs
        sage: get_vs(3, 1, 100)  # doctest: +NORMALIZE_WHITESPACE
        Computing the structural invariants ...  done.
        Computing big m ... Computing the finite product for p < 100 ...  done.
        done: we use big_m = 52 .
        Building indices ... done: there are 6 summands.
        -------------------
        For p + 3ZZ in frozenset({1})
        the product of 1 / (1 - p^{-1}) is between
        NaN
        and
        NaN
        (Obtained:  -1  correct decimal digits)
        -------------------
        For p + 3ZZ in frozenset({2})
        the product of 1 / (1 - p^{-1}) is between
        NaN
        and
        NaN
        (Obtained:  -1  correct decimal digits)
        ((frozenset({1}), frozenset({2})), ((NaN, NaN), (NaN, NaN)))

    """
    start = timer()
    #  Compute the structural invariants:
    if verbose >= 2:
        sys.stdout.write("Computing the structural invariants ... ")

    structure = ComponentStructure(q)
    # (theSGTuple, theClassTuple, nb_classes, theExponent,
    # phi_q, characterGroup, invertibles, invariantCharacters) = structure

    if verbose >= 2:
        print(" done.")

    #  Getting bigM:
    if verbose >= 2:
        sys.stdout.write("Computing big m ... ")

    allowed_primes = set(prime_factors(structure.phi_q))
    cte = structure.nb_class * prod([1 + 2 / (p - 1) for p in allowed_primes])
    big_m = 10
    while (float(log(cte) + log(1 + big_p / (big_m * s)) - s * big_m * log(big_p)
                 + (nb_decimals + 1) * log(10)) > 0):
        big_m = big_m + 1
    #  Initial computations:
    if verbose >= 2:
        sys.stdout.write("Computing the finite product for p < " + str(big_p) + " ... ")

    prec = ceil(nb_decimals * log(10) / log(2) + 10 + ceil(big_m))
    R = RealIntervalField(prec)
    #  Empty initial products are allowed:
    euler_prod_ini = tuple([R(1 / prod(
        flatten([1, [R(1 - 1 / p**s) for p in filter(lambda w: (w in Primes())
                and (w % q in structure.the_Class_tuple[i]), range(2, big_p))]])))
        for i in range(0, structure.nb_class)])
    if verbose >= 2:
        print(" done.")
    log_err = R(cte * (1 + big_p / (big_m * s)) / big_p**(s * big_m))

    if verbose >= 2:
        print("done: we use big_m =", big_m, ".")

    #  Build the set of indices di:
    if verbose >= 2:
        sys.stdout.write("Building indices ... ")

    my_indices = [m for m in filter(lambda w:
                                    set(prime_factors(w)).issubset(allowed_primes),
                                    range(1, big_m))]
    #  1 is indeed in my_indices.
    CAKm = structure.get_CA_Km(my_indices=my_indices)
    if verbose >= 2:
        print("done: there are", len(my_indices), "summands.")
    vs_approx = [R(0)] * structure.nb_class
    for m in my_indices:
        # print(q, m, s, bigP, prec)|-|
        aux = structure.get_gamma(m, s, big_p, prec)
        for ind_a in range(0, structure.nb_class):
            for ind_k in range(0, structure.nb_class):
                vs_approx[ind_a] += aux[ind_k] * CAKm[ind_a, ind_k, m] / m
    #  We now have to get the Euler products:
    eulerProds = tuple([(R(euler_prod_ini[i] * exp(vs_approx[i] - log_err)).lower(),
                       R(euler_prod_ini[i] * exp(vs_approx[i] + log_err)).upper())
                        for i in range(0, structure.nb_class)])

    if verbose >= 2:
        for i in range(0, structure.nb_class):
            nb_digits = nb_common_digits(eulerProds[i][1], eulerProds[i][0])  # type: ignore
            print("-------------------")
            print("For p + " + str(q) + "ZZ in", structure.the_Class_tuple[i])
            print("the product of 1 / (1 - p^{-" + str(s) + "}) is between")
            print(eulerProds[i][0])  # type: ignore
            print("and")
            print(eulerProds[i][1])  # type: ignore
            if with_laTeX == 1:
                print("LaTeX format:")
                how_many = min(nb_decimals, nb_digits)
                print(laTeX_for_number(eulerProds[i][0], how_many, 10))  # type: ignore
            print("(Obtained: ", nb_digits, " correct decimal digits)")

    end = timer()

    if verbose == 1:
        print("Time taken:", end - start, "seconds.")
    # print(my_indices)
    if verbose == -1:
        return([big_p, structure.phi_q, len(my_indices), structure.nb_class, big_m, end - start,
                - floor(log(eulerProds[0][1] - eulerProds[0][0]) / log(10))])  # type: ignore
    else:
        return structure.the_Class_tuple, eulerProds


def get_euler_products(q, s, f_init, h_init, nb_decimals=100, big_p=300, verbose=2, with_laTeX=0):
    r"""Summary for get_euler_products
    Computes an approximate value of the list ()
    for A in the lattice-invariant classes, ``LatticeInvariantClasses``.
    Careful! bigP may increase in the proof!
    We compute directly what happens for primes < big_p.

    to do

    assert F[0] = H[0] = 1

    GetEulerProds(3, 1, 1-x^2, 100)
    (q, s, Finit, Hinit, nbdecimals, bigP= 00, Verbose=2, WithLaTeX=0)
    (q, s, F, H, nbdecimals, bigP = 100, Verbose = 2, WithLaTeX = 0)


    INPUT:

    - ``q`` -- [type]
        [description]

    - ``s`` -- [type]
        [description]

    - ``f_init`` -- [type]
        [description]

    - ``h_init`` -- [type]
        [description]

    - ``nb_decimals`` -- [type]
        [description]

    - ``big_p`` : int, optional
        [description], by default 100

    - ``verbose`` : int, optional
        [description], by default 2

    - ``with_laTeX`` : int, optional
        [description], by default 0

    OUTPUT:

    [tuple]
        [return tuple of (the_Class_tuple, euler_prods)]

    EXCEPTIONS:
        ValueError   ('non convergent product')
        ValueError("f_init[0] and h_init[0] must be equal to 1")

    EXAMPLES::

        sage: from euler_product.lattice_invariant_euler_products import get_euler_products
        sage: get_euler_products(3, 1, 1-x^2, 1, 100, verbose=0)  # doctest: +NORMALIZE_WHITESPACE
        ((frozenset({1}), frozenset({2})),
        ((0.9671040753637981066150556834173635260473412207450092130719978569438733967843271277395717230016746853806050215621235810749643636399725665325875376146914709362753787689855429317947529895445140974345,
        0.9671040753637981066150556834173635260473412207450092130719978569438733967843271277395717230016746853806050215621235810749643636399725665325875376146914709362753787689855429317947529895445140974477),
        (0.7071813747951674302088659938984504109243584468119496848353517677901518159831128643782536704398941052120208041311403202957250160794697319584608281454011743387515885835706146696365506658500107821105,
        0.7071813747951674302088659938984504109243584468119496848353517677901518159831128643782536704398941052120208041311403202957250160794697319584608281454011743387515885835706146696365506658500107821228)))
        sage: get_euler_products(8, 1, 1-2*x-7*x^2-4*x^3, 1-2*x+x^2, 110, 50, 0, 0)  # doctest: +NORMALIZE_WHITESPACE
        ((frozenset({1}), frozenset({3}), frozenset({5}), frozenset({7})),
        ((0.95694534785160118343696705727389182875317497729139147890543260424601701644488885948144051203907950843120693625020812276103159906625162986027780621959661661058242880638494733,
        0.95694534785160118343696705727389182875317497729139147890543260424601701644488885948144051203907950843120698715116910955482476415178562075800144085574035025387260631491020912),
        (-1.1744975454490034177155999666151209160615842743206296883125115013500947685793723674107330894468690657084850070249226110372932226333018435422450675119284554631204742715594780,
        -1.1744975454490034177155999666151209160615842743206296883125115013500947685793723674107330894468690657084850694977072727918933849576550337569751438204088603874573384648809813),
        (0.41303522430115813462933064104403132614973953141027630692477294245200590806656802414019031361757979375573715775187373750238649707664458852273214341861159234461302771883822062,
        0.41303522430115813462933064104403132614973953141027630692477294245200590806656802414019031361757979375573717972166689137936668012188844661293239520403935555328934054712625207),
        (0.73440531676250743047848787987972815324702590578219161532691471695315326661191699963756307237940693554033848218351377181727690073086294638978098135674230697193936653008013195,
        0.73440531676250743047848787987972815324702590578219161532691471695315326661191699963756307237940693554033852124733175620467822747169360181588443153429051853637559530089233888)))
    """
    start = timer()
    # assert F[0] = H[0] = 1
    #  Compute the structural invariants:
    if verbose >= 2:
        sys.stdout.write("Computing the structural invariants ... ")

    structure = ComponentStructure(q=q)
    #  (theSGTuple, theClassTuple, nb_classes, theExponent,
    #  phi_q, characterGroup, invertibles, invariantCharacters) = structure

    if verbose >= 2:
        print(" done.")

    #  A small precision is enough:
    R0 = RealField(30)
    R0X = R0['x']
    (x,) = R0X._first_ngens(1)
    F0, H0 = R0X(f_init), R0X(h_init)
    if Integer(H0[0]) != Integer(F0[0]):
        raise ValueError("f_init[0] and h_init[0] must be equal to 1")
    my_delta = (F0 - H0).valuation()
    if my_delta * s <= 1:
      raise ValueError('non convergent product')
    #  Get my_beta, myDelta and big_p:
    my_beta = max(2, get_beta(F0 / F0[0]), get_beta(H0 / H0[0]))

    if verbose >= 2:
        print("We have Delta =", my_delta, "and beta =", my_beta)

    #  Getting bigM, prec and bigP:
    big_p = max(big_p, 2 * my_beta)
    cte = 4 * structure.nb_class**2 * (F0.degree() + H0.degree()) * (s + big_p)
    big_m = big_p + 10

    while (float(log(cte) + (big_m + 1) * log(big_p**s / my_beta) - (nb_decimals + 1) * log(10)) < 0):
        big_m = big_m + 10

    #  The coefficients CA(K,m,F/H) may increase like beta^m,
    #  This spoils the accuracy and has to be recovered:
    prec = ceil(nb_decimals * log(10) / log(2) + 10) + ceil(float(big_m * log(my_beta) / log(2)))

    if verbose >= 2:
        print("We use big_m =", big_m, ", big_p =", big_p, "and working prec =", prec)
    #  The precision has changed! Change the ring:
    R = RealIntervalField(prec)
    RF = RealField(prec + 1)
    log_err = R(cte * (my_beta / (big_p**s))**(big_m + 1))
    RX = R['x']
    (x,) = RX._first_ngens(1)
    F, H = RX(f_init / Integer(F0[0])), RX(h_init / Integer(H0[0]))

    #  Initial computations:
    if verbose >= 2:
        sys.stdout.write("Computing the finite products for p < " + str(big_p) + " ... ")
    #  Empty initial products are allowed:
    #  print([ p for p in range(2, big_p)])
    #  print([F(1/Integer(p)**s) for p in range(2, big_p)])
    #  print([H(1/Integer(p)**s) for p in range(2, big_p)])
    #  print([F(1/p**s)/H(1/p**s) for p in range(2, big_p)])
    #  print([R(F(1/p**s)/H(1/p**s)) for p in Primes()])
    #  prod_list = [1]
    #  for i in range(0, structure.nb_class):
    #      prod_list.append([R(F(1/Integer(p)**s) / H(1/Integer(p)**s))
    #                                       for p in filter(lambda w: (w in Primes())
    #                                       and (w%q in structure.the_Class_tuple[i]),
    #                                       range(2, big_p))])
    #  eulerProdIni = tuple(prod(flatten(prod_list)))

    eulerProdIni = tuple(prod(flatten([1, [R(F(1 / Integer(p)**s) / H(1 / Integer(p)**s))
                                           for p in filter(lambda w: (w in Primes())
                                                           and (w % q in structure.the_Class_tuple[i]),
                                           range(2, big_p))]]))
                         for i in range(0, structure.nb_class))

    if verbose >= 2:
        print(" done.")

    #  Compute CA(K, m, F/H):
    if verbose >= 1:
        sys.stdout.write("Computing C_A(K, m, F/H) ...\n")

    my_indices = [i for i in range(1, big_m + 1)]
    CAKmF_sur_H = structure.get_CA_Km_F_sur_H(my_indices, F.list(), H.list())
    logZs_approx = vector([R(0)] * structure.nb_class)

    #  Most of time is spent here.
    #  The main loop in m:
    for mm in range(my_delta, big_m + 1):
        aux = structure.get_gamma(mm, s, big_p, prec)
        for ind_a in range(0, structure.nb_class):
            for ind_k in range(0, structure.nb_class):
                logZs_approx[ind_a] += aux[ind_k] * CAKmF_sur_H[ind_a, ind_k, mm] / mm
    #  End of the main loop in m
    #  We now have to complete the Euler products:
    eulerProds = tuple([(R(eulerProdIni[i] * exp(logZs_approx[i] - log_err)).lower(),
                        R(eulerProdIni[i] * exp(logZs_approx[i] + log_err)).upper())
                        for i in range(0, structure.nb_class)])

    if verbose >= 2:
        for i in range(0, structure.nb_class):
            nb_digits = nb_common_digits(eulerProds[i][1], eulerProds[i][0])
            print("-------------------")
            print("For p + " + str(q) + " ZZ in ", structure.the_Class_tuple[i])
            print("For F(x) =", f_init)
            print("and H(x) =", h_init)
            if s == 1:
                print("the product of F(1 / p) / H( 1 / p) is between")
            else:
                print("the product of F(1 / p^" + str(s) + ")/H(1 / p^" + str(s) + ") is between")
            print(eulerProds[i][0])
            print("and")
            print(eulerProds[i][1])
            if with_laTeX == 1:
                print("LaTeX format:")
                how_many = min(nb_decimals, nb_digits)
                print(laTeX_for_number(eulerProds[i][0], how_many, 10))
            print("(Obtained: ", nb_digits, " correct decimal digits)")
    end = timer()
    if verbose >= 1:
        print("Time taken: ", end - start, "seconds.")
    #  print my_indices
    return structure.the_Class_tuple, eulerProds


def table_performance(min_q, max_q, nb_decimals=100, big_p=300):
    """Summary for table_performance
    Table Makers

    INPUT:

    - ``min_q`` -- [type]
        [description]

    - ``max_q`` -- [type]
        [description]

    - ``nb_decimals`` -- int, optional
        [description], by default 100

    - ``big_p`` -- int, optional
        [description], by default 300

    EXAMPLES:

        sage: from euler_product.lattice_invariant_euler_products import table_performance
        sage: table_performance(10, 30)  # random
        11 102 digits for the first product
        12 102 digits for the first product
        13 102 digits for the first product
        15 102 digits for the first product
        16 102 digits for the first product
        17 102 digits for the first product
        19 102 digits for the first product
        20 102 digits for the first product
        21 102 digits for the first product
        23 102 digits for the first product
        24 102 digits for the first product
        25 102 digits for the first product
        27 102 digits for the first product
        28 102 digits for the first product
        29 102 digits for the first product
        11& 10& 2& 8& 4& 21& 4 \\
        12& 4& 1& 5& 4& 21& 1 \\
        13& 12& 2& 10& 6& 21& 5 \\
        15& 8& 1& 5& 6& 21& 2 \\
        16& 8& 1& 5& 6& 21& 2 \\
        17& 16& 1& 5& 5& 21& 4 \\
        19& 18& 2& 10& 6& 21& 8 \\
        20& 8& 1& 5& 6& 21& 2 \\
        21& 12& 2& 10& 8& 21& 6 \\
        23& 22& 2& 6& 4& 21& 6 \\
        24& 8& 1& 5& 8& 21& 2 \\
        25& 20& 2& 8& 6& 21& 8 \\
        27& 18& 2& 10& 6& 21& 8 \\
        28& 12& 2& 10& 8& 21& 6 \\
        29& 28& 2& 7& 6& 21& 9 \\
    """
    ref_time = 0.1  # approximate time is s for q = 3
    res = {}
    for q in range(min_q, max_q + 1):
        if (q % 2 == 0) and (q % 4 != 0):
            pass
        else:
            sys.stdout.write(str(q) + " ")
            sys.stdout.flush()
            aux = get_vs(q, 2, nb_decimals, big_p, -1)
            sys.stdout.write(str(aux[6]) + " digits for the first product\n")  # type: ignore
            aux[5] = ceil(aux[5] * 1000 / ref_time)  # type: ignore
            res[q] = aux

    for q in range(min_q, max_q + 1):
        if q in res:
            str_res = str(q) + "& " + str(res[q][1])
            str_res = str_res + "& " + str(len(prime_divisors(res[q][1])))
            for i in range(2, 5):
                str_res = str_res + "& " + str(res[q][i])
            str_res = str_res + "& " + str(ceil(res[q][5] / 1000)) + " \\\\"
            print(str_res)
    return


def get_vs_checker(q, s, borne=10000):
    """summary for get_vs_checker
    Checking Engines

    INPUT:

    - ``q`` -- [type]
        [description]

    - ``s`` -- [type]
        [description]

    - ``borne`` -- int, optional
        boundary of computation, by default 10000

    EXAMPLES::

        sage: from euler_product.lattice_invariant_euler_products import get_vs_checker
        sage: get_vs_checker(8, 1)
        -------------------
        For p mod  8  in  frozenset({1})
        the product of 1/(1-p^{- 1 }) is about 1.3126486598837095
        -------------------
        For p mod  8  in  frozenset({3})
        the product of 1/(1-p^{- 1 }) is about 2.217318868012348
        -------------------
        For p mod  8  in  frozenset({5})
        the product of 1/(1-p^{- 1 }) is about 1.7913334693038976
        -------------------
        For p mod  8  in  frozenset({7})
        the product of 1/(1-p^{- 1 }) is about 1.5751024072482105

    """
    #  Computes an approximate value of the list (zeta(s; q, A))
    #  for A in the lattice-invariant classes.
    structure = ComponentStructure(q)
    #  (theSGTuple, theClassTuple, nb_class, theExponent,
    #  phi_q, characterGroup, invertibles, invariantCharacters) = structure
    vs_approx = [1 / prod([1.0 - 1 / p**s
                           for p in filter(lambda w: (w in Primes()) and (w % q in structure.the_Class_tuple[i]),
                                           range(2, borne))])
                 for i in range(0, structure.nb_class)]
    for i in range(0, structure.nb_class):
        print("-------------------")
        print("For p mod ", q, " in ", structure.the_Class_tuple[i])
        print("the product of 1/(1-p^{-", s, "}) is about", vs_approx[i])
