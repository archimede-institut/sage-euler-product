"""
Utils for Euler Product

utils_euler_product.py defines  class.

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

    sage: from euler_product import lattice_ivariant_euler_product as euler_p
    
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
import euler_product.lattice_invariant_euler_products 


from sage.symbolic.function import GinacFunction, BuiltinFunction
from sage.symbolic.symbols import register_symbol, symbol_table
from sage.functions.other import Function_ceil
from sage.rings.arith import prime_divisors


def table_performance(min_q, max_q, nb_decimals=100, big_p=300):
    """_summary_

    Parameters
    ----------
    min_q : int_
        minumum of q parameter
    max_q : integer
        maximum of parameter
    nb_decimals : int, optional
        decimal precision_, by default 100
    big_p : int, optional
        define a big number of p, by default 300
    """
    ref_time = 0.1 # approximate time is s for q = 3
    res = {}
    for q in range(min_q, max_q + 1):
        if (q%2 == 0) and (q%4 !=0):
            pass
        else:
            sys.stdout.write(str(q) + " ")
            sys.stdout.flush()
            aux = GetVs(q, 2, nb_decimals, big_p, -1)
            sys.stdout.write(str(aux[6]) + " digits for the first product\n")
            aux[5] = ceil(aux[5] * 1000/ref_time)
            res[q] = aux
            
    for q in range(min_q, max_q + 1):
        if q in res:
            str_res = str(q) + "& " + str(res[q][1])
            str_res = str_res + "& " + str(len(prime_divisors(res[q][1]))) 
            for i in range(2, 5):
                str_res = str_res + "& " + str(res[q][i])
                
            str_res = str_res + "& " + str(ceil(res[q][5]/1000)) + " \\\\"
            print(str_res)
    return