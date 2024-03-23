.. _index:


Tutorial Euler Product
======================

Introduction and principles
---------------------------

The main object of this software is to compute in a very fast manner Euler products with rational local factors over primes in some collections of arithmetic progressions modulo some fixed modulus :math:`q`. The computations are done with interval arithmetic, the results are certified and given in the form :code:`(lower_bound, upper_bound)`.
An Euler product is a product over all the (positive integer) primes, possibly in some subsequence. 
Here are two examples  :math:`E_1=\prod_{p\equiv 3,5[7]}(1-1/p^2)`, where the product is taken over every prime :math:`p` congruent to 3 or 5 modulo 7, and :math:`E_2 = \prod_{p\equiv 1[8]}\bigl(1-\frac{4}{p}\bigr)\bigl(\frac{p+1}{p-1}\bigr)^2`, the so-called Shanks's constant, where the product is taken over every prime congruent to 1 modulo 8.

During the computations, Euler products over rational functions such as :math:`E_2` are inferred from simpler Euler products of the shape 

.. math::
   \prod_{p\in\mathcal{A}\mod q}(1-1/p^s),

by following the 2021 paper of Ettahri, Surel and Ramaré,  where

 * :math:`\mathcal{A}` is some subset of :math:`G=(\mathbb{Z}/q\mathbb{Z})^\times`. The subset :math:`\mathcal{A}` has to be the union of "lattice invariant classes", as described below.
 * :math:`q` is a positive integer, the so-called "modulus". We have :math:`q=7` for :math:`E_1`.
 * :math:`s` is a real parameter that is strictly positif and \*in this example\* strictly larger than 1. A typical choice is :math:`s=2`.

.. seealso::

   The mathematical proof of this software is taken from the paper \*\* Fast multi-precision computation of some Euler products \*\* by Salma Ettahri, Olivier Ramaré and Léon Surel, published in 2021, in volume 90 of \*Mathematics of Computations\*, pages 2247 to 2265. 

In case :math:`q = 1`, the notion of lattice invariant classes is trivial. Let us start by describing this case.

Euler Product over every primes
-------------------------------

.. code-block:: default
     
     from euler_product.lattice_invariant_euler_products import get_euler_products
     get_euler_products(1, 2.1 , 1-x^2, 1+x^3, 103, 20, verbose = 0, with_Latex = 0, digits_offset = 10)

This computes the Euler product :math:`\prod_{p\ge2}\frac{1-1/p^{2s}}{1+1/p^{3s}}` where :math:`s=2.1` with potentially 103 correct digits and by computing directly the Euler product for all the primes less than :math:`P=20`. This value of :math:`P` is 100 by default. The level of comments :code:`verbose` can be set to 0, 1 or 2. The additional parameter :code:`with_Latex` is either equal to 1 or not equal to 1, with an obvious meaning. If the output does not have enough correct digits, the user is asked to increase the value 103 to 110 for instance. We decided not to automate this behaviour.

The parameter :code:`digits_offset` is not used as of now.

On the effect of the choice of :math:`s`, notice that the two calls

.. code-block:: default
     
     get_euler_products(1, 1 , 1-x^4, 1, 103)
     get_euler_products(1, 2 , 1-x^2, 1, 103)

give the same answer, which is readily seen to be an approximation of :math:`1/\zeta(4)`, where :math:`\zeta` is the Riemann-zeta function. Recall that we have :math:`\zeta(4)=\pi^4/90`, a fact that we may use to check our code.

Lattice Invariant Classes modulo :math:`q`
------------------------------------------


- Definition of Lattice Invariant Classes

We subdivide the multiplicative group :math:`G=(\mathbb{Z}/q\mathbb{Z})^\times` in classes called `Lattice Invariant Classes`.
Two points are in the same class if and only if they generate the same subgroup modulo :math:`q`.

When :math:`q = 15` these classes are obtained by

.. code-block:: default

       LatticeInvariant(15)[1]
       (frozenset({1}), frozenset({4}), frozenset({11}), frozenset({14}), frozenset({8, 2}), frozenset({13, 7}))


.. code-block:: default
     
        from euler_product.utils_euler_product import LatticeInvariant
        LatticeInvariant(15)
        ((frozenset({1}), frozenset({1, 4}), frozenset({1, 11}), frozenset({1, 14}), frozenset({8, 1, 2, 4}), frozenset({1, 4, 13, 7})), 
         (frozenset({1}), frozenset({4}), frozenset({11}), frozenset({14}), frozenset({8, 2}), frozenset({13, 7})))

The output is a couple whose first element is the tuple of the monogenic subgroups of :math:`G=(\mathbb{Z}/q\mathbb{Z})^\times` and whose second element is the tuple of lattice invariant classes.

- Low level tools

.. code-block:: default
       
       from euler_product.utils_euler_product import ComponentStructure
       mystructure = ComponentStructure(3)

This class proposes several quantities. It is used by the high level function :code:`get_vs` and :code:`get_euler_products`, so the user does not have to worry about it. However the quantities computed may have interest.

 * :code:`mystructure.q`: the modulus :math:`q`.
 * :code:`mystructure.phi_q`: the value of the Euler phi-function at :math:`q`.
 * :code:`mystructure.the_exponent`: the exponent of the group :math:`G=(\mathbb{Z}/q\mathbb{Z})^\times`.
 * :code:`mystructure.invertibles`: the tuple of invertibles in :math:`(\mathbb{Z}/q\mathbb{Z})`, i.e. an enumeration of :math:`G=(\mathbb{Z}/q\mathbb{Z})^\times`.
 * :code:`mystructure.the_SG_tuple`: the tuple of the subgroups of :math:`G=(\mathbb{Z}/q\mathbb{Z})^\times` that are generated by a single element. Such subgroups are also called \*monogenic\* subgroups.
 * :code:`mystructure.the_Class_tuple`: the tuple of the lattice invariant classes.
 * :code:`mystructure.nb_class`: the number of lattice invariant classes.
 * :code:`mystructure.character_group`: the character group of :math:`G=(\mathbb{Z}/q\mathbb{Z})^\times`.
 * :code:`mystructure.invariant_characters`: for each monogenic subgroup in :code:`mystructure.the_SG_tuple`, the list of (the indices of) the characters that has this subgroup in its kernel. The order of :code:`mystructure.invariant_characters` is the same as the one in :code:`mystructure.the_SG_tuple`.
 * Some methods are also available.

Euler Product over primes in arithmetic progression
---------------------------------------------------

We start with the three data:

* A modulus :code:`q`:math:`\ge 1`.
* A rational fraction given in the form :math:`F(x)/H(x)` where :math:`F(x)` and :math:`H(x)` are two polynomials with real coefficients and such that :math:`F(0)=H(0)=1`.
* A parameter :code:`s`.
* A wanted precision :code:`nb_decimals`, given as a number of decimal digits.

We have access to the lattice invariant classes, as per the preceding paragraph. For each of these classes :math:`(\mathcal{A})`, we compute 

.. math::
   \prod_{p\in\mathcal{A}}\frac{F(1/p^s)}{H(1/p^s)}.

There is a condition for this product to converge absolutely: on writing :math:`F(x)-H(x)=x^\Delta T(x)` for a :math:`\Delta\ge1` and a polynomial :math:`T(x)`, we need that :math:`\Delta s >1`. We assume this condition to hold.

.. code-block:: default
     
     from euler_product.lattice_invariant_euler_products import get_euler_products
     get_euler_products(q, s, F(x) , H(x), nb_decimals, big_p = 300, verbose = 0, with_Latex = 0, digits_offset = 10)

answers a couple whose first component is the tuple of the lattice invariant classes :math:`(\mathcal{A})`, and second component is the tuple of the values :math:`\prod_{p\in\mathcal{A}}\frac{F(1/p^s)}{H(1/p^s)}`, for example

.. code-block:: default
     
     from euler_product.lattice_invariant_euler_products import get_euler_products
     result = get_euler_products(5, 1, 1-x^2 , 1+x^3, 100, 300, 0)
     
     result[0][0]
     frozenset({1})
     result[1][0]
     (0.9884028950453419692925625250954713121182210521345380891771586345550561301333511982564965807673436742857698303688419181730105231677449, 0.9884028950453419692925625250954713121182210521345380891771586345550561301333511982564965807673437490090286957966947966907374203853849),

which means that

.. math::
   0.&9884028950453419692925625250954713121182210521345380891771586345550561301333511982564965807673436742857698303688419181730105231677449 
   
   &\le \prod_{p\equiv 1[5]} \frac{1-1/p^2}{1+1/p^3}
   
   &\le 0.9884028950453419692925625250954713121182210521345380891771586345550561301333511982564965807673437490090286957966947966907374203853849

With :code:`verbose = 1` or :code:`verbose = 2`, the results are more explicitly written.

To compute the specific quantities :math:`\prod_{p\in\mathcal{A}}(1-1/p^s)^{-1}` where the rational fraction is thus fixed, we have a shortcut:

.. code:: default

   from euler_product.lattice_invariant_euler_products import get_vs
   get_vs(q, s, nb_decimals=100, big_p=100, verbose=2, with_laTeX=0, digits_offset=10)

The output is similar to the one of :code:`get_euler_products`, with the same effect of the available parameters. However, there is the additional possible value :code:`verbose = -1`. In that case the output takes the shape

.. code:: default

   [big_p, phi_q, r, nb_invariant_class, big_m, time_end - time_start, difference]

which is rather explicit. The parameter :code:`big_m` is introduced in the reference paper and :code:`r` is the number of values of :math:`m`, as per Eq. (5) of the reference paper, that are being used. The timing is given in seconds, and :code:`difference` is an approximation of the number of correct decimals given.

- Auxiliaries

We finally provide two auxiliary functions.

.. code:: default

   from euler_product.lattice_invariant_euler_products import table_performance
   table_performance(min_q, max_q, nb_decimals = 100, big_p = 300)

This gives some timing info for :code:`get_vs(q, 2, nb_decimals, big_p, -1)`. The output has a LaTeX format of an array, the columns being :code:`q`, :code:`phi_q`, :code:`nb_prime_factors_phi_q`, :code:`r`, :code:`nb_invariant_class`, :code:`big_m` and finally :code:`time_end - time_start` in seconds / 10. The meanings are the same as in :code:`get_vs`.

.. code:: default

   from euler_product.lattice_invariant_euler_products import get_vs_checker
   get_vs_checker(q, s, borne = 10000):

This is a simple sanity check. The output :code:`get_vs` displays the Euler products computed by :code:`get_vs`, except that these products are only approximated by the truncated Euler product up to :code:`borne`.
