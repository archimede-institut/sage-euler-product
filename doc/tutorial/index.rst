.. _index:


Tutorial Euler Product
======================

- Introduction and principles

The main object of this software is to compute in a very fast manner Euler products with rational local factors over primes in some collections of arithmetic progressions modulo some fixed modulus :math:`q`.
Here are two examples  :math:`E_1=\prod_{p\equiv 3,5[7]}(1-1/p^2)`, where the product is taken over every prime congruent to 3 or 5 modulo 7, and :math:`E_2 = \prod_{p\equiv 1[8]}\bigl(1-\frac{4}{p}\biggr)\bigl(\frac{p+1}{p-1}\bigr)^2`, the so-called Shanks's constant, where the product is taken over every prime congruent to 1 modulo 8.

Euler products over rational functions as :math:`E_2` is inferred from simpler Euler products of the shape 

.. math::
\prod_{p\in\mathcal{A}\mod q}(1-1/p^s),

where

 * :math:`\mathcal{A}` is some subset of :math:`G=(\mathbb{Z}/q\mathbb{Z})^\times`. The subset :math:`\mathcal{A}` has to be the union of "lattice invariant classes", as described below.
 * :math:`q` is a positive integer, the so-called "modulus". We have :math:`q=7` for :math:`E_1`.
 * :math:`s` is a real parameter that strictly larger than 1. A typical choice is :math:`s=2`.


Euler Product over every primes
-------------------------------

.. code-block:: default
     
     from euler_product.utils_euler_product import ComponentStructure
     structure = ComponentStructure(3)
     from euler_product.lattice_invariant_euler_products import get_euler_products
     get_euler_products(1, 1, 1-x^2, 1+x^3, 100, verbose=0)



Euler Product over primes in arithmetic progression
---------------------------------------------------


- Definition Lattice Invariant

We subdivide the multiplicative group :math:`G=(\mathbb{Z}/q\mathbb{Z})^\times` in classes called  `LatticeInvariantClasses`
When :math:`q = 15` these classes are obtained by

.. code-block:: default

       LatticeInvariant(15)[1]
       (frozenset({1}), frozenset({4}), frozenset({11}), frozenset({14}), frozenset({8, 2}), frozenset({13, 7}))


.. code-block:: default
     
        from euler_product.utils_euler_product import LatticeInvariant
        LatticeInvariant(15)
        ((frozenset({1}), frozenset({1, 4}), frozenset({1, 11}), frozenset({1, 14}), frozenset({8, 1, 2, 4}), frozenset({1, 4, 13, 7})), 
         (frozenset({1}), frozenset({4}), frozenset({11}), frozenset({14}), frozenset({8, 2}), frozenset({13, 7})))


Two points are in the same class if and only if they generate the same subgroup modulo :math:`q`



