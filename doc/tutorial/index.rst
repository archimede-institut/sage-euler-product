.. _index:


Tutorial Euler Product
======================


Euler Product over every primes
-------------------------------

.. code-block:: default
     
    import euler_product.util


Euler Product over primes in arithmetic progression
---------------------------------------------------


- Definition Lattice Invariant

We subdivide the multiplicative group :math:`G=(\mathbb{Z}/q\mathbb{Z})^\times` in classes called  `LatticeInvariantClasses`
When :math:`q = 15` these classes are obtained by

.. code-block:: default

       LatticeInvariant(15)[1]
       (frozenset({1}), frozenset({4}), frozenset({11}), frozenset({14}), frozenset({8, 2}), frozenset({13, 7}))


.. code-block:: default
     
        sage: from euler_product.utils_euler_product import LatticeInvariant
        sage: LatticeInvariant(15)
        ((frozenset({1}), frozenset({1, 4}), frozenset({1, 11}), frozenset({1, 14}), frozenset({8, 1, 2, 4}), frozenset({1, 4, 13, 7})), 
         (frozenset({1}), frozenset({4}), frozenset({11}), frozenset({14}), frozenset({8, 2}), frozenset({13, 7})))


  Two points are in the same class if and only if they generate the same subgroup modulo `q`


