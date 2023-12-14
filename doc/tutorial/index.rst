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

.. code-block:: default
     
        sage: from euler_product.utils_euler_product import LatticeInvariant
        sage: LatticeInvariant(15)
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
    


