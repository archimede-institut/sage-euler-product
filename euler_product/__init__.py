from sage.misc.lazy_import import lazy_import

lazy_import('euler_product.utils_euler_product', 
            ['LatticeInvariant', 'ComponentStructure', 'laTeX_for_number', 
             'nb_common_digits'])

lazy_import('euler_product.lattice_invariant_euler_product', 
            ['get_euler_products', 'table_performance'])


