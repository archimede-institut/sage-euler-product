[![Test](https://github.com/archimede-institut/sage-euler-product/actions/workflows/tests.yml/badge.svg)](https://github.com/archimede-institut/sage-euler-product/actions/workflows/tests.yml)
[![documentation](https://github.com/archimede-institut/sage-euler-product/actions/workflows/manual.yml/badge.svg)](https://github.com/archimede-institut/sage-euler-product/actions/workflows/manual.yml)
[![Lint](https://github.com/archimede-institut/sage-euler-product/actions/workflows/lint.yml/badge.svg)](https://github.com/archimede-institut/sage-euler-product/actions/workflows/lint.yml)
[ ![pages-build-deployment](https://github.com/archimede-institut/sage-euler-product/actions/workflows/pages/pages-build-deployment/badge.svg)](https://github.com/archimede-institut/sage-euler-product/actions/workflows/pages/pages-build-deployment.yml)

  
# Euler Product for SageMath


 **Computing Lattice Invariant Euler Products** 

The **sage-euler-product** package for SageMath adds functionality related to Number Theory. It is based on SageMath <https://www.sagemath.org>_ and relies heavily on:

- gmp or mpir for arbitrary precision arithmetic
- PARI/GP for number field computations


## Prerequisites

Installing sage-euler-product requires a working Sage installation. 

## Installation

The module is distributed on PyPI and is easily installed through the Python package manager pip. If you downloaded a binary from the SageMath website (including the Cygwin version running on Windows) or compiled from source, run the following command:

    $ sage -pip install sage-euler-product [--user]

The --user option is optional and allows to install the module in your user space (and does not require administrator rights).

If you use Debian or Ubuntu and you installed Sage through the operating system's package manager (that is, the package sagemath), run these two commands:

    $ source /usr/share/sagemath/bin/sage-env
    $ pip install sage-euler-product --user

If you use Arch Linux, you need to install from source (see next section).

Install and use source version
This section provides detailed instructions on how to download, modify and install the development version of  **sage-euler-product**. In all commands,

PIP has to be replaced by either pip, pip2, or sage -pip
PYTHON has to be replaced by either python, python2 or sage -python
If you are an Arch Linux user with the sagemath package installed, use PIP=pip2 and PYTHON=python2. If you downloaded SageMath as a tarball or installed it from source use PIP='sage -pip' and PYTHON='sage -python'.

You can install the latest development version in one line with:

    $ PIP install git+https://github.com/archimede-institut/sage-euler-product [--user]

As before, the --user option is optional and when specified will install the module in your user space.

You can also perform a two stage installation that will allow you to modify the source code. The first step is to clone the repository:

    $ git clone https://github.com/archimede-institut/sage-euler-product

The above command creates a repository sage-euler-product with the source code, documentation and miscellaneous files. You can then change to the directory thus created and install the surface dynamics module with:

    $ cd sage-euler-product
    $ PIP install . [--user]
	
Do not forget the . that refers to the current directory.

When you don't want to install the package or you are testing some modifications to the source code, a more convenient way of using  **sage-euler-prodct** is to do everything locally. To do so, you need to compile the module in place via:

	$ PYTHON setup.py build_ext --inplace
	
Once done, you can import the sage-euler-product module. To check that you are actually using the right module (i.e. the local one) you can do in a SageMath session:

	sage: import euler_product
	sage: euler_product.__path__        # random
	['/home/you/sage-euler-product/euler_product/']

The result of the command must correspond to the path of the repository created by the command git clone given above. The compilation step PYTHON setup.py build_ext has to be redone each time you modify a C or Cython source file (i.e. with .c, .h, .pxd or .pyx extension). In other words, it is not needed if you only modify or create Python files (i.e. .py files).

If you wish to install your custom version of sage-euler-product just use PIP as indicated before.

## Documentation

complete module documentation: https://archimede-institut.github.io/sage-euler-product/

## Check

After installing  **sage-euler-product**, check that it works by launching Sage and typing the following commands. You should get the same output as below.

	sage: from euler_product.all import *
	sage: from euler_product.lattice_invariant_euler_produ import get_euler_products
	sage: get_euler_products(3, 1, 1-x^2,1, 100)
	Computing the structural invariants ...  done.
	We have Delta  = 2 and beta = 2
	We use big_m = 310 , big_p = 300 and working prec = 653 .
	Computing the finite products for p < 300 ...  done.
	Computing C_A(K, m, F/H) ... -------------------
	For p+3ZZ in  frozenset({1})
	For F(x) = -x^2 + 1
	and H(x) = 1
	the product of F(1/p)/H(1/p) is between
	0.9671040753637981066150556834173635260473412207450092130719978569438733967843271277395717230016746853806050215621235810749643636399725665325875376146914709362753787689855429317947529895445140974344
	and
	0.9671040753637981066150556834173635260473412207450092130719978569438733967843271277395717230016746853806050215621235810749643636399725665325875376146914709362753787689855429317947529895445140974475
	(Obtained:  193  correct decimal digits)
	-------------------
	For p+3ZZ in  frozenset({2})
	For F(x) = -x^2 + 1
	and H(x) = 1
	the product of F(1/p)/H(1/p) is between
	0.7071813747951674302088659938984504109243584468119496848353517677901518159831128643782536704398941052120208041311403202957250160794697319584608281454011743387515885835706146696365506658500107821107
	and
	0.7071813747951674302088659938984504109243584468119496848353517677901518159831128643782536704398941052120208041311403202957250160794697319584608281454011743387515885835706146696365506658500107821228
	(Obtained:  193  correct decimal digits)
	Time taken:  1.920718120993115 seconds.
	((frozenset({1}), frozenset({2})),
	 ((0.9671040753637981066150556834173635260473412207450092130719978569438733967843271277395717230016746853806050215621235810749643636399725665325875376146914709362753787689855429317947529895445140974344,
	   0.9671040753637981066150556834173635260473412207450092130719978569438733967843271277395717230016746853806050215621235810749643636399725665325875376146914709362753787689855429317947529895445140974475),
	  (0.7071813747951674302088659938984504109243584468119496848353517677901518159831128643782536704398941052120208041311403202957250160794697319584608281454011743387515885835706146696365506658500107821107,
	   0.7071813747951674302088659938984504109243584468119496848353517677901518159831128643782536704398941052120208041311403202957250160794697319584608281454011743387515885835706146696365506658500107821228)))






https://github.com/archimede-institut/sage-euler-product
Assuming you have the program git on your computer, you can install the development version with the command:

	$ sage -pip install git+https://github.com/archimede-institut/sage-euler-product [--user]





## Authors

Olivier Ramar\'e: see https://ramare-olivier.github.io/Maths/mcom3630.pdf for complete Mathematical references

Dominique Benielli: maintainerDeveloppement Cell, Institut Archimède Aix-Marseille Université


## How to cite this project

If you have used this project for please cite us as described on our zenodo site.

## Versions

The first release of sage-euler-product will appear soon as a sagemath spkg.

. 
