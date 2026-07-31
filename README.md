[![Linux Ubuntu Tests](https://github.com/MQCPack/mqcPack/actions/workflows/linux.yml/badge.svg)](https://github.com/MQCPack/mqcPack/actions/workflows/linux.yml)
[![MacOS Tests](https://github.com/MQCPack/mqcPack/actions/workflows/macos.yml/badge.svg)](https://github.com/MQCPack/mqcPack/actions/workflows/macos.yml)
<pre>
   **************************************************************************
   **************************************************************************
   **                                                                      **
   **   /$$      /$$ /$$$$$$  /$$$$$$ /$$$$$$$                  /$$        ** 
   **  | $$$    /$$$/$$__  $$/$$__  $| $$__  $$                | $$        **
   **  | $$$$  /$$$| $$  \ $| $$  \__| $$  \ $$/$$$$$$  /$$$$$$| $$   /$$  ** 
   **  | $$ $$/$$ $| $$  | $| $$     | $$$$$$$|____  $$/$$_____| $$  /$$/  **
   **  | $$  $$$| $| $$  | $| $$     | $$____/ /$$$$$$| $$     | $$$$$$/   **
   **  | $$\  $ | $| $$/$$ $| $$    $| $$     /$$__  $| $$     | $$_  $$   **
   **  | $$ \/  | $|  $$$$$$|  $$$$$$| $$    |  $$$$$$|  $$$$$$| $$ \  $$  **
   **  |__/     |__/\____ $$$\______/|__/     \_______/\_______|__/  \__/  **
   **          \__/                                                        **
   **                                                                      **
   **                 The Merced Quantum Chemistry Package                 **
   **                              (MQCPack)                               **
   **                            Version 26.7.3                            **
   **                             July 30, 2026                            **
   **                                                                      **
   **                                                                      **
   ** Written By:                                                          **
   **    Lee M. Thompson, Xianghai Sheng, Andrew D. Mahler, Dave Mullally  **
   **    and Hrant P. Hratchian                                            **
   **                                                                      **
   **                                                                      **
   **************************************************************************
   **************************************************************************
</pre>                                                                      

The MQC package is a library containing objects and procedures to streamline develoments 
of quantum chemical methodologies. MQC is built using object-orientated Fortran and has 
three key features.  

1. Algebra

MQC has a number of mathematical objects that are designed to hide away and automate 
manipulation of arrays. Arrays can be packed for efficient memory use and used for 
operations completely transparently to the user. Furthermore, there is no need to type 
arrays, as this can be manipulated on the fly. Arrays carry their own procedures, and 
use underlying lapack routines for efficiency. See the documentation for more 
information, but as an example: 

	scalar = vector%pop()   
	
removes the value from the end of a vector and returns it to scalar,  

	call matrix%eigensystem(metric,eigenvalues,rightvectors)  
	
computes the generalized eigenequation of the matrix given the metric matrix and returns 
the eigenvalues and the right eigenvectors (you can return the left eigenvectors as well 
if you wish, or just the left eigenvectors, or neither set of eigenvectors). MQC can 
deal with integer, real and complex arrays. 

2. Quantum chemistry 

MQC has objects that built over algebra objects designed specifically for quantum 
chemistry calculations. The most developed set of objects are the EST objects. This can 
work out whether you have a restricted, unrestricted or general spin framework and store 
the spin blocks as makes sense, it can automatically determine how to multiply different 
objects or carry out manipulations. Procedure interfaces allow you to use these objects 
just like normal matrices or take advantage of their higher-level organization.  

3. Interfacing with QM codes

A key feature of MQC is that it makes interacting with a range of computational 
chemistry codes entirely seamless. By having interfaces with different code APIs, MQC is 
able to fill the quantum chemistry objects and pass them back out in a manner that is 
completely transparent to the user. Thus, regardless regardles of what your favorite 
software is, your program will work without any modifications. 
 
# Getting Started

## Prerequisites
  + Required Software:
	 - autotools (only for developers)
	 - BLAS
	 - LAPACK library
         - gcc v8 or newer
         - gfortran v8 or newer, or pgfortran v19 or newer
  + Optional Software:
	 - gauopen
	 - Gaussian 16 or higher

## Installing

The interactive installer remains available:

```sh
./mqc_install
```

Enable MatrixFile support with the external GauOpen source directory and the
currently supported 8-byte integer ABI. The build compiles the unmodified
`qcmatrix.F` and `qcmatrixio.F` files directly with isolated legacy-source
flags. Configure identifies the supported public or frontier GauOpen API from
the selected sources:

```sh
./configure --prefix=/path/to/install \
  --with-gauopen=/path/to/gauopen \
  --with-gauopen-integer-bytes=8
make MQC_LAPACK="-llapack" MQC_BLAS="-lblas"
make check MQC_LAPACK="-llapack" MQC_BLAS="-lblas"
make install MQC_LAPACK="-llapack" MQC_BLAS="-lblas"
```

The default `--with-gauopen-api=auto` can be replaced by
`--with-gauopen-api=public` or `--with-gauopen-api=frontier` after verifying
the source interface. Public GauOpen supports numerical MatrixFile records but
does not provide character-record routines; MQCPack character FAF operations
therefore terminate with an explicit unsupported-feature diagnostic. Frontier
GauOpen provides the `TypeA` and character-buffer interface used for character
scalar and fixed-width-vector FAF records. To prevent positive `TypeA` records
from being silently treated as integers, the public profile also rejects
frontier version-3 FAF files.

Use `--without-gauopen` for a FormChk-only configuration; that path does not
need the GauOpen sources and produces a library that can link and use
the `MQC_Gaussian` formatted-checkpoint interface. MatrixFile types remain
available at compile time for shared source compatibility, but a MatrixFile
operation reports that GauOpen support was not configured and terminates.
A Git checkout may need
`autoreconf --install` once to generate `configure`; release source archives
normally include generated Autotools files.

# Licence
  This project is licensed under the Apache License - see the LICENSE file for details
