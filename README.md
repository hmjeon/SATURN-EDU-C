<a name="top"></a>

# SATURN-DISP4-pub

[![License](https://img.shields.io/badge/license-GNU%20GeneraL%20Public%20License%20v3,%20GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html) <br>
[![# of Language](https://img.shields.io/github/languages/count/hmjeon/SATURN-DISP4-pub.svg)](https://github.com/hmjeon/SATURN-DISP4-pub)
[![Language](https://img.shields.io/github/languages/top/hmjeon/SATURN-DISP4-pub.svg)](https://github.com/hmjeon/SATURN-DISP4-pub)
[![Tag](https://img.shields.io/github/tag/hmjeon/SATURN-DISP4-pub.svg)](https://github.com/hmjeon/SATURN-DISP4-pub/tags)
[![Last Commit](https://img.shields.io/github/last-commit/hmjeon/SATURN-DISP4-pub.svg)](https://github.com/hmjeon/SATURN-DISP4-pub)

**SATURN-DISP4-pub** is the part of the SATURN package and open-source written C/C++ that enables the finite element analysis under the plane stress condition with 4-node element.</br>

## Clone the repository
```git clone https://github.com/hmjeon/SATURN-DISP4-pub.git```</br>

---

[Requirements](#Requirements-to-compile-from-source) | [Features](#Features) | [Element types](#Element-types) | [Copyrights](#copyrights)

---

## Requirements to compile from source
* C/C++ compiler: Visual Studio or gcc compiler - SATURN-DISP4-pub
* Fortran compiler: [Intel Parallel Studio XE](https://software.intel.com/en-us/fortran-compilers) 2016, 2017 or 2018 or [PGI Fortran](https://www.pgroup.com/)
* We provide MakeFile which is a simple way to organize code compilation.

## Features
* Problem generator for several benchmark problems
* Supported for two storage methods; Skyline and [CRS](https://en.wikipedia.org/wiki/Sparse_matrix) (Compressed Row Storage) formats
* High performance linear equation solver; [PARDISO](https://www.pardiso-project.org/), [MUMPS](mumps.enseeiht.fr), [cuSPARSE](https://developer.nvidia.com/cusparse)
* Post-processing for [TecPlot](https://www.tecplot.com/) and [ParaView](https://www.paraview.org/)</br>
* MATLAB and Python scripts for plotting the element stress</br>
* Free and open source ([GNU General Public License, version 3.0](https://www.gnu.org/licenses/gpl-3.0.en.html/))</br>

## Element types
* MURCURY - MITC3, MITC4, MITC6, MITC9, MITC3P, MITC4P, MITC3E, MITC3PE
* SATURN - DISP3, DISP4, DISP4, DISP9, DISP3E
* JUPYTER - TRUSS2D, TRUSS3D, EBEAM2D, EBEAM3D, FRAME2D, FRAME3D
* All elements for linear, geometrical nonlinear , hyperelastic material, and GPU-enabled analyses
* These packages currently are private but it is avaiable upon request

## Copyrights
*Author*
* Dr. Hyungmin Jun - [E-mail](mailto:hyungminjun@outlook.com)</br>
* Web - [link](https://hyungminjun.com)</br>

*License*
* SATURN-DISP4-pub is an open-source software distributed under the [GPL license, version 3](https://www.gnu.org/licenses/gpl-3.0.en.html/)</br>
* © 2019 Hyungmin Jun ALL RIGHTS RESERVED</br>

Anyone who is interest to use, to develop or to contribute to SATURN, MURCURY and JUPYTER is always welcome!

Go to [Top](#top)
