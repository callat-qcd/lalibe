LALIBE
(Lattice Livermore Berkeley)

**Developers:**  
Arjun Gambhir  
David Brantley  
Jason Chang  
Henry Monge-Camacho  
Pavlos Vranas  
André Walker-Loud

[Copywrite](#copyright-notice)

This software links against Chroma (https://github.com/JeffersonLab/chroma) and performs various computations related to two-point functions and three point functions.  In particular
* [Feynman-Hellmann propagators and correlation functions](#feynman-hellmann)
* [three point functions](#three-point-functions)
* [HDF5 for correlation functions](#hdf5)
* [Hierarchical Probing](#hierarchical-probing)


## Feynman-Hellmann
This code constructs the Feynman-Hellmann propagators and correlation functions as described in [arXiv:1612.06963](https://arxiv.org/abs/1612.06963) and used to compute the nucleon axial charge in [arXiv:1805.12130](https://arxiv.org/abs/1805.12130).  

If this routine is used, please acknowledge the two publications
- *On the Feynman-Hellmann Theorem in Quantum Field Theory and the Calculation of Matrix Elements*  
  C. Bouchard, C. C. Chang, T. Kurth, K. Orginos and A. Walker-Loud  
  [Phys. Rev. D96 (2017) no.1, 014504](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.96.014504) [[arXiv:1612.06963](https://arxiv.org/abs/1612.06963)]
- *Simulating the weak death of the neutron in a femtoscale universe with near-Exascale computing*  
  E. Berkowitz, M.A. Clark, A. Gambhir, K. McElvain, A. Nicholson, E. Rinaldi, P. Vranas, A. Walker-Loud, C. C. Chang, B. Joó, T. Kurth, K. Orginos  
  [Supercomputing 2018 55:1-9](http://dl.acm.org/citation.cfm?id=3291656.3291730) [[arXiv:1810.01609](https://arxiv.org/abs/1810.01609)]


## Three Point Functions
This code constructs three point correlation functions of the proton following the conventional fixed source-sink time separation method.  The sequential source routine is heavily inspired by that in Chroma.  We modified the routine: 
* To use the local operators in the Dirac-Pauli basis ([arXiv:hep-lat/0506029](https://arxiv.org/abs/hep-lat/0506029))
* We re-wrote the sequential sink routine to support a time-dense sink, such that a few could be added together to make a coherent sequential sink ([arXiv:1001.3620](https://arxiv.org/abs/1001.3620)) in a more economical way, allowing for one to select specific times slices for each src-propagator rather than reconstructing the sink each time on the fly, minimizing the time spent constructing these coherent sinks.  

The formfactor routine is an edited version of that in Chroma allowing for
* Writing correlation functions in [HDF5](#hdf5)
* Writing full 4-dimensional correlation functions
* Improved support for the momentum at the current insertion

If this routine is used, please acknowledge the publication
* *coming soon!*


## HDF5
T. Kurth originally added support for HDF5 I/O routines in QDP++, [arXiv:1501.06992](https://arxiv.org/abs/1501.06992) with help from Bálint Joó.  (At the same time, Andrew Pochinsky added similar support in Qlua.)  Arjun extended this to inline measurements that can write propagators, gauge fields, etc. He also converted writing of any correlation function in lalibe to have an hdf5 option.

If this routine is used, please acknowledge the publications
* *High-Performance I/O: HDF5 for Lattice QCD*  
  T. Kurth, A. Pochinsky, A. Sarje, S. Syritsyn, A. Walker-Loud  
  [PoS LATTICE2014 (2015) 045](https://pos.sissa.it/214/045) [[arXiv:1501.06992](https://arxiv.org/abs/1501.06992)]
* *coming someday to Chroma!*

## Monte Carlo Method for Estimating Traces
This is a collection of measurements that generate random ZN sources and solve propagators off of them to be used for either computing disconnected diagrams or to create "stochastic propagators". The stochastic propagators can be used for the Stochastic Feynman-Hellmann Method or for stochastically computing three-point functions with the sequential source method. 

If these random sources are used, cite:

M. F. Hutchinson, *A stochastic estimator of the trace of the influence matrix for Laplacian
smoothing splines*, J. Commun. Statist. Simula., 19 (1990), pp. 433–450.

Additionally, if they are used for the Stochastic Feynman-Hellman method, cite:

A. S. Gambhir, E. Berkowitz,D. Brantley, C. C. Chang, M. A. Clark, T. Kurth, C. Monahan,
A. Nicholson, P. Vranas, A. Walker-Loud, *The Stochastic Feynman-Hellmann Method*, PoS
LATTICE 2018, 126 (2018) [arXiv:1905.03355 [hep-lat]].

## Hierarchical Probing
Andreas Stathopoulos 
binaryRecursiveColoring.cc
This file was originally written by Andreas Stathopoulos and adopted for use in lalibe by Arjun Singh Gambhir.

Copyright (c) 2015, College of William & Mary
All rights reserved.

Andreas Stathopoulos, Jesse Laeuchli, and Kostas Orginos, *Hierarchical probing for estimating the
trace of the matrix inverse on toroidal lattices*, 2013, SIAM J. Sci. Comput., 35(5), S299–S322, arXiv:1302.4018 [hep-lat]

Hierarchical Probing is an advanced method for computing the trace of a matrix inverse. It can be used for disconnected diagrams, the Stochastic Feynman-Hellmann Method, and stochastic three-point functions. For intermediate to heavy pion masses, it gives substantially lower statistical uncertainty than random noise techniques. If hierarchical probing is used, cite the paper referenced above. 


## Copyright Notice

lalibe Copyright (c) 2018, The Regents of the University of California, through Lawrence Berkeley National Laboratory, and Lawrence Livermore National Security, LLC, for the operation of Lawrence Livermore National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved.
If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Intellectual Property Office at IPO@lbl.gov.
NOTICE. This Software was developed under funding from the U.S. Department of Energy and the U.S. Government consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, prepare derivative works, and perform publicly and display publicly, and to permit other to do so.
