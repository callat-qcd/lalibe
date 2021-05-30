/*! 
 *  Functions that do baryon spin and color contractions.
 *  Color is nested inside of spin.
 *  Authors:
 *  Ben Hoerz
 *  Andre Walker-Loud
 */

#include "spin_basis.h"
#include "util/ferm/diractodr.h"


namespace Chroma 
{ 

  void rotate_to_Dirac_Basis(LatticePropagator & quark_to_be_rotated)
  {
    //I am lazy so I copy Robert to rotate this stuff...
    SpinMatrix U = DiracToDRMat();
    quark_to_be_rotated = adj(U)*quark_to_be_rotated*U ;

  }

  void rotate_from_Dirac_Basis(LatticePropagator & quark_to_be_rotated)
  {
    SpinMatrix U = DiracToDRMat();
    quark_to_be_rotated = U*quark_to_be_rotated*adj(U);
  }

} // End namespace Chroma
