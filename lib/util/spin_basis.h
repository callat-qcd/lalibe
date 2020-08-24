/*! 
 */

#include "chromabase.h"

#ifndef __spin_basis_h__
#define __spin_basis_h__

namespace Chroma 
{
  void rotate_to_Dirac_Basis(LatticePropagator & quark_to_be_rotated);

  void rotate_from_Dirac_Basis(LatticePropagator & quark_to_be_rotated);
}

#endif
