// -*- C++ -*-
/*! 
 *  \brief Construct the traceless Field Strength Tensor.
 *  If there are questions or bugs, blame David Brantley.
 */

#ifndef __tr_less_fields_h__
#define __tr_less_fields_h__

#include "chromabase.h"

namespace Chroma 
{
  /* Construct the field strength tensor using the definition in Zweig et. al. (1983)*/
  // This program reads in a gauge field, and constructs the resulting traceless
  // field strength tensor. The output is then returned as a 2d array of lattice color
  // matrices.
  //         **** COMPUTES BOTH THE UPPER AND LOWER TRIANGULAR PORTIONS ****

  multi1d<LatticeColorMatrix> trLessFieldST(const multi1d<LatticeColorMatrix>& u);
  
}

#endif
