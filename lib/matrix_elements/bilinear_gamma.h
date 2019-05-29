// -*- C++ -*-
// $Id: bilinear_gamma.h,  2018-03-1 16:13:30 Henry M
/*! \file
 *  \brief Create a bilinear LatticeSpinMatrix
 */

#ifndef __bilinear_gamma_h__
#define __bilinear_gamma_h__

namespace Chroma 
{ 
  //! Construct the "bilinear"
  /*!
   * \ingroup bilinear
   *
   * Arguments:
   *  \param b		 (bilinear ID)
   */

  void Bilinear_Gamma(std::string present_current, LatticePropagator& out_quark_src, LatticePropagator& quark_src, const multi1d<LatticeColorMatrix>& u);
}

#endif
