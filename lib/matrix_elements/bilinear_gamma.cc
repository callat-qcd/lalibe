// -*- C++ -*-
// $Id: bilinear_gamma.h,  2018-03-1 16:13:30 Henry M
/*! \file
 *  \brief Create a bilinear LatticeSpinMatrix
 */

#include "chromabase.h"
//This is needed for the chromomag operator.
#include "chromomag_seqsource_w.h"


namespace Chroma 
{ 
  //! Construct the "bilinear"
  /*!
   * \ingroup bilinear
   *
   * Arguments:
   *  \param b		 (bilinear ID)
   */
  void Bilinear_Gamma(std::string present_current, LatticePropagator& out_quark_src, LatticePropagator& quark_src, const multi1d<LatticeColorMatrix>& u){
    START_CODE();
	if (present_current == "S")
	   out_quark_src =  Gamma(0)* quark_src;
    else if (present_current == "P")
       out_quark_src =  Gamma(15)* quark_src;
    else if (present_current == "A1")
       out_quark_src =  Gamma(14)* quark_src;
    else if (present_current == "A2")
       out_quark_src =  Gamma(13)* -quark_src;
    else if (present_current == "A3")
       out_quark_src =  Gamma(11)* quark_src;
    else if (present_current == "A4")
       out_quark_src =  Gamma(7)* -quark_src;
    else if (present_current == "V1")
       out_quark_src =  Gamma(1)* quark_src;
    else if (present_current == "V2")
       out_quark_src =  Gamma(2)* quark_src;
    else if (present_current == "V3")
       out_quark_src =  Gamma(4)* quark_src;
    else if (present_current == "V4")
       out_quark_src =  Gamma(8)* quark_src;
    else if (present_current == "T12")
       out_quark_src =  Gamma(3)* quark_src;
    else if (present_current == "T13")
       out_quark_src =  Gamma(5)* quark_src;
    else if (present_current == "T14")
       out_quark_src =  Gamma(9)* quark_src;
    else if (present_current == "T23")
       out_quark_src =  Gamma(6)* quark_src;
    else if (present_current == "T24")
       out_quark_src =  Gamma(10)* quark_src;
    else if (present_current == "T34")
       out_quark_src =  Gamma(12)* quark_src;
    else if (present_current == "CHROMO_MAG")
	  out_quark_src =  chromoMagneticSeqSource(quark_src,u);
    else
    {
		QDPIO::cerr << present_current << ": NOT DEFINED YET " << std::endl;
    	QDP_abort(1);
    }
    END_CODE();
  }

} // ENd namespace Chroma
