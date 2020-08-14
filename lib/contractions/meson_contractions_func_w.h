/*! 
 *  Functions that do baryon spin and color contractions.
 *  Color is nested inside of spin.
 *  Authors:
 *  Arjun Gambhir
 */

#ifndef __meson_contractions_func_h__
#define __meson_contractions_func_h__

#include "../momentum/lalibe_sftmom.h"

namespace Chroma
{

  void anti_quark_FH_prop(LatticePropagator & FH_quark_prop,
                          LatticePropagator & FH_antiquark_prop,
                          std::string& cur);



  void FH_I_one_Iz_pm_one_contract(LatticePropagator & quark_1,
                             LatticePropagator & quark_2,
                             LatticeComplex & contracted,
                             bool& is_FH_antiquark,
                             std::string& cur);

  void I_one_Iz_pm_one_contract(LatticePropagator & quark_1,
                             LatticePropagator & quark_2,
                             LatticeComplex & contracted);

  void write_correlator(bool full_correlator,
			std::string meson_name,
#ifdef BUILD_HDF5
			std::string path,
			HDF5Writer & h5writer,
			HDF5Base::writemode & h5mode,
#endif
			int t_0,
			int Nt,
			multi1d<int> & source_coords,
			LalibeSftMom & FT,
			LatticeComplex & meson);

}

#endif
