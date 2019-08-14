/*! 
 *  Functions that do baryon spin and color contractions.
 *  Authors:
 *  Arjun Gambhir
 */

#ifndef __baryon_contractions_func_h__
#define __baryon_contractions_func_h__

#include "../momentum/lalibe_sftmom.h"

#include <vector>
#include <tuple>
#include <map>
#include <utility>

namespace Chroma 
{

  void rotate_to_Dirac_Basis(LatticePropagator & quark_to_be_rotated);
 
  void color_contraction(const LatticeColorMatrix & quark_1,
			 const LatticeColorMatrix & quark_2,
			 const LatticeColorMatrix & quark_3,
			 LatticeComplex & spin_color_contracted_thing);

  std::tuple<char,char,char> get_flavor_code(const std::string& baryon_name);

  std::vector<std::string> get_spin_components(const std::string& baryon_name);

  void do_contraction(const LatticePropagator & quark_1,
		      const LatticePropagator & quark_2,
		      const LatticePropagator & quark_3,
		      const std::string& baryon_name,
		      const std::string& spin,
		      LatticeComplex & baryon_contracted_thing);

  void write_correlator(bool full_correlator,
			bool antiperiodic,
			std::string baryon_name,
			std::string spin,
#ifdef BUILD_HDF5
			std::string path,
			HDF5Writer & h5writer,
			HDF5Base::writemode & h5mode,
#endif
			int t_0,
			int Nt,
			multi1d<int> & source_coords,
			LalibeSftMom & FT,
			LatticeComplex & baryon);

}

#endif
