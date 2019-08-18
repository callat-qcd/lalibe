/*! 
 *  Functions that do baryon spin and color contractions.
 *  Color is nested inside of spin.
 *  Authors:
 *  Arjun Gambhir
 */

#ifndef __proton_contractions_func_h__
#define __proton_contractions_func_h__

#include "../momentum/lalibe_sftmom.h"

namespace Chroma 
{
namespace LegacyProton {

  void rotate_to_Dirac_Basis(LatticePropagator & quark_to_be_rotated);
  
  void get_spin_wavefunctions(multi2d<int> & src_spins,       
			multi2d<int> & snk_spins,       
			multi1d<Real> & src_weights,  
			multi1d<Real> & snk_weights, 
			std::string baryon_name,
			std::string spin,
			int parity);
 
  void color_contraction(LatticeColorMatrix & quark_1,
			 LatticeColorMatrix & quark_2,
			 LatticeColorMatrix & quark_3,
			 LatticeComplex & spin_color_contracted_thing);

  void spin_contraction(LatticePropagator & quark_1,
                        LatticePropagator & quark_2,
                        LatticePropagator & quark_3,
			multi2d<int> & src_spins,       //Indices of length (N_src, 3)
			multi2d<int> & snk_spins,       //Indices of length (N_snk, 3)
			multi1d<Real> & src_weights,    //Index of length (N_src) 
			multi1d<Real> & snk_weights,    //Index of length (N_snk)
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
}

#endif
