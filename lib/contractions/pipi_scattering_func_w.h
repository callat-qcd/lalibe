#ifndef __pipi_scattering_h__
#define __pipi_scattering_h__

namespace Chroma
{

	//void pion_correlator(LatticeComplex& correlator, const LatticePropagator& quark_prop_1, const LatticePropagator& quark_prop_2);

	void pipi_correlator(multi3d<DComplex>& correlator, const LatticePropagator& quark_prop_1, const LatticePropagator& quark_prop_2, const LatticePropagator& quark_prop_3, const LatticePropagator& quark_prop_4, const multi3d<int>& mom_list, const multi2d<int>& origin_list, const int t0, const int j_decay);

} // End namespace Chroma

#endif

