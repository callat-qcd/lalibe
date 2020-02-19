#ifndef __pipi_scattering_h__
#define __pipi_scattering_h__

namespace Chroma
{

	void pion_correlator(multi1d<DComplex>& correlator, const LatticePropagator& quark_prop_1, const LatticePropagator& quark_prop_2, const multi2d<int>& mom_list, const int t0, const int j_decay);


void pipi_correlator(multi1d<DComplex>& correlator, const LatticePropagator& quark_prop_1, const LatticePropagator& quark_prop_2, const LatticePropagator& quark_prop_3, const LatticePropagator& quark_prop_4, const multi3d<int>& mom_list, const int t0, const int j_decay);

} // End namespace Chroma

#endif

