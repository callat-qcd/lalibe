#ifndef __pipi_scattering_h__
#define __pipi_scattering_h__

namespace Chroma
{

	namespace CorrelatorType
	{

		typedef std::tuple<int, int, int> momenta;
		typedef std::pair<momenta, momenta> momenta_pair;
		typedef std::map<const momenta_pair, multi1d<DComplex>> Correlator;

	} // namespace CorrelatorType

	void pipi_correlator(CorrelatorType::Correlator& correlator_out, const LatticePropagator& quark_prop_1, const LatticePropagator& quark_prop_2, const LatticePropagator& quark_prop_3, const LatticePropagator& quark_prop_4, const multi2d<int>& origin_list, const int p2max, const int ptot2max, const int t0, const int j_decay);
	void kpi_correlator(CorrelatorType::Correlator& correlator_out, const LatticePropagator& quark_prop_1, const LatticePropagator& quark_prop_2, const LatticePropagator& quark_prop_3, const LatticePropagator& quark_prop_4, const multi2d<int>& origin_list, const int p2max, const int ptot2max, const int t0, const int j_decay);

} // End namespace Chroma

#endif
