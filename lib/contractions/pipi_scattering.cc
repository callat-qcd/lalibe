/*!
 *  Functions that perform the contraction of quarks into the two-pion system.
 *  Based on qqbar_w.cc in Chroma
 *  Color is nested inside of spin.
 *  Authors:
 *  Haobo Yan
 *  Andre Walker-Loud
 *  Ben Horz
 */

#include "chromabase.h"
#include "pipi_scattering.h"
#include "../../../chroma/lib/meas/hadron/qqbar_w.cc"

namespace Chroma
{

	void pion_correlator(multi1d<DComplex>& correlator, const LatticePropagator& quark_prop_1, const LatticePropagator& quark_prop_2, const multi2d<int>& mom_list, const int t0, const int j_decay)
	{
		int G5 = Ns * Ns - 1;
		multi2d<DPropagator> single_pion; // function of (p, t)

		// Construct SftMom object using the argument 'mom_list' with the following constructor:
		// Constructor about origin, with a list of momenta 
		// SftMom(const multi2d<int> & moms, int j_decay = -1);
		SftMom phases(mom_list, j_decay);
		int Nt = phases.numSubsets(); // Length of lattice in j_decay direction by definition, but I understand it as the lattice number in time
		
		// Implement the momentum projection part using qqbar
		compute_qqbar(single_pion, quark_prop_1, quark_prop_2, phases, t0);

		// Loop over momenta
		for (int mom_num = 0; mom_num < phases.numMom(); ++mom_num)
			for (int t = 0; t < Nt; ++t)
			{
				// Contract colors and spins with gamma matrices
				correlator[mom_num] = correlator[mom_num] + trace(single_pion[mom_num][t] * Gamma(G5)); // I'm not sure if I can use += for this type
				// Fix the origin

			}
	}


	void pipi_correlator(multi1d<DComplex>& correlator, const LatticePropagator& quark_prop_1, const LatticePropagator& quark_prop_2, const LatticePropagator& quark_prop_3, const LatticePropagator& quark_prop_4, const multi3d<int>& mom_list, const int t0, const int j_decay)
	{
		int G5 = Ns * Ns - 1;
		multi2d<DPropagator> pipi1, pipi2;

		SftMom phases1(mom_list[0], j_decay);
		SftMom phases2(mom_list[1], j_decay);
		int Nt = phases1.numSubsets();

		compute_qqbar(pipi1, quark_prop_1, quark_prop_2, phases1, t0);
		compute_qqbar(pipi2, quark_prop_3, quark_prop_4, phases2, t0);

		multi1d<DComplex> correlator1(0), correlator2(0), correlator3(0), correlator4(0);

		for (int mom_num = 0; mom_num < phases1.numMom(); ++mom_num)
			for (int t = 0; t < Nt; ++t)
			{
				// Combination 1
				correlator1[mom_num] = correlator1[mom_num] + trace(pipi1[mom_num][t] * Gamma(G5)) * trace(pipi2[mom_num][t] * Gamma(G5));

				// Combination 2
				DColorMatrix cm_A, cm_B;
				DPropagator A, B;
				A = pipi1[mom_num][t] * Gamma(G5);
				B = pipi2[mom_num][t] * Gamma(G5);

				for (int s1 = 0; s1 < Ns; ++s1)
					for (int s2 = 0; s2 < Ns; ++s2)
					{
						cm_A = peekSpin(A, s1, s2);
						cm_B = peekSpin(B, s2, s1);
						for (int c1 = 0; c1 < Nc; ++c1)
							for (int c2 = 0; c2 < Nc; ++c1)
								correlator2[mom_num] = correlator2[mom_num] - peekColor(cm_A, c1, c2) * peekColor(cm_B, c1, c2); // I use - because of the minus sign picked up by Grossman numbers
					}

				// Combination 3
				A = pipi1[mom_num][t] * Gamma(G5);
				B = pipi2[mom_num][t] * Gamma(G5);

				for (int s1 = 0; s1 < Ns; ++s1)
					for (int s2 = 0; s2 < Ns; ++s2)
					{
						cm_A = peekSpin(A, s1, s2);
						cm_B = peekSpin(B, s2, s1);
						for (int c1 = 0; c1 < Nc; ++c1)
							for (int c2 = 0; c2 < Nc; ++c1)
								correlator3[mom_num] = correlator3[mom_num] + peekColor(cm_A, c1, c2) * peekColor(cm_B, c1, c2);
					}

				// Combination 4
				correlator4[mom_num] = correlator4[mom_num] + trace(pipi1[mom_num][t] * Gamma(G5)) * trace(pipi2[mom_num][t] * Gamma(G5));

				// Total correlator
				correlator[mom_num] = correlator1[mom_num] + correlator2[mom_num] + correlator3[mom_num] + correlator4[mom_num];
				// Fix the origin

			}


	}


} // End namespace Chroma
