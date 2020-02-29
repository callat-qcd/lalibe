/*!
 *  Functions that perform the contraction of quarks into the two-pion system.
 *  Based on qqbar_w.cc in Chroma
 *  Color is nested inside of spin.
 *  Authors:
 *  Haobo Yan
 *  Andre Walker-Loud
 *  Ben Horz
 */
 // Mathematically, the correlation function is
 //                  _____
 //                  \      -i(px+p'x')                   +     +
 // Cππ(p,p',t)   = >    e           < 0 | π(x)π(x')π(y')π(y) | 0 >
 //                  /
 //                  -----
 //                   x,x'

 // In this program, we implement the calculation with the following formula
 // C(p,p',t) = 
 //    tr(Q1(p,y,y) * g5)*tr(P1(p',y',y') * g5)
 //  - tr(Q2(p,y,y') * g5 * P2(p',y',y) * g5)
 //  - tr(Q3(p,y',y) * g5 * P3(p',y,y') * g5)
 //  + tr(Q4(p,y',y') * g5)*tr(P4(p',y,y) * g5)

 //S_u(x, y)  = quark_prop_1
 //S_d(x, y)  = quark_prop_2
 //S_u(x, y') = quark_prop_3
 //S_d(x, y') = quark_prop_4

#include "chromabase.h"
#include "pipi_scattering_func_w.h"
#include "../../../chroma/lib/meas/hadron/qqbar_w.cc"

namespace Chroma
{

	void pipi_correlator(CorrelatorType::Correlator& correlator_out, const LatticePropagator& quark_prop_1, const LatticePropagator& quark_prop_2, const LatticePropagator& quark_prop_3, const LatticePropagator& quark_prop_4, const multi2d<int>& origin_list, const int p2max, const int t0, const int j_decay)
	{
		int G5 = Ns * Ns - 1;

		// Construct SftMoms
		SftMom phases(p2max, false, j_decay);
		int Nt = phases.numSubsets();

		// Judge if the origin is the same
		bool same_origin = true;
		for (int o_comp = 0; o_comp < origin_list.size2(); ++o_comp)
			if (origin_list[0][o_comp] != origin_list[1][o_comp])
			{
				same_origin = false;
				break;
			}

		if (same_origin) // y=y'
		{
			multi2d<DPropagator> Q;

			compute_qqbar(Q, quark_prop_2, quark_prop_1, phases, t0);

			DComplex origin_fix;
			Double orgin_phases;
			multi1d<int> mom_comp1, mom_comp2;

			// Loop over all ps, p's and time
			for (int mom_num1 = 0; mom_num1 < phases.numMom(); ++mom_num1)
				for (int mom_num2 = 0; mom_num2 < phases.numMom(); ++mom_num2)
					for (int t = 0; t < Nt; ++t)
					{
						mom_comp1 = phases.numToMom(mom_num1);
						mom_comp2 = phases.numToMom(mom_num2);

						// Fix the origin
						orgin_phases = 0;

						for (int p_comp = 0; p_comp < 3; ++p_comp)
							orgin_phases += mom_comp1[p_comp] * origin_list[0][p_comp] + mom_comp2[p_comp] * origin_list[1][p_comp];
						origin_fix = cmplx(cos(orgin_phases), sin(orgin_phases));

						correlator_out[std::make_pair(std::make_tuple(mom_comp1[0], mom_comp1[1], mom_comp1[2]), std::make_tuple(mom_comp2[0], mom_comp2[1], mom_comp2[2]))][t] = 2 * (trace(Q[mom_num1][t] * Gamma(G5)) * trace(Q[mom_num2][t] * Gamma(G5)) - trace(Q[mom_num1][t] * Gamma(G5) * Q[mom_num2][t] * Gamma(G5))) * origin_fix;

					}
		}
		else // The general case
		{
			multi2d<DPropagator> Q1, P1, Q2, P2;

			compute_qqbar(Q1, quark_prop_2, quark_prop_1, phases, t0);
			compute_qqbar(P1, quark_prop_4, quark_prop_3, phases, t0);
			compute_qqbar(Q2, quark_prop_2, quark_prop_3, phases, t0);
			compute_qqbar(P2, quark_prop_4, quark_prop_1, phases, t0);

			//compute_qqbar(Q3, quark_prop_4, quark_prop_1, phases1, t0);
			//compute_qqbar(P3, quark_prop_2, quark_prop_3, phases2, t0);
			//compute_qqbar(Q4, quark_prop_4, quark_prop_3, phases1, t0);
			//compute_qqbar(P4, quark_prop_2, quark_prop_1, phases2, t0);

			// The above Q3 to P4 can be deduced from the first four. So really we only need to compute four of the eight qqbar blocks
			//Q3 = P2;
			//P3 = Q2;
			//Q4 = P1;
			//P4 = Q1;

			DComplex origin_fix;
			Double orgin_phases;
			multi1d<int> mom_comp1, mom_comp2;

			for (int mom_num1 = 0; mom_num1 < phases.numMom(); ++mom_num1)
				for (int mom_num2 = 0; mom_num2 < phases.numMom(); ++mom_num2)
					for (int t = 0; t < Nt; ++t)
					{
						mom_comp1 = phases.numToMom(mom_num1);
						mom_comp2 = phases.numToMom(mom_num2);

						// Fix the origin
						orgin_phases = 0;

						for (int p_comp = 0; p_comp < 3; ++p_comp)
							orgin_phases += mom_comp1[p_comp] * origin_list[0][p_comp] + mom_comp2[p_comp] * origin_list[1][p_comp];
						origin_fix = cmplx(cos(orgin_phases), sin(orgin_phases));

						correlator_out[std::make_pair(std::make_tuple(mom_comp1[0], mom_comp1[1], mom_comp1[2]), std::make_tuple(mom_comp2[0], mom_comp2[1], mom_comp2[2]))][t] = trace(Q1[mom_num1][t] * Gamma(G5)) * trace(P1[mom_num2][t] * Gamma(G5)) - trace(Q2[mom_num1][t] * Gamma(G5) * P2[mom_num2][t] * Gamma(G5)) - trace(P2[mom_num1][t] * Gamma(G5) * Q2[mom_num2][t] * Gamma(G5)) + trace(P1[mom_num1][t] * Gamma(G5)) * trace(Q1[mom_num2][t] * Gamma(G5)) * origin_fix;

					}

		}

	}


} // End namespace Chroma

