/*!
 *  Functions that perform the contraction of quarks into the two-pion system.
 *  Based on qqbar_w.cc in Chroma
 *  Color is nested inside of spin.
 *  Authors:
 *  Haobo Yan
 *  Ben Horz
 *  Andre Walker-Loud
 */

// Mathematically, the correlation function is
//                  _____
//                  \      -i(px+p'x')                   +     +
// Cππ(p,p',t)   = >    e           < 0 | π(x)π(x')π(y')π(y) | 0 >
//                  /
//                  -----
//                   x,x'

// We will always use y' = y
// In this program, we implement the calculation with the following formula
// pipi and KK
// C(p,p',t) =
//    tr(Q1(p,y,y) * g5)*tr(P1(p',y',y') * g5)
//  - tr(Q2(p,y,y') * g5 * P2(p',y',y) * g5)
//  - tr(Q3(p,y',y) * g5 * P3(p',y,y') * g5)
//  + tr(Q4(p,y',y') * g5)*tr(P4(p',y,y) * g5)

// piK
// C(p,p',t) =
//    tr(Q1(p,y,y) * g5)*tr(P1(p',y',y') * g5)
//  - tr(Q2(p,y,y') * g5 * P2(p',y',y) * g5)

// Propagator notations
// quark_prop_1 = this will be the propagator
// quark_prop_2 = this will be the anti-quark
// For pi-K, we will specify all 4 quarks
// quark_prop_3 = 2nd quark
// quark_prop_4 = 2nd anti-quark

// Chroma includes
#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/qqbar_w.h"
// Lalibe includes
#include "pipi_scattering_func_w.h"

namespace Chroma
{
    void pipi_correlator(CorrelatorType::Correlator& correlator_out,
                         const LatticePropagator& quark_prop_1, const LatticePropagator& a_quark_prop_2,
                         const multi1d<int>& origin, const int p2max, const int ptot2max,
                         const int t0, const int j_decay, const int diagram)
    {
        int G5 = Ns * Ns - 1;
        // Construct SftMoms
        SftMom phases(p2max, false, j_decay);
        int Nt = phases.numSubsets();

        multi1d<DComplex> tmp_multi1d(Nt); // Use this to give correlator a Initialization
        tmp_multi1d=zero;
        multi2d<DPropagator> Q(phases.numMom(), Nt);
        Q = zero;

        QDPIO::cout << "DEBUG: compute_qqbar(Q, quark_prop_1, a_quark_prop_2, phases, t0)" << std::endl;
        compute_qqbar(Q, quark_prop_1, a_quark_prop_2, phases, t0);
        QDPIO::cout << "DEBUG: success qqbar" << std::endl;

        DComplex origin_fix;
        Double origin_phases;
        multi1d<int> mom_comp1, mom_comp2;

        // Loop over all ps, p's and time
        QDPIO::cout << "DEBUG: loop over mom_num1" << std::endl;
        for (int mom_num1 = 0; mom_num1 < phases.numMom(); ++mom_num1){
            QDPIO::cout << "DEBUG:     loop over mom_num2" << std::endl;
            for (int mom_num2 = 0; mom_num2 < phases.numMom(); ++mom_num2){
                mom_comp1 = phases.numToMom(mom_num1);
                mom_comp2 = phases.numToMom(mom_num2);

                // Select momenta pairs that are less than equal to p_tot^2_max
                int ptot2 = 0;
                for (int p_comp = 0; p_comp < 3; ++p_comp){
                    ptot2 += (mom_comp1[p_comp] + mom_comp2[p_comp]) * (mom_comp1[p_comp] + mom_comp2[p_comp]);
                }
                if (ptot2 <= ptot2max){
                    // Fix the origin
                    origin_phases = 0;

                    for (int p_comp = 0; p_comp < 3; ++p_comp){
                        origin_phases -= (mom_comp1[p_comp] + mom_comp2[p_comp]) * (origin[p_comp]) * 2 * M_PI / Layout::lattSize()[p_comp];
                    }

                    origin_fix = cmplx(cos(origin_phases), sin(origin_phases));

                    std::pair< std::tuple<int, int, int>, std::tuple<int, int, int> > pq_pair;
                    pq_pair = std::make_pair(std::make_tuple(mom_comp1[0], mom_comp1[1], mom_comp1[2]), std::make_tuple(mom_comp2[0], mom_comp2[1], mom_comp2[2]));

                    correlator_out[pq_pair] = tmp_multi1d;

                    if (diagram == 1 || diagram == 4){
                        for (int t = 0; t < Nt; ++t){
                            correlator_out[pq_pair][t] = (
                                                          trace(Q[mom_num1][t] * Gamma(G5))
                                                          * trace(Q[mom_num2][t] * Gamma(G5))
                                                          ) * origin_fix;
                        }
                    } else if (diagram == 2 || diagram == 3){
                        for (int t = 0; t < Nt; ++t){
                            correlator_out[pq_pair][t] = -(
                                                           trace(Q[mom_num1][t] * Gamma(G5) * Q[mom_num2][t] * Gamma(G5))
                                                           ) * origin_fix;
                        }
                    } else {
                        for (int t = 0; t < Nt; ++t){
                            correlator_out[pq_pair][t] = 2 * (
                                                              trace(Q[mom_num1][t] * Gamma(G5)) * trace(Q[mom_num2][t] * Gamma(G5))
                                                              - trace(Q[mom_num1][t] * Gamma(G5) * Q[mom_num2][t] * Gamma(G5))
                                                              ) * origin_fix;
                        }
                    }
                }
            }
            QDPIO::cout << "DEBUG:     mom_num2 finish" << std::endl;
        }
        QDPIO::cout << "DEBUG: mom_num1 finish" << std::endl;
    }

    void pik_correlator(CorrelatorType::Correlator& correlator_out,
                        const LatticePropagator& quark_prop_1, const LatticePropagator& a_quark_prop_2,
                        const LatticePropagator& quark_prop_3, const LatticePropagator& a_quark_prop_4,
                        const multi1d<int>& origin, const int p2max, const int ptot2max,
                        const int t0, const int j_decay, const int diagram)
    {
        int G5 = Ns * Ns - 1;

        // Construct SftMoms
        SftMom phases(p2max, false, j_decay);
        int Nt = phases.numSubsets();

        multi1d<DComplex> tmp_multi1d(Nt);
        tmp_multi1d=zero;

        multi2d<DPropagator> Q1(phases.numMom(), Nt), P1(phases.numMom(), Nt);
        Q1 = zero;
        P1 = zero;

        compute_qqbar(Q1, quark_prop_1, a_quark_prop_2, phases, t0);
        compute_qqbar(P1, quark_prop_3, a_quark_prop_4, phases, t0);

        DComplex origin_fix;
        Double origin_phases;
        multi1d<int> mom_comp1, mom_comp2;

        for (int mom_num1 = 0; mom_num1 < phases.numMom(); ++mom_num1)
            for (int mom_num2 = 0; mom_num2 < phases.numMom(); ++mom_num2)
                {
                    mom_comp1 = phases.numToMom(mom_num1);
                    mom_comp2 = phases.numToMom(mom_num2);

                    int ptot2 = 0;
                    for (int p_comp = 0; p_comp < 3; ++p_comp)
                        ptot2 += (mom_comp1[p_comp] + mom_comp2[p_comp]) * (mom_comp1[p_comp] + mom_comp2[p_comp]);
                    if (ptot2 <= ptot2max)
                        {
                            origin_phases = 0;
                            for (int p_comp = 0; p_comp < 3; ++p_comp)
                                {
                                    origin_phases -= (mom_comp1[p_comp] + mom_comp2[p_comp]) * (origin[p_comp]) * 2 * M_PI / Layout::lattSize()[p_comp];
                                }
                            origin_fix = cmplx(cos(origin_phases), sin(origin_phases));

                            std::pair< std::tuple<int, int, int>, std::tuple<int, int, int> > pq_pair;
                            pq_pair = std::make_pair(std::make_tuple(mom_comp1[0], mom_comp1[1], mom_comp1[2]), std::make_tuple(mom_comp2[0], mom_comp2[1], mom_comp2[2]));

                            correlator_out[pq_pair] = tmp_multi1d;

                            if (diagram == 1)
                                {
                                    for (int t = 0; t < Nt; ++t)
                                        correlator_out[pq_pair][t] =
                                            trace(Q1[mom_num1][t] * Gamma(G5))
                                            * trace(P1[mom_num2][t] * Gamma(G5))
                                            * origin_fix;
                                }
                            else if (diagram == 2)
                                {
                                    for (int t = 0; t < Nt; ++t)
                                        correlator_out[pq_pair][t] =
                                            - trace(Q1[mom_num1][t] * Gamma(G5) * P1[mom_num2][t] * Gamma(G5))
                                            * origin_fix;
                                }
                            else
                                {
                                    for (int t = 0; t < Nt; ++t)
                                        correlator_out[pq_pair][t] = (
                                                                      trace(Q1[mom_num1][t] * Gamma(G5)) * trace(P1[mom_num2][t] * Gamma(G5))
                                                                      - trace(Q1[mom_num1][t] * Gamma(G5) * P1[mom_num2][t] * Gamma(G5))
                                                                      ) * origin_fix;
                                }
                        }
                }
    }

} // End namespace Chroma
