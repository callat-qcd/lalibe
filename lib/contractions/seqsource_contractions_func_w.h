/*
 * Helper routines for constructing sequential sources.
 * David Brantley
 * IMPORTANT NOTE: The quarks must be in the DIRAC BASIS for these routines to work.
 */


#ifndef __seqsource_contractions_func_w_h__
#define __seqsource_contractions_func_w_h__

//CHROMA STUFF
#include "chromabase.h"

//LALIBE STUFF
#include "baryon_seqsource_w.h"
#include "util/ferm/diractodr.h"



namespace Chroma
{
    void rotate_to_Dirac_Basis(LatticePropagator & quark_to_be_rotated);

    void rotate_from_Dirac_Basis(LatticePropagator & quark_to_be_rotated);

    SpinMatrix protonDiquarkSpin(int parity);

    SpinMatrix protonSpinUp(int parity);

    SpinMatrix protonSpinDn(int parity);

    LatticeColorMatrix dblEpsContract1( LatticeColorMatrix& cmat1,
                                        LatticeColorMatrix& cmat2
                                        );

    LatticeColorMatrix dblEpsContract2( LatticeColorMatrix& cmat1,
                                        LatticeColorMatrix& cmat2
                                        );

    LatticeColorMatrix dblEpsContract3( LatticeColorMatrix& cmat1,
                                        LatticeColorMatrix& cmat2
                                        );

    LatticePropagator contract1(LatticePropagator& quark_1,
                                LatticePropagator& quark_2
                                );

    LatticePropagator contract2(LatticePropagator& quark_1,
                                LatticePropagator& quark_2
                                );

    LatticePropagator contract3(LatticePropagator& quark_1,
                                LatticePropagator& quark_2
                                );

    LatticePropagator contract4(LatticePropagator& quark_1,
                                LatticePropagator& quark_2
                                );

    LatticeComplex singlePhase( multi1d<int>& mom,
                                multi1d<int>& origin_off,
                                int j_decay
                                );

    void baryonTimeOrder(   LatticeComplex& timeorder,
                            int t_source,
                            multi1d<int>& BC,
                            int j_decay,
                            int parity=0
                            );

    void projectTimeSlice( LatticePropagator& source_prop,
                                        int t_sink,
                                        int j_decay
                                        );

    LatticePropagator projectBaryonSeqSource(
                                            LatticePropagator& seq_source,
                                            multi1d<int>& mom,
                                            multi1d<int>& origin_off,
                                            int t_sink,
                                            int j_decay,
                                            multi1d<int>& BC,
                                            int parity=0,
                                            bool t_all=false
                                            );

}  // end namespace Chroma

#endif
