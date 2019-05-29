/*
 * Routines to construct baryon sequential sources.
 * David Brantley
   11-Feb-2019 AWL adding support for all time-slice seqsrc
 */

//CHROMA STUFF
#include "chromabase.h"

//LALIBE STUFF
#include "baryon_seqsource_w.h"
#include "util/ferm/diractodr.h"
#include "seqsource_contractions_func_w.h"


/*
Construct the baryon sequential sources.
Takes in:

    LatticePropagator quark flavors
    std::string       source_spin
    std::string       sink_spin
    multi1d<int>      sink_mom    - 3-momentum of sink.
    multi1d<int>      origin_off  - Source point of the propagators. This must be the same for all quarks.
    multi1d<int>      BC          - Boundary conditions for all propagators, again must be the same.
    int               t_sink      - Time to fix the sink.
    int               j_decay     - Decay direction.
    int               parity      - positive [0] or negative [1]
    bool              t_all       - do all time slices or slice t_sink

Returns a sequential propagator in DR basis suitable for inversion.
*/


namespace Chroma
{
    LatticePropagator ProtDtoD(
        LatticePropagator& up_quark,
        std::string & source_spin,
        std::string & sink_spin,
        multi1d<int>& sink_mom,
        multi1d<int>& origin_off,
        multi1d<int>& BC,
        int & t_sink,
        int& j_decay,
        int parity,
        bool t_all
    )
    {
        // quark_1 is the down, while quark_2 is the up quark.
        rotate_to_Dirac_Basis(up_quark);

        LatticePropagator tmp1 = up_quark;
        LatticePropagator tmp2 = up_quark;
        LatticePropagator tmp_seq_source;

        SpinMatrix sink_proj;
        SpinMatrix source_proj;
        //Parity is now passed into protonDiQuark and protonSpin functions.
        SpinMatrix diquark_proj = protonDiquarkSpin(parity);

        if(source_spin == "up")
        {
            if(sink_spin == "up")
            {
                source_proj = protonSpinUp(parity);
                sink_proj   = protonSpinUp(parity);
            }
            else if(sink_spin == "dn")
            {
                source_proj = protonSpinUp(parity);
                sink_proj   = protonSpinDn(parity);
            }
            else
            {
                QDPIO::cerr << "Sink-Source spin combination "
                            <<sink_spin<<"-"<<source_spin<<" unknown." <<std::endl;
                QDP_abort(1);
            }
        }
        else if(source_spin == "dn")
        {
            if(sink_spin == "up")
            {
                source_proj = protonSpinDn(parity);
                sink_proj   = protonSpinUp(parity);
            }
            else if(sink_spin == "dn")
            {
                source_proj = protonSpinDn(parity);
                sink_proj   = protonSpinDn(parity);
            }
            else
            {
                QDPIO::cerr << "Sink-Source spin combination "
                            <<sink_spin<<"-"<<source_spin<<" unknown." <<std::endl;
                QDP_abort(1);
            }
        }
        else
        {
            QDPIO::cerr << "Sink-Source spin combination "
                        <<sink_spin<<"-"<<source_spin<<" unknown." <<std::endl;
            QDP_abort(1);
        }

        tmp1 = tmp1*transpose(source_proj);
        tmp1 = -adj(diquark_proj)*tmp1;

        tmp2 = tmp2*adj(diquark_proj);
        tmp2 = sink_proj*tmp2;

        tmp_seq_source = contract1(tmp1,tmp2);

        tmp1 = up_quark;
        tmp2 = up_quark;

        tmp1 = tmp1*adj(diquark_proj);
        tmp1 = -adj(diquark_proj)*tmp1;

        tmp2 = sink_proj*tmp2;
        tmp2 = tmp2*transpose(source_proj);

        tmp_seq_source -= contract2(tmp1,tmp2);

        // Apply the g5 hermitian transformation before inversion.
        tmp_seq_source = transpose(tmp_seq_source);

        rotate_from_Dirac_Basis(tmp_seq_source);
        //    tmp_seq_source = Gamma(Nd*Nd-1)*conj(tmp_seq_source);
        tmp_seq_source = adj(Gamma(Nd*Nd-1)*tmp_seq_source*Gamma(Nd*Nd-1));

        tmp_seq_source = projectBaryonSeqSource(tmp_seq_source, sink_mom, origin_off, t_sink, j_decay, BC, parity, t_all);

        return tmp_seq_source;
    } // end ProtDtoD

    LatticePropagator ProtUtoU(
        LatticePropagator& up_quark,
        LatticePropagator& dn_quark,
        std::string & source_spin,
        std::string & sink_spin,
        multi1d<int>& sink_mom,
        multi1d<int>& origin_off,
        multi1d<int>& BC,
        int & t_sink,
        int& j_decay,
        int parity,
        bool t_all
    )
    {
        // quark_1 is the down, while quark_2 is the up quark.
        rotate_to_Dirac_Basis(up_quark);
        rotate_to_Dirac_Basis(dn_quark);

        LatticePropagator tmp1 = up_quark;
        LatticePropagator tmp2 = up_quark;
        LatticePropagator tmp_seq_source;

        SpinMatrix sink_proj;
        SpinMatrix source_proj;
        SpinMatrix diquark_proj = protonDiquarkSpin(parity);

        if(source_spin == "up")
        {
            if(sink_spin == "up")
            {
                source_proj = protonSpinUp(parity);
                sink_proj = protonSpinUp(parity);
            }
            else if(sink_spin == "dn")
            {
                source_proj = protonSpinUp(parity);
                sink_proj = protonSpinDn(parity);
            }
            else
            {
                QDPIO::cerr << "Sink-Source spin combination " <<sink_spin<<"-"<<source_spin<<" unknown." <<std::endl;
                QDP_abort(1);
            }
        }
        else if(source_spin == "dn")
        {
            if(sink_spin == "up")
            {
                source_proj = protonSpinDn(parity);
                sink_proj = protonSpinUp(parity);
            }
            else if(sink_spin == "dn")
            {
                source_proj = protonSpinDn(parity);
                sink_proj = protonSpinDn(parity);
            }
            else
            {
                QDPIO::cerr << "Sink-Source spin combination " <<sink_spin<<"-"<<source_spin<<" unknown." <<std::endl;
                QDP_abort(1);
            }
        }
        else
        {
            QDPIO::cerr << "Sink-Source spin combination " <<sink_spin<<"-"<<source_spin<<" unknown." <<std::endl;
            QDP_abort(1);
        }

        tmp1 = diquark_proj*dn_quark;
        tmp1 = -tmp1*diquark_proj; //Diquark. This will be in every contraction.

        tmp2 = sink_proj*up_quark;
        tmp2 = tmp2*transpose(source_proj);

        tmp_seq_source = -contract2(tmp1,tmp2); //First term.

        tmp2 = up_quark*transpose(source_proj);
        tmp2 = tmp2*source_proj;

        tmp2 = tmp2*transpose(sink_proj);
        tmp2 = transposeSpin(tmp2);

        tmp_seq_source += contract4(tmp2,tmp1); //Second term.

        tmp2 = transposeSpin(tmp1);
        tmp2 = contract3(up_quark,tmp2);
        tmp2 = source_proj*traceSpin(tmp2);
        tmp2 = transpose(sink_proj)*tmp2;

        tmp_seq_source -= tmp2; // Third term.

        tmp2 = sink_proj*up_quark;
        tmp2 = transpose(source_proj)*tmp2;
        tmp1 = transposeSpin(tmp1); // Final term so why not.

        tmp2 = contract4(tmp2,tmp1);
        tmp_seq_source += transposeSpin(tmp2); // Final term.

        // Apply the g5 hermitian transformation before inversion.
        tmp_seq_source = transpose(tmp_seq_source);

        rotate_from_Dirac_Basis(tmp_seq_source);
        //    tmp_seq_source = Gamma(Nd*Nd-1)*conj(tmp_seq_source);
        tmp_seq_source = adj(Gamma(Nd*Nd-1)*tmp_seq_source*Gamma(Nd*Nd-1));

        tmp_seq_source = projectBaryonSeqSource(tmp_seq_source, sink_mom, origin_off, t_sink, j_decay, BC, parity, t_all);

        return tmp_seq_source;
    }
}  // end namespace Chroma
