/*
 * Form factor code from chroma, modified to use lalibe's fourier transform routine.
 * Arjun Gambhir
 */

#ifndef __lalibe_baryon_seqsource_h__
#define __lalibe_baryon_seqsource_h__



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
                            );


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
                            );

}  // end namespace Chroma

#endif
