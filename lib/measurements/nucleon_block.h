/*
  Authors
  Andre Walker-Loud

  re-factoring of Arjun's port of latscat to lalibe
  We will have two routines now instead of 1
  - make blocks
  - contract linear combinations of blocks
*/

#ifndef __lalibe_Nucleon_Block_h__
#define __lalibe_Nucleon_Block_h__

//Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma
{
    namespace LalibeNucleonBlockEnv
    {
        extern const std::string name;
        bool registerAll();

        struct NucleonBlockParams
        {
            NucleonBlockParams();
            NucleonBlockParams(XMLReader& xml_in, const std::string& path);
            void writeXML(XMLWriter& xml_out, const std::string& path);

            unsigned long      frequency;

            // only same flavor bilinear operators are used
            struct NucleonBlock_t{
                unsigned int            fft_chunksize;         //originally the only parameter in FFTPar struct
                bool                    fft_tune;              //tune the fft?
                bool                    compute_locals;        // make blocks for local NN computation?
                multi1d< multi1d<int> > boosts;                //boosts
                bool                    in_dirac_basis;        //specifies if props are in dirac basis, this is false by default
                bool                    negative_parity;       //specifies if we will run contractions on the negative parity nucleons
                //                                               The momentumspace trucation is not used
                //int                     contractions_n_sq;     // truncate rel-momentum space to this value
            } nblockparam ;

            struct NamedObject_t
            {
                //Gauge field likely needed for latscat from this section.
                std::string    gauge_id;
                std::string    prop0_id;
                std::string    prop1_id;
                // nucleon block map
                std::string    block_map;

            } named_obj;
        };
        // This is used for boosts
        const std::string dirlist[4]={"x","y","z","t"};

        class InlineMeas : public AbsInlineMeasurement
            {
            public:
                ~InlineMeas() {}
                InlineMeas(const NucleonBlockParams& p) : params(p) {}
                InlineMeas(const InlineMeas& p)         : params(p.params) {}

                unsigned long getFrequency(void) const {return params.frequency;}

                //! Do the measurement
                void operator()(const unsigned long update_no, XMLWriter& xml_out);

            private:
                NucleonBlockParams params;
            };

    }

}

#endif
