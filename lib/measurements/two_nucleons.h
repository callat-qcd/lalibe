/*
  Authors
  Andre Walker-Loud

  re-factoring of Arjun's port of latscat to lalibe
  We will have two routines now instead of 1
  - make blocks
  - contract linear combinations of blocks
*/

#ifndef __lalibe_Two_Nucleons_h__
#define __lalibe_Two_Nucleons_h__

//Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma
{
    namespace LalibeTwoNucleonsEnv
    {
        extern const std::string name;
        bool registerAll();

        struct TwoNucleonsParams
        {
            TwoNucleonsParams();
            TwoNucleonsParams(XMLReader& xml_in, const std::string& path);
            void writeXML(XMLWriter& xml_out, const std::string& path);

            unsigned long      frequency;

            // only same flavor bilinear operators are used
            struct TwoNucleons_t
            {
                std::string             contractions_filename;  // filename of hdf5 file that contains contraction terms
                std::string             output_filename;        // output file
                int                     output_stripesize;      // output stripesize; default recommended
                multi1d<std::string>    parities;               // specifies list of parity strings
                multi1d<int>            origin;                 // X0 Y0 Z0 T0 of origin
                bool                    compute_locals;         // compute local sinks (still FFT separately)?
                bool                    compute_loc_o;          // compute local origin sinks?
                bool                    compute_proton;         // compute local sinks (still FFT separately)?
                multi1d< multi1d<int> > boosts;                 // boosts
                multi1d<std::string>    prop_list_delete;       // list of prop ids to detele associated blocks
            } twonucleonsparam ;

            struct NamedObject_t
            {
                std::string             gauge_id;
                multi1d<std::string>    prop0_list;
                multi1d<std::string>    prop1_list;
                multi1d<GroupXML_t>     nucleon_blocks;         // list of blocks and their weights
            } named_obj;
        };
        // This is used for boosts
        const std::string dirlist[4]={"x","y","z","t"};

        const std::vector<std::string> nonlocalterms = {"PP_SING0","PP_TRIPP","PP_TRIP0","PP_TRIPM", "PN_SING0","PN_TRIPP","PN_TRIP0","PN_TRIPM"};
        const std::vector<std::string> localterms0   = {"PP_SING0_loc0", "PN_TRIPP_loc0","PN_TRIP0_loc0","PN_TRIPM_loc0"};
        const std::vector<std::string> localterms1   = {"PP_SING0_loc1", "PN_TRIPP_loc1","PN_TRIP0_loc1","PN_TRIPM_loc1"};

        /*
        nonlocalterms = {"PP_SING0","PP_TRIPP","PP_TRIP0","PP_TRIPM", "PN_SING0","PN_TRIPP","PN_TRIP0","PN_TRIPM"};
        localterms0   = {"PP_SING0_loc0", "PN_TRIPP_loc0","PN_TRIP0_loc0","PN_TRIPM_loc0"};
        localterms1   = {"PP_SING0_loc1", "PN_TRIPP_loc1","PN_TRIP0_loc1","PN_TRIPM_loc1"};
        const std::string nonlocalterms[8]={
            "PP_SING0","PP_TRIPP","PP_TRIP0","PP_TRIPM",
            "PN_SING0","PN_TRIPP","PN_TRIP0","PN_TRIPM"
        };
        const std::string localterms0[4] = {"PP_SING0_loc0", "PN_TRIPP_loc0","PN_TRIP0_loc0","PN_TRIPM_loc0"};
        const std::string localterms1[4] = {"PP_SING0_loc1", "PN_TRIPP_loc1","PN_TRIP0_loc1","PN_TRIPM_loc1"};
        */

        class InlineMeas : public AbsInlineMeasurement
        {
            public:
                ~InlineMeas() {}
                InlineMeas(const TwoNucleonsParams& p) : params(p) {}
                InlineMeas(const InlineMeas& p)        : params(p.params) {}

                unsigned long getFrequency(void) const {return params.frequency;}

                //! Do the measurement
                void operator()(const unsigned long update_no, XMLWriter& xml_out);

            private:
                TwoNucleonsParams params;
        };

    }

}

#endif
