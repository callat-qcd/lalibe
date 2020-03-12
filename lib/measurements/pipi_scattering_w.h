/*
 * pipi scattering
 * Authors:
 * Arjun Gambhir
 * Andre Walker-Loud
 * Jason Chang
 * David Brantley
 * Ben Horz
 * Haobo Yan
 * Do pipi scattering and write out the two-point correlator in hdf5
 */


#ifndef __lalibe_pipi_scattering_h__
#define __lalibe_pipi_scattering_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma
{
    namespace LalibePipiScatteringEnv
    {
        extern const std::string name;
        bool registerAll();


        struct PipiParams
        {
            PipiParams();
            PipiParams(XMLReader& xml_in, const std::string& path);
            void writeXML(XMLWriter& xml_out, const std::string& path);

            unsigned long frequency;

            struct Param_t
            {

                bool output_full_correlator;          //If no momentum is specified, we output the full correlator.
                int p2_max;                           //max of momentum transfer squared, optional
                multi1d<std::string> particle_list;   //list of actual particles we gunna make from the contractions yo
#ifdef BUILD_HDF5
                std::string file_name;
                std::string obj_path;
#endif
            } param;

            struct NamedObject_t
            {
                //std::string  gauge_id;     //Grab the usual gauge field.
                std::string  quark_prop_1;
                bool is_prop_1;
                std::string  quark_prop_2;
                bool is_prop_2;
                std::string  quark_prop_3;
                bool is_prop_3;
                std::string  quark_prop_4;
                bool is_prop_4;

                //Above are various strings and bools that will identify the quarks and which quarks are present in the input.
                //This could just be done by checking if these strings are null (and ensuring they are with no input), but I like bools, so whatever.

            } named_obj;

        };


        class InlineMeas : public AbsInlineMeasurement
        {
        public:
            ~InlineMeas() {}
            InlineMeas(const PipiParams& p) : params(p) {}

            unsigned long getFrequency(void) const { return params.frequency; }

            void operator()(const unsigned long update_no,
                XMLWriter& xml_out);

        private:
            PipiParams params;
        };

    } // namespace LalibeMesonContractionsEnv

} // namespace Chroma

#endif
