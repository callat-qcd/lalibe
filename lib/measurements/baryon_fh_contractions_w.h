/*
 * Single baryon FH contractions
 * 
 * Andre Walker-Loud with Ben Horz's guidance - 19 Aug 2022
 * Do baryon FH contractions and write out the two-point correlator in hdf5
 */
#ifndef __lalibe_baryon_fh_contractions_h__
#define __lalibe_baryon_fh_contractions_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

#include <set>
#include <map>
#include <string>

namespace Chroma 
{ 
    namespace LalibeBaryonFHContractionsEnv
    {
        extern const std::string name;
        bool registerAll();

        struct FHBaryonParams 
        {
            FHBaryonParams();
            FHBaryonParams(XMLReader& xml_in, const std::string& path);
            void writeXML(XMLWriter& xml_out, const std::string& path);

            unsigned long frequency;

            struct Param_t
            {
                bool ng_parity;                       //Also write out negative parity partner (this can be used for doubling statistics)
                bool rotate_to_Dirac;                 //If the correlator is in DeGrand-Rossi basis, set this to true to rotate 
                                                      //to Dirac basis.
                bool is_antiperiodic;                 //Tracking anti-periodicity, if not specified, it's assumed to be true.
                bool output_full_correlator;          //If no momentum is specified, we output the full correlator.
                bool is_mom_max;                      //keeps track of which momentum mode we are using
                int p2_max;                           //max of momentum transfer squared, optional
                multi1d<multi1d<int>> mom_list;       //list of momenta insertions, optional
                multi2d<int> p_list;                  //momentum list the slow fourier transform needs
                std::set<std::pair<std::string, std::string>> particle_list;   //list of actual particles we gunna make from the contractions
                multi1d<std::string> currents;        //list of currents - should correspond to NamedObject FH props
                multi1d<std::string> flavors;         //list of flavors  - should correspond to NamedObject FH props
#ifdef BUILD_HDF5
                std::string file_name;
                std::string obj_path;
#endif
            } param;

            struct NamedObject_t
            {
                //std::string  gauge_id;     //Grab the usual gauge field.
                std::map<char, std::string>  quark_map;
                multi1d<std::string>         fh_quarks;

            } named_obj;

        };


        class InlineMeas : public AbsInlineMeasurement 
        {
        public:
            ~InlineMeas() {}
        InlineMeas(const FHBaryonParams& p) : params(p) {}

            unsigned long getFrequency(void) const {return params.frequency;}

            void operator()(const unsigned long update_no,
                            XMLWriter& xml_out); 

        private:
            FHBaryonParams params;
        };


        namespace {
            std::map<std::string, std::vector<std::string>> aliasMap = {
                { "octet",    {"proton", "lambda_z", "sigma_p", "xi_z"} },
                { "decuplet", {"delta_pp", "sigma_star_p", "xi_star_z", "omega_m"} },
                { "octet_iso",    {"proton", "neutron", "lambda_z", "lambda_to_sigma", "sigma_to_lambda", "sigma_p", "sigma_z", "sigma_m", "xi_z", "xi_m"} },
                { "decuplet_iso", {"delta_pp", "delta_p", "delta_z", "delta_m", "sigma_star_p", "sigma_star_z", "sigma_star_m", "xi_star_z", "xi_star_m", "omega_m"} },
            };
        }

    } // namespace LalibeBaryonContractionsEnv

} // namespace Chroma

#endif
