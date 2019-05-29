/*
 * Baryon contractions that intercept normal quark lines with an FH propagator.
 * We call these Feynman Hellman c2pt baryon contractions.
 * These assume initial and final interpolating states are the same (flavor conserving).
 * Authors:
 * Arjun Gambhir
 * Do FH baryon contractions and write out the two-point correlator in hdf5 or asci
 */
#ifndef __lalibe_flavor_conserving_fh_baryon_contractions_h__
#define __lalibe_flavor_conserving_fh_baryon_contractions_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma
{
  namespace LalibeFlavorConservingFHBaryonContractionsEnv
  {
    extern const std::string name;
    bool registerAll();


    struct FlavorConservingFHBaryonParams
    {
      FlavorConservingFHBaryonParams();
      FlavorConservingFHBaryonParams(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long frequency;

      struct Param_t
      {
	bool ng_parity;                       //Also write out negative parity partner (this can be used for doubling statistics)
	bool rotate_to_Dirac;                 //If the correlator is in DeGrand-Rossi basis, set this to true to rotate
					      //to Dirac basis.
	bool is_antiperiodic;                //Tracking anti-periodicity, if not specified, it's assumed to be true.
	bool output_full_correlator;         //If no momentum is specified, we output the full correlator.
	bool is_mom_max;                     //keeps track of which momentum mode we are using
	int p2_max;                          //max of momentum transfer squared, optional
	multi1d<multi1d<int>> mom_list;      //list of momenta insertions, optional
	multi2d<int> p_list;                 //momentum list the slow fourier transform needs
	multi1d<std::string> particle_list;  //list of actual particles we gunna make from the contractions yo
    std::string flavor;                  //flavor of the fh_prop - needed for determing wick contractions
    std::string current;                 //current used in FH prop
#ifdef BUILD_HDF5
	std::string file_name;
	std::string obj_path;
#endif
      } param;

      struct NamedObject_t
      {
	std::string  up_quark;
	bool is_up;
	std::string  down_quark;
	bool is_down;
	std::string  strange_quark;
	bool is_strange;
	std::string  charm_quark;
	bool is_charm;
        //Above are various strings and bools that will identify the quarks and which quarks are present in the input.
	//This could just be done by checking if these strings are null (and ensuring they are with no input), but I like bools, so whatever.
	std::string fh_quark;	//This is the fh_prop, it MUST be present for this code to run.

      } named_obj;

    };


    class InlineMeas : public AbsInlineMeasurement
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const FlavorConservingFHBaryonParams& p) : params(p) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out);

    private:
      FlavorConservingFHBaryonParams params;
    };

  } // namespace LalibeFlavorConservingFHBaryonContractionsEnv

} // namespace Chroma

#endif
