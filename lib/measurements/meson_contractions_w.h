/*
 * Single baryon contractions
 * Authors:
 * Arjun Gambhir
 * Andre Walker-Loud
 * Jason Chang
 * David Brantley
 * Do meson contractions and write out the two-point correlator in hdf5
 * Maybe we'll also support sdb, one day...
 */
#ifndef __lalibe_meson_contractions_h__
#define __lalibe_meson_contractions_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  namespace LalibeMesonContractionsEnv
  {
    extern const std::string name;
    bool registerAll();


    struct MesonParams
    {
      MesonParams();
      MesonParams(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long frequency;

      struct Param_t
      {

	bool rotate_to_Dirac;                 //If the correlator is in DeGrand-Rossi basis, set this to true to rotate
					      //to Dirac basis.
	bool output_full_correlator;          //If no momentum is specified, we output the full correlator.
	bool is_mom_max;                      //keeps track of which momentum mode we are using
	int p2_max;                           //max of momentum transfer squared, optional
	multi1d<multi1d<int>> mom_list;       //list of momenta insertions, optional
	multi2d<int> p_list;                  //momentum list the slow fourier transform needs
	multi1d<std::string> particle_list;   //list of actual particles we gunna make from the contractions yo
#ifdef BUILD_HDF5
	std::string file_name;
	std::string obj_path;
#endif
      } param;

      struct NamedObject_t
      {
	//std::string  gauge_id;     //Grab the usual gauge field.
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

      } named_obj;

    };


    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const MesonParams& p) : params(p) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    private:
      MesonParams params;
    };

  } // namespace LalibeMesonContractionsEnv

} // namespace Chroma

#endif
