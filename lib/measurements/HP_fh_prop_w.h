/*
Authors
Arjun Gambhir

INPUT
    Propagator
    List of Currents (spin, space, color, momentum)
    Information about scalar random noise 
    (should match with ZN inverter)
OUTPUT
    HP FH Propagator for each of the specified currents 
*/

#ifndef __lalibe_hp_fh_prop_w_h__
#define __lalibe_hp_fh_prop_w_h__

//Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  namespace LalibeHPFHPropagatorEnv
  {
    extern const std::string name;
    bool registerAll();


    struct HPFHParams 
    {
      HPFHParams();
      HPFHParams(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long      frequency;

      struct HPFHProp_t{
	multi1d<std::string> currents;        //list of currents
	bool is_mom_max;                      //keeps track of which momentum mode we are using
	int p2_max;                           //max of momentum transfer squared, optional
	multi1d<multi1d<int>> mom_list;       //list of momenta insertions, optional
	multi2d<int> p_list;                  //momentum list the slow fourier transform needs
	//Params below for the noise part
	//int vector_number;		      //old variable, before a noise loop was added in this measurement
	int starting_vector;		      //starting vector number, to simplify things, this should be 1 for a new gauge cfg
	int ending_vector;		      //ending vector number, to simplify things, this should be 1 for a new gauge cfg
	Seed ran_seed; 			      //seed value, keeping this the same for hp_fh and disco is crucial!
	int ZN;                               //the type of random noise
	bool delete_props;                    // Delete props after being summed, by default this is turned off.
      } hpfhparam ;

      struct NamedObject_t
      {
	std::string                     gauge_id;
	std::string                     src_prop_id;
	multi1d<std::string>            hp_prop_id;
	multi1d<std::string>            fh_prop_id;

      } named_obj;
    };

    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const HPFHParams& p) : params(p) {}
      InlineMeas(const InlineMeas& p) : params(p.params) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    private:
      HPFHParams params;
    };

  }

}

#endif
