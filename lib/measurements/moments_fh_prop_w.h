/*
Authors
Arjun Gambhir

FH Propagator Task
This computes a Feynman-Hellmann (FH) propagator for bi-linear currents with a z^2 insertion (z can be in any direction)
INPUT
    Propagator
    List of Currents (spin, space, color, momentum, direction  of z^2)
    Parameters for linear solver
OUTPUT
    FH Propagator for each of the specified currents
*/

#ifndef __lalibe_moments_fh_prop_w_h__
#define __lalibe_moments_fh_prop_w_h__

//Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  namespace LalibeMomentsFHPropagatorEnv
  {
    extern const std::string name;
    bool registerAll();


    struct MomentsFHParams 
    {
      MomentsFHParams();
      MomentsFHParams(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long      frequency;

      // only same flavor bilinear operators are used
      struct MomentsFHProp_t{
	multi1d<std::string> currents;        //list of currents
	bool is_mom_max;                      //keeps track of which momentum mode we are using
	int moments_direction;                //direction to apply z^2, for actual z the number should be 3
	//int p2_max;                           //max of momentum transfer squared, optional
	//multi1d<multi1d<int>> mom_list;       //list of momenta insertions, optional
	//multi2d<int> p_list;                  //momentum list the slow fourier transform needs
	// Moments FH computes momenta insertion in a different way from the slow fourier transform that is done in fh_prop.
	// This inputs are disabled here.
	ChromaProp_t prop_param;              //params for next lin solve
      } momentsfhparam ;

      struct NamedObject_t
      {
	std::string                     gauge_id;
	std::string                     src_prop_id;
	multi1d<std::string>            fh_prop_id;

      } named_obj;
    };

    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const MomentsFHParams& p) : params(p) {}
      InlineMeas(const InlineMeas& p) : params(p.params) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    private:
      MomentsFHParams params;
    };

  }

}

#endif
