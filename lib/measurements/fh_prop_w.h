/*
Authors
Arjun Gambhir
Andre Walker-Loud

FH Propagator Task
This computes a Feynman-Hellmann (FH) propagator for bi-linear currents
https://arxiv.org/abs/1612.06963
INPUT
    Propagator
    List of Currents (spin, space, color, momentum)
    Parameters for linear solver
OUTPUT
    FH Propagator for each of the specified currents
*/

#ifndef __lalibe_fh_prop_w_h__
#define __lalibe_fh_prop_w_h__

//Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace LalibeFHPropagatorEnv
  {
    extern const std::string name;
    bool registerAll();


    struct FHParams 
    {
      FHParams();
      FHParams(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long      frequency;

      // only same flavor bilinear operators are used
      struct FHProp_t{
	multi1d<std::string> currents;        //list of currents
	//int  t0;                            //t0 of input prop
	//int j_decay;                        //orthogonal direction of FT
	//Apparently I can read these things from a src prop.
	bool is_mom_max;                      //keeps track of which momentum mode we are using
	int p2_max;                           //max of momentum transfer squared, optional
	multi1d<multi1d<int>> mom_list;       //list of momenta insertions, optional
	multi2d<int> p_list;                  //momentum list the slow fourier transform needs
	ChromaProp_t prop_param;              //params for next lin solve
      } fhparam ;

      struct NamedObject_t
      {
	std::string                     gauge_id;
	std::string                     src_prop_id;
	multi1d<std::string>            fh_prop_id;

      } named_obj;
    };

    //! Compute a spacetime dependent sequential source
    /*! \ingroup inlinehadron */
    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const FHParams& p) : params(p) {}
      InlineMeas(const InlineMeas& p) : params(p.params) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      //! Do the measurement
      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    private:
      FHParams params;
    };

  }

}

#endif
