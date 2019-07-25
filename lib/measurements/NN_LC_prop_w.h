/*
Authors
Arjun Gambhir

First pass at porting the linear combo executable of latscat over to lalibe.
After initial prototype working, lots of optimizations are likely to follow.
Memory management when boosting needs to be looked into.
*/

#ifndef __lalibe_NN_LC_prop_w_h__
#define __lalibe_NN_LC_prop_w_h__

//Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace LalibeNucleonNucleonLinearComboPropagatorEnv
  {
    extern const std::string name;
    bool registerAll();


    struct NNLCPropParams 
    {
      NNLCPropParams();
      NNLCPropParams(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long      frequency;

      // only same flavor bilinear operators are used
      struct NNLCProp_t{
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
      InlineMeas(const NNLCPropParams& p) : params(p) {}
      InlineMeas(const InlineMeas& p) : params(p.params) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      //! Do the measurement
      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    private:
      NNLCPropParams params;
    };

  }

}

#endif
