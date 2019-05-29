/*
Arjun Gambhir

This computes a Feynman-Hellmann a fully diluted spin/color ZN (stochastic) propagator.
This is useful for disco calculations, but can also be used to intercept a quark line in a hadron to create a stochastic Feyman Hellman propagator.
This code will allow batches of ZN vectors to be inverted at once, but the user must be careful to name them uniquely based on the current "vector number".
INPUT
    Type of ZN noise.
    Random noise seed.
    Starting and ending vector number within the seed.
    Parameters for linear solver
OUTPUT
    ZN propagator, note the header info won't be right, but that's okay for disco and stochastic fh_prop
*/

#ifndef __lalibe_ZN_prop_w_h__
#define __lalibe_ZN_prop_w_h__

//Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace LalibeZNPropagatorEnv
  {
    extern const std::string name;
    bool registerAll();


    struct ZNParams 
    {
      ZNParams();
      ZNParams(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long      frequency;

      struct ZNProp_t{

	int starting_vector;		      //starting vector number, to simplify things, this should be 1 for a new gauge cfg
	int ending_vector;                    //ending vector number
	Seed ran_seed; 			      //seed value, keeping this the same for stoch_fh and disco is crucial!
	int ZN;                               //the type of random noise
	ChromaProp_t prop_param;              //params for next lin solve
      } znparam ;

      struct NamedObject_t
      {
	std::string                     gauge_id;
	multi1d<std::string>            zn_prop_id;

      } named_obj;
    };

    /*! \ingroup inlinehadron */
    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const ZNParams& p) : params(p) {}
      InlineMeas(const InlineMeas& p) : params(p.params) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      //! Do the measurement
      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    private:
      ZNParams params;
    };

  }

}

#endif
