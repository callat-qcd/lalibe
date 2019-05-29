/*
Authors
Arjun Gambhir

INPUT
    Propagator
    List of Currents (spin, space, color, momentum)
    Parameters for linear solver
    Information about scalar random noise 
    (should match between two bilinears to construct the four quark operator)
OUTPUT
    FH Propagator for each of the specified currents 
    (with random noise used as part of the second source) to be used in four quark contraction
*/

#ifndef __lalibe_stochastic_four_quark_fh_prop_w_h__
#define __lalibe_stochastic_four_quark_fh_prop_w_h__

//Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  namespace LalibeStochasticFourQuarkFHPropagatorEnv
  {
    extern const std::string name;
    bool registerAll();


    struct StochasticFourQuarkFHParams 
    {
      StochasticFourQuarkFHParams();
      StochasticFourQuarkFHParams(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long      frequency;

      // only same flavor bilinear operators are used
      struct StochasticFourQuarkFHProp_t{
	multi1d<std::string> currents;        //list of currents
	//int  t0;                            //t0 of input prop
	//int j_decay;                        //orthogonal direction of FT
	//Apparently I can read these things from a src prop.
	bool is_mom_max;                      //keeps track of which momentum mode we are using
	int p2_max;                           //max of momentum transfer squared, optional
	multi1d<multi1d<int>> mom_list;       //list of momenta insertions, optional
	multi2d<int> p_list;                  //momentum list the slow fourier transform needs
	ChromaProp_t prop_param;              //params for next lin solve
	//Params below for the noise part
	int vector_number;		      //starting vector number, to simplify things, this should be 1 for a new gauge cfg
	Seed ran_seed; 			      //seed value, keeping this the same for stoch_fh and disco is crucial!
	int ZN;                               //the type of random noise
	bool conjugate;                       //this is a switch for conjugating the noise vector
      } stochfourqfhparam ;

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
      InlineMeas(const StochasticFourQuarkFHParams& p) : params(p) {}
      InlineMeas(const InlineMeas& p) : params(p.params) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      //! Do the measurement
      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    private:
      StochasticFourQuarkFHParams params;
    };

  }

}

#endif
