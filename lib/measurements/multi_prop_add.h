/*
Authors
Arjun Gambhir
Andre Walker-Loud

Add multiple propagators
INPUT
    List of Propagators
    optional List of weights
    - if no weights are provide, do simple average
OUTPUT
    added propagator
*/

#ifndef __lalibe_multi_prop_add_h__
#define __lalibe_multi_prop_add_h__

//Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma
{
  /*! \ingroup inlinehadron */
  namespace LalibeMultiPropagatorAddEnv
  {
    extern const std::string name;
    bool registerAll();

    struct PropWeights
    {
      PropWeights();
      PropWeights(XMLReader& xml_in,  const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long      frequency;

      struct weights_t{
  	    multi1d<double> weights;        // optional list of weights
        bool have_weights;
	bool delete_props;                  // Delete props after being summed, by default this is turned off.
      } weights_lst ;

      struct NamedObject_t
      {
  	    multi1d<std::string> prop_ids;
        std::string result_prop;
      } named_obj;
    };

    //! Get list of props and weights for adding
    class InlineMeas : public AbsInlineMeasurement
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const PropWeights& p): params(p) {}
      InlineMeas(const InlineMeas&  p): params(p.params) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      //! Do the measurement
      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out);

    private:
      PropWeights params;
    };

  }

}

#endif
