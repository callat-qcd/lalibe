/*
Authors
Arjun Gambhir

INPUT
    Source
OUTPUT
    Propagator solved by gauge transporting other prop.
*/

#ifndef __lalibe_gauge_transported_prop_w_h__
#define __lalibe_gauge_transported_prop_w_h__

//Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  namespace LalibeGaugeTransportedPropagatorEnv
  {
    extern const std::string name;
    bool registerAll();


    struct GaugeTransportedParams 
    {
      GaugeTransportedParams();
      GaugeTransportedParams(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long      frequency;

      struct GaugeTransportedProp_t{
	ChromaProp_t prop_param;              //params for next lin solve
      } gtparam ;

      struct NamedObject_t
      {
	std::string                     gauge_id;
	std::string                     src_id;
	std::string                     prop_id;
	std::string                     shifted_prop_id;

      } named_obj;
    };

    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const GaugeTransportedParams& p) : params(p) {}
      InlineMeas(const InlineMeas& p) : params(p.params) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    private:
      GaugeTransportedParams params;
    };

  }

}

#endif
