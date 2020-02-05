/*
Authors
Arjun Gambhir

INPUT
    Source
OUTPUT
    Propagator solved by shifting spin/color indices of first solution.
*/

#ifndef __lalibe_shifted_spincolor_prop_w_h__
#define __lalibe_shifted_spincolor_prop_w_h__

//Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  namespace LalibeShiftedSpinColorPropagatorEnv
  {
    extern const std::string name;
    bool registerAll();


    struct ShiftedSpinColorParams 
    {
      ShiftedSpinColorParams();
      ShiftedSpinColorParams(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long      frequency;

      struct ShiftedSpinColorProp_t{
	ChromaProp_t prop_param;              //params for next lin solve
      } shiftedscparam ;

      struct NamedObject_t
      {
	std::string                     gauge_id;
	std::string                     src_id;
	std::string                     prop_id;

      } named_obj;
    };

    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const ShiftedSpinColorParams& p) : params(p) {}
      InlineMeas(const InlineMeas& p) : params(p.params) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    private:
      ShiftedSpinColorParams params;
    };

  }

}

#endif
