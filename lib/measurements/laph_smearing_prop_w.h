/*
Authors
Arjun Gambhir
Ben Horz

INPUT
    Propagator or Source
OUTPUT
    Propagator or Source that's projected in the laph smearing subspace.
*/

#ifndef __lalibe_laph_smearing_prop_w_h__
#define __lalibe_laph_smearing_prop_w_h__

//Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  namespace LalibeLaphSmearingPropagatorEnv
  {
    extern const std::string name;
    bool registerAll();


    struct LaphSmearingParams 
    {
      LaphSmearingParams();
      LaphSmearingParams(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long      frequency;

      struct LaphSmearingProp_t{
	std::string                     evec_file;
      } laphparam ;

      struct NamedObject_t
      {
	std::string                     gauge_id;
	std::string                     prop_id;
	std::string                     laph_prop_id;

      } named_obj;
    };

    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const LaphSmearingParams& p) : params(p) {}
      InlineMeas(const InlineMeas& p) : params(p.params) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    private:
      LaphSmearingParams params;
    };

  }

}

#endif
