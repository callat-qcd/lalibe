/*
Authors
Andre Walker-Loud

Add multiple sequential sinks
INPUT
    List of sinks
    value of t_sep
    value of j_decay
OUTPUT
    added sink
*/

#ifndef __inline_coherent_seqsource_h__
#define __inline_coherent_seqsource_h__

//Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma
{
    /*! \ingroup inlinehadron */
    namespace LalibeCoherentSeqsourceEnv
    {
        extern const std::string name;
        bool registerAll();

        struct SinkParams
        {
            SinkParams();
            SinkParams(XMLReader& xml_in,  const std::string& path);
            void writeXML(XMLWriter& xml_out, const std::string& path);

            unsigned long      frequency;

            struct SinkInfo_t
            {
  	            int j_decay;
                int t_sep;  // used to determine what time slice to use for each sink
            } ssparam ;

            struct NamedObject_t
            {
                //std::string  gauge_id;
  	            multi1d<std::string> sink_ids;
                std::string result_sink;
            } named_obj;
        };

        //! Get list of sinks for adding
        class InlineMeas : public AbsInlineMeasurement
        {
        public:
            ~InlineMeas() {}
            InlineMeas(const SinkParams& p): params(p) {}
            InlineMeas(const InlineMeas&  p): params(p.params) {}

            unsigned long getFrequency(void) const {return params.frequency;}

            //! Do the measurement
            void operator()(const unsigned long update_no, XMLWriter& xml_out);

        private:
            SinkParams params;
        };
    }
}

#endif
