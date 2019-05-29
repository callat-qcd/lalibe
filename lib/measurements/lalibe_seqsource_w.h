/*
Authors
David Brantley
Arjun Gambhir


INPUT
    Propagator
    List of hadrons to build seqprops.
    Sink time.
    Sink momentum.

OUTPUT
    Sequential source for a given baryon.
*/

#ifndef __lalibe_seqsource_w_h__
#define __lalibe_seqsource_w_h__

//Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma
{
    /*! \ingroup inlinehadron */
    namespace LalibeSeqSourceEnv
    {
        extern const std::string name;
        bool registerAll();


        struct SeqSourceParams
        {
            SeqSourceParams();
            SeqSourceParams(XMLReader& xml_in, const std::string& path);
            void writeXML(XMLWriter& xml_out, const std::string& path);

            unsigned long      frequency;

            // only same flavor bilinear operators are used
            struct SeqSource_t
            {
	            //Apparently I can read these things from a src prop.
	            std::string    particle;     //which particle.
	            std::string    flavor;       //which flavor insertion.
	            std::string    source_spin;  //Source spin state.
	            std::string    sink_spin;    //Sink spin state.
	            multi1d<int>   sink_mom;     //list of momenta insertions, optional
                int            t_0;          //source location
	            int            t_sink;       //which sink time to fix.
	            int            t_sep;        //User can specify tsep instead of t_sink.
                bool           t_all;        //if t_sink and t_sep are not entered, do all t
                bool           spin_zero_meson;  // Helper boolean.
            } ssparam ;

            PropSinkSmear_t    sink_header;


            struct NamedObject_t
            {
	            std::string  gauge_id;
	            std::string  up_quark;
	            bool         is_up;
	            std::string  down_quark;
	            bool         is_down;
	            std::string  strange_quark;
	            bool         is_strange;
	            std::string  charm_quark;
	            bool         is_charm;
                /*  Above are various strings and bools that will identify
                    the quarks and which quarks are present in the input.
	                This could just be done by checking if these strings are
                    null (and ensuring they are with no input),
                    but I like bools, so whatever.
                */
	            std::string  seqsource_id;
                int          num_props;
            } named_obj;
        };

    /*! \ingroup inlinehadron */
    class InlineMeas : public AbsInlineMeasurement
    {
        public:
            ~InlineMeas() {}
            InlineMeas(const SeqSourceParams& p) : params(p) {}
            InlineMeas(const InlineMeas& p) : params(p.params) {}

            unsigned long getFrequency(void) const {return params.frequency;}

            //! Do the measurement
            void operator()(const unsigned long update_no, XMLWriter& xml_out);

        private:
            SeqSourceParams params;
        };

    }

}

#endif
