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

// Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"

// Lalibe Stuff
#include "coherent_seqsource_w.h"
#include "../contractions/seqsource_contractions_func_w.h"
#include "../io/lalibe_qprop_io.h"
#include "lalibe_seqsource_w.h"

namespace Chroma
{
    namespace LalibeCoherentSeqsourceEnv
    {
        namespace
        {
            AbsInlineMeasurement* createMeasurement(XMLReader& xml_in,const std::string& path)
            {
                return new InlineMeas(SinkParams(xml_in, path));
            }
            //! Local registration flag
            bool registered = false;
        }
        const std::string name = "COHERENT_SEQSOURCE";

        //! Register all the factories
        bool registerAll()
        {
            bool success = true;
            if (! registered)
            {
                success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
                registered = true;
            }
            return success;
        }

        //! SinkInfo input
        void read(XMLReader& xml, const std::string& path, SinkParams::SinkInfo_t& input)
        {
            XMLReader inputtop(xml, path);
            read(inputtop, "j_decay", input.j_decay);
            read(inputtop, "t_sep", input.t_sep);
        }

        //! SinkInfo output
        void write(XMLWriter& xml, const std::string& path, const SinkParams::SinkInfo_t& input)
        {
            push(xml, path);
            write(xml, "j_decay", input.j_decay);
            write(xml, "t_sep", input.t_sep);
            pop(xml);
        }

        //! NamedObject input
        void read(XMLReader& xml, const std::string& path, SinkParams::NamedObject_t& input)
        {
            XMLReader inputtop(xml, path);
            //read(inputtop, "gauge_id", input.gauge_id);
            read(inputtop, "sink_ids", input.sink_ids);
            read(inputtop, "result_sink", input.result_sink);
        }

        //! NamedObject output
        void write(XMLWriter& xml, const std::string& path, const SinkParams::NamedObject_t& input)
        {
            push(xml, path);
            //write(xml, "gauge_id", input.gauge_id);
            write(xml, "sink_ids", input.sink_ids);
            write(xml, "result_sink", input.result_sink);
            pop(xml);
        }

        // Param stuff
        SinkParams::SinkParams()
        {
            frequency = 0;
        }

        SinkParams::SinkParams(XMLReader& xml_in, const std::string& path)
        {
            try
            {
                XMLReader paramtop(xml_in, path);
                if (paramtop.count("Frequency") == 1)
                    read(paramtop, "Frequency", frequency);
                else
                    frequency = 1;

                // Read SinkInfo
                read(paramtop, "SinkParams", ssparam);

                // Read in the NamedObject info
                read(paramtop, "NamedObject", named_obj);
            }
            catch(const std::string& e)
            {
                QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
                QDP_abort(1);
            }
        }

        void SinkParams::writeXML(XMLWriter& xml_out, const std::string& path)
        {
            push(xml_out, path);
            write(xml_out, "SinkParams", ssparam);
            write(xml_out, "NamedObject", named_obj);
            pop(xml_out);
        }

        // Function call
        void  InlineMeas::operator()(unsigned long update_no, XMLWriter& xml_out)
        {
            START_CODE();

            StopWatch snoop;
            snoop.reset();
            snoop.start();
            QDPIO::cout << "COHERENT_SEQSOURCE: start" << std::endl;

            int N_sinks = params.named_obj.sink_ids.size();
            QDPIO::cout << "COHERENT_SEQSOURCE: N_sinks = " << N_sinks << std::endl;
            QDPIO::cout << "                    t_sep   = "<<params.ssparam.t_sep << std::endl;
            int t_0;
            int Nt = QDP::Layout::lattSize()[params.ssparam.j_decay];
            int t_sink;

            // Instantiate xml stuff and LatticePropagator; Grab xml record from first prop being added
            XMLReader seqsource_file_xml, seqsource_record_xml;
            LalibeSeqSource_t seqsource_header;
            LatticePropagator multi_sink_result;
            LatticePropagator sink_to_add;
            try
            {
	            // Try the cast to see if this is a valid source
	            sink_to_add =
                    TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.sink_ids[0]);
	            TheNamedObjMap::Instance().get(params.named_obj.sink_ids[0]).getFileXML(seqsource_file_xml);
                TheNamedObjMap::Instance().get(params.named_obj.sink_ids[0]).getRecordXML(seqsource_record_xml);
            }
            catch (std::bad_cast)
            {
	            QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
	            QDP_abort(1);
            }
            catch (const std::string& e)
            {
	           QDPIO::cerr << name << ": error extracting source_header: " << e << std::endl;
	           QDP_abort(1);
            }
            // get t_0 and compute t_sink
            read(seqsource_record_xml, "/SequentialSource/SeqSource", seqsource_header);
            t_0 = seqsource_header.t_0;
            t_sink = t_0 + params.ssparam.t_sep;
            if ( t_sink >= Nt ) t_sink -= Nt; // pos parity, pos t_sep
            if ( t_sink <  0  ) t_sink += Nt; // neg parity, neg t_sep
            QDPIO::cout << "            sink 0: t_0 = "<<t_0<<", t_sink = "<<t_sink<< std::endl;
            // zero out all but t_sink
            projectTimeSlice(sink_to_add, t_sink, params.ssparam.j_decay);
            multi_sink_result = sink_to_add;

            // Start adding props to result
            for(int ni = 1; ni < N_sinks; ni++)
            {
                try
                {
                    sink_to_add = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.sink_ids[ni]);
                        TheNamedObjMap::Instance().get(params.named_obj.sink_ids[ni]).getFileXML(seqsource_file_xml);
                        TheNamedObjMap::Instance().get(params.named_obj.sink_ids[ni]).getRecordXML(seqsource_record_xml);
                    // ADD to result prop
                    read(seqsource_record_xml, "/SequentialSource/SeqSource", seqsource_header);
                    t_0 = seqsource_header.t_0;
                    t_sink = t_0 + params.ssparam.t_sep;
                    if ( t_sink >= Nt ) t_sink -= Nt; // pos parity, pos t_sep
                    if ( t_sink <  0  ) t_sink += Nt; // neg parity, neg t_sep
                    QDPIO::cout << "            sink "<<ni<<": t_0 = "<<t_0<<", t_sink = "<<t_sink<< std::endl;
                    // zero out all but t_sink
                    projectTimeSlice(sink_to_add, t_sink, params.ssparam.j_decay);
                    multi_sink_result += sink_to_add;
                }
                catch (std::bad_cast)
                {
    	            QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
    	            QDP_abort(1);
                }
                catch (const std::string& e)
                {
    	           QDPIO::cerr << name << ": error extracting source_header: " << e << std::endl;
    	           QDP_abort(1);
                }
            }

            /*
            *  Write the a source out to a named buffer
            */
            try
            {
                QDPIO::cout << "Attempt to store resulting coherent_seqsource" << std::endl;

                // Store the seqsource
                TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.result_sink);
                TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.result_sink) = multi_sink_result ;

                // For now, we copy the HACK in seqsource to create the xml to write, L1036-1090
                SeqSource_t seq_source_params;
                GroupXML_t seqsrc;
                std::string tmp;
                //NOTE:AWL I think the <SeqSource> should be removed
                tmp =      "\n  <SeqSource>\n    <SeqSourceType>"+seqsource_header.particle+"</SeqSourceType>\n";
                tmp +=     "    <flavor>"+seqsource_header.flavor+"</flavor>\n";
                tmp +=     "    <source_spin>"+seqsource_header.source_spin+"</source_spin>\n";
                tmp +=     "    <sink_spin>"+seqsource_header.source_spin+"</sink_spin>\n";
                tmp +=     "    <sink_mom>"+std::to_string(seqsource_header.sink_mom[0])+" "+std::to_string(seqsource_header.sink_mom[1])+" "+std::to_string(seqsource_header.sink_mom[2])+"</sink_mom>\n";
                if(seqsource_header.t_all)
                {
                    tmp += "    <t_all>true</t_all>\n";
                }
                else
                {
                    tmp += "    <t_all>false</t_all>\n";
                }
                tmp +=     "    <t_0>"+std::to_string(seqsource_header.t_0)+"</t_0>\n";
                tmp +=     "    <t_sink>"+std::to_string(seqsource_header.t_sink)+"</t_sink>\n";

                // NOTE: for t_sep - we want the input from the COHERENT_SEQSOURCE xml
                tmp +=     "    <t_sep>"+std::to_string(params.ssparam.t_sep)+"</t_sep>\n";
                tmp +=     "    <j_decay>"+std::to_string(params.ssparam.j_decay)+"</j_decay>\n";
                tmp +=     "  </SeqSource>\n";

                seqsrc.xml = tmp;
                seq_source_params.seqsrc = seqsrc;

                // Sequential source header + all other stuff in the xml
                SequentialSource_t seqsource_params_old;
                read(seqsource_record_xml, "/SequentialSource", seqsource_params_old);
                SequentialSource_t new_header;
                new_header.sink_header      = seqsource_params_old.sink_header;
                new_header.seqsource_header = seq_source_params;
                new_header.forward_props    = seqsource_params_old.forward_props;
                new_header.gauge_header     = seqsource_params_old.gauge_header;
                
                /*
                ForwardProp_t orig_header;
                //read(seqsource_record_xml, "/SequentialSource/SeqSourceSinkSmear", orig_header);
                read(seqsource_record_xml, "/SequentialSource", orig_header);
                PropSinkSmear_t    sink_header;
                sink_header   = orig_header.sink_header;
                ForwardProp_t forward_header;
                forward_header = orig_header.forward_props
                //forward_header.prop_header   = orig_header.prop_header;
                //forward_header.source_header = orig_header.source_header;

                SequentialSource_t new_header;
                new_header.sink_header      = sink_header;
                new_header.seqsource_header = seq_source_params;
                new_header.forward_props    = forward_header;
                new_header.gauge_header     = gauge_xml.printCurrentContext();
                read(seqsource_record_xml, "/SequentialSource/SeqSourceSinkSmear", orig_header);
                */
                XMLBufferWriter record_xml;
                write(record_xml, "SequentialSource", new_header);


                TheNamedObjMap::Instance().get(params.named_obj.result_sink).setFileXML(seqsource_file_xml);
                TheNamedObjMap::Instance().get(params.named_obj.result_sink).setRecordXML(record_xml);

                QDPIO::cout << "coherent sink sum successfully stored"  << std::endl;
            }
            catch (std::bad_cast)
            {
                QDPIO::cerr << name << ": dynamic cast error" << std::endl;
                QDP_abort(1);
            }
            catch (const std::string& e)
            {
                QDPIO::cerr << name << ": error storing seqsource: " << e << std::endl;
                QDP_abort(1);
            }

            snoop.stop();
            QDPIO::cout << LalibeCoherentSeqsourceEnv::name << ": total time = " << snoop.getTimeInSeconds() << " secs" << std::endl;
            QDPIO::cout << LalibeCoherentSeqsourceEnv::name<< ": ran successfully" << std::endl;
            END_CODE();
        }
    }// LalibeCoherentSeqsourceEnv
};
