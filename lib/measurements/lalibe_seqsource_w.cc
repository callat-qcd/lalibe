/*
Authors
David Brantley
Arjun Gambhir
11-Feb-2019 AWL adding support for all time-slice seqsrc

This code assumes positive t_sep only for positive parity
                  negative t_sep only for negative parity

INPUT
    Propagator
    List of hadrons to build seqprops.
    Sink time [optional]
    Sink momentum.

OUTPUT
    Sequential source for a given baryon.
*/

// Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/sinks/sink_smearing_factory.h"
#include "meas/sinks/sink_smearing_aggregate.h"
#include "util/info/unique_id.h"
#include "io/qprop_io.h"

// Lalibe Stuff
#include "../contractions/baryon_seqsource_w.h"
#include "lalibe_seqsource_w.h"
#include "../momentum/lalibe_sftmom.h"

namespace Chroma
{
    namespace LalibeSeqSourceEnv
    {
        namespace
        {
            AbsInlineMeasurement* createMeasurement(XMLReader& xml_in,
                const std::string& path)
            {
                return new InlineMeas(SeqSourceParams(xml_in, path));
            }
            //! Local registration flag
            bool registered = false;
        }
        const std::string name = "LALIBE_SEQSOURCE";

        //! Register all the factories
        bool registerAll()
        {
            bool success = true;
            if (! registered)
            {
                success &= QuarkSinkSmearingEnv::registerAll();
                success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
                registered = true;
            }
            return success;
        }

        void read(XMLReader& xml, const std::string& path, SeqSourceParams::SeqSource_t& par)
        {
            XMLReader paramtop(xml, path);
            read(paramtop, "particle" ,par.particle  ); //which particle.
            read(paramtop, "flavor"   ,par.flavor  ); //which flavor insertion.
            if (paramtop.count("t_sep") != 0)
            {
                read(paramtop, "t_sep" ,par.t_sep  ); //Read in relative tsep instead of an absolute t_snk.
                QDPIO::cout<<"Using relative source-sink separation of tsep: "<<par.t_sep<<std::endl;
                par.t_sink = -1; //This is set as a flag to make sure to use tsep.
                par.t_all  = false;
            }
            else if (paramtop.count("t_sink") != 0)
            {
                read(paramtop, "t_sink" ,par.t_sink  ); //which sink time to fix.
                QDPIO::cout<<"Using absolute timeslice for sink : " <<par.t_sink<<std::endl;
                par.t_all = false;
                // set t_sep=0 and compute later
                par.t_sep = 0;
            }
            else
            {
                par.t_all  = true;
                par.t_sink = -1;
                // set t_sep = 0, later we'll set to Nt to trigger contractions
                par.t_sep  = 0;
            }

            if (paramtop.count("sink_mom") != 0)
            {
                read(paramtop, "sink_mom" ,par.sink_mom);
                QDPIO::cout<<"Fixing momentum of the sink to be: "<<"px: "<<par.sink_mom[0]<<" py:"<<par.sink_mom[1]<<" pz: "<<par.sink_mom[2]<<std::endl;
            }
            else
            {
                QDPIO::cout << "Momentum of the sink must be specified, aborting.... "<<std::endl;
                QDP_abort(1);
            }
            par.spin_zero_meson = false;

            if ( par.particle == "piplus" )
            {
                par.spin_zero_meson = true;
                QDPIO::cout<<"Creating a spin zero meson sequential source. Ignoring spin tags. "<<std::endl;
            }
            else
            {
                read(paramtop, "source_spin" ,par.source_spin  ); //Source spin state.
                read(paramtop, "sink_spin" ,par.sink_spin  ); //Source spin state.
            }
        }// END read

        void write(XMLWriter& xml, const std::string& path, SeqSourceParams::SeqSource_t& par)
        {
            push(xml, path);

            write(xml, "particle" ,par.particle); //which particle
            write(xml, "flavor" ,par.flavor  ); //which flavor insertion.
            //ASG: This will write -1 if tsep is used instead.
            write(xml, "t_all",  "true" );
            write(xml, "t_0",    par.t_0);
            write(xml, "t_sink", par.t_sink ); //which sink time to fix.
            write(xml, "t_sep",  par.t_sep );

            if ( par.spin_zero_meson == false )
            {
                write(xml, "source_spin" ,par.source_spin  ); //Source spin state.
                write(xml, "sink_spin" ,par.sink_spin  ); //Sink spin state.
            }
            pop(xml);
        }// END write

        //! NamedObject input
        void read(XMLReader& xml, const std::string& path, SeqSourceParams::NamedObject_t& input)
        {
            input.num_props = 0;
            XMLReader inputtop(xml, path);
            read(inputtop, "gauge_id"     , input.gauge_id);
            //Logic to read quark propagators (whichever ones are present.)
            if (inputtop.count("up_quark") != 0)
            {
                read(inputtop, "up_quark" ,input.up_quark);
                QDPIO::cout<<"I found an up quark, here is its id: "<<input.up_quark<<std::endl;
                input.is_up = true;
                input.num_props += 1;
            }
            else
            {
                QDPIO::cout<<"I couldn't find an up quark, hope you don't need it for the inputted baryon contractions. "<<std::endl;
                input.is_up = false;
            }
            if (inputtop.count("down_quark") != 0)
            {
                read(inputtop, "down_quark" ,input.down_quark);
                QDPIO::cout<<"I found a down quark, here is its id: "<<input.down_quark<<std::endl;
                input.is_down = true;
                input.num_props += 1;
            }
            else
            {
                QDPIO::cout<<"I couldn't find a down quark, hope you don't need it for the inputted baryon contractions. "<<std::endl;
                input.is_down = false;
            }
            if (inputtop.count("strange_quark") != 0)
            {
                read(inputtop, "strange_quark" ,input.up_quark);
                QDPIO::cout<<"I found an strange quark, here is its id: "<<input.strange_quark<<std::endl;
                input.is_strange = true;
                input.num_props += 1;
            }
            else
            {
                QDPIO::cout<<"I couldn't find a strange quark, hope you don't need it for the inputted baryon contractions. "<<std::endl;
                input.is_strange = false;
            }
            if (inputtop.count("charm_quark") != 0)
            {
                read(inputtop, "charm_quark" ,input.charm_quark);
                QDPIO::cout<<"I found an charm quark, here is its id: "<<input.charm_quark<<std::endl;
                input.is_charm = true;
                input.num_props += 1;
            }
            else
            {
                QDPIO::cout<<"I couldn't find a charm quark, hope you don't need it for the inputted baryon contractions. "<<std::endl;
                input.is_charm = false;
            }
            read(inputtop, "seqsource_id" ,input.seqsource_id);
        } // END read NamedObject

        //! NamedObject output
        void write(XMLWriter& xml, const std::string& path, const SeqSourceParams::NamedObject_t& input)
        {
            push(xml, path);
            write(xml, "gauge_id", input.gauge_id    );

            if(input.is_up == true)
                write(xml, "up_quark" ,input.up_quark);
            if(input.is_down == true)
                write(xml, "down_quark" ,input.down_quark);
            if(input.is_strange == true)
                write(xml, "strange_quark" ,input.strange_quark);
            if(input.is_charm == true)
                write(xml, "charm_quark" ,input.charm_quark);

            write(xml, "seqsource_id" ,input.seqsource_id);
            pop(xml);
        }// END write NamedObject

        // Param stuff
        SeqSourceParams::SeqSourceParams()
        {
            frequency = 0;
        }

        SeqSourceParams::SeqSourceParams(XMLReader& xml_in, const std::string& path)
        {
            try
            {
                XMLReader paramtop(xml_in, path);
                if (paramtop.count("Frequency") == 1)
                    read(paramtop, "Frequency", frequency);
                else
                    frequency = 1;

                // Parameters for source construction
                read(paramtop, "SeqSourceParams", ssparam);

                // Read in the NamedObject info
                read(paramtop, "NamedObject", named_obj);
            }
            catch(const std::string& e)
            {
                QDPIO::cerr << __func__ << ": Caught Exception reading XML: "
                            << e << std::endl;
                QDP_abort(1);
            }
        }

        void SeqSourceParams::writeXML(XMLWriter& xml_out, const std::string& path)
        {
            push(xml_out, path);
            write(xml_out, "SeqSourceParams", ssparam);
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
            QDPIO::cout << "LALIBE_SEQSOURCE: start" << std::endl;

            // Test and grab a reference to the gauge field

            XMLBufferWriter gauge_xml;
            try
            {
                TheNamedObjMap::Instance().getData
                    <multi1d <LatticeColorMatrix> >(params.named_obj.gauge_id);
                TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
            }
            catch( std::bad_cast )
            {
                QDPIO::cerr << LalibeSeqSourceEnv::name
                    << ": caught dynamic cast error" << std::endl;
                QDP_abort(1);
            }
            catch (const std::string& e)
            {
                QDPIO::cerr << LalibeSeqSourceEnv::name
                    << ": map call failed: " << e << std::endl;
                QDP_abort(1);
            }
            const multi1d<LatticeColorMatrix>& u =
                TheNamedObjMap::Instance().getData
                    <multi1d <LatticeColorMatrix> >(params.named_obj.gauge_id);

            push(xml_out, "seqsource");
            write(xml_out, "update_no", update_no);

            // Write out the input
            params.writeXML(xml_out, "Input");

            // Write out the config header
            write(xml_out, "Config_info", gauge_xml);

            push(xml_out, "Seqsource_version");
            write(xml_out, "seq_version", 1);
            pop(xml_out);

            // Grab all the propagators and do some sanity checks...
            XMLReader up_prop_file_xml, up_prop_record_xml;
            LatticePropagator up_quark_propagator;
            XMLReader down_prop_file_xml, down_prop_record_xml;
            LatticePropagator down_quark_propagator;
            XMLReader strange_prop_file_xml, strange_prop_record_xml;
            LatticePropagator strange_quark_propagator;
            XMLReader charm_prop_file_xml, charm_prop_record_xml;
            LatticePropagator charm_quark_propagator;

            int num_props = params.named_obj.num_props;

            // Chroma sequential prop struct demands a multi1d<ForwardProp_t> structure containing all forward prop headers.
            multi1d<ForwardProp_t> forward_headers(num_props);

            // Bools for determining quark smearing type.
            multi1d<bool> smeared(num_props);

            //Need origin, j_decay, and t0 for fourier transform!
            //Need j_decay of bc to know what comes with a minus sign.
            int j_decay;
            int t_0;
            multi1d<int> bc(num_props);
            multi1d<int> origin(num_props);

            // For sanity checks.
            multi1d<int> j_vec(num_props);
            multi1d<multi1d<int>> origin_vec(num_props);

            if(params.named_obj.is_up == true)
            {
                QDPIO::cout << "Attempting to read up propagator" << std::endl;
                try
                {
                    up_quark_propagator =
                    TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.up_quark);
                    TheNamedObjMap::Instance().get(params.named_obj.up_quark).getFileXML(up_prop_file_xml);
                    TheNamedObjMap::Instance().get(params.named_obj.up_quark).getRecordXML(up_prop_record_xml);
                    //Get the origin  and j_decay for the FT, this assumes all quarks have the same origin.
                    //            MakeSourceProp_t  orig_header;

                    if (up_prop_record_xml.count("/Propagator") != 0)
                    {
                        Propagator_t orig_header;

                        QDPIO::cout<<"Up quark propagator is unsmeared, reading from Propagator tag..."<<std::endl;
                        read(up_prop_record_xml, "/Propagator", orig_header);

                        smeared[0] = false;

                        j_decay = orig_header.source_header.j_decay;

                        t_0 = orig_header.source_header.t_source;

                        origin = orig_header.source_header.getTSrce();

                        bc = getFermActBoundary(orig_header.prop_header.fermact);

                        // for future sanity check.
                        j_vec[0] = j_decay;

                        origin_vec[0] = origin;

                        // Save headers.
                        forward_headers[0].prop_header   = orig_header.prop_header;
                        forward_headers[0].source_header = orig_header.source_header;
                        forward_headers[0].gauge_header  = orig_header.gauge_header;

                    }
                    else if (up_prop_record_xml.count("/SinkSmear") != 0)
                    {
                        //QDPIO::cout<<"Please pass in an UNSMEARED propagator. Thank you."<<std::endl;
                        ForwardProp_t orig_header;

                        QDPIO::cout<<"Up quark propagator is smeared, reading from SinkSmear tag..."<<std::endl;
                        read(up_prop_record_xml, "/SinkSmear", orig_header);
                        smeared[0] = true;

                        //smear_header.sink_header   = orig_header.sink_header;
                        //smear_header.prop_header   = orig_header.prop_header;
                        //smear_header.source_header = orig_header.source_header;
                        //smear_header.gauge_header  = orig_header.gauge_header;

                        j_decay = orig_header.source_header.j_decay;
                        t_0 = orig_header.source_header.t_source;
                        origin = orig_header.source_header.getTSrce();

                        bc = getFermActBoundary(orig_header.prop_header.fermact);

                        // for future sanity check.
                        j_vec[0] = j_decay;
                        origin_vec[0] = origin;

                        // Save headers.

                        forward_headers[0].prop_header   = orig_header.prop_header;
                        forward_headers[0].source_header = orig_header.source_header;
                        forward_headers[0].gauge_header  = orig_header.gauge_header;

                        params.sink_header   = orig_header.sink_header;
                    }
                    else
                    {
                        QDPIO::cout<<"What type of weird propagator did you give me? I can't find the right tag to get j_decay, source location, and whatnot...exiting..."<<std::endl;
                    }

                }
                catch (std::bad_cast)
                {
                    QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
                    QDP_abort(1);
                }
                catch (const std::string& e)
                {
                    QDPIO::cerr << name << ": error reading src prop_header: "
                    << e << std::endl;
                    QDP_abort(1);
                }
            }

            if(params.named_obj.is_down == true)
            {
                QDPIO::cout << "Attempting to read down propagator" << std::endl;
                try
                {
                    down_quark_propagator = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.down_quark);
                    TheNamedObjMap::Instance().get(params.named_obj.down_quark).getFileXML(down_prop_file_xml);
                    TheNamedObjMap::Instance().get(params.named_obj.down_quark).getRecordXML(down_prop_record_xml);
                    //Get the origin  and j_decay for the FT, this assumes all quarks have the same origin.
                    //MakeSourceProp_t  orig_header;

                    if (down_prop_record_xml.count("/Propagator") != 0)
                    {
                        Propagator_t orig_header;

                        QDPIO::cout<<"Down quark propagator is unsmeared, reading from Propagator tag..."<<std::endl;
                        read(down_prop_record_xml, "/Propagator", orig_header);
                        smeared[1] = false;

                        j_decay = orig_header.source_header.j_decay;
                        t_0 = orig_header.source_header.t_source;
                        origin = orig_header.source_header.getTSrce();
                        bc = getFermActBoundary(orig_header.prop_header.fermact);

                        // for future sanity check.
                        j_vec[1] = j_decay;
                        origin_vec[1] = origin;

                        // Save headers.
                        forward_headers[1].prop_header   = orig_header.prop_header;
                        forward_headers[1].source_header = orig_header.source_header;
                        forward_headers[1].gauge_header  = orig_header.gauge_header;
                    }
                    else if (down_prop_record_xml.count("/SinkSmear") != 0)
                    {
                        //QDPIO::cout<<"Please pass in an UNSMEARED propagator. Thank you."<<std::endl;
                        ForwardProp_t orig_header;

                        QDPIO::cout<<"Down quark propagator is smeared, reading from SinkSmear tag..."<<std::endl;
                        read(down_prop_record_xml, "/SinkSmear", orig_header);
                        smeared[1] = true;

                        //smear_header.sink_header   = orig_header.sink_header;
                        //smear_header.prop_header   = orig_header.prop_header;
                        //smear_header.source_header = orig_header.source_header;
                        //smear_header.gauge_header  = orig_header.gauge_header;

                        j_decay = orig_header.source_header.j_decay;
                        t_0 = orig_header.source_header.t_source;
                        origin = orig_header.source_header.getTSrce();
                        bc = getFermActBoundary(orig_header.prop_header.fermact);

                        // for future sanity check.
                        j_vec[1] = j_decay;
                        origin_vec[1] = origin;

                        // Save headers.
                        forward_headers[1].prop_header   = orig_header.prop_header;
                        forward_headers[1].source_header = orig_header.source_header;
                        forward_headers[1].gauge_header  = orig_header.gauge_header;
                        params.sink_header   = orig_header.sink_header;
                    }
                    else
                    {
                        QDPIO::cout<<"What type of weird propagator did you give me? I can't find the right tag to get j_decay, source location, and whatnot...exiting..."<<std::endl;
                    }

                }
                catch (std::bad_cast)
                {
                    QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
                    QDP_abort(1);
                }
                catch (const std::string& e)
                {
                    QDPIO::cerr << name << ": error reading src prop_header: "
                    << e << std::endl;
                    QDP_abort(1);
                }
            }

            if(params.named_obj.is_strange == true)
            {
                QDPIO::cout << "Attempting to read strange propagator" << std::endl;
                try
                {
                    strange_quark_propagator = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.strange_quark);
                    TheNamedObjMap::Instance().get(params.named_obj.strange_quark).getFileXML(strange_prop_file_xml);
                    TheNamedObjMap::Instance().get(params.named_obj.strange_quark).getRecordXML(strange_prop_record_xml);
                    //Get the origin  and j_decay for the FT, this assumes all quarks have the same origin.
                    //            MakeSourceProp_t  orig_header;
                    ForwardProp_t orig_header;

                    if (strange_prop_record_xml.count("/Propagator") != 0)
                    {
                        QDPIO::cout<<"Strange quark propagator is unsmeared, reading from Propagator tag..."<<std::endl;
                        read(strange_prop_record_xml, "/Propagator", orig_header);
                        smeared[2] = false;

                        j_decay = orig_header.source_header.j_decay;
                        t_0 = orig_header.source_header.t_source;
                        origin = orig_header.source_header.getTSrce();
                        bc = getFermActBoundary(orig_header.prop_header.fermact);

                        // for future sanity check.
                        j_vec[2] = j_decay;
                        origin_vec[2] = origin;

                        // Save headers.
                        forward_headers[2].prop_header   = orig_header.prop_header;
                        forward_headers[2].source_header = orig_header.source_header;
                        forward_headers[2].gauge_header  = orig_header.gauge_header;
                    }
                    else if (strange_prop_record_xml.count("/SinkSmear") != 0)
                    {
                        QDPIO::cout<<"Please pass in an UNSMEARED propagator. Thank you."<<std::endl;
                        QDPIO::cout<<"Strange quark propagator is smeared, reading from SinkSmear tag..."<<std::endl;
                        read(strange_prop_record_xml, "/SinkSmear", orig_header);
                        smeared[2] = true;

                        j_decay = orig_header.source_header.j_decay;
                        t_0 = orig_header.source_header.t_source;
                        origin = orig_header.source_header.getTSrce();
                        bc = getFermActBoundary(orig_header.prop_header.fermact);

                        // for future sanity check.
                        j_vec[2] = j_decay;
                        origin_vec[2] = origin;

                        // Save headers.
                        forward_headers[2].prop_header   = orig_header.prop_header;
                        forward_headers[2].source_header = orig_header.source_header;
                        forward_headers[2].gauge_header  = orig_header.gauge_header;

                        //              smear_header.sink_header   = orig_header.sink_header;
                        //              smear_header.prop_header   = orig_header.prop_header;
                        //              smear_header.source_header = orig_header.source_header;
                        //              smear_header.gauge_header  = orig_header.gauge_header;

                        params.sink_header   = orig_header.sink_header;
                    }
                    else
                    {
                        QDPIO::cout<<"What type of weird propagator did you give me? I can't find the right tag to get j_decay, source location, and whatnot...exiting..."<<std::endl;
                    }
                }
                catch (std::bad_cast)
                {
                    QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
                    QDP_abort(1);
                }
                catch (const std::string& e)
                {
                    QDPIO::cerr << name << ": error reading src prop_header: "
                    << e << std::endl;
                    QDP_abort(1);
                }
            }

            if(params.named_obj.is_charm == true)
            {
                QDPIO::cout << "Attempting to read charm propagator" << std::endl;
                try
                {
                    charm_quark_propagator = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.charm_quark);
                    TheNamedObjMap::Instance().get(params.named_obj.charm_quark).getFileXML(charm_prop_file_xml);
                    TheNamedObjMap::Instance().get(params.named_obj.charm_quark).getRecordXML(charm_prop_record_xml);
                    //Get the origin  and j_decay for the FT, this assumes all quarks have the same origin.
                    //            MakeSourceProp_t  orig_header;

                    ForwardProp_t orig_header;
                    if (charm_prop_record_xml.count("/Propagator") != 0)
                    {
                        QDPIO::cout<<"Charm quark propagator is unsmeared, reading from Propagator tag..."<<std::endl;
                        read(charm_prop_record_xml, "/Propagator", orig_header);
                        smeared[3] = false;

                        j_decay = orig_header.source_header.j_decay;
                        t_0 = orig_header.source_header.t_source;
                        origin = orig_header.source_header.getTSrce();
                        bc = getFermActBoundary(orig_header.prop_header.fermact);

                        // for future sanity check.
                        j_vec[3] = j_decay;
                        origin_vec[3] = origin;

                        // Save headers.
                        forward_headers[3].prop_header   = orig_header.prop_header;
                        forward_headers[3].source_header = orig_header.source_header;
                        forward_headers[3].gauge_header  = orig_header.gauge_header;
                    }
                    else if (charm_prop_record_xml.count("/SinkSmear") != 0)
                    {
                        //QDPIO::cout<<"Please pass in an UNSMEARED propagator. Thank you."<<std::endl;
                        QDPIO::cout<<"Charm quark propagator is smeared, reading from SinkSmear tag..."<<std::endl;
                        read(charm_prop_record_xml, "/SinkSmear", orig_header);
                        smeared[3] = true;

                        j_decay = orig_header.source_header.j_decay;
                        t_0 = orig_header.source_header.t_source;
                        origin = orig_header.source_header.getTSrce();
                        bc = getFermActBoundary(orig_header.prop_header.fermact);

                        // for future sanity check.
                        j_vec[3] = j_decay;
                        origin_vec[3] = origin;

                        // Save headers.
                        forward_headers[3].prop_header   = orig_header.prop_header;
                        forward_headers[3].source_header = orig_header.source_header;
                        forward_headers[3].gauge_header  = orig_header.gauge_header;

                        //smear_header.sink_header   = orig_header.sink_header;
                        //smear_header.prop_header   = orig_header.prop_header;
                        //smear_header.source_header = orig_header.source_header;
                        //smear_header.gauge_header  = orig_header.gauge_header;

                        params.sink_header   = orig_header.sink_header;
                    }
                    else
                    {
                        QDPIO::cout<<"What type of weird propagator did you give me? I can't find the right tag to get j_decay, source location, and whatnot...exiting..."<<std::endl;
                    }
                }
                catch (std::bad_cast)
                {
                    QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
                    QDP_abort(1);
                }
                catch (const std::string& e)
                {
                    QDPIO::cerr << name << ": error reading src prop_header: "
                    << e << std::endl;
                    QDP_abort(1);
                }
            }

            // Do sanity checks to ensure all propagators are acceptable.
            for( int n = 0; n < num_props; n = n + 1 ) {
                if( j_vec[n] != j_decay || origin_vec[n][0] != origin[0] || origin_vec[n][1] != origin[1] || origin_vec[n][2] != origin[2] || origin_vec[n][3] != origin[3])
                {
                    QDPIO::cout<<"Propagator source or j_decay parameters do not match. Cannot continue..."<<std::endl;
                    QDP_abort(1);
                }
            }
            // Define Nt
            int Nt = QDP::Layout::lattSize()[j_decay];
            params.ssparam.t_0 = t_0;
            int parity = 0;
            if(params.ssparam.particle.substr( params.ssparam.particle.length() -3 ) == "_np")
                parity = 1;

            if(params.ssparam.t_all == true)
            {
                QDPIO::cout<< "t_all = " << params.ssparam.t_all << std::endl;
                QDPIO::cout<< "      will compute seqsource on all t" << std::endl;
                params.ssparam.t_sink = -1;
                params.ssparam.t_sep  = Nt;
            }
            else
            {
                //If relative tsep was used, t_sink is calculated accounting for the lattice boundary below.
                if(params.ssparam.t_sink == -1)
                {
                    params.ssparam.t_sink = t_0 + params.ssparam.t_sep;
                    if (params.ssparam.t_sink >= Nt)
                        params.ssparam.t_sink -= Nt; // Subtract Nt if we go over.
                    else if (params.ssparam.t_sink < 0)
                        params.ssparam.t_sink += Nt; // Add Nt if we are negative.
                }
                else
                {
                    if (parity == 0 )
                    {
                        params.ssparam.t_sep = params.ssparam.t_sink - t_0;
                    }
                    else if ( parity == 1 )
                    {
                        params.ssparam.t_sep = t_0 - params.ssparam.t_sink;
                    }
                    else
                    {
                        QDPIO::cout << "parity = "<< parity << "not recognized" << std::endl;
                        QDP_abort(1);
                    }
                }
                QDPIO::cout<<"Timeslice of sink is: "<<params.ssparam.t_sink<<std::endl;

                // A sanity check
                if (params.ssparam.t_sink < 0 || params.ssparam.t_sink >= Nt)
                {
                    QDPIO::cout << "Sink time coordinate greater than Nt = "<<Nt<< std::endl;
                    QDPIO::cout << "t_sink = " << params.ssparam.t_sink << std::endl;
                    QDPIO::cout << "Mod the source sink displacement to account for lattice boundary."<< std::endl;
                    QDP_abort(1);
                }
            }

            // Another sanity check, this time for smearing.
            bool smear_sink = smeared[0];

            for( int n = 0; n < num_props; n = n + 1 )
            {
                if( smear_sink != smeared[n] )
                {
                    QDPIO::cout << "Propagator smearing types do not match. All propagators must be equivalently smeared." << std::endl;
                    QDP_abort(1);
                }
            }

            // ***************    Now we really start working    ******************

            QDPIO::cout <<"Attempting to construct sequential source for the "
                        <<params.ssparam.particle<<" with flavor bilinear insertion "
                        <<params.ssparam.flavor<<std::endl;

            LatticePropagator seqsource;

            if(params.ssparam.particle == "proton")
            {
                //int parity = 0; // This is the normal proton, with neg_par turned off.
                // Execute some proton seqsource contractions.
                if(params.ssparam.flavor == "DD")
                {
                    QDPIO::cout <<"Starting "<<params.ssparam.particle<< " "
                                <<params.ssparam.flavor<<" sequential source construction."
                                <<std::endl;

                    QDPIO::cout <<params.ssparam.particle
                                <<" wave function source spin "<<params.ssparam.source_spin
                                <<" to sink spin "<<params.ssparam.sink_spin<<std::endl;

                    QDPIO::cout << "Sink momentum px: "<<params.ssparam.sink_mom[0]
                                <<" py: "<<params.ssparam.sink_mom[1]
                                <<" pz: "<<params.ssparam.sink_mom[2]
                                <<std::endl;

                    seqsource = ProtDtoD
                                (
                                    up_quark_propagator,
                                    params.ssparam.source_spin,
                                    params.ssparam.sink_spin,
                                    params.ssparam.sink_mom,
                                    origin,
                                    bc,
                                    params.ssparam.t_sink,
                                    j_decay,
                                    parity,
                                    params.ssparam.t_all
                                );
                }
                else if(params.ssparam.flavor == "UU")
                {
                    QDPIO::cout <<"Starting "<<params.ssparam.particle<< " "
                                <<params.ssparam.flavor<<" sequential source construction."
                                <<std::endl;
                    QDPIO::cout <<params.ssparam.particle
                                <<" wave function source spin "<<params.ssparam.source_spin
                                <<" to sink spin "<<params.ssparam.sink_spin
                                <<std::endl;

                    QDPIO::cout <<"Sink momentum px: "<<params.ssparam.sink_mom[0]
                                <<" py: "<<params.ssparam.sink_mom[1]
                                <<" pz: "<<params.ssparam.sink_mom[2]
                                <<std::endl;

                    seqsource = ProtUtoU
                                (
                                    up_quark_propagator,
                                    down_quark_propagator,
                                    params.ssparam.source_spin,
                                    params.ssparam.sink_spin,
                                    params.ssparam.sink_mom,
                                    origin,
                                    bc,
                                    params.ssparam.t_sink,
                                    j_decay,
                                    parity,
                                    params.ssparam.t_all
                                );
                }
                else if(params.ssparam.flavor == "UD")
                {
                    QDPIO::cout <<"Starting "<<params.ssparam.particle<< " "
                                <<params.ssparam.flavor<<" sequential source construction."
                                <<std::endl;
                    QDPIO::cout <<params.ssparam.particle<<" wave function source spin "
                                <<params.ssparam.source_spin<<" to sink spin "
                                <<params.ssparam.sink_spin
                                <<std::endl;

                    QDPIO::cout <<"Sink momentum px: "<<params.ssparam.sink_mom[0]
                                <<" py: "<<params.ssparam.sink_mom[1]
                                <<" pz: "<<params.ssparam.sink_mom[2]
                                <<std::endl;
                }
                else
                {
                    QDPIO::cout <<params.ssparam.particle<<" "<<params.ssparam.flavor
                                <<" flavor structure not allowed. Skipping..."
                                <<std::endl;
                }
            }
            else if(params.ssparam.particle == "proton_np")
            {
                //int parity = 1; // This is the time reversed proton, with neg_par turned on.
                // For proton_np, we call the same contractions as proton, but with the parity int switched to 1.
                if(params.ssparam.flavor == "DD")
                {
                    QDPIO::cout <<"Starting "<<params.ssparam.particle<< " "
                                <<params.ssparam.flavor<<" sequential source construction."
                                <<std::endl;

                    QDPIO::cout <<params.ssparam.particle<<" wave function source spin "
                                <<params.ssparam.source_spin
                                <<" to sink spin "<<params.ssparam.sink_spin
                                <<std::endl;

                    QDPIO::cout <<"Sink momentum px: "<<params.ssparam.sink_mom[0]
                                <<" py: "<<params.ssparam.sink_mom[1]
                                <<" pz: "<<params.ssparam.sink_mom[2]
                                <<std::endl;

                    seqsource = ProtDtoD
                                (
                                    up_quark_propagator,
                                    params.ssparam.source_spin,
                                    params.ssparam.sink_spin,
                                    params.ssparam.sink_mom,
                                    origin,
                                    bc,
                                    params.ssparam.t_sink,
                                    j_decay,
                                    parity,
                                    params.ssparam.t_all
                                );
                }
                else if(params.ssparam.flavor == "UU")
                {
                    QDPIO::cout <<"Starting "<<params.ssparam.particle<< " "
                                <<params.ssparam.flavor<<" sequential source construction."
                                <<std::endl;
                    QDPIO::cout <<params.ssparam.particle<<" wave function source spin "
                                <<params.ssparam.source_spin<<" to sink spin "
                                <<params.ssparam.sink_spin
                                <<std::endl;
                    QDPIO::cout <<"Sink momentum px: "<<params.ssparam.sink_mom[0]
                                <<" py: "<<params.ssparam.sink_mom[1]
                                <<" pz: "<<params.ssparam.sink_mom[2]
                                <<std::endl;

                    seqsource = ProtUtoU
                                (
                                    up_quark_propagator,
                                    down_quark_propagator,
                                    params.ssparam.source_spin,
                                    params.ssparam.sink_spin,
                                    params.ssparam.sink_mom,
                                    origin,
                                    bc,
                                    params.ssparam.t_sink,
                                    j_decay,
                                    parity,
                                    params.ssparam.t_all
                                );
                }
                else if(params.ssparam.flavor == "UD")
                {
                    QDPIO::cout <<"Starting "<<params.ssparam.particle<< " "
                                <<params.ssparam.flavor<<" sequential source construction."
                                <<std::endl;
                    QDPIO::cout <<params.ssparam.particle<<" wave function source spin "
                                <<params.ssparam.source_spin<<" to sink spin "
                                <<params.ssparam.sink_spin
                                <<std::endl;

                    QDPIO::cout <<"Sink momentum px: "<<params.ssparam.sink_mom[0]
                                <<" py: "<<params.ssparam.sink_mom[1]
                                <<" pz: "<<params.ssparam.sink_mom[2]
                                <<std::endl;
                }
                else
                {
                    QDPIO::cout <<params.ssparam.particle<<" "<<params.ssparam.flavor
                                <<" flavor structure not allowed. Skipping..."
                                <<std::endl;
                }
            }
            else if(params.ssparam.particle == "neutron")
            {
                // Execute some neutron seqsource contractions.
                if(params.ssparam.flavor == "DD")
                {
                    QDPIO::cout <<"Starting "<<params.ssparam.particle<< " "
                                <<params.ssparam.flavor<<" sequential source construction."
                                <<std::endl;
                    QDPIO::cout <<params.ssparam.particle<<" wave function source spin "
                                <<params.ssparam.source_spin<<" to sink spin "
                                <<params.ssparam.sink_spin
                                <<std::endl;

                    QDPIO::cout <<"Sink momentum px: "<<params.ssparam.sink_mom[0]
                                <<" py: "<<params.ssparam.sink_mom[1]
                                <<" pz: "<<params.ssparam.sink_mom[2]
                                <<std::endl;
                }
                else if(params.ssparam.flavor == "UU")
                {
                    QDPIO::cout <<"Starting "<<params.ssparam.particle<< " "
                                <<params.ssparam.flavor<<" sequential source construction."
                                <<std::endl;
                    QDPIO::cout <<params.ssparam.particle<<" wave function source spin "
                                <<params.ssparam.source_spin<<" to sink spin "
                                <<params.ssparam.sink_spin
                                <<std::endl;

                    QDPIO::cout <<"Sink momentum px: "<<params.ssparam.sink_mom[0]
                                <<" py: "<<params.ssparam.sink_mom[1]
                                <<" pz: "<<params.ssparam.sink_mom[2]
                                <<std::endl;
                }
                else if(params.ssparam.flavor == "UD")
                {
                    QDPIO::cout <<"Starting "<<params.ssparam.particle<< " "
                                <<params.ssparam.flavor<<" sequential source construction."
                                <<std::endl;
                    QDPIO::cout <<params.ssparam.particle<<" wave function source spin "
                                <<params.ssparam.source_spin<<" to sink spin "
                                <<params.ssparam.sink_spin
                                <<std::endl;

                    QDPIO::cout <<"Sink momentum px: "<<params.ssparam.sink_mom[0]
                                <<" py: "<<params.ssparam.sink_mom[1]
                                <<" pz: "<<params.ssparam.sink_mom[2]
                                <<std::endl;
                }
                else
                {
                    QDPIO::cout <<"Neutron "<<params.ssparam.flavor
                                <<" flavor structure not allowed. Skipping..."
                                <<std::endl;
                }
            }
            else if(params.ssparam.particle == "piplus")
            {
                // Execute some pion seqsource contractions.
                QDPIO::cout<<"Starting "<<params.ssparam.particle<< " "<<params.ssparam.flavor<<" sequential source construction."<<std::endl;
                QDPIO::cout<<params.ssparam.particle<<" wave function source spin "<<params.ssparam.source_spin<<" to sink spin "<<params.ssparam.sink_spin<<std::endl;

                QDPIO::cout<<"Sink momentum px: "<<params.ssparam.sink_mom[0]<<" py: "<<params.ssparam.sink_mom[1]<<" pz: "<<params.ssparam.sink_mom[2]<<std::endl;
                QDPIO::cout<<"You beat me to it. Don't have these ready yet."<<std::endl;
            }
            else
            {
                QDPIO::cout<<"Particle "<<params.ssparam.particle<<" unknown. Aborting."<<std::endl;
                QDP_abort(1);
            }

            // Sink smear the sequential source if necessary.
            QDPIO::cout << "Smearing the sequential source." << std::endl;

            // ************* Sink smear the propagators ***************.

            // Grab a copy of the sink smear xml and create sink smear object if props are smeared.
            if( smear_sink == true ){
                QDPIO::cout << "Propagators are smeared -> Smearing the sequential source. " << std::endl;
                QDPIO::cout << "Pulling smearing parameters from input propagators. " << std::endl;

                std::istringstream  xml_s(params.sink_header.sink.xml);
                XMLReader  sinktop(xml_s);
                QDPIO::cout << "Smearing parameters succesfully parsed. " << std::endl;
                QDPIO::cout << "Sink = " << params.sink_header.sink.id << std::endl;
                QDPIO::cout << "Smearing XML:" << std::endl;
                QDPIO::cout << params.sink_header.sink.xml << std::endl;

                Handle<QuarkSourceSink<LatticePropagator>> sinkSmearing(ThePropSinkSmearingFactory::Instance().createObject(params.sink_header.sink.id, sinktop, params.sink_header.sink.path, u));
                (*sinkSmearing)(seqsource);
                QDPIO::cout << "Sink successfully updated." << std::endl;

                for(int loop=0; loop < num_props; ++loop)
                {
                    forward_headers[loop].sink_header = params.sink_header;
                }
            }

            // HACKETY HACK HACK.
            // if we take in point propagators, eg. /Propagator, then we need to fake some XML stuff.
            if( smear_sink == false )
            {
                QDPIO::cout << "Propagators are UN-smeared. Faking XML stuff in the forward props. " << std::endl;
                params.sink_header.sink.xml = "<Sink><SinkType>POINT_SINK</SinkType><j_decay>"+std::to_string(j_decay)+"</j_decay></Sink>";
                params.sink_header.j_decay = j_decay;
                for(int loop=0; loop < num_props; ++loop)
                {
                    forward_headers[loop].sink_header = params.sink_header;
                }
            }

            // Output some stuff to the xml for verification.
            LalibeSftMom phases(0, true, j_decay);

            multi1d<Double> seqsource_corr = sumMulti(localNorm2(seqsource),
                          phases.getSet());

            push(xml_out,  "SeqSource_correlator");
            write(xml_out, "seqsource_corr", seqsource_corr);
            pop(xml_out);

            QDPIO::cout << "Attempt to store sequential source" << std::endl;
            try
            {
                XMLBufferWriter file_xml;
                push(file_xml,  "seqsource");
                write(file_xml, "id", uniqueId());  // NOTE: new ID form
                pop(file_xml);

                // Hackety hack hack the sequential source header.
                SeqSource_t seq_source_params;
                GroupXML_t seqsrc;
                std::string tmp;
                //NOTE:AWL I think the <SeqSource> should be removed
                tmp =      "\n  <SeqSource>\n    <SeqSourceType>"+params.ssparam.particle+"</SeqSourceType>\n";
                tmp +=     "    <flavor>"+params.ssparam.flavor+"</flavor>\n";

                if(params.ssparam.spin_zero_meson == false)
                {
                    tmp += "    <source_spin>"+params.ssparam.source_spin+"</source_spin>\n";
                    tmp += "    <sink_spin>"+params.ssparam.source_spin+"</sink_spin>\n";
                }
                else
                {
                    tmp += "    <source_spin>"+std::to_string(0)+"</source_spin>\n";
                    tmp += "    <sink_spin>"+std::to_string(0)+"</sink_spin>\n";
                }

                tmp +=     "    <sink_mom>"+std::to_string(params.ssparam.sink_mom[0])+" "+std::to_string(params.ssparam.sink_mom[1])+" "+std::to_string(params.ssparam.sink_mom[2])+"</sink_mom>\n";
                if(params.ssparam.t_all)
                {
                    tmp += "    <t_all>true</t_all>\n";
                }
                else
                {
                    tmp += "    <t_all>false</t_all>\n";
                }
                tmp +=     "    <t_0>"+std::to_string(params.ssparam.t_0)+"</t_0>\n";
                tmp +=     "    <t_sink>"+std::to_string(params.ssparam.t_sink)+"</t_sink>\n";
                tmp +=     "    <t_sep>"+std::to_string(params.ssparam.t_sep)+"</t_sep>\n";
                tmp +=     "    <j_decay>"+std::to_string(j_decay)+"</j_decay>\n";
                tmp +=     "  </SeqSource>\n";

                seqsrc.xml = tmp;
                seq_source_params.seqsrc = seqsrc;
                //        seq_source_params.sink_mom = params.ssparam.sink_mom;
                //        seq_source_params.t_sink = params.ssparam.t_sink;
                //        seq_source_params.j_decay = j_decay;

                // Sequential source header
                // Header composed of all forward prop headers
                SequentialSource_t new_header;
                new_header.sink_header      = params.sink_header;
                new_header.seqsource_header = seq_source_params;
                new_header.forward_props    = forward_headers;
                new_header.gauge_header     = gauge_xml.printCurrentContext();

                XMLBufferWriter record_xml;
                write(record_xml, "SequentialSource", new_header);

                // Store the seqsource
                TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.seqsource_id);
                TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.seqsource_id) = seqsource;
                TheNamedObjMap::Instance().get(params.named_obj.seqsource_id).setFileXML(file_xml);
                TheNamedObjMap::Instance().get(params.named_obj.seqsource_id).setRecordXML(record_xml);

                QDPIO::cout << "Sequential source successfully stored"  << std::endl;
            }
            catch (std::bad_cast)
            {
                QDPIO::cerr << name << ": dynamic cast error"<< std::endl;
                QDP_abort(1);
            }
            catch (const std::string& e)
            {
                QDPIO::cerr << name << ": error storing seqsource: " << e << std::endl;
                QDP_abort(1);
            }

            pop(xml_out);    // sequential source.

            snoop.stop();
            QDPIO::cout << LalibeSeqSourceEnv::name << ": total time = " << snoop.getTimeInSeconds() << " secs" << std::endl;
            QDPIO::cout << LalibeSeqSourceEnv::name<< ": ran successfully" << std::endl;
            END_CODE();
        }
    }// LalibeSeqSourceEnv
};
