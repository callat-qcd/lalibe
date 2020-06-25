/*
 * pipi scattering
 * Authors:
 * Arjun Gambhir
 * Andre Walker-Loud
 * Jason Chang
 * David Brantley
 * Ben Horz
 * Haobo Yan
 * Do pipi scattering and write out the two-point correlator in hdf5
 */

#include "pipi_scattering_w.h"
#include "../contractions/pipi_scattering_func_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/named_objmap.h"
#include "io/qprop_io.h"

namespace Chroma
{
    namespace LalibePipiScatteringEnv
    {
        namespace
        {
            AbsInlineMeasurement* createMeasurement(XMLReader& xml_in,
                const std::string& path)
            {
                return new InlineMeas(PipiParams(xml_in, path));
            }

            bool registered = false;
        }
        const std::string name = "PIPI_SCATTERING";

        bool registerAll()
        {
            bool success = true;
            if (!registered)
            {
                success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
                registered = true;
            }
            return success;
        }

        void read(XMLReader& xml, const std::string& path, PipiParams::Param_t& par)
        {

            XMLReader paramtop(xml, path);

#ifdef BUILD_HDF5
            read(paramtop, "h5_file_name", par.file_name);
            read(paramtop, "obj_path", par.obj_path);
            QDPIO::cout << "HDF5 found, writing to" << par.file_name << " to path " << par.obj_path << std::endl;
#endif

            //We set output_full_correlator to true if no momentum is specified.
            read(paramtop, "p2_max", par.p2_max);
            read(paramtop, "ptot2_max", par.ptot2_max);
            read(paramtop, "diagrams", par.diagrams);
            QDPIO::cout << "Reading momenta centered around the origin with a max of " << par.p2_max << "and total max of " << par.ptot2_max << std::endl;
        }

        void write(XMLWriter& xml, const std::string& path, const PipiParams::Param_t& par) {}

        void read(XMLReader& xml, const std::string& path, PipiParams::NamedObject_t& input)
        {

            XMLReader inputtop(xml, path);

            //Logic to read quark propagators (whichever ones are present.)

            if (inputtop.count("light_prop_1") != 0)
            {
                read(inputtop, "light_prop_1", input.light_prop_1);
                QDPIO::cout << "I found the first light quark propagator, here is its id: " << input.light_prop_1 << std::endl;
                input.have_light_prop_1 = true;
            }
            else
            {
                QDPIO::cout << "I couldn't find the first light quark, this is absurd, go change your xml. " << std::endl;
                input.have_light_prop_1 = false;
            }
            if (inputtop.count("light_prop_2") != 0)
            {
                read(inputtop, "light_prop_2", input.light_prop_2);
                QDPIO::cout << "I found the second light quark propagator, here is its id: " << input.light_prop_2 << std::endl;
                input.have_light_prop_2 = true;
            }
            else
            {
                QDPIO::cout << "I couldn't find the second light quark, we have only one origin for the inputted pipi scattering. " << std::endl;
                input.have_light_prop_2 = false;
            }
            if (inputtop.count("strange_prop_1") != 0)
            {
                read(inputtop, "strange_prop_1", input.strange_prop_1);
                QDPIO::cout << "I found the first strange quark propagator, here is its id: " << input.strange_prop_1 << std::endl;
                input.have_strange_prop_1 = true;
            }
            else
            {
                QDPIO::cout << "I couldn't find the first strange quark, no KK or pik scattering. " << std::endl;
                input.have_strange_prop_1 = false;
            }
            if (inputtop.count("strange_prop_2") != 0)
            {
                read(inputtop, "strange_prop_2", input.strange_prop_2);
                QDPIO::cout << "I found the second strange quark propagator, here is its id: " << input.strange_prop_2 << std::endl;
                input.have_strange_prop_2 = true;
            }
            else
            {
                input.have_strange_prop_2 = false;
                if (inputtop.count("strange_prop_1") != 0)
                    QDPIO::cout << "I couldn't find the second strange quark, we have only one origin for the inputted KK scattering. " << std::endl;
            }

        }

        void write(XMLWriter& xml, const std::string& path, const PipiParams::NamedObject_t& input) {}


        PipiParams::PipiParams()
        {
            frequency = 0;
        }

        PipiParams::PipiParams(XMLReader& xml_in, const std::string& path)
        {
            try
            {
                XMLReader paramtop(xml_in, path);

                if (paramtop.count("Frequency") == 1)
                    read(paramtop, "Frequency", frequency);
                else
                    frequency = 1;

                read(paramtop, "PipiParams", param);

                read(paramtop, "NamedObject", named_obj);

            }
            catch (const std::string& e)
            {
                QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
                QDP_abort(1);
            }
        }

        void
            PipiParams::writeXML(XMLWriter& xml_out, const std::string& path) {}


        void  InlineMeas::operator()(unsigned long update_no, XMLWriter& xml_out)
        {
            START_CODE();

            StopWatch snoop;
            snoop.reset();
            snoop.start();

            QDPIO::cout << "Pipi scattering starting..." << std::endl;

            //Grab all the propagators that are given.
            XMLReader light_prop_1_file_xml, light_prop_1_record_xml;
            LatticePropagator light_quark_propagator_1;
            XMLReader light_prop_2_file_xml, light_prop_2_record_xml;
            LatticePropagator light_quark_propagator_2;
            XMLReader strange_prop_1_file_xml, strange_prop_1_record_xml;
            LatticePropagator strange_quark_propagator_1;
            XMLReader strange_prop_2_file_xml, strange_prop_2_record_xml;
            LatticePropagator strange_quark_propagator_2;

            //Need origin, j_decay, and t0 for fourier transform!
            //Need j_decay of bc to know what comes with a minus sign.
            int j_decay;
            int t_0;
            multi2d<int> origin11(2,4), origin22(2,4), origin12(2,4);

            if (params.named_obj.have_light_prop_1 == true)
            {
                QDPIO::cout << "Attempting to read the first light propagator" << std::endl;
                try
                {
                    light_quark_propagator_1 = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.light_prop_1);
                    TheNamedObjMap::Instance().get(params.named_obj.light_prop_1).getFileXML(light_prop_1_file_xml);
                    TheNamedObjMap::Instance().get(params.named_obj.light_prop_1).getRecordXML(light_prop_1_record_xml);
                    //Get the origin  and j_decay for the FT, this assumes all quarks have the same origin.
                    MakeSourceProp_t  orig_header;
                    if (light_prop_1_record_xml.count("/Propagator") != 0)
                    {
                        QDPIO::cout << "The first light quark propagator is unsmeared, reading from Propagator tag..." << std::endl;
                        read(light_prop_1_record_xml, "/Propagator", orig_header);
                    }
                    else if (light_prop_1_record_xml.count("/SinkSmear") != 0)
                    {
                        QDPIO::cout << "The first light quark propagator is smeared, reading from SinkSmear tag..." << std::endl;
                        read(light_prop_1_record_xml, "/SinkSmear", orig_header);
                    }
                    else
                    {
                        QDPIO::cout << "What type of weird propagator did you give me? I can't find the right tag to get j_decay, source location, and whatnot...exiting..." << std::endl;
                    }

                    j_decay = orig_header.source_header.j_decay;
                    t_0 = orig_header.source_header.t_source;
                    origin11[0] = orig_header.source_header.getTSrce();
                    origin11[1] = orig_header.source_header.getTSrce();
                    origin12[0] = orig_header.source_header.getTSrce();
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

            if (params.named_obj.have_light_prop_2 == true)
            {
                QDPIO::cout << "Attempting to read the second light propagator" << std::endl;
                try
                {
                    light_quark_propagator_2 = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.light_prop_2);
                    TheNamedObjMap::Instance().get(params.named_obj.light_prop_2).getFileXML(light_prop_2_file_xml);
                    TheNamedObjMap::Instance().get(params.named_obj.light_prop_2).getRecordXML(light_prop_2_record_xml);
                    //Get the origin  and j_decay for the FT, this assumes all quarks have the same origin.
                    MakeSourceProp_t  orig_header;
                    if (light_prop_2_record_xml.count("/Propagator") != 0)
                    {
                        QDPIO::cout << "The second light quark propagator is unsmeared, reading from Propagator tag..." << std::endl;
                        read(light_prop_2_record_xml, "/Propagator", orig_header);
                    }
                    else if (light_prop_2_record_xml.count("/SinkSmear") != 0)
                    {
                        QDPIO::cout << "The second light quark propagator is smeared, reading from SinkSmear tag..." << std::endl;
                        read(light_prop_2_record_xml, "/SinkSmear", orig_header);
                    }
                    else
                    {
                        QDPIO::cout << "What type of weird propagator did you give me? I can't find the right tag to get j_decay, source location, and whatnot...exiting..." << std::endl;
                    }

                    j_decay = orig_header.source_header.j_decay;
                    t_0 = orig_header.source_header.t_source;
                    origin22[0] = orig_header.source_header.getTSrce();
                    origin22[1] = orig_header.source_header.getTSrce();
                    origin12[1] = orig_header.source_header.getTSrce();
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

            if (params.named_obj.have_strange_prop_1 == true)
            {
                QDPIO::cout << "Attempting to read the first strange propagator" << std::endl;
                try
                {
                    strange_quark_propagator_1 = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.strange_prop_1);
                    TheNamedObjMap::Instance().get(params.named_obj.strange_prop_1).getFileXML(strange_prop_1_file_xml);
                    TheNamedObjMap::Instance().get(params.named_obj.strange_prop_1).getRecordXML(strange_prop_1_record_xml);
                    //Get the origin  and j_decay for the FT, this assumes all quarks have the same origin.
                    MakeSourceProp_t  orig_header;
                    if (strange_prop_1_record_xml.count("/Propagator") != 0)
                    {
                        QDPIO::cout << "The first strange quark propagator is unsmeared, reading from Propagator tag..." << std::endl;
                        read(strange_prop_1_record_xml, "/Propagator", orig_header);
                    }
                    else if (strange_prop_1_record_xml.count("/SinkSmear") != 0)
                    {
                        QDPIO::cout << "The first strange quark propagator is smeared, reading from SinkSmear tag..." << std::endl;
                        read(strange_prop_1_record_xml, "/SinkSmear", orig_header);
                    }
                    else
                    {
                        QDPIO::cout << "What type of weird propagator did you give me? I can't find the right tag to get j_decay, source location, and whatnot...exiting..." << std::endl;
                    }

                    j_decay = orig_header.source_header.j_decay;
                    t_0 = orig_header.source_header.t_source;
                    if (orig_header.source_header.getTSrce() != origin11[0])
                    {
                        QDPIO::cout << "The origin of the first light and strange propagators are different, go change your propagators." << std::endl;
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

            if (params.named_obj.have_strange_prop_2 == true)
            {
                QDPIO::cout << "Attempting to read the second strange propagator" << std::endl;
                try
                {
                    strange_quark_propagator_2 = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.strange_prop_2);
                    TheNamedObjMap::Instance().get(params.named_obj.strange_prop_2).getFileXML(strange_prop_2_file_xml);
                    TheNamedObjMap::Instance().get(params.named_obj.strange_prop_2).getRecordXML(strange_prop_2_record_xml);
                    //Get the origin  and j_decay for the FT, this assumes all quarks have the same origin.
                    MakeSourceProp_t  orig_header;
                    if (strange_prop_2_record_xml.count("/Propagator") != 0)
                    {
                        QDPIO::cout << "The second strange quark propagator is unsmeared, reading from Propagator tag..." << std::endl;
                        read(strange_prop_2_record_xml, "/Propagator", orig_header);
                    }
                    else if (strange_prop_2_record_xml.count("/SinkSmear") != 0)
                    {
                        QDPIO::cout << "The second strange quark propagator is smeared, reading from SinkSmear tag..." << std::endl;
                        read(strange_prop_2_record_xml, "/SinkSmear", orig_header);
                    }
                    else
                    {
                        QDPIO::cout << "What type of weird propagator did you give me? I can't find the right tag to get j_decay, source location, and whatnot...exiting..." << std::endl;
                    }

                    j_decay = orig_header.source_header.j_decay;
                    t_0 = orig_header.source_header.t_source;
                    if (orig_header.source_header.getTSrce() != origin22[0])
                    {
                        QDPIO::cout << "The origin of the second light and strange are different, you idiot, go change your propagators." << std::endl;
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

            int Nt = Layout::lattSize()[j_decay];

#ifdef BUILD_HDF5
            //If we are writing with hdf5, the start up is done here.
            HDF5Writer h5out(params.param.file_name);
            //h5out.push(params.param.obj_path);
            HDF5Base::writemode wmode;
            wmode = HDF5Base::ate;
#endif

            //Next we call the physics code to do the pipi scattering.

            CorrelatorType::Correlator correlator_pipi11, correlator_pipi12, correlator_pipi22, correlator_kk11, correlator_kk12, correlator_kk22, correlator_pik11, correlator_pik12, correlator_pik22;

            std::map<std::string, CorrelatorType::Correlator> correlators;

            if (params.named_obj.have_light_prop_1 == true && params.named_obj.have_light_prop_2 == false && params.named_obj.have_strange_prop_1 == false && params.named_obj.have_strange_prop_2 == false)
            {
                pipi_correlator(correlators["pipi11"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay);
            }
            else if (params.named_obj.have_light_prop_1 == true && params.named_obj.have_light_prop_2 == true && params.named_obj.have_strange_prop_1 == false && params.named_obj.have_strange_prop_2 == false)
            {
                if (params.param.diagrams == 1)
                {
                    pipi_correlator_debug(correlators["pipi11_diagram0"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 0);
                    pipi_correlator_debug(correlators["pipi22_diagram0"], light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, origin22, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 0);
                    pipi_correlator_debug(correlators["pipi12_diagram0"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_2, light_quark_propagator_2, origin12, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 0);
                    pipi_correlator_debug(correlators["pipi11_diagram1"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 1);
                    pipi_correlator_debug(correlators["pipi22_diagram1"], light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, origin22, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 1);
                    pipi_correlator_debug(correlators["pipi12_diagram1"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_2, light_quark_propagator_2, origin12, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 1);
                    pipi_correlator_debug(correlators["pipi11_diagram2"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 2);
                    pipi_correlator_debug(correlators["pipi22_diagram2"], light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, origin22, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 2);
                    pipi_correlator_debug(correlators["pipi12_diagram2"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_2, light_quark_propagator_2, origin12, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 2);
                    pipi_correlator_debug(correlators["pipi11_diagram3"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 3);
                    pipi_correlator_debug(correlators["pipi22_diagram3"], light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, origin22, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 3);
                    pipi_correlator_debug(correlators["pipi12_diagram3"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_2, light_quark_propagator_2, origin12, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 3);
                    pipi_correlator_debug(correlators["pipi11_diagram4"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 4);
                    pipi_correlator_debug(correlators["pipi22_diagram4"], light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, origin22, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 4);
                    pipi_correlator_debug(correlators["pipi12_diagram4"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_2, light_quark_propagator_2, origin12, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 4);
                }
                else
                {
                    pipi_correlator(correlators["pipi11"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay);
                    pipi_correlator(correlators["pipi22"], light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, origin22, params.param.p2_max, params.param.ptot2_max, t_0, j_decay);
                    pipi_correlator(correlators["pipi12"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_2, light_quark_propagator_2, origin12, params.param.p2_max, params.param.ptot2_max, t_0, j_decay);
                }
            }
            else if (params.named_obj.have_light_prop_1 == true && params.named_obj.have_light_prop_2 == false && params.named_obj.have_strange_prop_1 == true && params.named_obj.have_strange_prop_2 == false)
            {
                pipi_correlator(correlators["pipi11"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay);

                pipi_correlator(correlators["kk11"], light_quark_propagator_1, strange_quark_propagator_1, light_quark_propagator_1, strange_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay);

                pik_correlator(correlators["pik11"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, strange_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay);
            }
            else if (params.named_obj.have_light_prop_1 == true && params.named_obj.have_light_prop_2 == true && params.named_obj.have_strange_prop_1 == true && params.named_obj.have_strange_prop_2 == false)
            {
                pipi_correlator(correlators["pipi11"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay);
                pipi_correlator(correlators["pipi22"], light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, origin22, params.param.p2_max, params.param.ptot2_max, t_0, j_decay);
                pipi_correlator(correlators["pipi12"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_2, light_quark_propagator_2, origin12, params.param.p2_max, params.param.ptot2_max, t_0, j_decay);

                pipi_correlator(correlators["kk11"], light_quark_propagator_1, strange_quark_propagator_1, light_quark_propagator_1, strange_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay);

                pik_correlator(correlators["pik11"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, strange_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay);
                pik_correlator(correlators["pik12"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_2, strange_quark_propagator_2, origin12, params.param.p2_max, params.param.ptot2_max, t_0, j_decay);
            }
            else if (params.named_obj.have_light_prop_1 == true && params.named_obj.have_light_prop_2 == true && params.named_obj.have_strange_prop_1 == true && params.named_obj.have_strange_prop_2 == true)
            {
                if (params.param.diagrams == 1)
                {
                    pipi_correlator_debug(correlators["pipi11_diagram0"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 0);
                    pipi_correlator_debug(correlators["pipi22_diagram0"], light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, origin22, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 0);
                    pipi_correlator_debug(correlators["pipi12_diagram0"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_2, light_quark_propagator_2, origin12, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 0);
                    pipi_correlator_debug(correlators["pipi11_diagram1"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 1);
                    pipi_correlator_debug(correlators["pipi22_diagram1"], light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, origin22, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 1);
                    pipi_correlator_debug(correlators["pipi12_diagram1"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_2, light_quark_propagator_2, origin12, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 1);
                    pipi_correlator_debug(correlators["pipi11_diagram2"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 2);
                    pipi_correlator_debug(correlators["pipi22_diagram2"], light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, origin22, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 2);
                    pipi_correlator_debug(correlators["pipi12_diagram2"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_2, light_quark_propagator_2, origin12, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 2);
                    pipi_correlator_debug(correlators["pipi11_diagram3"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 3);
                    pipi_correlator_debug(correlators["pipi22_diagram3"], light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, origin22, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 3);
                    pipi_correlator_debug(correlators["pipi12_diagram3"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_2, light_quark_propagator_2, origin12, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 3);
                    pipi_correlator_debug(correlators["pipi11_diagram4"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 4);
                    pipi_correlator_debug(correlators["pipi22_diagram4"], light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, origin22, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 4);
                    pipi_correlator_debug(correlators["pipi12_diagram4"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_2, light_quark_propagator_2, origin12, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 4);

                    pipi_correlator_debug(correlators["kk11_diagram0"], light_quark_propagator_1, strange_quark_propagator_1, light_quark_propagator_1, strange_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 0);
                    pipi_correlator_debug(correlators["kk22_diagram0"], light_quark_propagator_2, strange_quark_propagator_2, light_quark_propagator_2, strange_quark_propagator_2, origin22, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 0);
                    pipi_correlator_debug(correlators["kk12_diagram0"], light_quark_propagator_1, strange_quark_propagator_1, light_quark_propagator_2, strange_quark_propagator_2, origin12, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 0);
                    pipi_correlator_debug(correlators["kk11_diagram1"], light_quark_propagator_1, strange_quark_propagator_1, light_quark_propagator_1, strange_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 1);
                    pipi_correlator_debug(correlators["kk22_diagram1"], light_quark_propagator_2, strange_quark_propagator_2, light_quark_propagator_2, strange_quark_propagator_2, origin22, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 1);
                    pipi_correlator_debug(correlators["kk12_diagram1"], light_quark_propagator_1, strange_quark_propagator_1, light_quark_propagator_2, strange_quark_propagator_2, origin12, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 1);
                    pipi_correlator_debug(correlators["kk11_diagram2"], light_quark_propagator_1, strange_quark_propagator_1, light_quark_propagator_1, strange_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 2);
                    pipi_correlator_debug(correlators["kk22_diagram2"], light_quark_propagator_2, strange_quark_propagator_2, light_quark_propagator_2, strange_quark_propagator_2, origin22, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 2);
                    pipi_correlator_debug(correlators["kk12_diagram2"], light_quark_propagator_1, strange_quark_propagator_1, light_quark_propagator_2, strange_quark_propagator_2, origin12, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 2);
                    pipi_correlator_debug(correlators["kk11_diagram3"], light_quark_propagator_1, strange_quark_propagator_1, light_quark_propagator_1, strange_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 3);
                    pipi_correlator_debug(correlators["kk22_diagram3"], light_quark_propagator_2, strange_quark_propagator_2, light_quark_propagator_2, strange_quark_propagator_2, origin22, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 3);
                    pipi_correlator_debug(correlators["kk12_diagram3"], light_quark_propagator_1, strange_quark_propagator_1, light_quark_propagator_2, strange_quark_propagator_2, origin12, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 3);
                    pipi_correlator_debug(correlators["kk11_diagram4"], light_quark_propagator_1, strange_quark_propagator_1, light_quark_propagator_1, strange_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 4);
                    pipi_correlator_debug(correlators["kk22_diagram4"], light_quark_propagator_2, strange_quark_propagator_2, light_quark_propagator_2, strange_quark_propagator_2, origin22, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 4);
                    pipi_correlator_debug(correlators["kk12_diagram4"], light_quark_propagator_1, strange_quark_propagator_1, light_quark_propagator_2, strange_quark_propagator_2, origin12, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 4);

                    pik_correlator_debug(correlators["pik11_diagram0"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, strange_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 0);
                    pik_correlator_debug(correlators["pik12_diagram0"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_2, strange_quark_propagator_2, origin12, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 0);
                    pik_correlator_debug(correlators["pik22_diagram0"], light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, strange_quark_propagator_2, origin22, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 0);
                    pik_correlator_debug(correlators["pik11_diagram1"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, strange_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 1);
                    pik_correlator_debug(correlators["pik12_diagram1"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_2, strange_quark_propagator_2, origin12, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 1);
                    pik_correlator_debug(correlators["pik22_diagram1"], light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, strange_quark_propagator_2, origin22, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 1);
                    pik_correlator_debug(correlators["pik11_diagram2"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, strange_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 2);
                    pik_correlator_debug(correlators["pik12_diagram2"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_2, strange_quark_propagator_2, origin12, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 2);
                    pik_correlator_debug(correlators["pik22_diagram2"], light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, strange_quark_propagator_2, origin22, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 2);
                }
                else
                {
                    pipi_correlator(correlators["pipi11"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay);
                    pipi_correlator(correlators["pipi22"], light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, origin22, params.param.p2_max, params.param.ptot2_max, t_0, j_decay);
                    pipi_correlator(correlators["pipi12"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_2, light_quark_propagator_2, origin12, params.param.p2_max, params.param.ptot2_max, t_0, j_decay);

                    pipi_correlator(correlators["kk11"], light_quark_propagator_1, strange_quark_propagator_1, light_quark_propagator_1, strange_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay);
                    pipi_correlator(correlators["kk22"], light_quark_propagator_2, strange_quark_propagator_2, light_quark_propagator_2, strange_quark_propagator_2, origin22, params.param.p2_max, params.param.ptot2_max, t_0, j_decay);
                    pipi_correlator(correlators["kk12"], light_quark_propagator_1, strange_quark_propagator_1, light_quark_propagator_2, strange_quark_propagator_2, origin12, params.param.p2_max, params.param.ptot2_max, t_0, j_decay);

                    pik_correlator(correlators["pik11"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_1, strange_quark_propagator_1, origin11, params.param.p2_max, params.param.ptot2_max, t_0, j_decay);
                    pik_correlator(correlators["pik12"], light_quark_propagator_1, light_quark_propagator_1, light_quark_propagator_2, strange_quark_propagator_2, origin12, params.param.p2_max, params.param.ptot2_max, t_0, j_decay);
                    pik_correlator(correlators["pik22"], light_quark_propagator_2, light_quark_propagator_2, light_quark_propagator_2, strange_quark_propagator_2, origin22, params.param.p2_max, params.param.ptot2_max, t_0, j_decay);
                }
            }

            QDPIO::cout << "Calculation finished. Starting to write HDF5..." << std::endl;

            // Write out the correlator.

            // Temp variable for writing below.
            DComplex temp_element;
            std::string particles, collision_type;

            std::map<std::string, CorrelatorType::Correlator>::iterator correlator_iter;

            // Loop over all collisions, 9 for max.
            for (correlator_iter = correlators.begin(); correlator_iter != correlators.end(); correlator_iter++)
            {
                particles = correlator_iter->first;
                multi2d<int> origin(2,4);
                int number_start;
                if (particles[1] == 'k')
                {
                  number_start = 2;
                  collision_type = "kk";
                }
                else if (particles[2] == 'p')
                {
                  number_start = 4;
                  collision_type = "pipi";
                }
                else if (particles[2] == 'k')
                {
                  number_start = 3;
                  collision_type = "pik";
                }

                if (particles.substr(number_start, 2) == "11")
                  origin = origin11;
                else if (particles.substr(number_start, 2) == "22")
                  origin = origin22;
                else if (particles.substr(number_start, 2) == "12")
                  origin = origin12;

                if (params.param.diagrams == 1)
                {
                  std::string diagram_dir=particles.substr(number_start+3, 8);
                  collision_type += "/" + diagram_dir;
                }

#ifdef BUILD_HDF5
                // x and y stands for the 4-vectors of the two pions, respectivily
                std::string correlator_path = params.param.obj_path + "/" + collision_type + "/x_" + std::to_string(origin[0][0]) + "_" +  std::to_string(origin[0][1]) + "_" + std::to_string(origin[0][2]) + "_" + std::to_string(origin[0][3]) + "__y_" + std::to_string(origin[1][0]) + "_" +  std::to_string(origin[1][1]) + "_" + std::to_string(origin[1][2]) + "_" + std::to_string(origin[1][3]);
                h5out.push(correlator_path);
#else
                std::string correlator_path = collision_type + "_x_" + std::to_string(origin[0][0]) + "_" +  std::to_string(origin[0][1]) + "_" + std::to_string(origin[0][2]) + "_" + std::to_string(origin[0][3]) + "__y_" + std::to_string(origin[1][0]) + "_" +  std::to_string(origin[1][1]) + "_" + std::to_string(origin[1][2]) + "_" + std::to_string(origin[1][3]);
#endif

                std::map<CorrelatorType::momenta_pair, multi1d<DComplex>>::iterator iter;
                for (iter = correlator_iter->second.begin(); iter != correlator_iter->second.end(); iter++)
                {
                    //One more temp variable instanited inside loop (once again for writing.)
                    multi1d<DComplex> pipi_correlator_towrite;
                    pipi_correlator_towrite.resize(Nt);
                    std::tuple<int, int, int> momenta1, momenta2;
                    momenta1 = std::get<0>(iter->first);
                    momenta2 = std::get<1>(iter->first);

  #ifndef BUILD_HDF5
                    std::string correlator_path_mom = correlator_path + "_px" + std::to_string(std::get<0>(momenta1)) + "_py" + std::to_string(std::get<1>(momenta1)) + "_pz" + std::to_string(std::get<2>(momenta1)) + "_qx" + std::to_string(std::get<0>(momenta2)) + "_qy" + std::to_string(std::get<1>(momenta2)) + "_qz" + std::to_string(std::get<2>(momenta2));
                    TextFileWriter file_out(correlator_path_mom);
  #endif
                    for (int t = 0; t < Nt; t++)
                    {
                        temp_element = iter->second[t];
  #ifndef BUILD_HDF5
                        file_out << temp_element << "\n";
  #endif
                        pipi_correlator_towrite[t] = temp_element;
                    }
  #ifndef BUILD_HDF5
                    file_out.close();
  #else
                    //Change the name of string compred to 4d output so general correlator path is the same.
                    //I add a total momentum directory to make André happy
                    std::string correlator_path_tot_mom = "/ptotx" + std::to_string(std::get<0>(momenta1)+std::get<0>(momenta2)) + "_ptoty" + std::to_string(std::get<1>(momenta1)+std::get<1>(momenta2)) + "_ptotz" + std::to_string(std::get<2>(momenta1)+std::get<2>(momenta2));
                    std::string correlator_path_mom = correlator_path + correlator_path_tot_mom + "/px" + std::to_string(std::get<0>(momenta1)) + "_py" + std::to_string(std::get<1>(momenta1)) + "_pz" + std::to_string(std::get<2>(momenta1)) + "_qx" + std::to_string(std::get<0>(momenta2)) + "_qy" + std::to_string(std::get<1>(momenta2)) + "_qz" + std::to_string(std::get<2>(momenta2));
                    //Haobo: I don't know why should I push, but push is right!
                    h5out.push(correlator_path+correlator_path_tot_mom);

                    h5out.write(correlator_path_mom, pipi_correlator_towrite, wmode);
                    h5out.writeAttribute(correlator_path_mom, "is_shifted", 1, wmode);
                    h5out.cd("/");
  #endif
                }
            }

#ifdef BUILD_HDF5
            h5out.cd("/");
            h5out.close();
#endif
            snoop.stop();
            QDPIO::cout << name << ": total time = "
                << snoop.getTimeInSeconds()
                << " secs" << std::endl;

            QDPIO::cout << name << ": ran successfully" << std::endl;

            END_CODE();
        }

    }

}
