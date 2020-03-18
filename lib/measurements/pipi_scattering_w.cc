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

 // Deleted:
 //if (paramtop.count("rotate_to_Dirac") != 0)
 //else if (paramtop.count("mom_list") != 0)
 //
 //            write(xml, "rotate_to_Dirac", par.rotate_to_Dirac);
 //When write,             if (par.is_mom_max == true)
 //            SftMom ft(params.param.p2_max, j_decay);
//LalibeSftMom ft = params.param.is_mom_max ? LalibeSftMom(params.param.p2_max, origin, false, j_decay)
//: LalibeSftMom(params.param.p_list, origin, j_decay);
 //if (full_correlator == true)
// delete full_correlator


#include "pipi_scattering_w.h"
#include "../contractions/pipi_scattering_func_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
//#include "../momentum/lalibe_sftmom.h"
#include "meas/inline/io/named_objmap.h"
#include "io/qprop_io.h"
//#include "util/spin_basis.h"


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
            QDPIO::cout << "Reading momenta centered around the origin with a max of " << par.p2_max << "and total max of " << par.ptot2_max << std::endl;
            read(paramtop, "particle_list", par.particle_list);
        }

        void write(XMLWriter& xml, const std::string& path, const PipiParams::Param_t& par) {}

        void read(XMLReader& xml, const std::string& path, PipiParams::NamedObject_t& input)
        {

            XMLReader inputtop(xml, path);

            //read(inputtop, "gauge_id" , input.gauge_id);
            //Logic to read quark propagators (whichever ones are present.)

            if (inputtop.count("quark_prop_1") != 0)
            {
                read(inputtop, "quark_prop_1", input.quark_prop_1);
                QDPIO::cout << "I found the first quark propagator, here is its id: " << input.quark_prop_1 << std::endl;
                input.is_prop_1 = true;
            }
            else
            {
                QDPIO::cout << "I couldn't find the first quark, hope you don't need it for the inputted pipi scattering. " << std::endl;
                input.is_prop_1 = false;
            }
            if (inputtop.count("quark_prop_2") != 0)
            {
                read(inputtop, "quark_prop_2", input.quark_prop_2);
                QDPIO::cout << "I found the second quark propagator, here is its id: " << input.quark_prop_2 << std::endl;
                input.is_prop_2 = true;
            }
            else
            {
                QDPIO::cout << "I couldn't find the second quark, hope you don't need it for the inputted pipi scattering. " << std::endl;
                input.is_prop_2 = false;
            }
            if (inputtop.count("quark_prop_3") != 0)
            {
                read(inputtop, "quark_prop_3", input.quark_prop_3);
                QDPIO::cout << "I found the third quark propagator, here is its id: " << input.quark_prop_3 << std::endl;
                input.is_prop_3 = true;
            }
            else
            {
                QDPIO::cout << "I couldn't find the third quark, hope you don't need it for the inputted pipi scattering. " << std::endl;
                input.is_prop_3 = false;
            }
            if (inputtop.count("quark_prop_4") != 0)
            {
                read(inputtop, "quark_prop_4", input.quark_prop_4);
                QDPIO::cout << "I found the forth quark propagator, here is its id: " << input.quark_prop_4 << std::endl;
                input.is_prop_4 = true;
            }
            else
            {
                QDPIO::cout << "I couldn't find the forth quark, hope you don't need it for the inputted pipi scattering. " << std::endl;
                input.is_prop_4 = false;
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
            XMLReader prop_1_file_xml, prop_1_record_xml;
            LatticePropagator quark_propagator_1;
            XMLReader prop_2_file_xml, prop_2_record_xml;
            LatticePropagator quark_propagator_2;
            XMLReader prop_3_file_xml, prop_3_record_xml;
            LatticePropagator quark_propagator_3;
            XMLReader prop_4_file_xml, prop_4_record_xml;
            LatticePropagator quark_propagator_4;

            //Need origin, j_decay, and t0 for fourier transform!
            //Need j_decay of bc to know what comes with a minus sign.
            int j_decay;
            int t_0;
            multi2d<int> origin(2,4);

            if (params.named_obj.is_prop_1 == true)
            {
                QDPIO::cout << "Attempting to read the first propagator" << std::endl;
                try
                {
                    quark_propagator_1 = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.quark_prop_1);
                    TheNamedObjMap::Instance().get(params.named_obj.quark_prop_1).getFileXML(prop_1_file_xml);
                    TheNamedObjMap::Instance().get(params.named_obj.quark_prop_1).getRecordXML(prop_1_record_xml);
                    //Get the origin  and j_decay for the FT, this assumes all quarks have the same origin.
                    MakeSourceProp_t  orig_header;
                    if (prop_1_record_xml.count("/Propagator") != 0)
                    {
                        QDPIO::cout << "The first quark propagator is unsmeared, reading from Propagator tag..." << std::endl;
                        read(prop_1_record_xml, "/Propagator", orig_header);
                    }
                    else if (prop_1_record_xml.count("/SinkSmear") != 0)
                    {
                        QDPIO::cout << "The first quark propagator is smeared, reading from SinkSmear tag..." << std::endl;
                        read(prop_1_record_xml, "/SinkSmear", orig_header);
                    }
                    else
                    {
                        QDPIO::cout << "What type of weird propagator did you give me? I can't find the right tag to get j_decay, source location, and whatnot...exiting..." << std::endl;
                    }

                    j_decay = orig_header.source_header.j_decay;
                    t_0 = orig_header.source_header.t_source;
                    origin[0] = orig_header.source_header.getTSrce();

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

            if (params.named_obj.is_prop_2 == true)
            {
                QDPIO::cout << "Attempting to read the second propagator" << std::endl;
                try
                {
                    quark_propagator_2 = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.quark_prop_2);
                    TheNamedObjMap::Instance().get(params.named_obj.quark_prop_2).getFileXML(prop_2_file_xml);
                    TheNamedObjMap::Instance().get(params.named_obj.quark_prop_2).getRecordXML(prop_2_record_xml);
                    //Get the origin  and j_decay for the FT, this assumes all quarks have the same origin.
                    MakeSourceProp_t  orig_header;
                    if (prop_2_record_xml.count("/Propagator") != 0)
                    {
                        QDPIO::cout << "The second quark propagator is unsmeared, reading from Propagator tag..." << std::endl;
                        read(prop_2_record_xml, "/Propagator", orig_header);
                    }
                    else if (prop_2_record_xml.count("/SinkSmear") != 0)
                    {
                        QDPIO::cout << "The second quark propagator is smeared, reading from SinkSmear tag..." << std::endl;
                        read(prop_2_record_xml, "/SinkSmear", orig_header);
                    }
                    else
                    {
                        QDPIO::cout << "What type of weird propagator did you give me? I can't find the right tag to get j_decay, source location, and whatnot...exiting..." << std::endl;
                    }

                    j_decay = orig_header.source_header.j_decay;
                    t_0 = orig_header.source_header.t_source;
                    if (orig_header.source_header.getTSrce() != origin[0])
                    {
                        QDPIO::cout << "The origin of the first two propagators are different, you idiot, I'm gonna speak Chinese. 你个傻逼." << std::endl;
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

            if (params.named_obj.is_prop_3 == true)
            {
                QDPIO::cout << "Attempting to read the third propagator" << std::endl;
                try
                {
                    quark_propagator_3 = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.quark_prop_3);
                    TheNamedObjMap::Instance().get(params.named_obj.quark_prop_3).getFileXML(prop_3_file_xml);
                    TheNamedObjMap::Instance().get(params.named_obj.quark_prop_3).getRecordXML(prop_3_record_xml);
                    //Get the origin  and j_decay for the FT, this assumes all quarks have the same origin.
                    MakeSourceProp_t  orig_header;
                    if (prop_3_record_xml.count("/Propagator") != 0)
                    {
                        QDPIO::cout << "The third quark propagator is unsmeared, reading from Propagator tag..." << std::endl;
                        read(prop_3_record_xml, "/Propagator", orig_header);
                    }
                    else if (prop_3_record_xml.count("/SinkSmear") != 0)
                    {
                        QDPIO::cout << "The third quark propagator is smeared, reading from SinkSmear tag..." << std::endl;
                        read(prop_3_record_xml, "/SinkSmear", orig_header);
                    }
                    else
                    {
                        QDPIO::cout << "What type of weird propagator did you give me? I can't find the right tag to get j_decay, source location, and whatnot...exiting..." << std::endl;
                    }

                    j_decay = orig_header.source_header.j_decay;
                    t_0 = orig_header.source_header.t_source;
                    origin[1] = orig_header.source_header.getTSrce();

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

            if (params.named_obj.is_prop_4 == true)
            {
                QDPIO::cout << "Attempting to read the fourth propagator" << std::endl;
                try
                {
                    quark_propagator_4 = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.quark_prop_4);
                    TheNamedObjMap::Instance().get(params.named_obj.quark_prop_4).getFileXML(prop_4_file_xml);
                    TheNamedObjMap::Instance().get(params.named_obj.quark_prop_4).getRecordXML(prop_4_record_xml);
                    //Get the origin  and j_decay for the FT, this assumes all quarks have the same origin.
                    MakeSourceProp_t  orig_header;
                    if (prop_4_record_xml.count("/Propagator") != 0)
                    {
                        QDPIO::cout << "The forth quark propagator is unsmeared, reading from Propagator tag..." << std::endl;
                        read(prop_4_record_xml, "/Propagator", orig_header);
                    }
                    else if (prop_4_record_xml.count("/SinkSmear") != 0)
                    {
                        QDPIO::cout << "The forth quark propagator is smeared, reading from SinkSmear tag..." << std::endl;
                        read(prop_4_record_xml, "/SinkSmear", orig_header);
                    }
                    else
                    {
                        QDPIO::cout << "What type of weird propagator did you give me? I can't find the right tag to get j_decay, source location, and whatnot...exiting..." << std::endl;
                    }

                    j_decay = orig_header.source_header.j_decay;
                    t_0 = orig_header.source_header.t_source;
                    if (orig_header.source_header.getTSrce() != origin[1])
                    {
                        QDPIO::cout << "The origin of the last two propagators are different, you idot, I'm gonna speak Chinese. 你个傻逼." << std::endl;
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

            //Here's Nt, we need this.
            int Nt = Layout::lattSize()[j_decay];

#ifdef BUILD_HDF5
            //If we are writing with hdf5, the start up is done here.
            HDF5Writer h5out(params.param.file_name);
            //h5out.push(params.param.obj_path);
            HDF5Base::writemode wmode;
            wmode = HDF5Base::ate;
#endif

            //Next we do the contractions for the specified particles.
            //This is going to be a horrible set of if statements for now, may change this later.
            //Loop over list of particles.

            for (int particle_index = 0; particle_index < params.param.particle_list.size(); particle_index++)
            {

                if (params.param.particle_list[particle_index] == "piplus")
                {
                    QDPIO::cout << "Particle number " << (particle_index + 1) << " is the piplus." << std::endl;
                    QDPIO::cout << "Checking to make sure we have the correct quark propagators to compute the piplus." << std::endl;
                    //if (params.named_obj.is_prop_1 == true && params.named_obj.is_prop_2 == true && params.named_obj.is_prop_3 == true && params.named_obj.is_prop_4 == true)
                    {
                        QDPIO::cout << "Found all four quarks for the two pions. Starting calculation..." << std::endl;


                        CorrelatorType::Correlator correlator_out;

                        pipi_correlator(correlator_out, quark_propagator_1, quark_propagator_2, quark_propagator_3, quark_propagator_4, origin, params.param.p2_max, params.param.ptot2_max, t_0, j_decay);

                        QDPIO::cout << "Calculation finished. Starting to write HDF5..." << std::endl;

                        // Write out the correlator.

                        //Temp variable for writing below.
                        DComplex temp_element;

                        //Move the h5 pushing here, since all momentum keys will be written in the same general path.
#ifdef BUILD_HDF5
                        //std::string correlator_path = params.param.obj_path + "/" + "pipi" + "/x" + std::to_string(origin[0][0]) + "_y" + std::to_string(origin[0][1]) + "_z" + std::to_string(origin[0][2]) + "_t" + std::to_string(origin[0][3]) + "_xprime" + std::to_string(origin[1][0]) + "_yprime" + std::to_string(origin[1][1]) + "_zprime" + std::to_string(origin[1][2]) + "_tprime" + std::to_string(origin[1][3]);
                        // x and y stands for the 4-vectors of the two pions, respectivily
                        std::string correlator_path = params.param.obj_path + "/" + "pipi" + "/x_" + std::to_string(origin[0][0]) + "_" +  std::to_string(origin[0][1]) + "_" + std::to_string(origin[0][2]) + "_" + std::to_string(origin[0][3]) + "__y_" + std::to_string(origin[1][0]) + "_" +  std::to_string(origin[1][1]) + "_" + std::to_string(origin[1][2]) + "_" + std::to_string(origin[1][3]);

                        h5out.push(correlator_path);
#else
                        std::string pipi_string_name="pipi";
                        //std::string correlator_path = pipi_string_name + "_x" + std::to_string(origin[0][0]) + "_y" + std::to_string(origin[0][1]) + "_z" + std::to_string(origin[0][2]) + "_t" + std::to_string(origin[0][3]) + "_xprime" + std::to_string(origin[1][0]) + "_yprime" + std::to_string(origin[1][1]) + "_zprime" + std::to_string(origin[1][2]) + "_tprime" + std::to_string(origin[1][3]);
                        std::string correlator_path = pipi_string_name + "_x_" + std::to_string(origin[0][0]) + "_" +  std::to_string(origin[0][1]) + "_" + std::to_string(origin[0][2]) + "_" + std::to_string(origin[0][3]) + "__y_" + std::to_string(origin[1][0]) + "_" +  std::to_string(origin[1][1]) + "_" + std::to_string(origin[1][2]) + "_" + std::to_string(origin[1][3]);
#endif
                        std::map<CorrelatorType::momenta_pair, multi1d<DComplex>>::iterator iter;
                        for (iter = correlator_out.begin(); iter != correlator_out.end(); iter++)
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
                                int t_relative = t - t_0;
                                if (t_relative < 0)
                                    t_relative += Nt;
#ifndef BUILD_HDF5
                                file_out << temp_element << "\n";
#endif
                                pipi_correlator_towrite[t_relative] = temp_element;
                            }
#ifndef BUILD_HDF5
                            file_out.close();
#else
                            //Change the name of string compred to 4d output so general correlator path is the same.
                            //std::string correlator_path_mom = correlator_path + "/px" + std::to_string(std::get<0>(momenta1)) + "_py" + std::to_string(std::get<1>(momenta1)) + "_pz" + std::to_string(std::get<2>(momenta1)) + "_qx" + std::to_string(std::get<0>(momenta2)) + "_qy" + std::to_string(std::get<1>(momenta2)) + "_qz" + std::to_string(std::get<2>(momenta2));
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
                    //else
                      //  QDPIO::cout << "Sorry, I couldn't find all four quark. Skipping the pipi scattering..." << std::endl;
                }

            }

#ifdef BUILD_HDF5
            h5out.cd("/");
            h5out.close();
#endif

            //pop(xml_out);

            snoop.stop();
            QDPIO::cout << name << ": total time = "
                << snoop.getTimeInSeconds()
                << " secs" << std::endl;

            QDPIO::cout << name << ": ran successfully" << std::endl;

            END_CODE();
        }

    }

}
