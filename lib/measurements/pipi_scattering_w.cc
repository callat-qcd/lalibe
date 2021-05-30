/*
 * pipi scattering
 * Authors:
 * Haobo Yan
 * Ben Horz
 * Andre Walker-Loud
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
            //  xmlQuarkMapper map internal character identifiers to the expected
            //  XML tag
            std::vector<char> quark_flavs { 'u', 'd', 's' };
            std::map<char, std::string> xmlQuarkMapper = {
                {'u', "up_quark"},
                {'d', "down_quark"},
                {'s', "strange_quark"},
            };

            for (auto aFlav : quark_flavs)
            {
                if (inputtop.count(xmlQuarkMapper[aFlav]) != 0)
                {
                    std::string tmpName;
                    read(inputtop, xmlQuarkMapper[aFlav], tmpName);
                    input.quark_map[aFlav] = tmpName;
                    QDPIO::cout<<"I found a "<<aFlav<<" quark, here is its id: " <<input.quark_map[aFlav]<<std::endl;
                    }
                else
                {
                    QDPIO::cout<<"I couldn't find an "<<aFlav<<" quark, hope you don't need it for the inputted meson-meson contractions. "<<std::endl;
                    }
            }
        }
        // don't proliferate xml out that no one reads
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

            //Need origin, j_decay, and t0 for fourier transform!
            //Need j_decay of bc to know what comes with a minus sign.
            int j_decay;
            int t_0;
            multi1d<int> origin(4);
            // if we read more than one propagator, check that origin, t_0 and
            // j_decay match between them
            bool checkFlag = false;

            // check if we have down and strange props
            bool have_down    = params.named_obj.quark_map.count('d');
            bool have_strange = params.named_obj.quark_map.count('s');

            // map for propagator memory
            std::map<char, const LatticePropagator*> prop_map;
            for (auto& qIt : params.named_obj.quark_map) {
                try
                    {
                    prop_map[qIt.first] = &TheNamedObjMap::Instance().getData<LatticePropagator>(qIt.second);

                    XMLReader prop_file_xml, prop_record_xml;
                    TheNamedObjMap::Instance().get(qIt.second).getFileXML(prop_file_xml);
                    TheNamedObjMap::Instance().get(qIt.second).getRecordXML(prop_record_xml);
                    //Get the origin  and j_decay for the FT, this assumes all quarks have the same origin.
                    MakeSourceProp_t  orig_header;
                    if (prop_record_xml.count("/Propagator") != 0)
                        {
                            QDPIO::cout<< qIt.first <<" quark propagator is unsmeared, reading from Propagator tag..."<<std::endl;
                            read(prop_record_xml, "/Propagator", orig_header);
                        }
                    else if (prop_record_xml.count("/SinkSmear") != 0)
                        {
                            QDPIO::cout<< qIt.first <<" quark propagator is smeared, reading from SinkSmear tag..."<<std::endl;
                            read(prop_record_xml, "/SinkSmear", orig_header);
                        }
                    else
                        {
                            QDPIO::cout<<"What type of weird propagator did you give me? I can't find the right tag to get j_decay, source location, and whatnot...exiting..."<<std::endl;
                        }
                    // set info
                    if (!checkFlag) {
                        j_decay = orig_header.source_header.j_decay;
                        t_0     = orig_header.source_header.t_source;
                        origin  = orig_header.source_header.getTSrce();
                    }
                    // check info
                    else
                        {
                        if (j_decay != orig_header.source_header.j_decay) {
                            QDPIO::cerr << "j_decay doesn't match between propagators" << std::endl;
                            QDP_abort(1);
                            }
                        if (t_0 != orig_header.source_header.t_source) {
                            QDPIO::cerr << "t_source doesn't match between propagators" << std::endl;
                            QDP_abort(1);
                            }
                        for (unsigned int iC = 0; iC < origin.size(); ++iC) {
                            if (origin[iC] != orig_header.source_header.getTSrce()[iC]) {
                                QDPIO::cerr << "origin doesn't match between propagators" << std::endl;
                                QDP_abort(1);
                                }
                            }
                        } // check of origin, t0, j_decay between props
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
                checkFlag = true;
            }// end prop memory assign loop

            // if we do not have a down propagator, assign it's memory to the up propagator
            if (! have_down ){
                prop_map['d'] = prop_map['u'];
                }

            int Nt = Layout::lattSize()[j_decay];

#ifdef BUILD_HDF5
            //If we are writing with hdf5, the start up is done here.
            HDF5Writer h5out(params.param.file_name);
            //h5out.push(params.param.obj_path);
            HDF5Base::writemode wmode;
            wmode = HDF5Base::ate;
#endif
            // Next we call the physics code to do the pipi scattering
            // and collect the results in a map
            std::map<std::string, CorrelatorType::Correlator> correlators;

            // we always do pi+ pi+
            QDPIO::cout << "    pi+ pi+" << std::endl;
            pipi_correlator(correlators["pip_pip"],
                            *prop_map['u'], *prop_map['d'],
                            origin, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 0);
            if (params.param.diagrams == 1) {
                // do all 4 diagrams
                for (int dd=1; dd<=4; dd++){
                    pipi_correlator(correlators["pip_pip_d" + std::to_string(dd)],
                                    *prop_map['u'], *prop_map['d'],
                                    origin, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, dd);
                                }
                }
            // do we have a strange quark?
            if ( have_strange){
                // add pi+ K+
                QDPIO::cout << "    pi+ K+" << std::endl;
                pik_correlator(correlators["pip_kp"],
                                *prop_map['u'], *prop_map['d'], *prop_map['u'], *prop_map['s'],
                                origin, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 0);
                if (params.param.diagrams == 1) {
                    for (int dd=1; dd<=2; dd++){
                        pik_correlator(correlators["pip_kp_d" + std::to_string(dd)],
                                        *prop_map['u'], *prop_map['d'], *prop_map['u'], *prop_map['s'],
                                        origin, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, dd);
                                    }
                }

                // add K+ K+
                QDPIO::cout << "    K+ K+" << std::endl;
                pipi_correlator(correlators["kp_kp"],
                                *prop_map['u'], *prop_map['s'],
                                origin, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 0);
                if (params.param.diagrams == 1) {
                    for (int dd=1; dd<=4; dd++){
                        pipi_correlator(correlators["kp_kp_d" + std::to_string(dd)],
                                        *prop_map['u'], *prop_map['s'],
                                        origin, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, dd);
                                    }
                }

                // do we also have a down quark?  (breaking isospin)
                if (have_down){
                    // add K0 K0
                    QDPIO::cout << "    K0 K0" << std::endl;
                    pipi_correlator(correlators["k0_k0"],
                                    *prop_map['d'], *prop_map['s'],
                                    origin, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 0);
                    if (params.param.diagrams == 1) {
                        for (int dd=1; dd<=4; dd++){
                            pipi_correlator(correlators["k0_k0_d" + std::to_string(dd)],
                                            *prop_map['d'], *prop_map['s'],
                                            origin, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, dd);
                        }
                    }
                    // add pi+ K0Bar
                    QDPIO::cout << "    pi+ K0Bar" << std::endl;
                    pik_correlator(correlators["pip_k0b"],
                                    *prop_map['u'], *prop_map['d'], *prop_map['s'], *prop_map['d'],
                                    origin, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, 0);
                    if (params.param.diagrams == 1) {
                        for (int dd=1; dd<=2; dd++){
                            pik_correlator(correlators["pip_k0b_d"+ std::to_string(dd)],
                                            *prop_map['u'], *prop_map['d'], *prop_map['s'], *prop_map['d'],
                                            origin, params.param.p2_max, params.param.ptot2_max, t_0, j_decay, dd);
                        }
                    }
                }// end have_down
            }// end have_strange

            QDPIO::cout << "Calculation finished. Starting to write HDF5..." << std::endl;

            if (params.param.diagrams == 1)
                QDPIO::cout << "I will print separate diagrams" << std::endl;

            // Write out the correlator.

            // Temp variable for writing below.
            DComplex temp_element;
            std::string particles, collision_type;

            std::map<std::string, CorrelatorType::Correlator>::iterator correlator_iter;

            // Loop over all collisions, 9 for max.
            for (correlator_iter = correlators.begin(); correlator_iter != correlators.end(); correlator_iter++)
            {
#ifdef BUILD_HDF5
                // x and y stands for the 4-vectors of the two pions, respectivily
                std::string correlator_path = params.param.obj_path + "/" + correlator_iter->first + "/x" + std::to_string(origin[0]) + "_y" +  std::to_string(origin[1]) + "_z" + std::to_string(origin[2]) + "_t" + std::to_string(origin[3]);
                h5out.push(correlator_path);
#else
                std::string correlator_path = correlator_iter->first + "_x" + std::to_string(origin[0]) + "_y" +  std::to_string(origin[1]) + "_z" + std::to_string(origin[2]) + "_t" + std::to_string(origin[3]);
#endif
                // Now loop over momentum
                std::map<CorrelatorType::momenta_pair, multi1d<DComplex>>::iterator pq_iter;
                for (pq_iter = correlator_iter->second.begin(); pq_iter != correlator_iter->second.end(); pq_iter++)
                {
                    //One more temp variable instanited inside loop (once again for writing.)
                    multi1d<DComplex> pipi_correlator_towrite;
                    pipi_correlator_towrite.resize(Nt);
                    std::tuple<int, int, int> momenta1, momenta2;
                    momenta1 = std::get<0>(pq_iter->first);
                    momenta2 = std::get<1>(pq_iter->first);

#ifndef BUILD_HDF5
                    std::string correlator_path_mom = correlator_path + "_px" + std::to_string(std::get<0>(momenta1)) + "_py" + std::to_string(std::get<1>(momenta1)) + "_pz" + std::to_string(std::get<2>(momenta1)) + "_qx" + std::to_string(std::get<0>(momenta2)) + "_qy" + std::to_string(std::get<1>(momenta2)) + "_qz" + std::to_string(std::get<2>(momenta2));
                    TextFileWriter file_out(correlator_path_mom);
#endif
                    for (int t = 0; t < Nt; t++)
                    {
                        temp_element = pq_iter->second[t];
#ifndef BUILD_HDF5
                        file_out << temp_element << "\n";
#endif
                        pipi_correlator_towrite[t] = temp_element;
                    }
#ifndef BUILD_HDF5
                    file_out.close();
#else
                    //Change the name of string compared to 4d output so general correlator path is the same.
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
