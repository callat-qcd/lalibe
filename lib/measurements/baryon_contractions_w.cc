/*
 * Single baryon contractions
 * Do baryon contractions and write out the two-point correlator in hdf5
 * Maybe we'll also support sdb, one day...
 */


#include "baryon_contractions_w.h"
#include "../contractions/baryon_contractions_func_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "../momentum/lalibe_sftmom.h"
#include "meas/inline/io/named_objmap.h"
#include "io/qprop_io.h"
#include "util/spin_basis.h"

#include <set>


namespace Chroma
{
    namespace LalibeBaryonContractionsEnv
    {
        namespace
        {
            AbsInlineMeasurement* createMeasurement(XMLReader& xml_in,
                                                    const std::string& path)
            {
                return new InlineMeas(BaryonParams(xml_in, path));
            }

            bool registered = false;
        }
        const std::string name = "BARYON_CONTRACTIONS";

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



        void read(XMLReader& xml, const std::string& path, BaryonParams::Param_t& par)
        {
            XMLReader paramtop(xml, path);

            read(paramtop, "ng_parity", par.ng_parity);
#ifdef BUILD_HDF5
            read(paramtop, "h5_file_name", par.file_name);
            read(paramtop, "path", par.obj_path);
            QDPIO::cout<<"HDF5 found, writing to"<<par.file_name<<" to path "<<par.obj_path<<std::endl;
#endif
            //We set output_full_correlator to true if no momentum is specified.
            //read(paramtop, "output_full_correlator", par.output_full_correlator);
            if (paramtop.count("rotate_to_Dirac") != 0)
                {
                    read(paramtop, "rotate_to_Dirac" ,par.rotate_to_Dirac);
                    QDPIO::cout<<"Rotating to Dirac basis from DeGrand-Rossi basis option set to "<<par.rotate_to_Dirac<<std::endl;
                }
            else
                {
                    QDPIO::cout<<"No XML option was specified regarding basis of propagators. Assuming DeGrand-Rossi and converting appropriately..."<<std::endl;
                    par.rotate_to_Dirac = true;
                }
            if (paramtop.count("is_antiperiodic") != 0)
                {
                    read(paramtop, "is_antiperiodic" ,par.is_antiperiodic);
                    QDPIO::cout<<"Ther user says that antiperiodic flag is set to "<<par.is_antiperiodic<<std::endl;
                }
            else
                {
                    QDPIO::cout<<"No XML option was specified regarding antiperiodicity. Assuming the quarks come from an antiperiodic lattice."<<std::endl;
                    par.is_antiperiodic = true;
                }
            if (paramtop.count("p2_max") != 0)
                {
                    read(paramtop, "p2_max" ,par.p2_max);
                    par.output_full_correlator = false;
                    QDPIO::cout<<"Reading momenta centered around the origin with a max of "<<par.p2_max<<std::endl;
                    par.is_mom_max = true;
                }
            else if (paramtop.count("mom_list") != 0)
                {
                    read(paramtop, "mom_list" ,par.mom_list);
                    par.output_full_correlator = false;
                    QDPIO::cout<<"Using custom momentum list."<<std::endl;
                    par.is_mom_max = false;
                    //Assumes the length is 3 for the inner dimension.
                    for(int iter = 0; iter < par.mom_list.size(); iter++)
                        QDPIO::cout<<"Momentum px: "<<par.mom_list[iter][0]<<" py: "<<par.mom_list[iter][1]<<" pz: "<<par.mom_list[iter][2]<<std::endl;

                    //nested multi1d to multi2d
                    par.p_list.resize(par.mom_list.size(), Nd -1);
                    for(int mom = 0; mom < par.mom_list.size(); mom++)
                        {
                            par.p_list[mom][0] = par.mom_list[mom][0];
                            par.p_list[mom][1] = par.mom_list[mom][1];
                            par.p_list[mom][2] = par.mom_list[mom][2];
                        }
                }
            else
                {
                    QDPIO::cout << "No momentum specified, not doing any FTs and dumping full 4d-correlator."<<std::endl;
                    par.output_full_correlator = true;
                    //Below is only so SftMom does not crash...
                    par.is_mom_max = true;
                    par.p2_max = 0;
                }

            multi1d<std::string> tmpPartList;
            read(paramtop, "particle_list", tmpPartList);

            // now go parse the particle list
            par.particle_list.clear();

            // small inline function to add all spin components for a given baryon
            // name to the list
            // if negative parity is enabled, add those computations too
            auto addFlavor = [&par](const std::string& flavName) {
                auto spinComponents = get_spin_components(flavName);
                for (auto aS : spinComponents)
                    par.particle_list.insert(std::make_pair(flavName, aS));

                if (par.ng_parity) {
                    auto spinComponents = get_spin_components(flavName+"_np");
                    for (auto aS : spinComponents)
                        par.particle_list.insert(std::make_pair(flavName+"_np", aS));
                }
            };

            for(unsigned int iPart = 0; iPart < tmpPartList.size(); iPart++) {

                // go look if this 'particle' is an alias for a group of particles
                auto alIt = aliasMap.find(tmpPartList[iPart]);
                if (alIt == aliasMap.end())
                    addFlavor(tmpPartList[iPart]);
                else {
                    for (auto aFl : alIt->second)
                        addFlavor(aFl);
                }
            }

        }


        void write(XMLWriter& xml, const std::string& path, const BaryonParams::Param_t& par)
        {
            push(xml, path);

            write(xml, "ng_parity", par.ng_parity);
            write(xml, "rotate_to_Dirac", par.rotate_to_Dirac);
            write(xml, "is_antiperiodic", par.is_antiperiodic);
#ifdef BUILD_HDF5
            write(xml, "h5_file_name", par.file_name);
            write(xml, "path", par.obj_path);
#endif
            //write(xml, "output_full_correlator", par.output_full_correlator);
            if(par.is_mom_max == true)
                write(xml, "p2_max" ,par.p2_max);
            else
                write(xml, "mom_list" ,par.mom_list);
            //write(xml, "particle_list", par.particle_list);

            pop(xml);
        }

        void read(XMLReader& xml, const std::string& path, BaryonParams::NamedObject_t& input)
        {
            XMLReader inputtop(xml, path);

            //read(inputtop, "gauge_id" , input.gauge_id);
            //Logic to read quark propagators (whichever ones are present.)
            //  xmlQuarkMapper map internal character identifiers to the expected
            //  XML tag
            std::vector<char> quark_flavs { 'u', 'd', 's', 'c' };
            std::map<char, std::string> xmlQuarkMapper = {
                {'u', "up_quark"},
                {'d', "down_quark"},
                {'s', "strange_quark"},
                {'c', "charm_quark"}
            };

            for (auto aFlav : quark_flavs) {
                if (inputtop.count(xmlQuarkMapper[aFlav]) != 0)
                    {
                        std::string tmpName;
                        read(inputtop, xmlQuarkMapper[aFlav], tmpName);
                        input.quark_map[aFlav] = tmpName;
                        QDPIO::cout<<"I found a "<<aFlav<<" quark, here is its id: "<<input.quark_map[aFlav]<<std::endl;
                    }
                else
                    {
                        QDPIO::cout<<"I couldn't find an "<<aFlav<<" quark, hope you don't need it for the inputted baryon contractions. "<<std::endl;
                    }
            }
        }

        void write(XMLWriter& xml, const std::string& path, const BaryonParams::NamedObject_t& input)
        {
            push(xml, path);
            //write(xml, "gauge_id" , input.gauge_id);
            std::vector<char> quark_flavs { 'u', 'd', 's', 'c' };
            std::map<char, std::string> xmlQuarkMapper = {
                {'u', "up_quark"},
                {'d', "down_quark"},
                {'s', "strange_quark"},
                {'c', "charm_quark"}
            };

            for (auto aFlav : quark_flavs) {
                auto qIt = input.quark_map.find(aFlav);
                if (qIt != input.quark_map.end())
                    write(xml, xmlQuarkMapper[aFlav], qIt->second);
            }
            pop(xml);
        }


        BaryonParams::BaryonParams()
        {
            frequency = 0;
        }

        BaryonParams::BaryonParams(XMLReader& xml_in, const std::string& path)
        {
            try
                {
                    XMLReader paramtop(xml_in, path);

                    if (paramtop.count("Frequency") == 1)
                        read(paramtop, "Frequency", frequency);
                    else
                        frequency = 1;

                    read(paramtop, "BaryonParams", param);

                    read(paramtop, "NamedObject", named_obj);

                }
            catch(const std::string& e)
                {
                    QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
                    QDP_abort(1);
                }
        }


        void
        BaryonParams::writeXML(XMLWriter& xml_out, const std::string& path)
        {
            push(xml_out, path);

            write(xml_out, "BaryonParams", param);
            write(xml_out, "NamedObject", named_obj);

            pop(xml_out);
        }


        void  InlineMeas::operator()(unsigned long update_no, XMLWriter& xml_out)
        {
            START_CODE();

            StopWatch snoop;
            snoop.reset();
            snoop.start();

            QDPIO::cout<<"Baryon contractions starting..."<<std::endl;

            //Grab all the propagators that are given.
            std::vector<char> quark_flavs { 'u', 'd', 's', 'c' };
            std::map<char, LatticePropagator> prop_map;
            //Need origin, j_decay, and t0 for fourier transform!
            //Need j_decay of bc to know what comes with a minus sign.
            int j_decay;
            int t_0;
            multi1d<int> origin;

            // if we read more than one propagator, check that origin, t_0 and
            // j_decay match between them
            bool checkFlag = false;

            for (auto aFlav : quark_flavs) {
                auto qIt = params.named_obj.quark_map.find(aFlav);
                if (qIt != params.named_obj.quark_map.end()) {

                    QDPIO::cout << "Attempting to read "<<aFlav<<" propagator" << std::endl;

                    try
                        {
                            prop_map[aFlav] = TheNamedObjMap::Instance().getData<LatticePropagator>(qIt->second);

                            XMLReader prop_file_xml, prop_record_xml;
                            TheNamedObjMap::Instance().get(qIt->second).getFileXML(prop_file_xml);
                            TheNamedObjMap::Instance().get(qIt->second).getRecordXML(prop_record_xml);
                            //Get the origin  and j_decay for the FT, this assumes all quarks have the same origin.
                            MakeSourceProp_t  orig_header;
                            if (prop_record_xml.count("/Propagator") != 0)
                                {
                                    QDPIO::cout<<aFlav<<" quark propagator is unsmeared, reading from Propagator tag..."<<std::endl;
                                    read(prop_record_xml, "/Propagator", orig_header);
                                }
                            else if (prop_record_xml.count("/SinkSmear") != 0)
                                {
                                    QDPIO::cout<<aFlav<<" quark propagator is smeared, reading from SinkSmear tag..."<<std::endl;
                                    read(prop_record_xml, "/SinkSmear", orig_header);
                                }
                            else
                                {
                                    QDPIO::cout<<"What type of weird propagator did you give me? I can't find the right tag to get j_decay, source location, and whatnot...exiting..."<<std::endl;
                                }
                            // set info
                            if (!checkFlag) {
                                j_decay = orig_header.source_header.j_decay;
                                t_0 = orig_header.source_header.t_source;
                                origin = orig_header.source_header.getTSrce();
                            }
                            // check info
                            else {
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
                            //If we need to rotate, we do it now.
                            if(params.param.rotate_to_Dirac == true)
                                rotate_to_Dirac_Basis(prop_map[aFlav]);
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

                    // we've read the first header; from now on compare that remaining
                    // props have matching origin
                    checkFlag = true;
                } // check flavor
            } // flavor loop


            //Initialize FT stuff here, whether this is used or not below is another story...
            LalibeSftMom ft = params.param.is_mom_max ? LalibeSftMom(params.param.p2_max, origin, false, j_decay)
                : LalibeSftMom(params.param.p_list, origin, j_decay);

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
            //Loop over list of particles and check that all the necessary flavors are present.
            std::set<char> reqProps;
            for(auto aParticle : params.param.particle_list){
                auto flavCode = get_flavor_code(aParticle.first);
                reqProps.insert(std::get<0>(flavCode));
                reqProps.insert(std::get<1>(flavCode));
                reqProps.insert(std::get<2>(flavCode));
            }

            for (auto aProp : reqProps) {
                if (prop_map.find(aProp) == prop_map.end()) {
                    QDPIO::cerr << "Could not find required propagator for "<<aProp<<" quark"<<std::endl;
                    QDP_abort(1);
                }
            }

            //If flavor check has passed, now we loop through particles and do the contractions.
            for (auto aParticle : params.param.particle_list){
                QDPIO::cout<<"Starting "<<aParticle.first<<" "<<aParticle.second<<" contraction..."<<std::endl;

                auto flavCode = get_flavor_code(aParticle.first);

                LatticeComplex baryon = zero;

                do_contraction(prop_map[std::get<0>(flavCode)],
                               prop_map[std::get<1>(flavCode)],
                               prop_map[std::get<2>(flavCode)],
                               aParticle.first,
                               aParticle.second,
                               baryon);

                write_correlator(params.param.output_full_correlator, params.param.is_antiperiodic,
                                 aParticle.first, aParticle.second,
#ifdef BUILD_HDF5
                                 params.param.obj_path, h5out, wmode,
#endif
                                 t_0, Nt, origin, ft, baryon);

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
