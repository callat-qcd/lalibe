/*
Authors
Arjun Gambhir

INPUT
    Propagator or Source
OUTPUT
    Propagator or Source that's projected in the laph smearing subspace.
*/

// Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
#include "momentum/lalibe_sftmom.h"
#include "fermact.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "util/info/unique_id.h"
#include "util/ferm/transf.h"

// Lalibe Stuff
#include "laph_smearing_prop_w.h"

namespace Chroma
{
    namespace LalibeLaphSmearingPropagatorEnv
    {
        namespace
        {
            AbsInlineMeasurement* createMeasurement(XMLReader& xml_in,
                const std::string& path)
            {
                return new InlineMeas(LaphSmearingParams(xml_in, path));
            }
            //! Local registration flag
            bool registered = false;
        }
        const std::string name = "LAPH_SMEARING_PROPAGATOR";

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

        void read(XMLReader& xml, const std::string& path, LaphSmearingParams::LaphSmearingProp_t& par)
        {
            XMLReader paramtop(xml, path);
	    read(paramtop, "laph_evecs" ,par.evec_file  ); //list of currents
        }

        void write(XMLWriter& xml, const std::string& path, LaphSmearingParams::LaphSmearingProp_t& par)
        {
            push(xml, path);
	    write(xml, "laph_evecs" ,par.evec_file); //list of currents
        }

        //! NamedObject input
        void read(XMLReader& xml, const std::string& path, LaphSmearingParams::NamedObject_t& input)
        {
            XMLReader inputtop(xml, path);
            read(inputtop, "gauge_id"     , input.gauge_id);
            read(inputtop, "prop_id"  ,     input.prop_id);
	    read(inputtop, "laph_prop_id"  , input.laph_prop_id);
        }

        //! NamedObject output
        void write(XMLWriter& xml, const std::string& path, const LaphSmearingParams::NamedObject_t& input)
        {
            push(xml, path);
            write(xml, "gauge_id"     , input.gauge_id    );
            write(xml, "prop_id"  ,     input.prop_id     );
            write(xml, "laph_prop_id" , input.laph_prop_id);
            pop(xml);
        }

        // Param stuff
        LaphSmearingParams::LaphSmearingParams()
        {
            frequency = 0;
        }

        LaphSmearingParams::LaphSmearingParams(XMLReader& xml_in, const std::string& path)
        {
            try
            {
                XMLReader paramtop(xml_in, path);
                if (paramtop.count("Frequency") == 1)
                    read(paramtop, "Frequency", frequency);
                else
                    frequency = 1;

                // Parameters for source construction
                read(paramtop, "LaphSmearingParams", laphparam);

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

        void LaphSmearingParams::writeXML(XMLWriter& xml_out, const std::string& path)
        {
            push(xml_out, path);
            write(xml_out, "LaphSmearingParams", laphparam);
            write(xml_out, "NamedObject", named_obj);
            pop(xml_out);
        }

        void  InlineMeas::operator()(unsigned long update_no, XMLWriter& xml_out)
        {
            START_CODE();

            StopWatch snoop;
            snoop.reset();
            snoop.start();
            QDPIO::cout << "LAPH_SMEARING_PROPAGATOR: start" << std::endl;

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
                QDPIO::cerr << LalibeLaphSmearingPropagatorEnv::name
                    << ": caught dynamic cast error" << std::endl;
                QDP_abort(1);
            }
            catch (const std::string& e)
            {
                QDPIO::cerr << LalibeLaphSmearingPropagatorEnv::name
                    << ": map call failed: " << e << std::endl;
                QDP_abort(1);
            }
            const multi1d<LatticeColorMatrix>& u =
                TheNamedObjMap::Instance().getData
                    <multi1d <LatticeColorMatrix> >(params.named_obj.gauge_id);

	    //Either do a normal source or sequential source.
	    bool propagatorP = false;
	    bool makesourceP = false;

	    // Read quark propagator
            XMLReader prop_file_xml, prop_record_xml;
            LatticePropagator quark_propagator;
            QDPIO::cout << "Attempt to read forward propagator" << std::endl;
            try
            {
                quark_propagator = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop_id);
                TheNamedObjMap::Instance().get(params.named_obj.prop_id).getFileXML(prop_file_xml);
                TheNamedObjMap::Instance().get(params.named_obj.prop_id).getRecordXML(prop_record_xml);
		if (prop_record_xml.count("/Propagator") != 0)
		{
		    propagatorP = true;
		    MakeSourceProp_t  orig_header;
		    read(prop_record_xml, "/Propagator", orig_header);
		}
		else if (prop_record_xml.count("/MakeSource") != 0)
		{
		    makesourceP = true;
		    SequentialSource_t   orig_header;
		    read(prop_record_xml, "/MakeSource", orig_header);
		}
		else
		{
		  throw std::string("This is some stupid source I dont' know how to deal with.");
		}
            }
            catch (std::bad_cast)
            {
                QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
                QDP_abort(1);
            }
            catch (const std::string& e)
            {
                QDPIO::cerr << name << ": error reading prop_header: "
                    << e << std::endl;
                QDP_abort(1);
            }

	    QDPIO::cout << "Attempt to read laph evecs." << std::endl;
	    StopWatch read_snoop;
	    read_snoop.reset();
	    read_snoop.start();
	    multi1d<LatticeColorVector> evecs;
	    int evecs_length;
	    // Open evecs file
	    XMLReader evec_file_xml;
	    //QDPFileReader to(evec_file_xml,params.laphparam.evec_file,QDPIO_SERIAL);
	    QDPFileReader to(evec_file_xml,params.laphparam.evec_file,QDPIO_PARALLEL);
	    read(evec_file_xml, "nEv", evecs_length);
	    evecs.resize(evecs_length);
	    XMLReader evec_record_xml;
	    read(to, evec_record_xml, evecs);
	    read_snoop.stop();
	    QDPIO::cout << "Reading "<< evecs_length << " eigenvectors completed in " 
		<< read_snoop.getTimeInSeconds() << " seconds." << std::endl;

	    //LatticeFermion prop_ferm = zero;
	    LatticeFermion evec_ferm = zero;
	    LatticePropagator laph_prop = zero;
	    //LatticeFermion laph_ferm = zero;
	    
	    // Optimization: do the copying once, outside of the loop.
	    multi1d<LatticeFermion> prop_ferm;
	    multi1d<LatticeFermion> laph_ferm;
	    prop_ferm.resize(Ns*Nc);
	    laph_ferm.resize(Ns*Nc);
	    
	    for(int color_source(0);color_source<Nc;color_source++){
		for(int spin_source=0; spin_source < Ns; spin_source++){
		    PropToFerm(quark_propagator, prop_ferm[spin_source + color_source*Ns], color_source, spin_source);
		    PropToFerm(laph_prop, laph_ferm[spin_source + color_source*Ns], color_source, spin_source);
		}
	    }
	    
	    // Evecs loop
	    for(int vec_index = 0; vec_index < evecs.size(); vec_index++)
	    {
		  QDPIO::cout << "Projecting vector " << vec_index << std::endl;
		  for(int eigenvec_spin_source=0; eigenvec_spin_source < Ns; eigenvec_spin_source++){
		      // Brute force way of doing the spin identity, find a shortcut once this is verified.
		      CvToFerm(evecs[vec_index], evec_ferm, eigenvec_spin_source);
		      for(int color_source(0);color_source<Nc;color_source++){
			  for(int spin_source=0; spin_source < Ns; spin_source++){
			      //PropToFerm(quark_propagator, prop_ferm, color_source, spin_source);
			      //PropToFerm(laph_prop, laph_ferm, color_source, spin_source);
			      //laph_ferm += evec_ferm*innerProduct(evec_ferm, prop_ferm);
			      //FermToProp(laph_ferm, laph_prop, color_source, spin_source);
			      laph_ferm[spin_source + color_source*Ns] += 
				  evec_ferm*innerProduct(evec_ferm, prop_ferm[spin_source + color_source*Ns]);
			  }
		      }
		  }
	    }
	    
	    for(int color_source(0);color_source<Nc;color_source++){
		for(int spin_source=0; spin_source < Ns; spin_source++){
		    FermToProp(laph_ferm[spin_source + color_source*Ns], laph_prop, color_source, spin_source);
		}
	    }
	    LatticePropagator smearedProp = zero;
	    for (unsigned int srcSpin=0; srcSpin<4; ++srcSpin) {
	      for (unsigned int snkSpin=0; snkSpin<4; ++snkSpin) {
		LatticeColorMatrix tmp = peekSpin(quark_propagator, snkSpin, srcSpin);
		for (unsigned int srcCol=0; srcCol<3; ++srcCol) {
		  LatticeColorVector sinkVec = zero;
		    // peek source color
		  for (unsigned int snkCol=0; snkCol<3; ++snkCol)
		    pokeColor(sinkVec, peekColor(tmp, snkCol, srcCol), snkCol);
		    // smeared component of the propagator
		  LatticeColorVector smearVec = zero;
		  for (unsigned int iEv=0; iEv<evecs.size(); ++iEv) {
		    smearVec += innerProduct(evecs[iEv], sinkVec) * evecs[iEv];
		  } 
		    // poke source color
		  for (unsigned int snkCol=0; snkCol<3; ++snkCol)
		    pokeColor(tmp, peekColor(smearVec, snkCol), snkCol, srcCol);
		}
		pokeSpin(smearedProp, tmp, snkSpin, srcSpin);
	      }
	    }
	    LatticePropagator diff = smearedProp - laph_prop;
	    for(int color_source(0);color_source<Nc;color_source++){
		for(int spin_source=0; spin_source < Ns; spin_source++){
		    LatticeFermion temp = zero;
		    PropToFerm(diff, temp, color_source, spin_source);
		    QDPIO::cout<<" Diff Norm 2 spin: "<<spin_source<<" color: "<<color_source<<" "
			<<innerProduct(temp, temp)<<std::endl;
		}
	    }
	      
	    // Not even needed since xml record is not being modified in any way.
	    /*if (propagatorP)
	    {
		//MakeSourceProp_t  orig_header;
		Propagator_t  orig_header;
		read(prop_record_xml, "/Propagator", orig_header);
		write(record_xml, "Propagator", new_header);
	      }
	      else if (makesourceP)
	      {
		MakeSourceProp_t  orig_header;
		read(prop_record_xml, "/MakeSource", orig_header);
		write(record_xml, "/MakeSource", orig_header);
	      }*/

	    TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.laph_prop_id);
	    TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.laph_prop_id) = laph_prop;
	    TheNamedObjMap::Instance().get(params.named_obj.laph_prop_id).setFileXML(prop_file_xml);
	    TheNamedObjMap::Instance().get(params.named_obj.laph_prop_id).setRecordXML(prop_record_xml);

	    snoop.stop();
	    QDPIO::cout << LalibeLaphSmearingPropagatorEnv::name << ": total time = " << snoop.getTimeInSeconds() << " secs" << std::endl;
	    QDPIO::cout << LalibeLaphSmearingPropagatorEnv::name<< ": ran successfully" << std::endl;
	    END_CODE();

	}
    }// LalibeLaphSmearingPropagatorEnv
  };
