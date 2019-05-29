/*
Arjun Gambhir

This computes a Feynman-Hellmann a fully diluted spin/color ZN (stochastic) propagator.
This is useful for disco calculations, but can also be used to intercept a quark line in a hadron to create a stochastic Feyman Hellman propagator.
This code will allow batches of ZN vectors to be inverted at once, but the user must be careful to name them uniquely based on the current "vector number".
INPUT
    Type of ZN noise.
    Random noise seed.
    Starting and ending vector number within the seed.
    Parameters for linear solver
OUTPUT
    ZN propagator, note the header info won't be right, but that's okay for disco and stochastic fh_prop
*/

// Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
#include "fermact.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "util/info/unique_id.h"
//We need this to insert Fermions to a Prop, this is imperitive for dilution.
#include "util/ferm/transf.h"

// Lalibe Stuff
#include "ZN_prop_w.h"

namespace Chroma
{
    namespace LalibeZNPropagatorEnv
    {
        namespace
        {
            AbsInlineMeasurement* createMeasurement(XMLReader& xml_in,
                const std::string& path)
            {
                return new InlineMeas(ZNParams(xml_in, path));
            }
            //! Local registration flag
            bool registered = false;
        }
        const std::string name = "ZN_PROPAGATOR";

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

        void read(XMLReader& xml, const std::string& path, ZNParams::ZNProp_t& par)
        {
            XMLReader paramtop(xml, path);
            read(paramtop, "PropagatorParam" ,par.prop_param ); 
	    read(paramtop, "starting_vector" ,par.starting_vector ); 
	    read(paramtop, "ending_vector" ,par.ending_vector ); 
	    read(paramtop, "Seed" ,par.ran_seed ); 
	    read(paramtop, "ZN" ,par.ZN ); 
        }

        void write(XMLWriter& xml, const std::string& path, ZNParams::ZNProp_t& par)
        {
            push(xml, path);
            write(xml, "PropagatorParam" ,par.prop_param);
	    write(xml, "starting_vector" ,par.starting_vector); 
	    write(xml, "ending_vector" ,par.ending_vector); 
	    write(xml, "Seed" ,par.ran_seed);
	    write(xml, "ZN" ,par.ZN);
            pop(xml);

        }

        //! NamedObject input
        void read(XMLReader& xml, const std::string& path, ZNParams::NamedObject_t& input)
        {
            XMLReader inputtop(xml, path);
	    read(inputtop, "gauge_id"     , input.gauge_id);
            read(inputtop, "zn_prop_id"   , input.zn_prop_id);
        }

        //! NamedObject output
        void write(XMLWriter& xml, const std::string& path, const ZNParams::NamedObject_t& input)
        {
            push(xml, path);
	    write(xml, "gauge_id"     , input.gauge_id    );
            write(xml, "zn_prop_id"   , input.zn_prop_id);
            pop(xml);
        }

        // Param stuff
        ZNParams::ZNParams()
        {
            frequency = 0;
        }

        ZNParams::ZNParams(XMLReader& xml_in, const std::string& path)
        {
            try
            {
                XMLReader paramtop(xml_in, path);
                if (paramtop.count("Frequency") == 1)
                    read(paramtop, "Frequency", frequency);
                else
                    frequency = 1;

                // Parameters for source construction
                read(paramtop, "ZNParams", znparam);

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

        void ZNParams::writeXML(XMLWriter& xml_out, const std::string& path)
        {
            push(xml_out, path);
            write(xml_out, "ZNParams", znparam);
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
            QDPIO::cout << "ZN_PROPAGATOR: start" << std::endl;

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
		QDPIO::cerr << LalibeZNPropagatorEnv::name
		    << ": caught dynamic cast error" << std::endl;
		QDP_abort(1);
	    }
	    catch (const std::string& e)
	    {
		QDPIO::cerr << LalibeZNPropagatorEnv::name
		    << ": map call failed: " << e << std::endl;
		QDP_abort(1);
	    }
	    const multi1d<LatticeColorMatrix>& u =
		TheNamedObjMap::Instance().getData
		    <multi1d <LatticeColorMatrix> >(params.named_obj.gauge_id);

	    // Make an action and all other stuff needed for a solver.

	    typedef LatticeFermion T;
	    typedef multi1d<LatticeColorMatrix> P;
	    typedef multi1d<LatticeColorMatrix> Q;

	    std::istringstream xml_action(params.znparam.prop_param.fermact.xml);
	    XMLReader action_reader(xml_action);
	    Handle<FermionAction<T, P, Q>> action(TheFermionActionFactory::Instance().createObject(params.znparam.prop_param.fermact.id, action_reader, params.znparam.prop_param.fermact.path));
	    Handle<FermState<T, P, Q>> action_state(action->createState(u));
	    Handle<SystemSolver<LatticeFermion>> solver = action->qprop(action_state, params.znparam.prop_param.invParam);
	    QDPIO::cout<<"Our action and fermion state are doing A-okay so far."<<std::endl;
	    //We are going to manually do Nc*Ns inversions, looping over spin/color and creating and solved noise propagator which doesn't require header stuff.

	    //Here's where the actual noisy stuff happens.
	    multi1d <LatticeComplex> vectors;
	    vectors.resize(params.znparam.ending_vector - params.znparam.starting_vector + 1);
	    LatticeComplex vec ;
	    Seed ran_seed;
	    QDP::RNG::savern(ran_seed);

	    // Set the seed to desired value
	    QDP::RNG::setrn(params.znparam.ran_seed);

	    for(int current_vec = 0; current_vec < params.znparam.ending_vector; current_vec++)
	    {
	      int N  = params.znparam.ZN;
	      LatticeReal rnd1, theta;
	      // twopi defined in chroma/lib/chromabase.h
	      Real twopiN = Chroma::twopi / N; 
	      random(rnd1); 
	      theta = twopiN * floor(N*rnd1);
	      vec = cmplx(cos(theta),sin(theta));
	      if ((current_vec + 1) >= params.znparam.starting_vector)
	        vectors[current_vec - params.znparam.starting_vector + 1] = vec; 
	    }
	    
	    //restore the seed
	    QDP::RNG::setrn(ran_seed);


	    for(int vec_index = 0; vec_index < vectors.size(); vec_index++)
	    {
	       LatticePropagator noise_prop = zero;
	       //Now let's do some dilution...
	       for(int color_source(0);color_source<Nc;color_source++){
		 QDPIO::cout << "color_source = " << color_source << std::endl; 
		 
		 LatticeColorVector vec_srce = zero ;
		 //	      pokeColor(tmp,wave[k],color_source) ;
		 pokeColor(vec_srce,vectors[vec_index],color_source) ;
		 
		 for(int spin_source=0; spin_source < Ns; ++spin_source){
		   QDPIO::cout << "spin_source = " << spin_source << std::endl; 
		   
		   // Insert a ColorVector into spin index spin_source
		   // This only overwrites sections, so need to initialize first
		   // Also, do the solve here
		   LatticeFermion chi = zero;
		   LatticeFermion noise_soln = zero;
		   CvToFerm(vec_srce, chi, spin_source);
		   SystemSolverResults_t res = (*solver)(noise_soln, chi);
		   FermToProp(noise_soln, noise_prop, color_source, spin_source); 
		 }
	       }

	      //Fake some propagator info that isn't relevant for stochastic ones.
	      XMLBufferWriter file_xml;
	      push(file_xml, "propagator");
	      write(file_xml, "id", uniqueId());  // NOTE: new ID form
	      pop(file_xml);
	      XMLBufferWriter record_xml;
	      push(record_xml, "propagator");
	      write(record_xml, "id", uniqueId());  // NOTE: new ID form
	      pop(record_xml);

  	      // Pass the propagator to the Named Object Buffer.
              // No other info is passed here.
  	      std::string current_id = params.named_obj.zn_prop_id[vec_index];
  	      TheNamedObjMap::Instance().create<LatticePropagator>(current_id);
  	      TheNamedObjMap::Instance().getData<LatticePropagator>(current_id) = noise_prop;
	      TheNamedObjMap::Instance().get(current_id).setFileXML(file_xml);
	      TheNamedObjMap::Instance().get(current_id).setRecordXML(record_xml);
  	      QDPIO::cout<<"Passed noise vector: "<<current_id<<"to the Named Object Buffer."<<std::endl;
	    }
	    snoop.stop();
	    QDPIO::cout << LalibeZNPropagatorEnv::name << ": total time = " << snoop.getTimeInSeconds() << " secs" << std::endl;
	    QDPIO::cout << LalibeZNPropagatorEnv::name<< ": ran successfully" << std::endl;
	    END_CODE();

	}
    }// LalibeZNPropagatorEnv
  };
