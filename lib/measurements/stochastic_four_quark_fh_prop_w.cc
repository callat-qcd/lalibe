/*
Authors
Arjun Gambhir

INPUT
    Propagator
    List of Currents (spin, space, color, momentum)
    Parameters for linear solver
    Information about scalar random noise
    (should match between two bilinears to construct the four quark operator)
OUTPUT
    FH Propagator for each of the specified currents
    (with random noise used as part of the second source) to be used in four quark contraction
*/

// Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
#include "util/ft/sftmom.h"
#include "fermact.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "util/info/unique_id.h"

// Lalibe Stuff
#include "stochastic_four_quark_fh_prop_w.h"
#include "../matrix_elements/chromomag_seqsource_w.h"
#include "../matrix_elements/bilinear_gamma.h"

namespace Chroma
{
    namespace LalibeStochasticFourQuarkFHPropagatorEnv
    {
        namespace
        {
            AbsInlineMeasurement* createMeasurement(XMLReader& xml_in,
                const std::string& path)
            {
                return new InlineMeas(StochasticFourQuarkFHParams(xml_in, path));
            }
            //! Local registration flag
            bool registered = false;
        }
        const std::string name = "STOCHASTIC_FOUR_QUARK_FH_PROPAGATOR";

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

        void read(XMLReader& xml, const std::string& path, StochasticFourQuarkFHParams::StochasticFourQuarkFHProp_t& par)
        {
            XMLReader paramtop(xml, path);
            read(paramtop, "currents" ,par.currents  ); //list of currents
            /*read(paramtop, "t0"      ,par.t0     );   //t0 of input prop
            if (paramtop.count("j_decay") != 0)
                read(paramtop, "j_decay" ,par.j_decay  ); //orthogoal direction of FT
            else{
                par.j_decay = Nd - 1;
                QDPIO::cout << "j_decay not specified - setting default j_decay = " <<
                par.j_decay << std::endl;
                }*/
            read(paramtop, "PropagatorParam" ,par.prop_param ); //params for next lin solve
	    read(paramtop, "vector_number" ,par.vector_number );
	    read(paramtop, "Seed" ,par.ran_seed );
	    read(paramtop, "ZN" ,par.ZN );
	    read(paramtop, "conjugate" ,par.conjugate );
	    if (paramtop.count("p2_max") != 0)
	    {
  	      read(paramtop, "p2_max" ,par.p2_max);
  	      QDPIO::cout<<"Reading momenta centered around the origin with a max of "<<par.p2_max<<std::endl;
  	      par.is_mom_max = true;
	    }
	    else if (paramtop.count("mom_list") != 0)
	    {
		read(paramtop, "mom_list" ,par.mom_list);
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
      		QDPIO::cout << "No momentum specified, setting FT to zero-momentum transfer only. "<<std::endl;
      		par.is_mom_max = true;
      		par.p2_max = 0;
	    }
        }

        void write(XMLWriter& xml, const std::string& path, StochasticFourQuarkFHParams::StochasticFourQuarkFHProp_t& par)
        {
            push(xml, path);
            write(xml, "currents" ,par.currents); //list of currents
            //write(xml, "t0"      ,par.t0     ); //t0 of input prop
            //write(xml, "j_decay" ,par.j_decay); //orthogoal direction of FT
            write(xml, "PropagatorParam" ,par.prop_param); //params for next lin solve
	    write(xml, "vector_number" ,par.vector_number);
	    write(xml, "Seed" ,par.ran_seed);
	    write(xml, "ZN" ,par.ZN);
	    write(xml, "conjugate" ,par.conjugate);
	    if(par.is_mom_max == true)
	      write(xml, "p2_max" ,par.p2_max);
	    else
	      write(xml, "mom_list" ,par.mom_list);
            pop(xml);

        }

        //! NamedObject input
        void read(XMLReader& xml, const std::string& path, StochasticFourQuarkFHParams::NamedObject_t& input)
        {
            XMLReader inputtop(xml, path);
            read(inputtop, "gauge_id"     , input.gauge_id);
            read(inputtop, "src_prop_id"  , input.src_prop_id);
            read(inputtop, "fh_prop_id"   , input.fh_prop_id);
        }

        //! NamedObject output
        void write(XMLWriter& xml, const std::string& path, const StochasticFourQuarkFHParams::NamedObject_t& input)
        {
            push(xml, path);
            write(xml, "gauge_id"     , input.gauge_id    );
            write(xml, "src_prop_id"  , input.src_prop_id     );
            write(xml, "fh_prop_id"   , input.fh_prop_id);
            pop(xml);
        }

        // Param stuff
        StochasticFourQuarkFHParams::StochasticFourQuarkFHParams()
        {
            frequency = 0;
        }

        StochasticFourQuarkFHParams::StochasticFourQuarkFHParams(XMLReader& xml_in, const std::string& path)
        {
            try
            {
                XMLReader paramtop(xml_in, path);
                if (paramtop.count("Frequency") == 1)
                    read(paramtop, "Frequency", frequency);
                else
                    frequency = 1;

                // Parameters for source construction
                read(paramtop, "StochasticFourQuarkFHParams", stochfourqfhparam);

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

        void StochasticFourQuarkFHParams::writeXML(XMLWriter& xml_out, const std::string& path)
        {
            push(xml_out, path);
            write(xml_out, "StochasticFourQuarkFHParams", stochfourqfhparam);
            write(xml_out, "NamedObject", named_obj);
            pop(xml_out);
        }

        void  InlineMeas::operator()(unsigned long update_no, XMLWriter& xml_out)
        {
            START_CODE();

            StopWatch snoop;
            snoop.reset();
            snoop.start();
            QDPIO::cout << "STOCHASTIC_FOUR_QUARK_FH_PROPAGATOR: start" << std::endl;

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
                QDPIO::cerr << LalibeStochasticFourQuarkFHPropagatorEnv::name
                    << ": caught dynamic cast error" << std::endl;
                QDP_abort(1);
            }
            catch (const std::string& e)
            {
                QDPIO::cerr << LalibeStochasticFourQuarkFHPropagatorEnv::name
                    << ": map call failed: " << e << std::endl;
                QDP_abort(1);
            }
            const multi1d<LatticeColorMatrix>& u =
                TheNamedObjMap::Instance().getData
                    <multi1d <LatticeColorMatrix> >(params.named_obj.gauge_id);

            // Read "src" quark propagator
            XMLReader prop_file_xml, prop_record_xml;
            LatticePropagator quark_propagator;
	    int t0;
	    int j_decay;
	    //Need origin for fourier transform!
	    multi1d<int> origin;
	    //We need this stuff to call quarkprop, it's pretty dumb, but I haven't found a way around it...
            QDPIO::cout << "Attempt to read forward propagator" << std::endl;
            try
            {
                quark_propagator = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.src_prop_id);
                TheNamedObjMap::Instance().get(params.named_obj.src_prop_id).getFileXML(prop_file_xml);
                TheNamedObjMap::Instance().get(params.named_obj.src_prop_id).getRecordXML(prop_record_xml);
		//This all assumes the incoming propagating is coming from a makesource, otherwise we are in a ton of trouble ~_~.
		MakeSourceProp_t  orig_header;
		read(prop_record_xml, "/Propagator", orig_header);
		j_decay = orig_header.source_header.j_decay;
		t0      = orig_header.source_header.t_source;
		origin = orig_header.source_header.getTSrce();
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

	    // Initialize FT stuff based on which type of momentum mode was selected.
	    // This assumes j_decay is the orthogonal direction: typically Nd - 1.
	    SftMom ft = params.stochfourqfhparam.is_mom_max ? SftMom(params.stochfourqfhparam.p2_max, origin, false, j_decay)
          : SftMom(params.stochfourqfhparam.p_list, j_decay);
	    // Make an action and all other stuff needed for a solver.

	    typedef LatticeFermion T;
	    typedef multi1d<LatticeColorMatrix> P;
	    typedef multi1d<LatticeColorMatrix> Q;

	    std::istringstream xml_action(params.stochfourqfhparam.prop_param.fermact.xml);
	    XMLReader action_reader(xml_action);
	    Handle<FermionAction<T, P, Q>> action(TheFermionActionFactory::Instance().createObject(params.stochfourqfhparam.prop_param.fermact.id, action_reader, params.stochfourqfhparam.prop_param.fermact.path));
	    Handle<FermState<T, P, Q>> action_state(action->createState(u));
	    QDPIO::cout<<"Our action and fermion state are doing A-okay so far."<<std::endl;
	    //Handle<SystemSolver<LatticeFermion>> solver = action->qprop(action_state, params.stochfourqfhparam.prop_param.invParam);
	    //Above is for a single fermion, but we want to loop over spin/color and solve for the full propagator.

	    int ncg_had = 0; //This appears in the propagator task, I am just copying it here.

      	    LatticePropagator fh_prop_src = quark_propagator;
	    LatticePropagator fh_prop_solution;

	    //Here's where the actual noisy stuff happens.
	    LatticeComplex vec ;
	    Seed ran_seed;
	    QDP::RNG::savern(ran_seed);

	    // Set the seed to desired value
	    QDP::RNG::setrn(params.stochfourqfhparam.ran_seed);

	    for(int current_vec = 0; current_vec < params.stochfourqfhparam.vector_number; current_vec++)
	    {
	      int N  = params.stochfourqfhparam.ZN;
	      LatticeReal rnd1, theta;
	      // twopi defined in chroma/lib/chromabase.h
	      Real twopiN = Chroma::twopi / N;
	      random(rnd1);
	      theta = twopiN * floor(N*rnd1);
	      vec = cmplx(cos(theta),sin(theta));
	    }

	    //restore the seed
	    QDP::RNG::setrn(ran_seed);

	    //conjugate noise if needed
	    //This may have to be current dependent, but we'll change that later.
	    if(params.stochfourqfhparam.conjugate == true)
	      vec = conj(vec);

            QDPIO::cout << "STOCHASTIC_FOUR_QUARK_FH_PROPAGATOR: N_currents " << params.stochfourqfhparam.currents.size() << std::endl;
	    for(int current_index = 0; current_index < params.stochfourqfhparam.currents.size(); current_index++)
	    {
	      std::string present_current = params.stochfourqfhparam.currents[current_index];
              QDPIO::cout << "STOCHASTIC_FOUR_QUARK_FH_PROPAGATOR: current " << present_current << std::endl;
	      fh_prop_solution = zero;
              // WE SHOULD MAKE THIS A FACTORY

	      //SpinMatrix bilinear_current;
	      //bilinear_current=Bilinear_Gamma(present_current);
	      //fh_prop_src = bilinear_current*quark_propagator;
	      Bilinear_Gamma(present_current, fh_prop_src, quark_propagator, u);

	      //Multiply the new source, by the scalar random noise.
	      fh_prop_src = vec*fh_prop_src;

	      //Momentum loop
	      for(int mom = 0; mom < ft.numMom(); mom++)
	      {

		fh_prop_src = ft[mom]*fh_prop_src;
	        //Now, we do the actual solve.
		action->quarkProp(fh_prop_solution, xml_out, fh_prop_src, t0, j_decay, action_state,
                params.stochfourqfhparam.prop_param.invParam,
                params.stochfourqfhparam.prop_param.quarkSpinType,
                params.stochfourqfhparam.prop_param.obsvP, ncg_had);

		push(xml_out,"Relaxation_Iterations");
  	        write(xml_out, "ncg_had", ncg_had);
  	        pop(xml_out);

  	        QDPIO::cout << "Writing propagator info, cause why not?" << std::endl;
  	        XMLBufferWriter file_xml;
  	        push(file_xml, "propagator");
  	        write(file_xml, "id", uniqueId());  // NOTE: new ID form
  	        pop(file_xml);

  	        //If ths src is not from make source these thing is not going to work...
  	        XMLBufferWriter record_xml;
  	        MakeSourceProp_t  orig_header;
  	        read(prop_record_xml, "/Propagator", orig_header);
  	        Propagator_t  new_header;   // note, abandoning state_info
  	        new_header.prop_header   = params.stochfourqfhparam.prop_param;
  	        new_header.source_header = orig_header.source_header;
  	        new_header.gauge_header  = orig_header.gauge_header;
  	        write(record_xml, "Propagator", new_header);

  	        // Pass the propagator info to the Named Object Buffer.
                // Looping over currents and momenta, momenta is inner most index
                // current_id = flattened index running over both these indices
  	        std::string current_id = params.named_obj.fh_prop_id[ft.numMom()*current_index + mom];
  	        TheNamedObjMap::Instance().create<LatticePropagator>(current_id);
  	        TheNamedObjMap::Instance().getData<LatticePropagator>(current_id) = fh_prop_solution;
  	        TheNamedObjMap::Instance().get(current_id).setFileXML(file_xml);
  	        TheNamedObjMap::Instance().get(current_id).setRecordXML(record_xml);
  	        QDPIO::cout<<"YAAAY! We finished current: "<<current_id<<std::endl;
	      }
	    }
	    snoop.stop();
	    QDPIO::cout << LalibeStochasticFourQuarkFHPropagatorEnv::name << ": total time = " << snoop.getTimeInSeconds() << " secs" << std::endl;
	    QDPIO::cout << LalibeStochasticFourQuarkFHPropagatorEnv::name<< ": ran successfully" << std::endl;
	    END_CODE();

	}
    }// LalibeStochasticFourQuarkFHPropagatorEnv
  };
