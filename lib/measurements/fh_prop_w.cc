/*
Authors
Arjun Gambhir
Andre Walker-Loud

FH Propagator Task
This computes a Feynman-Hellmann (FH) propagator for bi-linear currents
https://arxiv.org/abs/1612.06963
INPUT
    Propagator
    List of Currents (spin, space, color, momentum)
    Parameters for linear solver
OUTPUT
    FH Propagator for each of the specified currents
*/

// Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
#include "fermact.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "util/info/unique_id.h"

// Lalibe Stuff
#include "../momentum/lalibe_sftmom.h"
#include "fh_prop_w.h"
#include "../matrix_elements/bilinear_gamma.h"

namespace Chroma
{
    namespace LalibeFHPropagatorEnv
    {
        namespace
        {
            AbsInlineMeasurement* createMeasurement(XMLReader& xml_in,
                const std::string& path)
            {
                return new InlineMeas(FHParams(xml_in, path));
            }
            //! Local registration flag
            bool registered = false;
        }
        const std::string name = "FH_PROPAGATOR";

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

        void read(XMLReader& xml, const std::string& path, FHParams::FHProp_t& par)
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

        void write(XMLWriter& xml, const std::string& path, FHParams::FHProp_t& par)
        {
            push(xml, path);
            write(xml, "currents" ,par.currents); //list of currents
            //write(xml, "t0"      ,par.t0     ); //t0 of input prop
            //write(xml, "j_decay" ,par.j_decay); //orthogoal direction of FT
            write(xml, "PropagatorParam" ,par.prop_param); //params for next lin solve
	    if(par.is_mom_max == true)
	      write(xml, "p2_max" ,par.p2_max);
	    else
	      write(xml, "mom_list" ,par.mom_list);
            pop(xml);

        }

        //! NamedObject input
        void read(XMLReader& xml, const std::string& path, FHParams::NamedObject_t& input)
        {
            XMLReader inputtop(xml, path);
            read(inputtop, "gauge_id"     , input.gauge_id);
            read(inputtop, "src_prop_id"  , input.src_prop_id);
            read(inputtop, "fh_prop_id"   , input.fh_prop_id);
        }

        //! NamedObject output
        void write(XMLWriter& xml, const std::string& path, const FHParams::NamedObject_t& input)
        {
            push(xml, path);
            write(xml, "gauge_id"     , input.gauge_id    );
            write(xml, "src_prop_id"  , input.src_prop_id     );
            write(xml, "fh_prop_id"   , input.fh_prop_id);
            pop(xml);
        }

        // Param stuff
        FHParams::FHParams()
        {
            frequency = 0;
        }

        FHParams::FHParams(XMLReader& xml_in, const std::string& path)
        {
            try
            {
                XMLReader paramtop(xml_in, path);
                if (paramtop.count("Frequency") == 1)
                    read(paramtop, "Frequency", frequency);
                else
                    frequency = 1;

                // Parameters for source construction
                read(paramtop, "FHParams", fhparam);

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

        void FHParams::writeXML(XMLWriter& xml_out, const std::string& path)
        {
            push(xml_out, path);
            write(xml_out, "FHParams", fhparam);
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
            QDPIO::cout << "FH_PROPAGATOR: start" << std::endl;

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
                QDPIO::cerr << LalibeFHPropagatorEnv::name
                    << ": caught dynamic cast error" << std::endl;
                QDP_abort(1);
            }
            catch (const std::string& e)
            {
                QDPIO::cerr << LalibeFHPropagatorEnv::name
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
	    LalibeSftMom ft = params.fhparam.is_mom_max ? LalibeSftMom(params.fhparam.p2_max, origin, false, j_decay)
              : LalibeSftMom(params.fhparam.p_list, origin, j_decay);
	    // Make an action and all other stuff needed for a solver.

	    typedef LatticeFermion T;
	    typedef multi1d<LatticeColorMatrix> P;
	    typedef multi1d<LatticeColorMatrix> Q;

	    std::istringstream xml_action(params.fhparam.prop_param.fermact.xml);
	    XMLReader action_reader(xml_action);
	    Handle<FermionAction<T, P, Q>> action(TheFermionActionFactory::Instance().createObject(params.fhparam.prop_param.fermact.id, action_reader, params.fhparam.prop_param.fermact.path));
	    Handle<FermState<T, P, Q>> action_state(action->createState(u));
	    QDPIO::cout<<"Our action and fermion state are doing A-okay so far."<<std::endl;
	    //Handle<SystemSolver<LatticeFermion>> solver = action->qprop(action_state, params.fhparam.prop_param.invParam);
	    //Above is for a single fermion, but we want to loop over spin/color and solve for the full propagator.

	    int ncg_had = 0; //This appears in the propagator task, I am just copying it here.

      	    LatticePropagator fh_prop_src = quark_propagator;
	    LatticePropagator fh_prop_solution;

            QDPIO::cout << "FH_PROPAGATOR: N_currents " << params.fhparam.currents.size() << std::endl;
	    for(int current_index = 0; current_index < params.fhparam.currents.size(); current_index++)
	    {
	      std::string present_current = params.fhparam.currents[current_index];
              QDPIO::cout << "FH_PROPAGATOR: current " << present_current << std::endl;
	      fh_prop_solution = zero;
              // WE SHOULD MAKE THIS A FACTORY
             Bilinear_Gamma(present_current, fh_prop_src, quark_propagator, u);

	      //Momentum loop
	      for(int mom = 0; mom < ft.numMom(); mom++)
	      {

		multi1d<int> momenta = ft.numToMom(mom);
		QDPIO::cout << "Injecting momentum - px: "<<std::to_string(momenta[0])<<" py: "+std::to_string(momenta[1])<<" pz: "+std::to_string(momenta[2])<<std::endl;
		fh_prop_src = ft[mom]*fh_prop_src;
	        //Now, we do the actual solve.
		action->quarkProp(fh_prop_solution, xml_out, fh_prop_src, t0, j_decay, action_state,
                params.fhparam.prop_param.invParam,
                params.fhparam.prop_param.quarkSpinType,
                params.fhparam.prop_param.obsvP, ncg_had);

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
  	        new_header.prop_header   = params.fhparam.prop_param;
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
	    QDPIO::cout << LalibeFHPropagatorEnv::name << ": total time = " << snoop.getTimeInSeconds() << " secs" << std::endl;
	    QDPIO::cout << LalibeFHPropagatorEnv::name<< ": ran successfully" << std::endl;
	    END_CODE();

	}
    }// LalibeFHPropagatorEnv
  };
