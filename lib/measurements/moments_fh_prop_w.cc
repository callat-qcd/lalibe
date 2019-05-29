/*
Authors
Arjun Gambhir

FH Propagator Task
This computes a Feynman-Hellmann (FH) propagator for bi-linear currents with a z^2 insertion (z can be in any direction)
INPUT
    Propagator
    List of Currents (spin, space, color, momentum, direction  of z^2)
    Parameters for linear solver
OUTPUT
    FH Propagator for each of the specified currents
*/

// Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
//#include "util/ft/sftmom.h"
#include "fermact.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "util/info/unique_id.h"

// Lalibe Stuff
#include "moments_fh_prop_w.h"
#include "../matrix_elements/chromomag_seqsource_w.h"

namespace Chroma
{
    namespace LalibeMomentsFHPropagatorEnv
    {
        namespace
        {
            AbsInlineMeasurement* createMeasurement(XMLReader& xml_in,
                const std::string& path)
            {
                return new InlineMeas(MomentsFHParams(xml_in, path));
            }
            //! Local registration flag
            bool registered = false;
        }
        const std::string name = "MOMENTS_FH_PROPAGATOR";

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

        void read(XMLReader& xml, const std::string& path, MomentsFHParams::MomentsFHProp_t& par)
        {
            XMLReader paramtop(xml, path);
            read(paramtop, "currents" ,par.currents  ); //list of currents
            read(paramtop, "PropagatorParam" ,par.prop_param ); //params for next lin solve
	    read(paramtop, "z_2_direction" ,par.moments_direction );
	    /*if (paramtop.count("p2_max") != 0)
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
      		QDPIO::cout << "No momentum specified, settings FT to zero-momentum transfer only. "<<std::endl;
      		par.is_mom_max = true;
      		par.p2_max = 0;
	    }*/
        }

        void write(XMLWriter& xml, const std::string& path, MomentsFHParams::MomentsFHProp_t& par)
        {
            push(xml, path);
            write(xml, "currents" ,par.currents); //list of currents
            write(xml, "PropagatorParam" ,par.prop_param); //params for next lin solve
	    write(xml, "z_2_direction" ,par.moments_direction);
	    /*if(par.is_mom_max == true)
	      write(xml, "p2_max" ,par.p2_max);
	    else
	      write(xml, "mom_list" ,par.mom_list);*/
            pop(xml);

        }

        //! NamedObject input
        void read(XMLReader& xml, const std::string& path, MomentsFHParams::NamedObject_t& input)
        {
            XMLReader inputtop(xml, path);
            read(inputtop, "gauge_id"     , input.gauge_id);
            read(inputtop, "src_prop_id"  , input.src_prop_id);
            read(inputtop, "fh_prop_id"   , input.fh_prop_id);
        }

        //! NamedObject output
        void write(XMLWriter& xml, const std::string& path, const MomentsFHParams::NamedObject_t& input)
        {
            push(xml, path);
            write(xml, "gauge_id"     , input.gauge_id    );
            write(xml, "src_prop_id"  , input.src_prop_id     );
            write(xml, "fh_prop_id"   , input.fh_prop_id);
            pop(xml);
        }

        // Param stuff
        MomentsFHParams::MomentsFHParams()
        {
            frequency = 0;
        }

        MomentsFHParams::MomentsFHParams(XMLReader& xml_in, const std::string& path)
        {
            try
            {
                XMLReader paramtop(xml_in, path);
                if (paramtop.count("Frequency") == 1)
                    read(paramtop, "Frequency", frequency);
                else
                    frequency = 1;

                // Parameters for source construction
                read(paramtop, "MomentsFHParams", momentsfhparam);

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

        void MomentsFHParams::writeXML(XMLWriter& xml_out, const std::string& path)
        {
            push(xml_out, path);
            write(xml_out, "MomentsFHParams", momentsfhparam);
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
            QDPIO::cout << "MOMENTS_FH_PROPAGATOR: start" << std::endl;

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
                QDPIO::cerr << LalibeMomentsFHPropagatorEnv::name
                    << ": caught dynamic cast error" << std::endl;
                QDP_abort(1);
            }
            catch (const std::string& e)
            {
                QDPIO::cerr << LalibeMomentsFHPropagatorEnv::name
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
	    int zsq_origin; //This is needed for the moments,, read below/
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
		zsq_origin = orig_header.source_header.getTSrce()[params.momentsfhparam.moments_direction];
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
	    //SftMom ft = params.momentsfhparam.is_mom_max ? SftMom(params.momentsfhparam.p2_max, false, j_decay)
          //: SftMom(params.momentsfhparam.p_list, j_decay);
            // No momentum stuff done here.
	    // Make an action and all other stuff needed for a solver.

	    typedef LatticeFermion T;
	    typedef multi1d<LatticeColorMatrix> P;
	    typedef multi1d<LatticeColorMatrix> Q;

	    std::istringstream xml_action(params.momentsfhparam.prop_param.fermact.xml);
	    XMLReader action_reader(xml_action);
	    Handle<FermionAction<T, P, Q>> action(TheFermionActionFactory::Instance().createObject(params.momentsfhparam.prop_param.fermact.id, action_reader, params.momentsfhparam.prop_param.fermact.path));
	    Handle<FermState<T, P, Q>> action_state(action->createState(u));
	    QDPIO::cout<<"Our action and fermion state are doing A-okay so far."<<std::endl;

	    int ncg_had = 0; //This appears in the propagator task, I am just copying it here.

	    // Multiply quark_propagator by z^2 (this is done before the inversion, ie: fh version of the moments method)
	    int z_L = Layout::lattSize()[params.momentsfhparam.moments_direction];
	    LatticeInteger moments_insertion = Layout::latticeCoordinate(params.momentsfhparam.moments_direction) - zsq_origin;
	    // Find the taxi-driver distance to the origin from each site. Since we square, the sign doesn't matter.
	    // The abs() function doesn't work for LatticeIntegers, darn....
	    // Instead I use a where function below, this is a little yucky.
	    //LatticeInteger abs_moments_insertion = abs(moments_insertion);
	    LatticeInteger abs_moments_insertion = where(moments_insertion < 0, -moments_insertion, moments_insertion);
            LatticeInteger moments_insertion_bc = where(abs_moments_insertion > z_L/2, z_L - abs_moments_insertion, abs_moments_insertion);
	    // Intercept usual code to multiply the source by z^2, since this commutes with all gammas.
	    // WARNING: This only works for zero-momentum transfer, for non-zero momentum
	    quark_propagator = moments_insertion_bc*moments_insertion_bc*quark_propagator;
	    LatticePropagator fh_prop_src = quark_propagator;
	    LatticePropagator fh_prop_solution;

            QDPIO::cout << "FH_PROPAGATOR: N_currents " << params.momentsfhparam.currents.size() << std::endl;
	    for(int current_index = 0; current_index < params.momentsfhparam.currents.size(); current_index++)
	    {
	      std::string present_current = params.momentsfhparam.currents[current_index];
              QDPIO::cout << "FH_PROPAGATOR: current " << present_current << std::endl;
	      fh_prop_solution = zero;
              // WE SHOULD MAKE THIS A FACTORY
	      if (present_current == "S")
		fh_prop_src =  Gamma(0)*quark_propagator;
              else if (present_current == "P")
                fh_prop_src =  Gamma(15)*quark_propagator;
              else if (present_current == "A3")
                fh_prop_src =  Gamma(11)*quark_propagator;
              else if (present_current == "V4")
                fh_prop_src =  Gamma(8)*quark_propagator;
              else if (present_current == "T12")
                fh_prop_src =  Gamma(3)*quark_propagator;
	      else if (present_current == "CHROMO_MAG")
		fh_prop_src =  chromoMagneticSeqSource(quark_propagator,u);
              else
              {
		QDPIO::cerr << present_current << ": NOT DEFINED YET " << std::endl;
		QDP_abort(1);
	      }

	      //As before, no momentum-transfer in this version.
	      //Momentum loop
	      //for(int mom = 0; mom < ft.numMom(); mom++)
	      //{

		//fh_prop_src = ft[mom]*fh_prop_src;
	        //Now, we do the actual solve.
		action->quarkProp(fh_prop_solution, xml_out, fh_prop_src, t0, j_decay, action_state,
                params.momentsfhparam.prop_param.invParam,
                params.momentsfhparam.prop_param.quarkSpinType,
                params.momentsfhparam.prop_param.obsvP, ncg_had);

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
  	        new_header.prop_header   = params.momentsfhparam.prop_param;
  	        new_header.source_header = orig_header.source_header;
  	        new_header.gauge_header  = orig_header.gauge_header;
  	        write(record_xml, "Propagator", new_header);

  	        // Pass the propagator info to the Named Object Buffer.
                // Looping over currents and momenta, momenta is inner most index
                // current_id = flattened index running over both these indices
  	        //std::string current_id = params.named_obj.fh_prop_id[ft.numMom()*current_index + mom];
  	        // Only one index here since momenta are not done in this routine.
		std::string current_id = params.named_obj.fh_prop_id[current_index];
  	        TheNamedObjMap::Instance().create<LatticePropagator>(current_id);
  	        TheNamedObjMap::Instance().getData<LatticePropagator>(current_id) = fh_prop_solution;
  	        TheNamedObjMap::Instance().get(current_id).setFileXML(file_xml);
  	        TheNamedObjMap::Instance().get(current_id).setRecordXML(record_xml);
  	        QDPIO::cout<<"YAAAY! We finished current: "<<current_id<<std::endl;
	      //}
	    }
	    snoop.stop();
	    QDPIO::cout << LalibeMomentsFHPropagatorEnv::name << ": total time = " << snoop.getTimeInSeconds() << " secs" << std::endl;
	    QDPIO::cout << LalibeMomentsFHPropagatorEnv::name<< ": ran successfully" << std::endl;
	    END_CODE();

	}
    }// LalibeMomentsFHPropagatorEnv
  };
