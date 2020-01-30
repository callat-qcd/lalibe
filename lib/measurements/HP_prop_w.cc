/*
Arjun Gambhir

This computes a Feynman-Hellmann a fully diluted spin/color hierarchical probing (HP) propagator (with an element-wise product with random noise).
This is useful for disco calculations, but can also be used to intercept a quark line in a hadron to create a stochastic Feyman Hellman propagator.
This code will allow batches of HP vectors to be inverted at once, but the user must be careful to name them uniquely based on the current "vector number".
INPUT
    Type of ZN noise to be used for element-wise Hadamard product.
    Random noise seed.
    Starting and ending vector number for hierarchical probing.
    Parameters for linear solver
OUTPUT
    HP propagator, note the header info won't be right, but that's okay for disco and stochastic fh_prop
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
#include "HP_prop_w.h"
#include "../numerics/binaryRecursiveColoring.h"

namespace Chroma
{
    namespace LalibeHPPropagatorEnv
    {
        namespace
        {
            AbsInlineMeasurement* createMeasurement(XMLReader& xml_in,
                const std::string& path)
            {
                return new InlineMeas(HPParams(xml_in, path));
            }
            //! Local registration flag
            bool registered = false;
        }
        const std::string name = "HP_PROPAGATOR";

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

        void read(XMLReader& xml, const std::string& path, HPParams::HPProp_t& par)
        {
            XMLReader paramtop(xml, path);
            read(paramtop, "PropagatorParam" ,par.prop_param ); 
	    read(paramtop, "starting_vector" ,par.starting_vector ); 
	    read(paramtop, "ending_vector" ,par.ending_vector ); 
	    read(paramtop, "Seed" ,par.ran_seed ); 
	    read(paramtop, "ZN" ,par.ZN ); 
        }

        void write(XMLWriter& xml, const std::string& path, HPParams::HPProp_t& par)
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
        void read(XMLReader& xml, const std::string& path, HPParams::NamedObject_t& input)
        {
            XMLReader inputtop(xml, path);
	    read(inputtop, "gauge_id"     , input.gauge_id);
            read(inputtop, "hp_prop_id"   , input.hp_prop_id);
        }

        //! NamedObject output
        void write(XMLWriter& xml, const std::string& path, const HPParams::NamedObject_t& input)
        {
            push(xml, path);
	    write(xml, "gauge_id"     , input.gauge_id    );
            write(xml, "hp_prop_id"   , input.hp_prop_id);
            pop(xml);
        }

        // Param stuff
        HPParams::HPParams()
        {
            frequency = 0;
        }

        HPParams::HPParams(XMLReader& xml_in, const std::string& path)
        {
            try
            {
                XMLReader paramtop(xml_in, path);
                if (paramtop.count("Frequency") == 1)
                    read(paramtop, "Frequency", frequency);
                else
                    frequency = 1;

                // Parameters for source construction
                read(paramtop, "HPParams", hpparam);

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

        void HPParams::writeXML(XMLWriter& xml_out, const std::string& path)
        {
            push(xml_out, path);
            write(xml_out, "HPParams", hpparam);
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
            QDPIO::cout << "HP_PROPAGATOR: start" << std::endl;

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
		QDPIO::cerr << LalibeHPPropagatorEnv::name
		    << ": caught dynamic cast error" << std::endl;
		QDP_abort(1);
	    }
	    catch (const std::string& e)
	    {
		QDPIO::cerr << LalibeHPPropagatorEnv::name
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

	    std::istringstream xml_action(params.hpparam.prop_param.fermact.xml);
	    XMLReader action_reader(xml_action);
	    Handle<FermionAction<T, P, Q>> action(TheFermionActionFactory::Instance().createObject(params.hpparam.prop_param.fermact.id, action_reader, params.hpparam.prop_param.fermact.path));
	    Handle<FermState<T, P, Q>> action_state(action->createState(u));
	    Handle<SystemSolver<LatticeFermion>> solver = action->qprop(action_state, params.hpparam.prop_param.invParam);
	    QDPIO::cout<<"Our action and fermion state are doing A-okay so far."<<std::endl;
	    //We are going to manually do Nc*Ns inversions, looping over spin/color and creating and solved noise propagator which doesn't require header stuff.

	    //Here's where the actual noisy stuff happens.
	    LatticeComplex noise_vec ;
	    Seed ran_seed;
	    QDP::RNG::savern(ran_seed);

	    // Set the seed to desired value
	    QDP::RNG::setrn(params.hpparam.ran_seed);

	    int ZN  = params.hpparam.ZN;
	    LatticeReal rnd1, theta;
	    // twopi defined in chroma/lib/chromabase.h
	    Real twopiN = Chroma::twopi / ZN; 
	    random(rnd1); 
	    theta = twopiN * floor(ZN*rnd1);
	    noise_vec = cmplx(cos(theta),sin(theta));
	    
	    //restore the seed
	    QDP::RNG::setrn(ran_seed);
	    
	    //Here is where the HP stuff is created; adopted from Andrea's main function. 
	    struct meshVars mesh;
	    unsigned int i, N, Hpsize, sample;
	    unsigned int *perm, *Hperm;
	    //double *RHS;

	    mesh.d = Nd;
	    mesh.ptsPerDim = (unsigned int *)malloc(sizeof(unsigned int)*mesh.d);
	    N=1;
	    for (i=0;i<mesh.d; i++){
	      //printf("dim(%u): ",i);
	      //scanf("%u", &mesh.ptsPerDim[i]);
	      //This logic is taken from Bit Twiddling Hacks by Sean Anderson, September 5, 2010
	      unsigned int v = Layout::lattSize()[i]; // compute the next highest power of 2 of 32-bit v
              v--;
              v |= v >> 1;
              v |= v >> 2;
              v |= v >> 4;
              v |= v >> 8;
              v |= v >> 16;
              v++;
	      //mesh.ptsPerDim[i] = Layout::lattSize()[i];
	      //Extend HP to beyond purely powers of two.
	      mesh.ptsPerDim[i] = v;
	      N *= mesh.ptsPerDim[i];
	    }
	    
	    /* Set up the hierarchical data structs */
	    hierOrderSetUp(&mesh);

	    /* Find the row permutation once for each point in local mesh*/
	    /* Here it's done for all nodes N. But can be done only for local ones */
	    perm = (unsigned int *)malloc(sizeof(unsigned int)*N);
	    hierPerm(&mesh, perm, N);
	    
	    /* Now with the perm obtained we do not need the mesh any more */
	    freeMeshVars(&mesh);

	    /* Create the column permutation of the Hadamard vectors. Do it once */
	    /* Create only for as many columns as you need. For example 1024 or 4096 */
	    /* CAUTION: N is the global size (total spatial dimension of the mesh) */
	    /* CAUTION: N must be a power of 2 */
	    //Hpsize = 1024;
	    Hpsize = params.hpparam.ending_vector;
	    Hperm = (unsigned int *)malloc(sizeof(unsigned int)*Hpsize);
	    hadaColPerm(N, Hperm, Hpsize);
	    /**********************************************************************/
	    /* At this point we are ready to create the permuted Hadamard vectors */
	    /* What is needed is perm and Hperm */
	    /**********************************************************************/
	    
	    //Object to hold HP vectors we want to invert
	    multi1d <LatticeInteger> vectors;
	    vectors.resize(params.hpparam.ending_vector - params.hpparam.starting_vector + 1);
	    
	    //Back to code adopted from Andrea's main function.
	    //for (sample=0; sample<Hpsize; sample++) {
	       ///* I create all N here, but it can be done on local rows only */
	       //for (i=0;i<N;i++) 
		  //RHS[i] = (double) Hada_element(perm[i], Hperm[sample]);

	       ///* use RHS in trace computation */
	    //}
	   
	    int Nx = Layout::lattSize()[0];
	    int Ny = Layout::lattSize()[1];
	    int Nz = Layout::lattSize()[2];
	    int Nt = Layout::lattSize()[3];
	    for(int HP_index = params.hpparam.starting_vector; HP_index < params.hpparam.ending_vector + 1; HP_index++)
	    {
	      QDPIO::cout<<"Building HP vector number "<<HP_index<<std::endl;
	      for(int x = 0; x < Nx; x++)
		for(int y = 0; y < Ny; y++)
		  for(int z = 0; z < Nz; z++)
		    for(int t = 0; t < Nt; t++)
		    {
		      multi1d<int> chroma_coords;
		      chroma_coords.resize(Nd);
		      chroma_coords[0] = x; chroma_coords[1] = y; chroma_coords[2] = z; chroma_coords[3] = t;
		      int i = Layout::linearSiteIndex(chroma_coords);
		      int element = Hada_element(perm[i], (Hperm[HP_index - 1]));
		      Integer chroma_element = element;
		      pokeSite(vectors[HP_index - params.hpparam.starting_vector], chroma_element, chroma_coords);
		    }
	    }

	    //For debugging HP vectors.
	    ComplexD Trace = 0.0;
	    for(int vec_index = 0; vec_index < vectors.size(); vec_index++)
	    {
	       QDPIO::cout<<"Inverting hierarchnical probing vector number "<<vec_index+1<<std::endl;
	       LatticePropagator noise_prop = zero;
	       //Now let's do some dilution...
	       for(int color_source(0);color_source<Nc;color_source++){
		 QDPIO::cout << "color_source = " << color_source << std::endl; 
		 
		 LatticeColorVector vec_srce = zero ;
		 LatticeComplex HP_vec = noise_vec*vectors[vec_index];
		 pokeColor(vec_srce,HP_vec,color_source) ;
		 
		 for(int spin_source=0; spin_source < Ns; ++spin_source){
		   QDPIO::cout << "spin_source = " << spin_source << std::endl; 
		   
		   // Insert a ColorVector into spin index spin_source
		   // This only overwrites sections, so need to initialize first
		   // Also, do the solve here
		   LatticeFermion chi = zero;
		   LatticeFermion noise_soln = zero;
		   CvToFerm(vec_srce, chi, spin_source);
		   SystemSolverResults_t res = (*solver)(noise_soln, chi);
		   Trace += innerProduct(noise_soln, chi);
		   FermToProp(noise_soln, noise_prop, color_source, spin_source); 
		 }
	       }
	      
	      QDPIO::cout<<"HP Trace "<<vec_index+1<<" - "<<(Trace/(vec_index+1))<<std::endl;
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
  	      std::string current_id = params.named_obj.hp_prop_id[vec_index];
  	      TheNamedObjMap::Instance().create<LatticePropagator>(current_id);
  	      TheNamedObjMap::Instance().getData<LatticePropagator>(current_id) = noise_prop;
	      TheNamedObjMap::Instance().get(current_id).setFileXML(file_xml);
	      TheNamedObjMap::Instance().get(current_id).setRecordXML(record_xml);
  	      QDPIO::cout<<"Passed noise vector: "<<current_id<<"to the Named Object Buffer."<<std::endl;
	    }
	    snoop.stop();
	    QDPIO::cout << LalibeHPPropagatorEnv::name << ": total time = " << snoop.getTimeInSeconds() << " secs" << std::endl;
	    QDPIO::cout << LalibeHPPropagatorEnv::name<< ": ran successfully" << std::endl;
	    END_CODE();

	}
    }// LalibeHPPropagatorEnv
  };
