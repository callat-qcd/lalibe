/*
Authors
Arjun Gambhir

INPUT
    Propagator
    List of Currents (spin, space, color, momentum)
    Parameters for linear solver
    Information about scalar random noise
    (should match with ZN inverter)
OUTPUT
    HP FH Propagator for each of the specified currents
*/

// Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
#include "../momentum/lalibe_sftmom.h"
#include "fermact.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "util/info/unique_id.h"
//We need this to insert Fermions to a Prop, this is imperitive for dilution.
#include "util/ferm/transf.h"

// Lalibe Stuff
#include "HP_fh_prop_w.h"
#include "../matrix_elements/bilinear_gamma.h"
#include "numerics/binaryRecursiveColoring.h"

namespace Chroma
{
    namespace LalibeHPFHPropagatorEnv
    {
        namespace
        {
            AbsInlineMeasurement* createMeasurement(XMLReader& xml_in,
                const std::string& path)
            {
                return new InlineMeas(HPFHParams(xml_in, path));
            }
            //! Local registration flag
            bool registered = false;
        }
        const std::string name = "HP_FH_PROPAGATOR";

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

        void read(XMLReader& xml, const std::string& path, HPFHParams::HPFHProp_t& par)
        {
            XMLReader paramtop(xml, path);
            read(paramtop, "currents" ,par.currents  ); //list of currents
	    //read(paramtop, "vector_number" ,par.vector_number );
	    read(paramtop, "starting_vector" ,par.starting_vector );
	    read(paramtop, "ending_vector" ,par.ending_vector );
	    read(paramtop, "Seed" ,par.ran_seed );
	    read(paramtop, "ZN" ,par.ZN );
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
	    if (paramtop.count("delete_props") != 0)
	    {
		read(paramtop, "delete_props", par.delete_props); //
		QDPIO::cout<<"User specified to turn delete_props to "<<par.delete_props<<std::endl;
	    }
	    else
	    {
		par.delete_props = false;
		QDPIO::cout<<"By default, no props will be deleted after sum "<<std::endl;
	    }
        }

        void write(XMLWriter& xml, const std::string& path, HPFHParams::HPFHProp_t& par)
        {
            push(xml, path);
            write(xml, "currents" ,par.currents); //list of currents
	    //write(xml, "vector_number" ,par.vector_number);
	    write(xml, "starting_vector" ,par.starting_vector);
	    write(xml, "ending_vector" ,par.ending_vector);
	    write(xml, "Seed" ,par.ran_seed);
	    write(xml, "ZN" ,par.ZN);
	    if(par.is_mom_max == true)
	      write(xml, "p2_max" ,par.p2_max);
	    else
	      write(xml, "mom_list" ,par.mom_list);
	    write(xml, "delete_props" ,par.delete_props);
            pop(xml);

        }

        //! NamedObject input
        void read(XMLReader& xml, const std::string& path, HPFHParams::NamedObject_t& input)
        {
            XMLReader inputtop(xml, path);
            read(inputtop, "gauge_id"     , input.gauge_id);
            read(inputtop, "src_prop_id"  , input.src_prop_id);
	    read(inputtop, "hp_prop_id"  , input.hp_prop_id);
            read(inputtop, "fh_prop_id"   , input.fh_prop_id);
        }

        //! NamedObject output
        void write(XMLWriter& xml, const std::string& path, const HPFHParams::NamedObject_t& input)
        {
            push(xml, path);
            write(xml, "gauge_id"     , input.gauge_id    );
            write(xml, "src_prop_id"  , input.src_prop_id     );
	    write(xml, "hp_prop_id"  , input.hp_prop_id     );
            write(xml, "fh_prop_id"   , input.fh_prop_id);
            pop(xml);
        }

        // Param stuff
        HPFHParams::HPFHParams()
        {
            frequency = 0;
        }

        HPFHParams::HPFHParams(XMLReader& xml_in, const std::string& path)
        {
            try
            {
                XMLReader paramtop(xml_in, path);
                if (paramtop.count("Frequency") == 1)
                    read(paramtop, "Frequency", frequency);
                else
                    frequency = 1;

                // Parameters for source construction
                read(paramtop, "HPFHParams", hpfhparam);

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

        void HPFHParams::writeXML(XMLWriter& xml_out, const std::string& path)
        {
            push(xml_out, path);
            write(xml_out, "HPFHParams", hpfhparam);
            write(xml_out, "NamedObject", named_obj);
            pop(xml_out);
        }

        void  InlineMeas::operator()(unsigned long update_no, XMLWriter& xml_out)
        {
            START_CODE();

            StopWatch snoop;
            snoop.reset();
            snoop.start();
            QDPIO::cout << "HP_FH_PROPAGATOR: start" << std::endl;

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
                QDPIO::cerr << LalibeHPFHPropagatorEnv::name
                    << ": caught dynamic cast error" << std::endl;
                QDP_abort(1);
            }
            catch (const std::string& e)
            {
                QDPIO::cerr << LalibeHPFHPropagatorEnv::name
                    << ": map call failed: " << e << std::endl;
                QDP_abort(1);
            }
            const multi1d<LatticeColorMatrix>& u =
                TheNamedObjMap::Instance().getData
                    <multi1d <LatticeColorMatrix> >(params.named_obj.gauge_id);

	    //Either do a normal source or sequential source.
	    bool propagatorP = false;
	    bool seqsourceP = false;

	    // Read "src" quark propagator
            XMLReader prop_file_xml, prop_record_xml;
            LatticePropagator quark_propagator;
	    int t0;
	    int j_decay;
	    //Need origin for fourier transform!
	    multi1d<int> origin;
            QDPIO::cout << "Attempt to read forward propagator" << std::endl;
            try
            {
                quark_propagator = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.src_prop_id);
                TheNamedObjMap::Instance().get(params.named_obj.src_prop_id).getFileXML(prop_file_xml);
                TheNamedObjMap::Instance().get(params.named_obj.src_prop_id).getRecordXML(prop_record_xml);
		if (prop_record_xml.count("/Propagator") != 0)
		{
		propagatorP = true;
		MakeSourceProp_t  orig_header;
		read(prop_record_xml, "/Propagator", orig_header);
		j_decay = orig_header.source_header.j_decay;
		t0      = orig_header.source_header.t_source;
		origin = orig_header.source_header.getTSrce();
		}
		else if (prop_record_xml.count("/SequentialSource") != 0)
		{
		  seqsourceP = true;
		  SequentialSource_t   orig_header;
		  read(prop_record_xml, "/SequentialSource", orig_header);

		  j_decay = orig_header.seqsource_header.j_decay;
		  t0      = orig_header.seqsource_header.t_sink;
		  origin = orig_header.forward_props[0].source_header.getTSrce();
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
                QDPIO::cerr << name << ": error reading src prop_header: "
                    << e << std::endl;
                QDP_abort(1);
            }


	    // Initialize FT stuff based on which type of momentum mode was selected.
	    // This assumes j_decay is the orthogonal direction: typically Nd - 1.
	    LalibeSftMom ft = params.hpfhparam.is_mom_max ? LalibeSftMom(params.hpfhparam.p2_max, origin, false, j_decay)
              : LalibeSftMom(params.hpfhparam.p_list, origin, j_decay);
	    // Make an action and all other stuff needed for a solver.

	    LatticePropagator fh_prop_src = quark_propagator;
	    LatticePropagator hp_src = zero;
	    LatticePropagator hp_fh_prop = zero;

	    //Here's where the actual noisy stuff happens.
	    LatticeComplex noise_vec ;
	    Seed ran_seed;
	    QDP::RNG::savern(ran_seed);

	    // Set the seed to desired value
	    QDP::RNG::setrn(params.hpfhparam.ran_seed);

	    int ZN  = params.hpfhparam.ZN;
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
	    Hpsize = params.hpfhparam.ending_vector;
	    Hperm = (unsigned int *)malloc(sizeof(unsigned int)*Hpsize);
	    hadaColPerm(N, Hperm, Hpsize);
	    /**********************************************************************/
	    /* At this point we are ready to create the permuted Hadamard vectors */
	    /* What is needed is perm and Hperm */
	    /**********************************************************************/

	    //Object to hold HP vectors we want to invert
	    multi1d <LatticeInteger> vectors;
	    vectors.resize(params.hpfhparam.ending_vector - params.hpfhparam.starting_vector + 1);

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
	    for(int HP_index = params.hpfhparam.starting_vector; HP_index < params.hpfhparam.ending_vector + 1; HP_index++)
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
		      pokeSite(vectors[HP_index - params.hpfhparam.starting_vector], chroma_element, chroma_coords);
		    }

	    //This is the old way of doing the noise bit.
	    //Now let's do some dilution and turn this noise into a proper quark!
	    /*for(int color_source(0);color_source<Nc;color_source++){
	      QDPIO::cout << "color_source = " << color_source << std::endl;

	      LatticeColorVector vec_srce = zero ;
	      //	      pokeColor(tmp,wave[k],color_source) ;
	      pokeColor(vec_srce,vec,color_source) ;

	      for(int spin_source=0; spin_source < Ns; ++spin_source){
		QDPIO::cout << "spin_source = " << spin_source << std::endl;

		// Insert a ColorVector into spin index spin_source
		// This only overwrites sections, so need to initialize first
		LatticeFermion chi = zero;
		CvToFerm(vec_srce, chi, spin_source);
		FermToProp(chi, noise_src, color_source, spin_source);
	      }
	    }*/

	    //In order to the sum over noise vectors inside here instead of post-processing, we need accumulators.
	    multi1d <LatticePropagator> accumulated_props;
	    accumulated_props.resize(params.named_obj.fh_prop_id.size());
	    //Zero out the accumulators, just to be safe from dirty memory.
	    for(int accumulation_index = 0; accumulation_index < accumulated_props.size(); accumulation_index++)
	      accumulated_props[accumulation_index] = zero;

	    ComplexD Trace = 0.0;
	    //Noise loop
	    for(int vec_index = 0; vec_index < vectors.size(); vec_index++)
	    {
              // Read "hp" quark propagator, now this is done inside a loop over hp vectors
	      XMLReader hp_prop_file_xml, hp_prop_record_xml;
	      LatticePropagator hp_quark_propagator;
	      QDPIO::cout << "Attempt to read HP propagator number " << vec_index << std::endl;
	      try
	      {
		  hp_quark_propagator = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.hp_prop_id[vec_index]);
		  TheNamedObjMap::Instance().get(params.named_obj.hp_prop_id[vec_index]).getFileXML(hp_prop_file_xml);
		  TheNamedObjMap::Instance().get(params.named_obj.hp_prop_id[vec_index]).getRecordXML(hp_prop_record_xml);
	      }
              catch (std::bad_cast)
	      {
		  QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
		  QDP_abort(1);
	      }
	      catch (const std::string& e)
	      {
		  QDPIO::cerr << name << ": error reading hp prop_header: "
		      << e << std::endl;
		  QDP_abort(1);
	      }
	      for(int current_index = 0; current_index < params.hpfhparam.currents.size(); current_index++)
	      {
		std::string present_current = params.hpfhparam.currents[current_index];
		QDPIO::cout << "HP_FH_PROPAGATOR: current " << present_current << std::endl;
		// WE SHOULD MAKE THIS A FACTORY
		Bilinear_Gamma(present_current, fh_prop_src, quark_propagator, u);
		for(int mom = 0; mom < ft.numMom(); mom++)
		{

		  multi1d<int> momenta = ft.numToMom(mom);
		  QDPIO::cout << "Injecting momentum - px: "<<std::to_string(momenta[0])<<" py: "+std::to_string(momenta[1])<<" pz: "+std::to_string(momenta[2])<<std::endl;
		  fh_prop_src = ft[mom]*fh_prop_src;
		  //Below is the thing that matters.
		  //The hp_prop is tied with the source_prop and hp_src to make the hp FH prop.
		  //hp_fh_prop = hp_quark_propagator*innerProduct(hp_src, fh_prop_src);
		  //Instead of doing the above, I will manually specify indices.
		  for(int fh_color_source(0);fh_color_source<Nc;fh_color_source++){
		    //Nested spin/color loop so our dilution vectors hit every propagator component.
		    //QDPIO::cout<<"Extracing Feynman Hellman color "<<fh_color_source<<std::endl;
		    for(int fh_spin_source=0; fh_spin_source < Ns; fh_spin_source++){
		      LatticeFermion fh_ferm = zero;
		      //QDPIO::cout<<"Tying Feynman Hellman spin "<<fh_spin_source<<" to diluted hp."<<std::endl;
		      PropToFerm(fh_prop_src, fh_ferm, fh_color_source, fh_spin_source);
		      LatticeFermion hp_fh_ferm = zero;
		      for(int color_source(0);color_source<Nc;color_source++){
			LatticeColorVector vec_srce = zero ;
			LatticeComplex HP_vec = noise_vec*vectors[vec_index];
			pokeColor(vec_srce,HP_vec,color_source) ;
			for(int spin_source=0; spin_source < Ns; spin_source++){
			  LatticeFermion chi = zero;
			  CvToFerm(vec_srce, chi, spin_source);
			  LatticeFermion hp_ferm = zero;
			  PropToFerm(hp_quark_propagator, hp_ferm, color_source, spin_source);
			  hp_fh_ferm += hp_ferm*innerProduct(chi, fh_ferm);
 			  //Only compute trace once.
 			  if(current_index == 0 && mom == 0 && fh_color_source == 0 && fh_spin_source == 0)
			    Trace += innerProduct(hp_ferm, chi);
			}
		      }
		      FermToProp(hp_fh_ferm, hp_fh_prop, fh_color_source, fh_spin_source);
		    }
		  }
		  //Below we accumulate.
		  std::string current_id = params.named_obj.fh_prop_id[current_index*ft.numMom() + mom];
		  QDPIO::cout<<"Adding outer product to prop with id "<<current_id<<std::endl;
		  accumulated_props[current_index*ft.numMom() + mom] += hp_fh_prop;
		}
	      }
	      //Print debugger
	      QDPIO::cout<<"HP Trace "<<vec_index+1<<" - "<<(Trace/(vec_index+1))<<std::endl;
	      if (params.hpfhparam.delete_props)
	      {
		// Deleting the object.
		TheNamedObjMap::Instance().erase(params.named_obj.hp_prop_id[vec_index]);
	      }
	    }

	    //Divide by the total number of noise_vecs to create an average, then ship out result.
	    for(int accumulation_index = 0; accumulation_index < accumulated_props.size(); accumulation_index++)
	    {
	      accumulated_props[accumulation_index] /= vectors.size();
	      push(xml_out,"Relaxation_Iterations");
	      //ncg_had not needed since no inversion is happening
	      //write(xml_out, "ncg_had", ncg_had);
	      pop(xml_out);

	      QDPIO::cout << "Writing propagator info, cause why not?" << std::endl;
	      XMLBufferWriter file_xml;
	      push(file_xml, "propagator");
	      write(file_xml, "id", uniqueId());  // NOTE: new ID form
	      pop(file_xml);

	      //If ths src is not from make source these thing is not going to work...
	      XMLBufferWriter record_xml;
	      if (propagatorP)
	      {
		//MakeSourceProp_t  orig_header;
		Propagator_t  orig_header;
		read(prop_record_xml, "/Propagator", orig_header);
		//Note below:
		//Normally we copy relevant header info and use new prop_params, but this task doesn't do any inversions...
		//This requies changing types from Propagator to MakeSourceProp.
		//Propagator_t  new_header;   // note, abandoning state_info
		//Therefore, here the src_header is just completely copied over.
		//new_header.prop_header   = orig_header.params.hpfhparam.prop_param;
		//new_header.source_header = orig_header.source_header;
		//new_header.gauge_header  = orig_header.gauge_header;
		//MakeSourceProp_t  new_header = orig_header;
		Propagator_t  new_header = orig_header;
		write(record_xml, "Propagator", new_header);
	      }
	      else if (seqsourceP)
	      {
		SequentialSource_t  orig_header;
		read(prop_record_xml, "/SequentialSource", orig_header);

		SequentialProp_t  new_header;   // note, abandoning state_info
		//Dig into first quark and copy its chroma_prop xml stuff.
		new_header.seqprop_header   = orig_header.forward_props[0].prop_header;
		new_header.sink_header      = orig_header.sink_header;
		new_header.seqsource_header = orig_header.seqsource_header;
		new_header.forward_props    = orig_header.forward_props;
		new_header.gauge_header     = orig_header.gauge_header;
		write(record_xml, "SequentialProp", new_header);
	      }

	      // Pass the propagator info to the Named Object Buffer.
	      // Looping over currents and momenta, momenta is inner most index
	      // current_id = flattened index running over both these indices
	      //std::string current_id = params.named_obj.fh_prop_id[current_index*vectors.size()*ft.numMom() + vec_index*ft.numMom() + mom];
	      //Above is old index ordering, doesn't make much sense now.
	      std::string current_id = params.named_obj.fh_prop_id[accumulation_index];
	      TheNamedObjMap::Instance().create<LatticePropagator>(current_id);
	      TheNamedObjMap::Instance().getData<LatticePropagator>(current_id) = accumulated_props[accumulation_index];
	      TheNamedObjMap::Instance().get(current_id).setFileXML(file_xml);
	      TheNamedObjMap::Instance().get(current_id).setRecordXML(record_xml);
	      QDPIO::cout<<"YAAAY! We finished current: "<<current_id<<std::endl;
	    }


	    snoop.stop();
	    QDPIO::cout << LalibeHPFHPropagatorEnv::name << ": total time = " << snoop.getTimeInSeconds() << " secs" << std::endl;
	    QDPIO::cout << LalibeHPFHPropagatorEnv::name<< ": ran successfully" << std::endl;
	    END_CODE();

	}
    }// LalibeHPFHPropagatorEnv
  };
