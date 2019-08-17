/*
Authors
Arjun Gambhir

First pass at porting the linear combo executable of latscat over to lalibe.
After initial prototype working, lots of optimizations are likely to follow.
Memory management when boosting needs to be looked into.
*/

// Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
#include "fermact.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "util/info/unique_id.h"
#include "io/xml_group_reader.h"

// Lalibe Stuff
#include "../momentum/lalibe_sftmom.h"
#include "NN_LC_prop_w.h"
#include "../matrix_elements/bilinear_gamma.h"
//TODO: bilinear_gamma above shouldn't be needed once boiler plate is removed.

// Latscat Stuff
#ifndef CUFFT
#include "../NN/fourier_cpu.h"
#else
#include "../NN/fourier_cuda.h"
#endif
#include "../NN/NN_LC_w.h"

namespace Chroma
{
  namespace LalibeNucleonNucleonLinearComboPropagatorEnv
  {
    namespace
    {
        AbsInlineMeasurement* createMeasurement(XMLReader& xml_in,
            const std::string& path)
        {
            return new InlineMeas(NNLCPropParams(xml_in, path));
        }
        //! Local registration flag
        bool registered = false;
    }
    const std::string name = "NUCLEON_NUCLEON_LINEAR_COMBO_PROPAGATOR";

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

    void read(XMLReader& xml, const std::string& path, NNLCPropParams::NNLCProp_t& par)
    {
      XMLReader paramtop(xml, path);

      //Thorsten's sink types read here.
      //TODO: Change this to completely use chroma's sink smearing.
      XMLReader xml_tmp(paramtop, "sink");
      
      XMLBufferWriter xml_out;
      push(xml_out,"sink");
      xml_out << xml_tmp;
      pop(xml_out);
      
      XMLReader xml_in(xml_out);
      par.sink_xml = readXMLGroup(xml_in, "/sink", "combo"); //Thorsten's sink smearing

      //FH prop params are below here.
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
      //Where latscat params start.
      read(paramtop, "contractions_filename", par.contractions_filename); //hdf5 file containing contractions
      if (paramtop.count("contractions_n_sq") != 0)
        read(paramtop, "contractions_n_sq", par.contractions_n_sq); //FIXME Description needed here.
      else
        par.contractions_n_sq = -1;
      if (paramtop.count("fft_chunksize") != 0)
        read(paramtop, "fft_chunksize", par.fft_chunksize); //originally the only parameter in FFTPar struct
      else
        par.fft_chunksize = 0;
      if (paramtop.count("fft_tune") != 0)
        read(paramtop, "fft_tune", par.fft_tune); //tune fft?
      else
        par.fft_tune = false;
      if (paramtop.count("boosts") != 0)
        read(paramtop, "boosts", par.boosts); //boosts
      else 
      {
        par.boosts.resize(Nd - 1);
        for (int p = 0; p < (Nd - 1); p++)
          par.boosts[p] = 0; //The default is no boost, ie populate with zeros.
      }
      read(paramtop, "output_filename", par.output_filename); //output file
      if (paramtop.count("output_stripesize") != 0)
        read(paramtop, "output_stripesize", par.output_stripesize); //output stripesize; default recommended
      else
        par.output_stripesize = -1;
      if (paramtop.count("dirac_basis") != 0)
        read(paramtop, "dirac_basis", par.dirac_basis); //specifies props in dirac basis, this is false by default
      else
        par.dirac_basis = false;
    }

    void write(XMLWriter& xml, const std::string& path, NNLCPropParams::NNLCProp_t& par)
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
      //Where latscat params start.
      write(xml, "contractions_filename", par.contractions_filename); //hdf5 file containing contractions
      write(xml, "contractions_n_sq", par.contractions_n_sq);         //FIXME Comment needed here as well.
      write(xml, "fft_chunksize", par.fft_chunksize);                 //originally the only parameter in FFTPar struct
      write(xml, "fft_tune", par.fft_chunksize);                      //tune the fft?
      write(xml, "boosts", par.boosts);                               //boosts
      write(xml, "output_filename", par.output_filename);             //output filename
      write(xml, "output_stripesize", par.output_stripesize);         //output stripesize; default recommended
      write(xml, "dirac_basis", par.dirac_basis);                     //specifies props in dirac basis, this is false by default
      //Sink stuff is a GroupXML, we are not writing it.
      //write(xml, "sink_type", par.sink_xml);                          //sink info
      pop(xml);
    }

    //! NamedObject input
    void read(XMLReader& xml, const std::string& path, NNLCPropParams::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);
      read(inputtop, "gauge_id"     , input.gauge_id);
      read(inputtop, "src_prop_id"  , input.src_prop_id);
      read(inputtop, "fh_prop_id"   , input.fh_prop_id);
    }

    //! NamedObject output
    void write(XMLWriter& xml, const std::string& path, const NNLCPropParams::NamedObject_t& input)
    {
      push(xml, path);
      write(xml, "gauge_id"     , input.gauge_id    );
      write(xml, "src_prop_id"  , input.src_prop_id     );
      write(xml, "fh_prop_id"   , input.fh_prop_id);
      pop(xml);
    }

    // Param stuff
    NNLCPropParams::NNLCPropParams()
    {
        frequency = 0;
    }

    NNLCPropParams::NNLCPropParams(XMLReader& xml_in, const std::string& path)
    {
      try
      {
        XMLReader paramtop(xml_in, path);
        if (paramtop.count("Frequency") == 1)
          read(paramtop, "Frequency", frequency);
        else
          frequency = 1;

        // Parameters for source construction
        read(paramtop, "NNLCPropParams", nnlcparam);

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

    void NNLCPropParams::writeXML(XMLWriter& xml_out, const std::string& path)
    {
      push(xml_out, path);
      write(xml_out, "NNLCPropParams", nnlcparam);
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
      QDPIO::cout << "NUCLEON_NUCLEON_LINEAR_COMBO_PROPAGATOR: start" << std::endl;

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
        QDPIO::cerr << LalibeNucleonNucleonLinearComboPropagatorEnv::name
            << ": caught dynamic cast error" << std::endl;
        QDP_abort(1);
      }
      catch (const std::string& e)
      {
        QDPIO::cerr << LalibeNucleonNucleonLinearComboPropagatorEnv::name
            << ": map call failed: " << e << std::endl;
        QDP_abort(1);
      }
      const multi1d<LatticeColorMatrix>& u =
        TheNamedObjMap::Instance().getData
        <multi1d <LatticeColorMatrix> >(params.named_obj.gauge_id);

      //Code from latscat ported below, with comments for missing or adjusted features.
      //Chroma initialize stuff was here, not needed in lalibe.
      //Latscat has an option for its own singnal handlers that goes here, not being ported.
#ifdef DEBUG
      QDPIO::cout << "Warning, DEBUG mode enabled!" << std::endl;
#endif
      //Stuff about LatticePars is read next, this is already read in our main function.
      //Layout creation/setup is done next, chroma handles this for us.
      
      //set up fourier stuff:   
      QDPIO::cout << "Setting up Communicators for FFT..." << std::flush;
      /*Fourier fft(latpars.tDir);
      Fourier fftblock(latpars.tDir);*/
      //I am going to assume that tDir is the same as j_decay and manually set it for now.
      //How j_decay is determined/set may change later.
      int j_decay = Nd - 1; 
      //This code only works when we link against FFTW. 
      //If there is a need to make FFTW optional for this measurement, this can be changed.
      //After fh prop boiler plate is removed, this directive can be moved up to where latscat code starts.
#ifdef BUILD_FFTW
      Fourier fft(j_decay);
      Fourier fftblock(j_decay);
#ifdef PROFILE
      fft.print_flops(true);
      fftblock.print_flops(true);
#endif
      QDPIO::cout << "done!" << std::endl;
      //ContractionPars contpars;
      //This struct is not needed anymore either.
      if(params.nnlcparam.contractions_n_sq >= 0)
        QDPIO::cout << "Truncating output files to n_sq <= " << params.nnlcparam.contractions_n_sq <<  "!" << std::endl;
      initTopologies(params.nnlcparam.contractions_filename, params.nnlcparam.contractions_n_sq, j_decay); 
      //tDir changed to j_decay again in the above line.
     
      //Number of configs reading was done here, but I am going to avoid that altogether.
      //We will read the config the way all other lalibe measurements do.

      //QDPIO::cout << "Initializing sourcepars to hold the boosts..." << std::flush;
      //SourcePars sourcepars;
      //We aren't actually going to use SourcePars since the only param read in is boost.
      if(params.nnlcparam.fft_chunksize !=0 ) 
        QDPIO::cout << "Using chunksize " << params.nnlcparam.fft_chunksize << " for the Baryon-Block FFT!" << std::endl;

      //Do FFT tuning if it's enabled; comment below comes from latscat.
      // This should probably be put closer to the other fft specification.  Is there a reason it cannot?
      if(params.nnlcparam.fft_tune)
      {
        QDPIO::cout << "Tuning FFT for better performance..." << std::flush;
        StopWatch swatch_fftune;
        swatch_fftune.reset();
        swatch_fftune.start();
        fftblock.tune(sizeof(HalfBaryonblock),true);
        swatch_fftune.stop();
        QDPIO::cout << "done! Time " << swatch_fftune.getTimeInSeconds() << std::endl;
      }

      //Parsing sinks code goes here in the original latscat. 
      //I am going to let the sink construction measurement handle that.
      //TODO: Use chroma's smearing to do this. Already done in our c3pt code, will do this on second pass.
      QDPIO::cout << "Parsing sinks..." << std::endl;
      multi1d<sink*>   sinks(params.nnlcparam.sink_xml.size());
      multi1d<Complex> weights(params.nnlcparam.sink_xml.size());
      for(int s=0; s<params.nnlcparam.sink_xml.size(); s++){
        // QDPIO::cout << "Parsing sink " << s << " ... " << std::flush;
        parseSink(params.nnlcparam.sink_xml[s],sinks[s],weights[s],fft,j_decay);
        //Last param above used to be latpars.tDir.
        // QDPIO::cout << sinks[s]->type() << std::endl;
      }
      QDPIO::cout << "Found " << sinks.size() << " sinks." << std::endl; 

#else
      QDPIO::cout << "This measurement only works if we have linked against FFTW. Please rebuild." << std::endl;
#endif
     
      //Everything beyond this is from the old file...will eventually be deleted.
      // Read "src" quark propagator
      XMLReader prop_file_xml, prop_record_xml;
      LatticePropagator quark_propagator;
      int t0;
      //Because of latscat code, I have commented out the declaration of j_decay.
      //int j_decay;
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
      LalibeSftMom ft = params.nnlcparam.is_mom_max ? LalibeSftMom(params.nnlcparam.p2_max, origin, false, j_decay)
          : LalibeSftMom(params.nnlcparam.p_list, origin, j_decay);
      // Make an action and all other stuff needed for a solver.

      typedef LatticeFermion T;
      typedef multi1d<LatticeColorMatrix> P;
      typedef multi1d<LatticeColorMatrix> Q;


      std::istringstream xml_action(params.nnlcparam.prop_param.fermact.xml);
      XMLReader action_reader(xml_action);
      Handle<FermionAction<T, P, Q>> action(TheFermionActionFactory::Instance().createObject(params.nnlcparam.prop_param.fermact.id, action_reader, params.nnlcparam.prop_param.fermact.path));
      Handle<FermState<T, P, Q>> action_state(action->createState(u));
      QDPIO::cout<<"Our action and fermion state are doing A-okay so far."<<std::endl;
      //Handle<SystemSolver<LatticeFermion>> solver = action->qprop(action_state, params.nnlcparam.prop_param.invParam);
      //Above is for a single fermion, but we want to loop over spin/color and solve for the full propagator.

      int ncg_had = 0; //This appears in the propagator task, I am just copying it here.

      LatticePropagator fh_prop_src = quark_propagator;
      LatticePropagator fh_prop_solution;

      QDPIO::cout << "NUCLEON_NUCLEON_LINEAR_COMBO_PROPAGATOR: N_currents " << params.nnlcparam.currents.size() << std::endl;
      for(int current_index = 0; current_index < params.nnlcparam.currents.size(); current_index++)
      {
        std::string present_current = params.nnlcparam.currents[current_index];
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
              params.nnlcparam.prop_param.invParam,
              params.nnlcparam.prop_param.quarkSpinType,
              params.nnlcparam.prop_param.obsvP, ncg_had);
          
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
          new_header.prop_header   = params.nnlcparam.prop_param;
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
      QDPIO::cout << LalibeNucleonNucleonLinearComboPropagatorEnv::name << ": total time = " << snoop.getTimeInSeconds() << " secs" << std::endl;
      QDPIO::cout << LalibeNucleonNucleonLinearComboPropagatorEnv::name<< ": ran successfully" << std::endl;
      END_CODE();
    
    }
  }// LalibeNucleonNucleonLinearComboPropagatorEnv
};
