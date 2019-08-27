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
#include "util/ferm/paulitodr.h"

// Lalibe Stuff
#include "../momentum/lalibe_sftmom.h"
#include "NN_LC_prop_w.h"
#include "../matrix_elements/bilinear_gamma.h"
//TODO: bilinear_gamma above shouldn't be needed once boiler plate is removed.
#include "../contractions/baryon_contractions_func_w.h"

// Latscat Stuff
#ifndef CUFFT
#include "../NN/fourier_cpu.h"
#else
#include "../NN/fourier_cuda.h"
#endif
#include "../NN/NN_LC_w.h"
#include "../NN/spinstuff.h"
#include "../NN/transform.h"
#include "../NN/checkpointstuff.h"
#include "../NN/momstuff.h"

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
      //par.sink_xml = readXMLGroup(xml_in, "/sink", "combo"); //Thorsten's sink smearing
      //Since sink_xml is a multi1d of XMLGroups, we need to array instead.
      par.sink_xml = readXMLArrayGroup(xml_in, "/sink", "sink_type"); //Thorsten's sink smearing
      QDPIO::cout<<par.sink_xml.size()<<" sinks successfully read."<<std::endl;

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
      read(paramtop, "boosts", par.boosts); //boosts
      read(paramtop, "output_filename", par.output_filename); //output file
      if (paramtop.count("output_stripesize") != 0)
        read(paramtop, "output_stripesize", par.output_stripesize); //output stripesize; default recommended
      else
        par.output_stripesize = -1;
      if (paramtop.count("is_dirac_basis") != 0)
        read(paramtop, "is_dirac_basis", par.dirac_basis); //specifies props in dirac basis, this is false by default
      else
        par.dirac_basis = false;
    }

    void write(XMLWriter& xml, const std::string& path, NNLCPropParams::NNLCProp_t& par)
    {
      push(xml, path);
      write(xml, "contractions_filename", par.contractions_filename); //hdf5 file containing contractions
      write(xml, "contractions_n_sq", par.contractions_n_sq);         //FIXME Comment needed here as well.
      write(xml, "fft_chunksize", par.fft_chunksize);                 //originally the only parameter in FFTPar struct
      write(xml, "fft_tune", par.fft_chunksize);                      //tune the fft?
      write(xml, "boosts", par.boosts);                               //boosts
      write(xml, "output_filename", par.output_filename);             //output filename
      write(xml, "output_stripesize", par.output_stripesize);         //output stripesize; default recommended
      write(xml, "is_dirac_basis", par.dirac_basis);                     //specifies props in dirac basis, this is false by default
      //Sink stuff is a GroupXML, we are not writing it.
      //write(xml, "sink_type", par.sink_xml);                          //sink info
      pop(xml);
    }

    //! NamedObject input
    void read(XMLReader& xml, const std::string& path, NNLCPropParams::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);
      read(inputtop, "gauge_id"     , input.gauge_id);
      read(inputtop, "prop0_id"     , input.prop0_id);
      read(inputtop, "prop1_id"     , input.prop1_id);
    }

    //! NamedObject output
    void write(XMLWriter& xml, const std::string& path, const NNLCPropParams::NamedObject_t& input)
    {
      push(xml, path);
      write(xml, "gauge_id"     , input.gauge_id    );
      write(xml, "prop0_id"     , input.prop0_id);
      write(xml, "prop1_id"     , input.prop1_id);
      pop(xml);
    }

    //This is a function from latscat that should just be local to this translation unit. 
    //I'll just keep it here.
    std::string boost_string(const multi1d<int> &boost){
        std::string boostdir("boost_");
        for(unsigned int d=0; d<(Nd-1); d++){
            if(boost[d]>=0) boostdir+="p";
            else boostdir+="m";
            boostdir+=dirlist[d]+std::to_string(abs(boost[d]));
        }
        if(boost.size() > (Nd-1) && boost[Nd-1] !=0){
            if(boost[Nd-1]>=0) boostdir+="p";
            else boostdir+="m";
            boostdir+=dirlist[Nd-1]+std::to_string(abs(boost[Nd-1]));
        }
        return boostdir;
    }
    
    //Similar category for this function. This comes from utils.cc, I don't want to take all the other stuff from there.
    //This is an isolated copy.
    const std::string getTimestamp() {
      time_t     now = time(0);
      struct tm  tstruct;
      char       buf[80];
      tstruct = *localtime(&now);
      strftime(buf, sizeof(buf), "%Y-%m-%d-%X", &tstruct);
      std::string result(buf);
      result.erase(std::remove(result.begin(), result.end(), ':'), result.end());
    
      return result;
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

      //For now, this measurement is only going to crunch stuff if we have built with hdf5.
#ifdef BUILD_HDF5
      
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

      multi1d<int> pos0(Nd), pos1(Nd), disp(Nd);
      //latscat reads a bunch of prop stuff from XML after this, I am going to directly extract it instead.
      //Read/set up prop 0.
      XMLReader prop0_file_xml, prop0_record_xml;
      LatticePropagator uprop_p1;
      //From earlier latscat code I have set j_decay to Nd -1, so not reading it here.
      QDPIO::cout << "Attempt to read propagator 0" << std::endl;
      try
      {
        uprop_p1 = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop0_id);
        TheNamedObjMap::Instance().get(params.named_obj.prop0_id).getFileXML(prop0_file_xml);
        TheNamedObjMap::Instance().get(params.named_obj.prop0_id).getRecordXML(prop0_record_xml);
        //This all assumes the incoming propagator is coming from a makesource, otherwise we are in a ton of trouble ~_~.
        MakeSourceProp_t  orig_header;
        read(prop0_record_xml, "/Propagator", orig_header);
        pos0 = orig_header.source_header.getTSrce();
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
      //Now read/set up prop 1.
      XMLReader prop1_file_xml, prop1_record_xml;
      LatticePropagator uprop_p2;
      //From earlier latscat code I have set j_decay to Nd -1, so not reading it here.
      QDPIO::cout << "Attempt to read propagator 0" << std::endl;
      try
      {
        uprop_p2 = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop1_id);
        TheNamedObjMap::Instance().get(params.named_obj.prop1_id).getFileXML(prop1_file_xml);
        TheNamedObjMap::Instance().get(params.named_obj.prop1_id).getRecordXML(prop1_record_xml);
        //This all assumes the incoming propagator is coming from a makesource, otherwise we are in a ton of trouble ~_~.
        MakeSourceProp_t  orig_header;
        read(prop1_record_xml, "/Propagator", orig_header);
        pos1 = orig_header.source_header.getTSrce();
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
      //Now back to latscat logic of calculating displacements and stuff.
      QDPIO::cout << "    Propagator 0: " << params.named_obj.prop0_id << std::endl;
      QDPIO::cout << "        Location: "; for(int i = 0; i<Nd; i++){ QDPIO::cout << pos0[i] << " " ;}; QDPIO::cout << std::endl;
      QDPIO::cout << "    Propagator 1: " << params.named_obj.prop1_id << std::endl;
      QDPIO::cout << "        Location: "; for(int i = 0; i<Nd; i++){ QDPIO::cout << pos1[i] << " " ;}; QDPIO::cout << std::endl;

      disp = pos1 - pos0;
      QDPIO::cout << "    Displacement: "; for(int i = 0; i<Nd; i++){ QDPIO::cout << disp[i] << " " ;}; QDPIO::cout << std::endl;
      for(unsigned int d=0; d<Nd; d++){
          QDPIO::cout << disp[d] << "%" << Layout::lattSize()[d] ;
          disp[d] = (disp[d] + Layout::lattSize()[d]) % Layout::lattSize()[d];
          QDPIO::cout << " --> " << disp[d] ;
          if(disp[d] > Layout::lattSize()[d]/2){ disp[d] = disp[d] - Layout::lattSize()[d]; }
          QDPIO::cout << " --> " << disp[d] << std::endl;
      }
      QDPIO::cout << "    Real Displacement: "; for(int i = 0; i<Nd; i++){ QDPIO::cout << disp[i] << " " ;}; QDPIO::cout << std::endl;
      
    // //displacements
    std::string displacedir;
    for(unsigned int d=0; d<Nd; d++){
        if(disp[d] >= 0) displacedir+="p";
        else displacedir+="m";
        displacedir+=dirlist[d]+std::to_string(abs(disp[d]));

    }
    
    // TODO: it would be better to make these more object-oriented.
    //relevant spin projectors, in DP rep
    //These are comments from latscat.
    std::map<std::string,HalfSpinMatrix> projectors;
    projectors["SING0"]=getHalfProjector("SING0");
    projectors["TRIPP"]=getHalfProjector("TRIPP1");
    projectors["TRIP0"]=getHalfProjector("TRIP0");
    projectors["TRIPM"]=getHalfProjector("TRIPM1");

    //In this (lalibe) version, we have already grabbed the gauge field from the named object store.
    //We are just going to do a single set of sink refreshes with this gauge field.
    QDPIO::cout << "Refreshing sinks for our gauge field..." << std::endl;
    for(int s=0; s<sinks.size(); s++){
        if( "GAUGE_INV_GAUSS" == sinks[s]->type()){
            dynamic_cast<gauge_gauss_sink*>(sinks[s])->set_gauge(u);
            QDPIO::cout << "    Sink " << s << " refreshed." << std::endl;
        }
    }
    QDPIO::cout << "... sinks refreshed." << std::endl;

    //Checkpointing isn't suppose to be supported, I include the bare minimum here...
    checkpoint chk(params.nnlcparam.output_filename+".NN_w.chk",params.nnlcparam.output_stripesize);
    
    QDPIO::cout << "Creating timers..." << std::flush;
    // TODO: make timings more meaningful?
    //timings:
    StopWatch swatch_pdir, swatch_ploc, swatch_pdisp, swatch_bblock, swatch_io_read, swatch_io_write;
    swatch_io_read.reset();
    swatch_io_write.reset();
    swatch_pdir.reset();
    swatch_ploc.reset();
    swatch_pdisp.reset();
    swatch_bblock.reset();
    QDPIO::cout << "done!" << std::endl;


    multi1d<int> protonpos1(Nd), protonpos2(Nd);
    protonpos1 = pos0;  // I know.
    protonpos2 = pos1;  // I'm sorry.
    
    //set up map
    //Timeshiftmap tshiftmap(protonpos1[latpars.tDir],latpars.tDir,latpars.tLength);
    //I do the appropriate lalibe substitutions here.
    Timeshiftmap tshiftmap(protonpos1[j_decay],j_decay,Layout::lattSize()[j_decay]);
    QDPIO::cout << "Skipping source setup --- inversions already accomplished." << std::endl;
   
    //If the propagators aren't specified to already be in the Dirac basis, we rotate them now.
    //This is done the lalibe way, instead of passing a spin_basis string, like in latscat.
    if(params.nnlcparam.dirac_basis == false)
    {
      QDPIO::cout << "Rotating the propagators to Dirac basis." << std::endl;
      rotate_to_Dirac_Basis(uprop_p1);
      rotate_to_Dirac_Basis(uprop_p2);
    }
    
    //Necessary Fields:
    multi1d<BaryOp> Nup=get_local_MA_single(0,"DP");
    LatticeComplex tmplatcomp_P, tmplatcomp_P_34, token;
    std::map<std::string,LatticeHalfSpinMatrix> tmplatmats, tmplatmats_34;
    for(unsigned int s=0; s<12; s++){
        tmplatmats[contterms[s]]=LatticeHalfSpinMatrix();
        tmplatmats_34[contterms[s]]=LatticeHalfSpinMatrix();
    }
    
    // Get g5 in Dirac-Pauli (chiral) basis from Degrand-Rossi basis.
    SpinMatrixD gamma5= adj(PauliToDRMat()) * Gamma(15) * PauliToDRMat();
    
    multi1d<LatticePropagator> prop_0(sinks.size());
    multi1d<LatticePropagator> prop_0_34(sinks.size());
    multi1d<LatticePropagator> prop_1(sinks.size());
    multi1d<LatticePropagator> prop_1_34(sinks.size());
    
    QDPIO::cout << "Sinking propagator 0:" << std::endl;
    for(unsigned int i=0; i<sinks.size(); i++){
        prop_0[i] = sinks[i]->operator()(uprop_p1);
        prop_0_34[i] = gamma5 * prop_0[i] * gamma5;
    }

    QDPIO::cout << "Sinking propagator 1:" << std::endl;
    for(unsigned int i=0; i<sinks.size(); i++){
        prop_1[i] = sinks[i]->operator()(uprop_p2);
        prop_1_34[i] = gamma5 * prop_1[i] * gamma5;
    }

    /********************************************************
    *                                                       *
    *   LOOP OVER BOOSTS                                    *
    *                                                       *
    ********************************************************/
    for(unsigned int b=0; b < params.nnlcparam.boosts.size(); b++){
        multi1d<int> boost = params.nnlcparam.boosts[b];
        std::string boostdir=boost_string(boost);

        momentum mom(params.nnlcparam.boosts[b]);

        QDPIO::cout << mom << "    " << boostdir << std::endl;

        chk.create_directory(boostdir);
        chk.close();

        LatticeComplex phases=get_phases(mom,j_decay,+1);

        swatch_pdisp.start();
        QDPIO::cout << "Computing contractions..." << std::endl;
        // if displacement between sources of uprop_p{1,2} is nonzero:
        QDPIO::cout << "    Doing displaced only." << std::endl;
        QDPIO::cout << "    For now, if you want local, just pass the same propagator twice." << std::endl;
        contract(tmplatcomp_P, tmplatmats, prop_0, prop_1, Nup[0].get_gamma(0), phases, fftblock, false, weights);
        QDPIO::cout << "    Contractions done!" << std::endl;

        QDPIO::cout << "Computing contractions with negative parity blocks ..." << std::endl;
        contract(tmplatcomp_P_34, tmplatmats_34, prop_0_34, prop_1_34, Nup[0].get_gamma(0), phases, fftblock, false, weights);
        QDPIO::cout << "    done!" << std::endl;


        QDPIO::cout << "Writing!" << std::endl;

        chk.open();
        chk.set_consistency(false);

        //QDPIO::cout << "Skipping the single proton." << std::endl;
        /*******************************************
         * SINGLE PROTON
         *******************************************/
        if(b==0){
            QDPIO::cout << "Single proton." << std::endl;
            //positive parity proton
            token=zero;
            token+=/* sourcefact * tmpprefact * */ tshiftmap(tmplatcomp_P,true);
            swatch_io_write.start();
            chk.set_parameter("proton1",token);
            swatch_io_write.stop();

            //negative parity proton
            token=zero;
            token+=/* sourcefact * tmpprefact * */ tshiftmap(tmplatcomp_P_34,true);
            swatch_io_write.start();
            chk.set_parameter("proton1_34",token);
            swatch_io_write.stop();
        }
    
        //TODO: This print comes from latscat, but it also serves as a placeholder to extend this measurement and make it
        //more general in a second or third pass.
        QDPIO::cout << "Skipping the local sources." << std::endl;
        
        for(std::map<std::string,LatticeHalfSpinMatrix>::iterator it=tmplatmats.begin(); it!=tmplatmats.end(); ++it){
            std::string idstring=it->first;
            if(idstring.find("loc")!=std::string::npos) continue;

            size_t firstpos=idstring.find("_");
            size_t secondpos=idstring.find(firstpos);
            std::string conttype=idstring.substr(0,firstpos);
            std::transform(conttype.begin(),conttype.end(),conttype.begin(),::tolower);

            std::string srcspin=idstring.substr(firstpos+1,secondpos);
            std::string srcspinval(&srcspin[srcspin.size()-1]);
            srcspin=srcspin.substr(0,srcspin.size()-1);

            for(std::map<std::string,HalfSpinMatrix>::iterator innerit=projectors.begin(); innerit!=projectors.end(); ++innerit){
                std::string inneridstring=innerit->first;
                std::string snkspin=inneridstring.substr(0,inneridstring.size()-1);
                if(snkspin!=srcspin) continue;
                std::string snkspinval(&inneridstring[inneridstring.size()-1]);

                std::string corrname=conttype+"corr"+"_"+srcspin+"_"+snkspinval+"_"+srcspinval+"_"+displacedir; //s[0]; // 0 = mu hardcoded for this version.
                token=zero;
                // TODO src always is 0 for this version!!
                // if(src==0) token=zero;
                // else{
                //     swatch_io_read.start();
                //     chk.get_parameter(boostdir+"/"+corrname,token);
                //     swatch_io_read.stop();
                // }
                token+=/* sourcefact * tmpprefact * */ tshiftmap(LatticeComplex(trace((innerit->second)*(it->second))));
                swatch_io_write.start();
                chk.set_parameter(boostdir+"/"+corrname,token);
                swatch_io_write.stop();
            }
        }
        
        // other-parity (34 entries)
        for(std::map<std::string,LatticeHalfSpinMatrix>::iterator it=tmplatmats_34.begin(); it!=tmplatmats_34.end(); ++it){
            std::string idstring=it->first;
            if(idstring.find("loc")!=std::string::npos) continue;

            size_t firstpos=idstring.find("_");
            size_t secondpos=idstring.find(firstpos);
            std::string conttype=idstring.substr(0,firstpos);
            std::transform(conttype.begin(),conttype.end(),conttype.begin(),::tolower);

            std::string srcspin=idstring.substr(firstpos+1,secondpos);
            std::string srcspinval(&srcspin[srcspin.size()-1]);
            srcspin=srcspin.substr(0,srcspin.size()-1);

            for(std::map<std::string,HalfSpinMatrix>::iterator innerit=projectors.begin(); innerit!=projectors.end(); ++innerit){
                std::string inneridstring=innerit->first;
                std::string snkspin=inneridstring.substr(0,inneridstring.size()-1);
                if(snkspin!=srcspin) continue;
                std::string snkspinval(&inneridstring[inneridstring.size()-1]);

                std::string corrname=conttype+"corr"+"_"+srcspin+"_"+snkspinval+"_"+srcspinval+"_"+displacedir; //s[0]; // 0 = mu hardcoded for this version.
                token=zero;
                // TODO src always is 0 for this version!!
                // if(src==0) token=zero;
                // else{
                //     swatch_io_read.start();
                //     chk.get_parameter(boostdir+"/"+corrname,token);
                //     swatch_io_read.stop();
                // }
                token+=/* sourcefact * tmpprefact * */ tshiftmap(LatticeComplex(trace((innerit->second)*(it->second))));
                swatch_io_write.start();
                chk.set_parameter(boostdir+"/"+corrname+"_34",token);
                swatch_io_write.stop();
            }
        }
        
        //state
        chk.set_parameter("mucurrent",static_cast<unsigned int>(0+1)); // 0 = mu hardcoded for this version.
        // if((mu+1)<sourcepars.displacements.nrows()) chk.set_consistency(true);
        chk.set_consistency(true);  // can hardcode THIS too, because mu is always 0!
        chk.close();
        swatch_pdisp.stop();

/*#ifndef NO_SIGNAL_HANDLERS
        sighand::exit_code_on_abort();
#endif*/
        //This signal handling stuff is not necessary. I don't want the baggage of Thorsten's specialized XML readexs.

        QDPIO::cout << "NN-corr-disp: time=" << swatch_pdisp.getTimeInSeconds() << std::endl;
        swatch_pdisp.reset();
        // if(src!=0) QDPIO::cout << "NN-corr-io-read: time=" << swatch_io_read.getTimeInSeconds() << std::endl;
        QDPIO::cout << "NN-corr-io-write: time=" << swatch_io_write.getTimeInSeconds() << std::endl;
        swatch_io_read.reset();
        swatch_io_write.reset();
    
    } //This is manually here to close loop while porting.
    /********************************************************
    *                                                       *
    *   END LOOP OVER BOOSTS                                *
    *                                                       *
    ********************************************************/

    //move the checkpoint file:
    std::string timestamp=getTimestamp();
    rename(std::string(params.nnlcparam.output_filename+".NN_w.chk").c_str(),std::string(params.nnlcparam.output_filename).c_str());

    //This is where the end of the loop over configurations went in latscat, obviously don't need that here.
    
    //clear baryon blocks:
    clearTopologies();
    
    //This is where latscat's swatch_everything stops. I am going to use the timer that's printed at the end.

#else
      QDPIO::cout << "This measurement only works if we have linked against FFTW. Please rebuild." << std::endl;
#endif

#else
      QDPIO::cout << "This measurement only works if we have enabled HDF5. Please rebuild." << std::endl;
#endif
     
      snoop.stop();
      QDPIO::cout << LalibeNucleonNucleonLinearComboPropagatorEnv::name << ": total time = " << snoop.getTimeInSeconds() << " secs" << std::endl;
      QDPIO::cout << LalibeNucleonNucleonLinearComboPropagatorEnv::name<< ": ran successfully" << std::endl;
      END_CODE();
    
    }
  }// LalibeNucleonNucleonLinearComboPropagatorEnv
};
