/*
 * Single baryon contractions
 * Authors:
 * David Brantley
 * Arjun Gambhir
 * Andre Walker-Loud
 * Jason Chang
 * Do meson contractions and write out the two-point correlator in hdf5
 * Maybe we'll also support sdb, one day...
 */


#include "meson_contractions_w.h"
#include "../contractions/meson_contractions_func_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "../momentum/lalibe_sftmom.h"
#include "meas/inline/io/named_objmap.h"
#include "io/qprop_io.h"
#include "util/spin_basis.h"


namespace Chroma
{
  namespace LalibeMesonContractionsEnv
  {
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in,
					      const std::string& path)
      {
	return new InlineMeas(MesonParams(xml_in, path));
      }

      bool registered = false;
    }
    const std::string name = "MESON_CONTRACTIONS";

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

    void read(XMLReader& xml, const std::string& path, MesonParams::Param_t& par)
    {

      XMLReader paramtop(xml, path);

#ifdef BUILD_HDF5
      read(paramtop, "h5_file_name", par.file_name);
      read(paramtop, "obj_path", par.obj_path);
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

      read(paramtop, "particle_list", par.particle_list);
    }

    void write(XMLWriter& xml, const std::string& path, const MesonParams::Param_t& par)
    {
      push(xml, path);

      write(xml, "rotate_to_Dirac", par.rotate_to_Dirac);
#ifdef BUILD_HDF5
      write(xml, "h5_file_name", par.file_name);
      write(xml, "obj_path", par.obj_path);
#endif
      //write(xml, "output_full_correlator", par.output_full_correlator);
      if(par.is_mom_max == true)
	write(xml, "p2_max" ,par.p2_max);
      else
	write(xml, "mom_list" ,par.mom_list);
      write(xml, "particle_list", par.particle_list);

      pop(xml);
    }

    void read(XMLReader& xml, const std::string& path, MesonParams::NamedObject_t& input)
    {

      XMLReader inputtop(xml, path);

      //read(inputtop, "gauge_id" , input.gauge_id);
      //Logic to read quark propagators (whichever ones are present.)
      if (inputtop.count("up_quark") != 0)
      {
	read(inputtop, "up_quark" ,input.up_quark);
	QDPIO::cout<<"I found an up quark, here is its id: "<<input.up_quark<<std::endl;
	input.is_up = true;
      }
      else
      {
	  QDPIO::cout<<"I couldn't find an up quark, hope you don't need it for the inputted meson contractions. "<<std::endl;
	  input.is_up = false;
      }
      if (inputtop.count("down_quark") != 0)
      {
	read(inputtop, "down_quark" ,input.down_quark);
	QDPIO::cout<<"I found a down quark, here is its id: "<<input.down_quark<<std::endl;
	input.is_down = true;
      }
      else
      {
	  QDPIO::cout<<"I couldn't find a down quark, hope you don't need it for the inputted meson contractions. "<<std::endl;
	  input.is_down = false;
      }
      if (inputtop.count("strange_quark") != 0)
      {
	read(inputtop, "strange_quark" ,input.strange_quark);
	QDPIO::cout<<"I found an strange quark, here is its id: "<<input.strange_quark<<std::endl;
	input.is_strange = true;
      }
      else
      {
	  QDPIO::cout<<"I couldn't find a strange quark, hope you don't need it for the inputted meson contractions. "<<std::endl;
	  input.is_strange = false;
      }
      if (inputtop.count("charm_quark") != 0)
      {
	read(inputtop, "charm_quark" ,input.charm_quark);
	QDPIO::cout<<"I found an charm quark, here is its id: "<<input.charm_quark<<std::endl;
	input.is_charm = true;
      }
      else
      {
	  QDPIO::cout<<"I couldn't find a charm quark, hope you don't need it for the inputted meson contractions. "<<std::endl;
	  input.is_charm = false;
      }
    }

    void write(XMLWriter& xml, const std::string& path, const MesonParams::NamedObject_t& input)
    {
      push(xml, path);
      //write(xml, "gauge_id" , input.gauge_id);
      if(input.is_up == true)
	write(xml, "up_quark" ,input.up_quark);
      if(input.is_down == true)
	write(xml, "down_quark" ,input.down_quark);
      if(input.is_strange == true)
	write(xml, "strange_quark" ,input.strange_quark);
      if(input.is_charm == true)
	write(xml, "charm_quark" ,input.charm_quark);
      pop(xml);
    }


    MesonParams::MesonParams()
    {
      frequency = 0;
    }

    MesonParams::MesonParams(XMLReader& xml_in, const std::string& path)
    {
      try
      {
	XMLReader paramtop(xml_in, path);

	if (paramtop.count("Frequency") == 1)
	  read(paramtop, "Frequency", frequency);
	else
	  frequency = 1;

	read(paramtop, "MesonParams", param);

	read(paramtop, "NamedObject", named_obj);

      }
      catch(const std::string& e)
      {
	QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
	QDP_abort(1);
      }
    }

    void
    MesonParams::writeXML(XMLWriter& xml_out, const std::string& path)
    {
      push(xml_out, path);

      write(xml_out, "MesonParams", param);
      write(xml_out, "NamedObject", named_obj);

      pop(xml_out);
    }


    void  InlineMeas::operator()(unsigned long update_no, XMLWriter& xml_out)
    {
      START_CODE();

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      QDPIO::cout<<"Meson contractions starting..."<<std::endl;

      //Grab all the propagators that are given.
      XMLReader up_prop_file_xml, up_prop_record_xml;
      LatticePropagator up_quark_propagator;
      XMLReader down_prop_file_xml, down_prop_record_xml;
      LatticePropagator down_quark_propagator;
      XMLReader strange_prop_file_xml, strange_prop_record_xml;
      LatticePropagator strange_quark_propagator;
      XMLReader charm_prop_file_xml, charm_prop_record_xml;
      LatticePropagator charm_quark_propagator;
      //Need origin, j_decay, and t0 for fourier transform!
      //Need j_decay of bc to know what comes with a minus sign.
      int j_decay;
      int t_0;
      multi1d<int> origin;

      if(params.named_obj.is_up == true)
      {
	QDPIO::cout << "Attempting to read up propagator" << std::endl;
	try
	{
	    up_quark_propagator = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.up_quark);
	    TheNamedObjMap::Instance().get(params.named_obj.up_quark).getFileXML(up_prop_file_xml);
	    TheNamedObjMap::Instance().get(params.named_obj.up_quark).getRecordXML(up_prop_record_xml);
	    //Get the origin  and j_decay for the FT, this assumes all quarks have the same origin.
	    MakeSourceProp_t  orig_header;
	    if (up_prop_record_xml.count("/Propagator") != 0)
	    {
	      QDPIO::cout<<"Up quark propagator is unsmeared, reading from Propagator tag..."<<std::endl;
	      read(up_prop_record_xml, "/Propagator", orig_header);
	    }
	    else if (up_prop_record_xml.count("/SinkSmear") != 0)
	    {
	      QDPIO::cout<<"Up quark propagator is smeared, reading from SinkSmear tag..."<<std::endl;
	      read(up_prop_record_xml, "/SinkSmear", orig_header);
	    }
	    else
	    {
	      QDPIO::cout<<"What type of weird propagator did you give me? I can't find the right tag to get j_decay, source location, and whatnot...exiting..."<<std::endl;
	    }

	    j_decay = orig_header.source_header.j_decay;
	    t_0 = orig_header.source_header.t_source;
	    origin = orig_header.source_header.getTSrce();
	    //If we need to rotate, we do it now.
	    if(params.param.rotate_to_Dirac == true)
	      rotate_to_Dirac_Basis(up_quark_propagator);
            QDPIO::cout << name << ": Read up quark" << std::endl;

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
      }

      if(params.named_obj.is_down == true)
      {
	QDPIO::cout << "Attempting to read down propagator" << std::endl;
	try
	{
	    down_quark_propagator = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.down_quark);
	    TheNamedObjMap::Instance().get(params.named_obj.down_quark).getFileXML(down_prop_file_xml);
	    TheNamedObjMap::Instance().get(params.named_obj.down_quark).getRecordXML(down_prop_record_xml);
	    //Get the origin  and j_decay for the FT, this assumes all quarks have the same origin.

	    MakeSourceProp_t  orig_header;
	    if (down_prop_record_xml.count("/Propagator") != 0)
	    {
	      QDPIO::cout<<"Down quark propagator is unsmeared, reading from Propagator tag..."<<std::endl;
	      read(down_prop_record_xml, "/Propagator", orig_header);
	    }
	    else if (down_prop_record_xml.count("/SinkSmear") != 0)
	    {
	      QDPIO::cout<<"Down quark propagator is smeared, reading from SinkSmear tag..."<<std::endl;
	      read(down_prop_record_xml, "/SinkSmear", orig_header);
	    }
	    else
	    {
	      QDPIO::cout<<"What type of weird propagator did you give me? I can't find the right tag to get j_decay, source location, and whatnot...exiting..."<<std::endl;
	    }

	    j_decay = orig_header.source_header.j_decay;
	    origin = orig_header.source_header.getTSrce();
	    //If we need to rotate, we do it now.
	    if(params.param.rotate_to_Dirac == true)
	      rotate_to_Dirac_Basis(down_quark_propagator);
            QDPIO::cout << name << ": Read down quark" << std::endl;

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
      }

      if(params.named_obj.is_strange == true)
      {
	QDPIO::cout << "Attempting to read strange propagator" << std::endl;
	try
	{
	    strange_quark_propagator = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.strange_quark);
	    TheNamedObjMap::Instance().get(params.named_obj.strange_quark).getFileXML(strange_prop_file_xml);
	    TheNamedObjMap::Instance().get(params.named_obj.strange_quark).getRecordXML(strange_prop_record_xml);
	    //Get the origin  and j_decay for the FT, this assumes all quarks have the same origin.

	    MakeSourceProp_t  orig_header;
	    if (strange_prop_record_xml.count("/Propagator") != 0)
	    {
	      QDPIO::cout<<"Strange quark propagator is unsmeared, reading from Propagator tag..."<<std::endl;
	      read(strange_prop_record_xml, "/Propagator", orig_header);
	    }
	    else if (strange_prop_record_xml.count("/SinkSmear") != 0)
	    {
	      QDPIO::cout<<"Strange quark propagator is smeared, reading from SinkSmear tag..."<<std::endl;
	      read(strange_prop_record_xml, "/SinkSmear", orig_header);
	    }
	    else
	    {
	      QDPIO::cout<<"What type of weird propagator did you give me? I can't find the right tag to get j_decay, source location, and whatnot...exiting..."<<std::endl;
	    }

	    j_decay = orig_header.source_header.j_decay;
	    origin = orig_header.source_header.getTSrce();
	    //If we need to rotate, we do it now.
	    if(params.param.rotate_to_Dirac == true)
	      rotate_to_Dirac_Basis(strange_quark_propagator);
            QDPIO::cout << name << ": Read strange quark" << std::endl;

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
      }

      if(params.named_obj.is_charm == true)
      {
	QDPIO::cout << "Attempting to read charm propagator" << std::endl;
	try
	{
	    charm_quark_propagator = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.charm_quark);
	    TheNamedObjMap::Instance().get(params.named_obj.charm_quark).getFileXML(charm_prop_file_xml);
	    TheNamedObjMap::Instance().get(params.named_obj.charm_quark).getRecordXML(charm_prop_record_xml);
	    //Get the origin  and j_decay for the FT, this assumes all quarks have the same origin.
	    MakeSourceProp_t  orig_header;

	    if (charm_prop_record_xml.count("/Propagator") != 0)
	    {
	      QDPIO::cout<<"Charm quark propagator is unsmeared, reading from Propagator tag..."<<std::endl;
	      read(charm_prop_record_xml, "/Propagator", orig_header);
	    }
	    else if (charm_prop_record_xml.count("/SinkSmear") != 0)
	    {
	      QDPIO::cout<<"Charm quark propagator is smeared, reading from SinkSmear tag..."<<std::endl;
	      read(charm_prop_record_xml, "/SinkSmear", orig_header);
	    }
	    else
	    {
	      QDPIO::cout<<"What type of weird propagator did you give me? I can't find the right tag to get j_decay, source location, and whatnot...exiting..."<<std::endl;
	    }

	    j_decay = orig_header.source_header.j_decay;
	    origin = orig_header.source_header.getTSrce();
	    //If we need to rotate, we do it now.
	    if(params.param.rotate_to_Dirac == true)
	      rotate_to_Dirac_Basis(charm_quark_propagator);
            QDPIO::cout << name << ": Read charm quark" << std::endl;

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
      }

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
      //This is going to be a horrible set of if statements for now, may change this later.
      //Loop over list of particles.

      // Declare lists to check
      std::vector<std::string> pi_octet{"piplus", "octet", "octet_iso"};
      std::vector<std::string> kp_octet{"kplus",  "octet", "octet_iso"};
      std::vector<std::string> ss_octet{"ss_conn","octet", "octet_iso"};

      std::vector<std::string> pz_octet{"pizero_conn", "octet_iso"};
      std::vector<std::string> kz_octet{"kzero",       "octet_iso"};
      std::vector<std::string> uu_octet{"uu_conn",     "octet_iso"};
      std::vector<std::string> dd_octet{"dd_conn",     "octet_iso"};
      std::vector<std::string> ee_octet{"e8_conn",     "octet_iso"};

      for(int particle_index = 0; particle_index < params.param.particle_list.size(); particle_index++){
          // PiPlus
          if( std::find(std::begin(pi_octet), std::end(pi_octet), params.param.particle_list[particle_index]) != std::end(pi_octet)){
              QDPIO::cout<<"Particle number "<<(particle_index+1)<<" is the piplus."<<std::endl;
              QDPIO::cout<<"Checking to make sure we have the correct quark propagators to compute the piplus."<<std::endl;
              if(params.named_obj.is_up == true && params.named_obj.is_down == true){
                  QDPIO::cout<<"Found up and down quarks for the piplus. Starting calculation..."<<std::endl;

                  LatticeComplex piplus = zero;
                  I_one_Iz_pm_one_contract(up_quark_propagator, down_quark_propagator, piplus);
                  // Write out the correlator.
                  write_correlator(params.param.output_full_correlator,
                                   "piplus",
#ifdef BUILD_HDF5
                                   params.param.obj_path, h5out, wmode,
#endif
                                   t_0, Nt, origin, ft, piplus);
              }
              else{
                  QDPIO::cerr<<"You requested PiPlus but did not provide UP and DOWN quarks - aborting" << std::endl;
                  QDP_abort(1);
              }
          }
          // KPlus
          if( std::find(std::begin(kp_octet), std::end(kp_octet), params.param.particle_list[particle_index]) != std::end(kp_octet)){
              QDPIO::cout<<"Particle number "<<(particle_index+1)<<" is the kplus."<<std::endl;
              QDPIO::cout<<"Checking to make sure we have the correct quark propagators to compute the kplus."<<std::endl;
              if(params.named_obj.is_up == true && params.named_obj.is_strange == true){
                  QDPIO::cout<<"Found up and strage quarks for the kplus. Starting calculation..."<<std::endl;

                  LatticeComplex kplus = zero;
                  I_one_Iz_pm_one_contract(up_quark_propagator, strange_quark_propagator, kplus);
                  // Write out the correlator.
                  write_correlator(params.param.output_full_correlator,
                                   "kplus",
#ifdef BUILD_HDF5
                                   params.param.obj_path, h5out, wmode,
#endif
                                   t_0, Nt, origin, ft, kplus);
              }
              else{
                  QDPIO::cerr<<"You requested KPlus but did not provide UP and STRANGE quarks - aborting" << std::endl;
                  QDP_abort(1);
              }
          }
          // Sbar S
          if( std::find(std::begin(ss_octet), std::end(ss_octet), params.param.particle_list[particle_index]) != std::end(ss_octet)){
              QDPIO::cout<<"Particle number "<<(particle_index+1)<<" is the ss_conn."<<std::endl;
              QDPIO::cout<<"Checking to make sure we have the correct quark propagators to compute the ss_conn."<<std::endl;
              if(params.named_obj.is_strange == true){
                  QDPIO::cout<<"Found strange quark for the ss_conn. Starting calculation..."<<std::endl;

                  LatticeComplex ss_conn = zero;
                  I_one_Iz_pm_one_contract(strange_quark_propagator, strange_quark_propagator, ss_conn);
                  // Write out the correlator.
                  write_correlator(params.param.output_full_correlator,
                                   "ss_conn",
#ifdef BUILD_HDF5
                                   params.param.obj_path, h5out, wmode,
#endif
                                   t_0, Nt, origin, ft, ss_conn);
              }
              else{
                  QDPIO::cerr<<"You requested ss_conn but did not provide STRANGE quarks - aborting" << std::endl;
                  QDP_abort(1);
              }
          }
          // Ubar U
          if( std::find(std::begin(uu_octet), std::end(uu_octet), params.param.particle_list[particle_index]) != std::end(uu_octet)){
              QDPIO::cout<<"Particle number "<<(particle_index+1)<<" is the uu_conn."<<std::endl;
              QDPIO::cout<<"Checking to make sure we have the correct quark propagators to compute the uu_conn."<<std::endl;
              if(params.named_obj.is_up == true){
                  QDPIO::cout<<"Found up quark for the uu_conn. Starting calculation..."<<std::endl;

                  LatticeComplex uu_conn = zero;
                  I_one_Iz_pm_one_contract(up_quark_propagator, up_quark_propagator, uu_conn);
                  // Write out the correlator.
                  write_correlator(params.param.output_full_correlator,
                                   "uu_conn",
#ifdef BUILD_HDF5
                                   params.param.obj_path, h5out, wmode,
#endif
                                   t_0, Nt, origin, ft, uu_conn);
              }
              else{
                  QDPIO::cerr<<"You requested uu_conn but did not provide UP quarks - aborting" << std::endl;
                  QDP_abort(1);
              }
          }
          // Dbar D
          if( std::find(std::begin(dd_octet), std::end(dd_octet), params.param.particle_list[particle_index]) != std::end(dd_octet)){
              QDPIO::cout<<"Particle number "<<(particle_index+1)<<" is the dd_conn."<<std::endl;
              QDPIO::cout<<"Checking to make sure we have the correct quark propagators to compute the dd_conn."<<std::endl;
              if(params.named_obj.is_down == true){
                  QDPIO::cout<<"Found down quark for the dd_conn. Starting calculation..."<<std::endl;

                  LatticeComplex dd_conn = zero;
                  I_one_Iz_pm_one_contract(down_quark_propagator, down_quark_propagator, dd_conn);
                  // Write out the correlator.
                  write_correlator(params.param.output_full_correlator,
                                   "dd_conn",
#ifdef BUILD_HDF5
                                   params.param.obj_path, h5out, wmode,
#endif
                                   t_0, Nt, origin, ft, dd_conn);
              }
              else{
                  QDPIO::cerr<<"You requested dd_conn but did not provide DOWN quarks - aborting" << std::endl;
                  QDP_abort(1);
              }
          }
          // Pi Zero
          if( std::find(std::begin(pz_octet), std::end(pz_octet), params.param.particle_list[particle_index]) != std::end(pz_octet)){
              QDPIO::cout<<"Particle number "<<(particle_index+1)<<" is the pizero_conn."<<std::endl;
              QDPIO::cout<<"Checking to make sure we have the correct quark propagators to compute the pizero_conn."<<std::endl;
              if(params.named_obj.is_up == true && params.named_obj.is_down == true){
                  QDPIO::cout<<"Found up and down quarks for the pizero_conn. Starting calculation..."<<std::endl;

                  LatticeComplex uu = zero;
                  LatticeComplex dd = zero;
                  I_one_Iz_pm_one_contract(up_quark_propagator,   up_quark_propagator, uu);
                  I_one_Iz_pm_one_contract(down_quark_propagator, down_quark_propagator, dd);
                  LatticeComplex pizero_conn = 0.5 * uu + 0.5 * dd;
                  // Write out the correlator.
                  write_correlator(params.param.output_full_correlator,
                                   "pizero_conn",
#ifdef BUILD_HDF5
                                   params.param.obj_path, h5out, wmode,
#endif
                                   t_0, Nt, origin, ft, pizero_conn);
              }
              else{
                  QDPIO::cerr<<"You requested pizero_conn but did not provide UP and DOWN quarks - aborting" << std::endl;
                  QDP_abort(1);
              }
          }
          // K Zero
          if( std::find(std::begin(kz_octet), std::end(kz_octet), params.param.particle_list[particle_index]) != std::end(kz_octet)){
              QDPIO::cout<<"Particle number "<<(particle_index+1)<<" is the kzero."<<std::endl;
              QDPIO::cout<<"Checking to make sure we have the correct quark propagators to compute the kzero."<<std::endl;
              if(params.named_obj.is_down == true && params.named_obj.is_strange == true){
                  QDPIO::cout<<"Found down and strage quarks for the kplus. Starting calculation..."<<std::endl;

                  LatticeComplex kzero = zero;
                  I_one_Iz_pm_one_contract(down_quark_propagator, strange_quark_propagator, kzero);
                  // Write out the correlator.
                  write_correlator(params.param.output_full_correlator,
                                   "kzero",
#ifdef BUILD_HDF5
                                   params.param.obj_path, h5out, wmode,
#endif
                                   t_0, Nt, origin, ft, kzero);
              }
              else{
                  QDPIO::cerr<<"You requested KZero but did not provide DOWN and STRANGE quarks - aborting" << std::endl;
                  QDP_abort(1);
              }
          }
          // Eta 8 connected
          if( std::find(std::begin(ee_octet), std::end(ee_octet), params.param.particle_list[particle_index]) != std::end(ee_octet)){
              QDPIO::cout<<"Particle number "<<(particle_index+1)<<" is the e8_conn."<<std::endl;
              QDPIO::cout<<"Checking to make sure we have the correct quark propagators to compute the e8_conn."<<std::endl;
              if(params.named_obj.is_up == true && params.named_obj.is_down == true && params.named_obj.is_strange == true){
                  QDPIO::cout<<"Found up and down and strange quarks for the e8_conn. Starting calculation..."<<std::endl;

                  LatticeComplex uu_e = zero;
                  LatticeComplex dd_e = zero;
                  LatticeComplex ss_e = zero;
                  I_one_Iz_pm_one_contract(up_quark_propagator,      up_quark_propagator, uu_e);
                  I_one_Iz_pm_one_contract(down_quark_propagator,    down_quark_propagator, dd_e);
                  I_one_Iz_pm_one_contract(strange_quark_propagator, strange_quark_propagator, ss_e);
                  // 1/6 UU + 1/6 DD + 4/6 SS
                  LatticeComplex ee_conn = 0.16666666666666667 * uu_e + 0.16666666666666667 * dd_e + 4 * 0.16666666666666667 * ss_e;
                  // Write out the correlator.
                  write_correlator(params.param.output_full_correlator,
                                   "e8_conn",
#ifdef BUILD_HDF5
                                   params.param.obj_path, h5out, wmode,
#endif
                                   t_0, Nt, origin, ft, ee_conn);
              }
              else{
                  QDPIO::cerr<<"You requested e8_conn but did not provide UP and DOWN and STRANGE quarks - aborting" << std::endl;
                  QDP_abort(1);
              }
          }

          // Other specific mesons
          if(params.param.particle_list[particle_index] == "piminus"){
              QDPIO::cout<<"Particle number "<<(particle_index+1)<<" is the piminus."<<std::endl;
              QDPIO::cout<<"Checking to make sure we have the correct quark propagators to compute the piminus."<<std::endl;
              if(params.named_obj.is_up == true && params.named_obj.is_down == true){
                  QDPIO::cout<<"Found up and down quarks for the piminus. Starting calculation..."<<std::endl;

                  LatticeComplex piminus = zero;
                  I_one_Iz_pm_one_contract(down_quark_propagator, up_quark_propagator, piminus);
                  // Write out the correlator.
                  write_correlator(params.param.output_full_correlator,
                                   "piminus",
#ifdef BUILD_HDF5
                                   params.param.obj_path, h5out, wmode,
#endif
                                   t_0, Nt, origin, ft, piminus);
              }
              else{
                  QDPIO::cerr<<"You requested PiMinus but did not provide UP and DOWN quarks - aborting" << std::endl;
                  QDP_abort(1);
              }
          }
          
          if(params.param.particle_list[particle_index] == "kminus"){
              QDPIO::cout<<"Particle number "<<(particle_index+1)<<" is the kminus."<<std::endl;
              QDPIO::cout<<"Checking to make sure we have the correct quark propagators to compute the kminus."<<std::endl;
              if(params.named_obj.is_up == true && params.named_obj.is_strange == true){
                  QDPIO::cout<<"Found up and strange quarks for the kminus. Starting calculation..."<<std::endl;

                  LatticeComplex kminus = zero;
                  I_one_Iz_pm_one_contract(strange_quark_propagator,up_quark_propagator, kminus);
                  // Write out the correlator.
                  write_correlator(params.param.output_full_correlator,
                                   "kminus",
#ifdef BUILD_HDF5
                                   params.param.obj_path, h5out, wmode,
#endif
                                   t_0, Nt, origin, ft, kminus);
              }
              else{
                  QDPIO::cerr<<"You requested KMinus but did not provide UP and STRANGE quarks - aborting" << std::endl;
                  QDP_abort(1);
              }
          }
          
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
