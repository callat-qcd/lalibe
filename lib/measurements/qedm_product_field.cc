#include "qedm_product_field.h"
#include "meas/inline/abs_inline_measurement_factory.h"
//#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
#include "util/info/proginfo.h"
//#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"
//#include "../external/external_field.h"
//#include "meas/glue/mesplq.h"
#include "util/gauge/unit_check.h"

using namespace std;

namespace Chroma 
{ 
  namespace LalibeQEDMProductFieldEnv
  {
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, const std::string& path) {
	return new InlineProductField(InlineProductFieldParams(xml_in, path));
      }
      
      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "QEDM_PRODUCT_FIELD";
    
    //! Register all the factories
    bool registerAll() {
      bool success = true; 
      if (! registered)
	{
	  success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	  registered = true;
	}
      return success;
    }
  }


  //! Reader for parameters
  void read(XMLReader& xml, const string& path, InlineProductFieldParams::Param_t& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "quark_type", param.quark_type);
    if (param.quark_type != "up" && param.quark_type != "down" && param.quark_type != "strange") {
      QDPIO::cerr << "Illegal quark_type: " << param.quark_type << endl;
      QDP_abort(1);
    }
    read(paramtop, "external_field", param.ext_field_filename);
    read(paramtop, "coupling", param.coupling);

  }


  //! Writer for parameters
  void write(XMLWriter& xml, const string& path, const InlineProductFieldParams::Param_t& param)
  {
    push(xml, path);

    write(xml, "quark_type", param.quark_type);
    write(xml, "external_field", param.ext_field_filename);
    write(xml, "coupling", param.coupling);
    //write(xml, "debug", param.debug);
    //write(xml, "debugOutfile", param.debug_outfile);
    pop(xml);
  }


  //! Gauge input
  void read(XMLReader& xml, const string& path, InlineProductFieldParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    //read(inputtop, "external_id", input.external_id);
    read(inputtop, "gauge_in_id", input.gauge_in_id);
    read(inputtop, "gauge_out_id", input.gauge_out_id);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlineProductFieldParams::NamedObject_t& input)
  {
    push(xml, path);

    //write(xml, "external_id", input.external_id);
    write(xml, "gauge_in_id", input.gauge_in_id);
    write(xml, "gauge_out_id", input.gauge_out_id);

    pop(xml);
  }


  // Param stuff
  InlineProductFieldParams::InlineProductFieldParams() { frequency = 0; }

  InlineProductFieldParams::InlineProductFieldParams(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;

      // Parameters for source construction
      read(paramtop, "Param", param);

      // Read in the input/output configuration info
      read(paramtop, "NamedObject", named_obj);

      // Possible alternate XML file pattern
      //if (paramtop.count("xml_file") != 0) 
      //{
      //  read(paramtop, "xml_file", xml_file);
      //}
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }


  void InlineProductFieldParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    Chroma::write(xml_out, "Param", param);
    Chroma::write(xml_out, "NamedObject", named_obj);
    //QDP::write(xml_out, "xml_file", xml_file);

    pop(xml_out);
  }


  // Function call
  void InlineProductField::operator()(unsigned long update_no, XMLWriter& xml_out) 
  {
    // If xml file not empty, then use alternate
    //if (params.xml_file != "")
    //{
    //  string xml_file = makeXMLFileName(params.xml_file, update_no);

    //  push(xml_out, "product_field");
    //  write(xml_out, "update_no", update_no);
    //  write(xml_out, "xml_file", xml_file);
    //  pop(xml_out);

    //  XMLFileWriter xml(xml_file);
    //  func(update_no, xml);
    //}
    //else
    {
      func(update_no, xml_out);
    }
  }


  // Real work done here
  void InlineProductField::func(unsigned long update_no, XMLWriter& xml_out) 
  {
    START_CODE();

    StopWatch snoop;
    snoop.reset();
    snoop.start();

    // Test and grab a reference to the gauge field
    XMLBufferWriter gauge_xml;
    try
    {
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_in_id);
      TheNamedObjMap::Instance().get(params.named_obj.gauge_in_id).getRecordXML(gauge_xml);
    }
    catch( std::bad_cast ) 
    {
      QDPIO::cerr << LalibeQEDMProductFieldEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << LalibeQEDMProductFieldEnv::name << ": map call failed: " << e 
		  << endl;
      QDP_abort(1);
    }
    const multi1d<LatticeColorMatrix>& u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_in_id);

    // Read in the external field
    StopWatch swatch;
    swatch.reset(); swatch.start();
    multi1d<LatticeReal> ext_field = importExternalField(params.param.ext_field_filename);
    swatch.stop();

    push(xml_out, "product_field");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << LalibeQEDMProductFieldEnv::name << ": modification by external U(1) field" << endl;

    QDPIO::cout << endl << "     Gauge group: SU(" << Nc << ")" << endl;

    QDPIO::cout << "     volume: " << Layout::lattSize()[0];
    for (int i=1; i<Nd; ++i) {
      QDPIO::cout << " x " << Layout::lattSize()[i];
    }
    QDPIO::cout << endl;

    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.write(xml_out, "Input");

    // Write out the config info
    write(xml_out, "Config_info", gauge_xml);

    push(xml_out, "Output_version");
    write(xml_out, "out_version", 1);
    pop(xml_out);

    // First calculate some gauge invariant observables just for info.
    //MesPlq(xml_out, "Input_Gauge_Observables", u);

    multi1d<LatticeColorMatrix> unew(Nd);

    Real q_charge;
    if (params.param.quark_type == "up")
      q_charge = Real(2./3.);
    else if (params.param.quark_type == "down" || params.param.quark_type == "strange")
      q_charge = Real(-1./3.);
    else {
      QDPIO::cerr << LalibeQEDMProductFieldEnv::name << ": Unknown quark type" << params.param.quark_type <<endl;
      QDP_abort(1);
    }

    // Multiply SU(3) by compactified U(1)
    Real coupling = params.param.coupling;
    Complex I=cmplx(0.0,Real(1.0));
    for(int mu = 0; mu < Nd; mu++) {
      unew[mu]=u[mu] * ( cos(q_charge*coupling*ext_field[mu]) + I * sin(q_charge*coupling*ext_field[mu]) );
    }
    
    QDPIO::cout << "Gauge field modification completed." <<  endl;

    // Write out initial and final fields in human readable format for the lads
    //multi1d<int> nrow(Nd);
    //nrow= Layout::lattSize();

    //if(params.param.debug) {
    //  PrintGaugeField(u, "in_" + params.param.debug_outfile, nrow);
    //  PrintGaugeField(unew, "out_" + params.param.debug_outfile, nrow);
    //}

    // Write out new gauge field measurements.
    unitarityCheck(unew);
    //MesPlq(xml_out, "Output_Gauge_Observables", unew);

    // Now store the configuration to a memory object (see lib/meas/inline/gfix/inline_coulgauge.cc)
    
    XMLBufferWriter out_file_xml, out_record_xml;
    push(out_file_xml, "gauge");
    write(out_file_xml, "id", int(0));
    pop(out_file_xml);
    out_record_xml << gauge_xml;  // THIS HEADER WILL BE INCORRECT: NEEDS FIXING !!!
    
    // Store the gauge field
    TheNamedObjMap::Instance().create< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_out_id);
    TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_out_id) = unew;
    TheNamedObjMap::Instance().get(params.named_obj.gauge_out_id).setFileXML(out_file_xml);
    TheNamedObjMap::Instance().get(params.named_obj.gauge_out_id).setRecordXML(out_record_xml);
    
    pop(xml_out);  // external_field

    snoop.stop();
    QDPIO::cout << LalibeQEDMProductFieldEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs"<<endl
		<< "  IO time: " << swatch.getTimeInSeconds()<<" secs"<<endl;

    QDPIO::cout << LalibeQEDMProductFieldEnv::name << ": ran successfully" << endl;

    END_CODE();
  } 					 


  multi1d<LatticeReal> InlineProductField::importExternalField(const std::string& ext_field_filename)
  {

    if (Nd!=4) { QDP_error_exit("Must have Nd=4."); }

    multi1d<int> nrow = Layout::lattSize(); 
    multi1d<LatticeReal> ext_field(Nd);

    // Read in the U(1) gauge field from a file and process
    multi1d<int> currentSite(Nd);
    double elem;
    int nu;

    ifstream infile( ext_field_filename.c_str(), ios::in | ios::binary );

    if ( infile ) {
      for(currentSite(3)=0; currentSite(3)<nrow(3); ++currentSite(3))
      for(currentSite(0)=0; currentSite(0)<nrow(0); ++currentSite(0))
      for(currentSite(1)=0; currentSite(1)<nrow(1); ++currentSite(1))
      for(currentSite(2)=0; currentSite(2)<nrow(2); ++currentSite(2)) {
        for(int mu=0; mu<4; ++mu) {
          nu=(mu+3)%4; // Converts mike-world (t,x,y,z) to chroma-world (x,y,z,t)
          infile.read( (char*)(&elem), sizeof(double) );
          pokeSite(ext_field[nu],Real(elem),currentSite);
        }
      }
      infile.close();
    } else {
      QDP_error_exit("File doesn't exist.", ext_field_filename.c_str() );
    }

    return ext_field;
  }

};
