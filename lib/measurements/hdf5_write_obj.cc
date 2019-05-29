// Oct 9, 2017
/*! Arjun Singh Gambhir
 * What up my Lattice peeps?
 * Inline tasks to write NamedObject stuff to disk with h5.
 * cc file that sets things up...
 */

//Protect everything with a preprocessing directive.
#ifdef BUILD_HDF5

#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/named_objmap.h"
#include "io/enum_io/enum_qdpvolfmt_io.h"

//LALIBE stuff
#include "hdf5_write_obj.h"
#include "../io/hdf5_write_obj_funcmap.h"

namespace Chroma 
{ 
  namespace LalibeHDF5WriteNamedObjEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMeas(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;

      const std::string name = "HDF5_WRITE_NAMED_OBJECT";
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// Datatype writer
	success &= HDF5WriteObjCallMapEnv::registerAll();

	// Inline measurement
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);

	registered = true;
      }
      return success;
    }


    //! Object buffer
    void write(XMLWriter& xml, const std::string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "object_id", input.object_id);
      write(xml, "object_type", input.object_type);

      pop(xml);
    }

    //! File output
    void write(XMLWriter& xml, const std::string& path, const Params::File_t& input)
    {
      push(xml, path);

      write(xml, "file_name", input.file_name);
      write(xml, "path", input.path);
      write(xml, "obj_name", input.obj_name);
      //write(xml, "writemode", input.enum_wmode);

      pop(xml);
    }


    //! Object buffer
    void read(XMLReader& xml, const std::string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "object_id", input.object_id);
      read(inputtop, "object_type", input.object_type);
    }

    //! File output
    void read(XMLReader& xml, const std::string& path, Params::File_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "file_name", input.file_name);
      read(inputtop, "path", input.path);
      read(inputtop, "obj_name", input.obj_name);
      //read(inputtop, "writemode", input.enum_wmode);

    }


    // Param stuff
    Params::Params() { frequency = 0; }

    Params::Params(XMLReader& xml_in, const std::string& path) 
    {
      try 
      {
	XMLReader paramtop(xml_in, path);

	if (paramtop.count("Frequency") == 1)
	  read(paramtop, "Frequency", frequency);
	else
	  frequency = 1;

	// Parameters for source construction
	read(paramtop, "NamedObject", named_obj);

	// Read in the destination
	read(paramtop, "File", file);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << __func__ << ": I am really sorry, but I didn't undestand something about your XML. Fix it, run again, and have a nice day! :) " << e << std::endl;
	QDP_abort(1);
      }
    }


    void
    Params::writeXML(XMLWriter& xml_out, const std::string& path) 
    {
      push(xml_out, path);
    
      // Parameters for source construction
      write(xml_out, "NamedObject", named_obj);

      // Write out the destination
      write(xml_out, "File", file);

      pop(xml_out);
    }


    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      START_CODE();

      push(xml_out, "hdf5_write_named_obj");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": object writer" << std::endl;
      StopWatch swatch;

      try
      {
	swatch.reset();

	// Write the object
	swatch.start();
	HDF5Base::writemode wmode;
	//Writer hardcoded to append-to-end instead of truncation.
	wmode = HDF5Base::ate;
	//wmode = HDF5Base::trunc;
	//std::string mode = params.file.enum_wmode;
	//Possible broadcasting, but it looks like it's not required.
	/*QDPIO::cout<<"Write mode is: "<<mode<<std::endl;
	if (params.file.enum_wmode.compare("ate") == 0)
	{
	  wmode = HDF5Base::ate;
	}
	else if (params.file.enum_wmode.compare("trunc") == 0)
	  wmode = HDF5Base::trunc;
	else
	  QDPIO::cerr << __func__ << ": The writemode you have selected doesn't exist. Try either ate or trunc." << std::endl;
	  QDP_abort(1);*/
	HDF5WriteObjCallMapEnv::TheHDF5WriteObjFuncMap::Instance().callFunction(params.named_obj.object_type,
									      params.named_obj.object_id,
									      params.file.file_name,
									      params.file.obj_name,
									      params.file.path, wmode);
	swatch.stop();

	QDPIO::cout << "Object successfully written: time= " 
		    << swatch.getTimeInSeconds() 
		    << " secs" << std::endl;
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << name << ": cast error" 
		    << std::endl;
	QDP_abort(1);
      }
      catch (const std::string& e) 
      {
	QDPIO::cerr << name << ": error message: " << e 
		    << std::endl;
	QDP_abort(1);
      }
    
      QDPIO::cout << name << ": ran successfully" << std::endl;

      pop(xml_out);  // write_named_obj

      END_CODE();
    } 

  }

}

#endif
