/*! Arjun Singh Gambhir
 * Inline tasks to write NamedObject stuff to disk with h5 and then free the object from memory.
 */

//Protect everything with a preprocessing directive.
#ifdef BUILD_HDF5

#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/named_objmap.h"

//LALIBE stuff
#include "hdf5_write_erase_obj.h"

namespace Chroma 
{ 
  namespace LalibeHDF5WriteEraseNamedObjEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMeas(LalibeHDF5WriteNamedObjEnv::Params(xml_in, path));
      }

      bool registered = false;

      const std::string name = "HDF5_WRITE_ERASE_NAMED_OBJECT";
    }

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


    // Func
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      START_CODE();

      push(xml_out, "hdf5_write_erase_named_obj");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": HDF5 Object Writer" << std::endl;
      StopWatch swatch;

      QDPIO::cout << "Attempting to write then delete an object named " << params.named_obj.object_id << std::endl;
      write(xml_out, "object_id", params.named_obj.object_id);
      try
      {
	swatch.reset();

	// Calling the hdf5 writer generically.
	LalibeHDF5WriteNamedObjEnv::InlineMeas hdf5_writer(params);
	hdf5_writer(update_no, xml_out);
	// Deleting the object.
	TheNamedObjMap::Instance().erase(params.named_obj.object_id);
	QDPIO::cout << "Object is gone" << std::endl;
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
      pop(xml_out);  

      END_CODE();
    } 

  }

}

#endif
