// Oct 9, 2017
/*! Arjun Singh Gambhir
 * What up my Lattice peeps?
 * Inline tasks to write NamedObject stuff to disk with h5.
 * The guts of the code are in my funcmap file, but hey, someone's gotta register this stuff into the factory. :)
 */

//Protect everything with a preprocessing directive.
#ifdef BUILD_HDF5

#ifndef __lalibe_hdf5_write_obj_h__
#define __lalibe_hdf5_write_obj_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlineio */
  namespace LalibeHDF5WriteNamedObjEnv 
  {
    bool registerAll();

    //! Parameter structure
    /*! \ingroup inlineio */
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long frequency;

      struct NamedObject_t
      {
	std::string   object_id;
	std::string   object_type;
      } named_obj;

      struct File_t
      {
	std::string   file_name;
	std::string   path;
	std::string   obj_name;
	//std::string   enum_wmode;
      } file;
    };

    //! Inline writing of memory objects
    /*! \ingroup inlineio */
    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const Params& p) : params(p) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      //! Do the writing
      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    private:
      Params params;
    };

  }

}

#endif

#endif
