// -*- C++ -*-
/*! This file reads an object written from the totally snazzy hdf5 writer and puts in named object storage.
 *  If you find bugs, blame Arjun Singh Gambhir
 */

//Protect everything with a preprocessing directive.
#ifdef BUILD_HDF5

#ifndef __lalibe_hdf5_read_obj_h__
#define __lalibe_hdf5_read_obj_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlineio */
  namespace LalibeHDF5ReadNamedObjEnv 
  {
    bool registerAll();

    //! Parameter structure
    /*! \ingroup inlineio */
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);

      unsigned long frequency;

      struct NamedObject_t
      {
	std::string   object_id;
	std::string   object_type;
      };

      struct File_t
      {
	std::string   file_name;
	std::string   path;
	std::string   obj_name;
      };

      File_t file;
      NamedObject_t named_obj;
      GroupXML_t    named_obj_xml;  /*!< Holds standard named objects */
    };

    //! Inline reading of hdf5 objects
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

  } // namespace LalibeHDF5ReadNamedObjEnv 

} // namespace Chroma

#endif

#endif
