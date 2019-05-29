// -*- C++ -*-
/*! Arjun Singh Gambhir
 *  The meat of the hdf5 writing is being done here.
 */

//Protect everything with a preprocessing directive.
#ifdef BUILD_HDF5

#ifndef __hdf5_write_obj_funcmap_h__
#define __hdf5_write_obj_funcmap_h__

#include "singleton.h"
#include "funcmap.h"
#include "chromabase.h"

namespace Chroma
{

  //! Write object function std::map
  /*! \ingroup inlineio */
  namespace HDF5WriteObjCallMapEnv
  { 
    struct DumbDisambiguator {};

    //! Write object function std::map
    /*! \ingroup inlineio */
    typedef SingletonHolder< 
      FunctionMap<DumbDisambiguator,
		  void,
		  std::string,
		  TYPELIST_5(const std::string&,
			     const std::string&, 
			     const std::string&, const std::string&, HDF5Base::writemode),
		  void (*)(const std::string& buffer_id,
			   const std::string& filename, 
			   const std::string& obj_name, const std::string& path, HDF5Base::writemode wmode),
		  StringFunctionMapError> >
    TheHDF5WriteObjFuncMap;

    bool registerAll();
  }

} // end namespace Chroma


#endif

#endif
