/*! Arjun Singh Gambhir
 * Inline tasks to write NamedObject stuff to disk with h5 and then free the object from memory.
 */

//Protect everything with a preprocessing directive.
#ifdef BUILD_HDF5

#ifndef __lalibe_hdf5_write_erase_obj_h__
#define __lalibe_hdf5_write_erase_obj_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

//LALIBE stuff
#include "hdf5_write_obj.h"

namespace Chroma 
{ 
  namespace LalibeHDF5WriteEraseNamedObjEnv 
  {
    bool registerAll();

    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const LalibeHDF5WriteNamedObjEnv::Params& p) : params(p) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    private:
      LalibeHDF5WriteNamedObjEnv::Params params;
    };

  }

}

#endif

#endif
