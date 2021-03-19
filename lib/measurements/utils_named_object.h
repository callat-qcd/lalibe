/*
    Utilities for reading xml info from named object map objects

*/

#ifndef __lalibe_Utils_Named_Object_h__
#define __lalibe_Utils_Named_Object_h__

//Chroma Stuff
#include "chromabase.h"

namespace Chroma
{
    namespace LalibeUtilsNambedObjEnv
    {
        extern const std::string name;

        // Read a Propagator position
        multi1d<int> get_prop_position(const std::string& named_obj_id);

    }
}

#endif
