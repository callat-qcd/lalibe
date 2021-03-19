/*
    Utilities for reading xml info from named object map objects

*/

//Chroma Stuff
#include "chromabase.h"
#include "io/qprop_io.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma
{
    namespace LalibeUtilsNambedObjEnv
    {
        std::string name = "UTILS_NAMED_OBJECT";

        // Read a Propagator position
        multi1d<int> get_prop_position(const std::string& named_obj_id){
            try{
                XMLReader prop_record_xml;
                TheNamedObjMap::Instance().get(named_obj_id).getRecordXML(prop_record_xml);
                MakeSourceProp_t  orig_header;
                if (prop_record_xml.count("/Propagator") != 0){
                    read(prop_record_xml, "/Propagator", orig_header);
                }
                // Or if we pass a smeared propagator
                else if (prop_record_xml.count("/SinkSmear") != 0){
                    read(prop_record_xml, "/SinkSmear", orig_header);
                }

                multi1d<int> pos0 = orig_header.source_header.getTSrce();

                return pos0;
            }
            catch (const std::string& e) {
                QDPIO::cerr << LalibeUtilsNambedObjEnv::name << ": ERROR reading src prop_header: "
                            << e << std::endl;
                QDP_abort(1);
            }
        }
    }
}
