/*
    Utilities for reading xml info from named object map objects

*/

//Chroma Stuff
#include "chromabase.h"
#include "io/qprop_io.h"
/*
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
#include "util/info/unique_id.h"
#include "io/xml_group_reader.h"
*/

namespace Chroma
{
    namespace LalibeUtilsNambedObjEnv
    {
        // Read a Propagator position
        multi1d<int> get_prop_position(const std::string& named_obj_id){
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
            pos0 = orig_header.source_header.getTSrce();

            return pos0
        }
    }
}
