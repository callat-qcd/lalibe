//Arjun Singh Gambhir
//Real authors: Ben Hoerz, John Bulava, Colin Morningstar, ??
//Code directly ported from chroma_laph, more acknowledgements coming later.
#include "gauge_configuration_handler.h"
#include "chroma.h"
#include "xml_help.h"

namespace Chroma {
  namespace LaphEnv {

 // **********************************************************

GaugeConfigurationHandler::GaugeConfigurationHandler()
     : gauge_info(0), cfg(0) {}


GaugeConfigurationHandler::GaugeConfigurationHandler(
                      const GaugeConfigurationInfo& gaugeinfo,
                      const std::string& gauge_id)
     : gauge_info(0), cfg(0)
{
 set_info(gaugeinfo,gauge_id);
}


void GaugeConfigurationHandler::setInfo(const GaugeConfigurationInfo& gaugeinfo,
                                        const std::string& gauge_id)
{
 clear();
 set_info(gaugeinfo,gauge_id);
}

void GaugeConfigurationHandler::set_info(const GaugeConfigurationInfo& gaugeinfo, 
                                         const std::string& gauge_id)
{
 objmap_gauge_id=gauge_id;
 if (objmap_gauge_id.empty()) objmap_gauge_id="default_gauge_field";
 try{
    gauge_info = new GaugeConfigurationInfo(gaugeinfo);}
 catch(...){
    QDPIO::cerr << "problem allocating GaugeConfigurationInfo"<<std::endl;
    QDP_abort(1);}    
// QDPIO::cout << "GaugeConfig info set in GaugeConfigurationHandler"<<endl;
}

void GaugeConfigurationHandler::setData()
{
 check_info_set("setData");

// string gauge_id = gauge_info->getGaugeId();
// if (gauge_id.empty()){
//    QDPIO::cerr << "empty gauge_id in GaugeConfigurationHandler" << endl;
//    QDP_abort(1);}

            //Assign the cfg from the named object map
 XMLBufferWriter gauge_xml_buff;
 try{
    TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(objmap_gauge_id);
    TheNamedObjMap::Instance().get(objmap_gauge_id).getRecordXML(gauge_xml_buff);
    }
 catch( std::bad_cast ){
    QDPIO::cerr << __func__ << ": caught dynamic cast error" << std::endl;
    QDP_abort(1);}
 catch (const std::string& err){
    QDPIO::cerr << __func__ << ": map call failed: " << err << std::endl;
    QDP_abort(1);}

 cfg = &(TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(objmap_gauge_id));

}


GaugeConfigurationHandler::~GaugeConfigurationHandler()
{
 clear();
}
    
void GaugeConfigurationHandler::clear()
{
 try {delete gauge_info;} catch(...) {QDP_abort(1);}
 objmap_gauge_id.clear();
 gauge_info=0;
 cfg=0;
}

const multi1d<LatticeColorMatrix>& GaugeConfigurationHandler::getData() 
{
 check_info_set("getData");
 if (!isDataSet()) setData();
 return *cfg;
}

const GaugeConfigurationInfo& GaugeConfigurationHandler::getGaugeConfigurationInfo() const
{
 check_info_set("getGaugeConfigurationInfo");
 return *gauge_info;
}

void GaugeConfigurationHandler::check_info_set(const std::string& name) const
{
 if (!isInfoSet()){
    QDPIO::cerr << "error in GaugeConfigurationHandler:"<<std::endl;
    QDPIO::cerr << "  must setInfo before calling "<<name<<std::endl;
    QDP_abort(1);}
}

std::string GaugeConfigurationHandler::getFullRecordXML() const
{
 check_info_set("getFullRecordXML");
 std::string gauge_xml;

        //Test reference to the cfg in the named object map
 XMLBufferWriter gauge_header_buff;
 try{   
    TheNamedObjMap::Instance().get(objmap_gauge_id).getRecordXML(gauge_header_buff);
    }
 catch( std::bad_cast ){
    QDPIO::cerr << "caught dynamic cast error in GaugeConfigurationInfo" << std::endl;
    throw(std::string("error"));
    }
 catch (const std::string& err){
    QDPIO::cerr << "TheNamedObjMap call failed in GaugeConfigurationInfo: " << err << std::endl;
    throw(std::string("error"));
    }
 XMLReader xmlg0(gauge_header_buff);   
 XMLReader xmlg(xmlg0,"/");         // required in case of XMLReader bug

    // make the gauge_header string (the configuration could be
    // generated from Monte Carlo updating, or could be an
    // initial configuration from a cold or hot start)

 if (xml_tag_count(xmlg,"MCControl")>0){   // from Monte Carlo generated config
    try{
       XMLReader xmlg1(xmlg,".//MCControl");
       XMLReader xmlg2(xmlg,".//HMCTrj");
       std::ostringstream oss;
       oss << "<GaugeConfigHeader>"<<std::endl;
       xmlg1.print(oss);
       xmlg2.print(oss);
       oss << "</GaugeConfigHeader>";
       gauge_xml = oss.str();}
    catch(const std::string& err){
       QDPIO::cerr << "Could not create gauge_header in GaugeConfigurationInfo"<<std::endl;
       throw(std::string("error"));}}
 else{
    std::ostringstream temp_strm;                 // from a dummy config
    temp_strm << "<GaugeConfigHeader>"<<std::endl;
    xmlg.printCurrentContext(temp_strm);
    temp_strm << "</GaugeConfigHeader>";
    gauge_xml = temp_strm.str();}

 return gauge_xml;
}


// **********************************************************************
  }
}

