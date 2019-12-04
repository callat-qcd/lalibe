//Arjun Singh Gambhir
//Real authors: Ben Hoerz, John Bulava, Colin Morningstar, ??
//Code directly ported from chroma_laph, more acknowledgements coming later.
//The definitions for the GaugeConfiguration class

extern "C" {
   #include "lime_reader.h"
   }

#include "gauge_configuration_info.h"
#include "chroma.h"
#include "xml_help.h"
#include <sstream>

using namespace std;

namespace Chroma {
  namespace LaphEnv {
  

// *************************************************************

   //  Constructor expect XML of the form
   //     <GaugeConfigurationInfo>
   //          .....
   //     </GaugeConfigurationInfo>
   //
   //  Use a tag <gauge_id>....</gauge_id> to set from TheNamedObjMap.
   //  Use a tag <getFromFile> ....</getFromFile> to set from a config file.
   //  Otherwise, all tags must be explicitly given.


GaugeConfigurationInfo::GaugeConfigurationInfo(XmlReader& xml_in)
{
 string gauge_xml;
 set_info(xml_in,gauge_xml);
}


GaugeConfigurationInfo::GaugeConfigurationInfo(XmlReader& xml_in,
                                               string& gauge_xml)
{
 set_info(xml_in,gauge_xml,true);
}


void GaugeConfigurationInfo::set_info(XmlReader& xml_in,
                                      string& gauge_xml, bool read_gauge_xml)
{
 xml_tag_assert(xml_in,"GaugeConfigurationInfo");
 XmlReader xmlr(xml_in,"./descendant-or-self::GaugeConfigurationInfo");
 int n_id = xml_tag_count(xmlr,"gauge_id");
 int n_fn = xml_tag_count(xmlr,"getFromFile");

 if ((n_id==1)&&(n_fn==0)){
    string gauge_id;
    xmlread(xmlr, "gauge_id", gauge_id, "GaugeConfigurationInfo");
    gauge_id=tidyString(gauge_id);
    if (gauge_id.empty()){
       xml_cerr(xml_in,"empty gauge_id in GaugeConfigurationInfo");
       xmlreadfail(xmlr,"GaugeConfigurationInfo");}
    set_from_named_obj_map(gauge_id,gauge_xml);
    }
 else if ((n_fn==1)&&(n_id==0)){
    string config_file;
    xmlread(xmlr, "getFromFile", config_file, "GaugeConfigurationInfo");
    config_file=tidyString(config_file);
    if (config_file.empty()){
       xml_cerr(xml_in,"empty file name in GaugeConfigurationInfo");
       xmlreadfail(xmlr,"GaugeConfigurationInfo");}
    set_from_file_name(config_file,gauge_xml);
    }
 else if ((n_fn==0)&&(n_id==0)){
    set_from_xml(xmlr);
    if (read_gauge_xml && config_type!="CERN" && file_name != "none") 
       gauge_xml=getFileXML();
    }
 else{
    xml_cerr(xml_in,"Bad XML input to GaugeConfigurationInfo");
    xmlreadfail(xml_in,"GaugeConfigurationInfo");}

 checkDimensions();
}



void GaugeConfigurationInfo::set_from_named_obj_map(const string& gauge_id,
                                                    string& gauge_xml)
{
 XMLBufferWriter gauge_header_buff;
 try{   
    TheNamedObjMap::Instance().get(gauge_id).getRecordXML(gauge_header_buff);
    }
 catch( std::bad_cast ){
    QDPIO::cerr << "caught dynamic cast error in GaugeConfigurationInfo" << endl;
    throw(string("error"));
    }
 catch (const string& err){
    QDPIO::cerr << "TheNamedObjMap call failed in GaugeConfigurationInfo: " 
         << err << endl;
    throw(string("error"));
    }
 XmlReader xmlg(gauge_header_buff);   
 gauge_xml=gauge_header_buff.str();
 try{set_from_file_xml(xmlg);}
 catch(...){
    QDPIO::cerr << "could not set GaugeConfigurationInfo from file xml"<<endl;
    QDPIO::cerr << "gauge_xml = "<<gauge_xml<<endl;
    QDP_abort(1);}
}


void GaugeConfigurationInfo::set_from_file_name(const string& configfile,
                                                string& gauge_xml)
{
 try { readGaugeXML(configfile,gauge_xml);}
 catch(...){
     QDPIO::cerr << "could not read file XML in GaugeConfigurationInfo"<<endl;
     QDP_abort(1);}
 stringstream oss;
 oss << gauge_xml;
 XmlReader xmlg(oss);
 try{set_from_file_xml(xmlg);}
 catch(...){
    QDPIO::cerr << "could not set GaugeConfigurationInfo from file xml"<<endl;
    QDPIO::cerr << "gauge_xml = "<<gauge_xml<<endl;
    QDP_abort(1);}
 file_name=tidyString(configfile);  // instead of from file_xml
}


string GaugeConfigurationInfo::getFileXML() const
{
 string gauge_xml;
 try { readGaugeXML(file_name,gauge_xml);}
 catch(...){
     QDPIO::cerr << "could not read file XML in GaugeConfigurationInfo"<<endl;
     QDPIO::cerr << "check that filename stored in configuration file is not problem"<<endl;
     throw(string("error"));}
 return gauge_xml;
}


void GaugeConfigurationInfo::set_from_file_xml(XmlReader& xmlg)
{
 if (xmlg.count(".//StartUpdateNum") != 0){
    xmlread(xmlg,"StartUpdateNum",traj_num,"GaugeConfigurationInfo");
    xmlread(xmlg,"HMCTrj/MDIntegrator/t_dir",time_dir,"GaugeConfigurationInfo");
    xmlread(xmlg,"cfg_type",config_type,"GaugeConfigurationInfo");
//    xmlread(xmlg,"cfg_file",file_name,"GaugeConfigurationInfo");
    xmlread(xmlg,"HMCTrj/nrow",extents,"GaugeConfigurationInfo");}
 else{
    traj_num = 1000;          // the case of a random or unit gauge config
    time_dir = 3;
    config_type = "dummy";
    file_name = "none";
    extents=QDP::Layout::lattSize();}
 
 config_type=tidyString(config_type);
// file_name=tidyString(file_name);
 number_dir=extents.size();
 if ((time_dir<0)||(time_dir>=number_dir)){
    xmlreadfail(xmlg,"time_dir invalid in GaugeConfigurationInfo");}
 time_extent=extents[time_dir];
}


void GaugeConfigurationInfo::set_from_xml(XmlReader& xmlg)
{
 xmlread(xmlg,"HMCTrajectoryNumber",traj_num,"GaugeConfigurationInfo");
 xmlread(xmlg,"TimeDir",time_dir,"GaugeConfigurationInfo");
 xmlread(xmlg,"ConfigType",config_type,"GaugeConfigurationInfo");
 xmlread(xmlg,"FileName",file_name,"GaugeConfigurationInfo");
 xmlread(xmlg,"TimeExtent",time_extent,"GaugeConfigurationInfo");
 xmlread(xmlg,"NumberOfDir",number_dir,"GaugeConfigurationInfo");
 xmlread(xmlg,"LatticeExtents",extents,"GaugeConfigurationInfo");
 file_name=tidyString(file_name);
}

void GaugeConfigurationInfo::checkDimensions()
{
 if ((number_dir != 4)||(time_dir<0)||(time_dir>=number_dir)){
    QDPIO::cerr << "Invalid dimensions in GaugeConfigurationInfo"<<endl;
    QDP_abort(1);}
 if (time_extent != extents[time_dir]){
    QDPIO::cerr << "Invalid time extent in GaugeConfigurationInfo"<<endl;
    QDP_abort(1);}

#if (QDP_ND == 4)     // check lattice dimensions against layout

 for (int dir4d=0;dir4d<QDP::Nd;dir4d++){
    if (extents[dir4d]!=QDP::Layout::lattSize()[dir4d]){
       QDPIO::cerr << "GaugeConfigurationInfo lattice size does not match"
                   <<" the current layout"<<endl;
       QDPIO::cerr << "Layout: "<<Layout::lattSize()[0]<<" x "
           <<Layout::lattSize()[1]<<" x "<<Layout::lattSize()[2]<<" x "
           <<Layout::lattSize()[3]<<endl;
       QDPIO::cerr << "GaugeConfigurationInfo: "<<extents[0]<<" x "
           <<extents[1]<<" x "<<extents[2]<<" x " <<extents[3]<<endl;
       QDP_abort(1);}}

#elif (QDP_ND == 3)      // check 3 spatial dimensions against layout

 multi1d<int> extents3d = QDP::Layout::lattSize();     // 3-d extents
 int dir3d=0;
 for (int dir4d=0;dir4d<=QDP::Nd;dir4d++){
    if (dir4d!=time_dir){
       if (extents[dir4d]!=extents3d[dir3d]){
          QDPIO::cerr << "Dimensions of 3-d lattice do not match"
                      <<" the 3 spatial dimensions of the 4-d lattice"<<endl;
          QDPIO::cerr << "3-d lattice dimensions: "<<extents3d[0]<<" x "
                      << extents3d[1]<<" x "<<extents3d[2]<<endl;
          QDPIO::cerr << "4-d lattice dimensions: "<<extents[0]<<" x "
                      << extents[1]<<" x "<<extents[2]<<" x "
                      << extents[3]<<" time direction = "<<time_dir<<endl;
          QDP_abort(1);}
       dir3d++;}}

#endif
}

// *******************************************************************

    // copy constructor
    
GaugeConfigurationInfo::GaugeConfigurationInfo(
                const GaugeConfigurationInfo& rhs) 
         : file_name(rhs.file_name), config_type(rhs.config_type), 
           traj_num(rhs.traj_num), time_dir(rhs.time_dir), 
           time_extent(rhs.time_extent), extents(rhs.extents),
           number_dir(rhs.number_dir) {}


GaugeConfigurationInfo& GaugeConfigurationInfo::operator=(
               const GaugeConfigurationInfo& rhs)
{
 file_name = rhs.file_name;
 traj_num = rhs.traj_num;
 config_type = rhs.config_type;
 time_dir = rhs.time_dir;
 time_extent = rhs.time_extent;
 number_dir = rhs.number_dir;
 extents = rhs.extents;
 return *this;
}

void GaugeConfigurationInfo::checkEqual(const GaugeConfigurationInfo& rhs) const
{
 if  ((tidyFileName(file_name)!=tidyFileName(rhs.file_name))
    ||(traj_num!=rhs.traj_num)
    ||(config_type!=rhs.config_type)
    ||(time_dir!=rhs.time_dir)
    ||(time_extent!=rhs.time_extent)
    ||(number_dir!=rhs.number_dir)
    ||(extents!=rhs.extents)){
    QDPIO::cerr << "GaugeConfigurationInfo checkEqual failed"<<endl;
    QDPIO::cerr << "LHS:"<<endl<<output()<<endl<<"RHS:"<<endl<<rhs.output()<<endl;
    throw string("GaugeConfigurationInfo checkEqual failed");}
}


bool GaugeConfigurationInfo::operator==(const GaugeConfigurationInfo& rhs) const
{
 return ((tidyFileName(file_name)==tidyFileName(rhs.file_name))
       &&(traj_num==rhs.traj_num)
       &&(config_type==rhs.config_type)
       &&(time_dir==rhs.time_dir)
       &&(time_extent==rhs.time_extent)
       &&(number_dir==rhs.number_dir)
       &&(extents==rhs.extents));
}


string GaugeConfigurationInfo::output(int indent) const
{
 string pad(3*indent,' ');
 ostringstream oss;
 oss << pad << "<GaugeConfigurationInfo>"<<endl;
 oss << pad << "   <FileName>" << file_name << "</FileName>"<<endl;
 oss << pad << "   <ConfigType>" << config_type << "</ConfigType>"<<endl;
 oss << pad << "   <HMCTrajectoryNumber>" << traj_num << "</HMCTrajectoryNumber>"<<endl;
 oss << pad << "   <TimeDir>" << time_dir << "</TimeDir>"<<endl;
 oss << pad << "   <TimeExtent>" << time_extent << "</TimeExtent>"<<endl;
 oss << pad << "   <NumberOfDir>" << number_dir << "</NumberOfDir>"<<endl;
 oss << pad << "   <LatticeExtents>";
 if (extents.size()>0) oss << extents[0];
 for (int k=1;k<extents.size();k++) oss << " "<<extents[k];
 oss << "</LatticeExtents>"<<endl;
 oss << pad << "</GaugeConfigurationInfo>"<<endl;
 return oss.str();
}

void GaugeConfigurationInfo::output(XmlWriter& xmlout) const
{
 push(xmlout,"GaugeConfigurationInfo");
 write(xmlout,"FileName",file_name);
 write(xmlout,"ConfigType",config_type);
 write(xmlout,"HMCTrajectoryNumber",traj_num);
 write(xmlout,"TimeDir",time_dir);
 write(xmlout,"TimeExtent",time_extent);
 write(xmlout,"NumberOfDir",number_dir);
 write(xmlout,"LatticeExtents",extents);
 pop(xmlout);
}


// *********************************************************

   // This routine reads the file xml content of a gauge configuration
   // without reading the binary data or checking sizes.

   
void GaugeConfigurationInfo::readGaugeXML(
                  const string& config_file_name,
                  string& gauge_xml) const
{

 gauge_xml.clear();
 n_uint64_t nbytes, read_bytes;
 int status;
 bool flag=false;
 FILE *fp=0;
 LimeReader *reader=0;

 if (Layout::primaryNode()){

  try{
   
   fp = DCAPL(fopen)(config_file_name.c_str(), "r");
   if (fp==(FILE*)NULL){    
      throw(string("Unable to open file "+config_file_name+" for reading"));}
   reader = limeCreateReader(fp);
   if (reader==(LimeReader*)NULL){ 
      throw(string("Unable to open LimeReader"));}

      // loop over records to find the non-private file-xml
     
   while ( (status = limeReaderNextRecord(reader)) != LIME_EOF ){
    
      if ( status != LIME_SUCCESS ) { 
         throw(string("limeReaderNextRecord returned bad status"));}

      nbytes = limeReaderBytes(reader);
      char *lime_type = limeReaderType(reader);
      //size_t bytes_pad = limeReaderPadBytes(reader);
      //int MB_flag = limeReaderMBFlag(reader);
      //int ME_flag = limeReaderMEFlag(reader);

      //QDPIO::cout << endl<<endl;
      //QDPIO::cout << "Type:           "<< lime_type << endl;
      //QDPIO::cout << "Data Length:    "<< size_t(nbytes) << endl;
      //QDPIO::cout << "Padding Length: "<< bytes_pad <<endl;
      //QDPIO::cout << "MB flag:        "<< MB_flag <<endl;
      //QDPIO::cout << "ME flag:        "<< ME_flag <<endl;

        // read file header information
        
      if ((strstr(lime_type,"record-xml")!=0) && (strstr(lime_type,"private")==0)){
         //QDPIO::cout << "This is the file xml!!"<<endl;
         char* strbuf = new char[nbytes]; 
         read_bytes = nbytes;
         status = limeReaderReadData((void *)strbuf, &read_bytes, reader);
         if ((status<0 && status!=LIME_EOR) || (read_bytes!=nbytes)){ 
            throw(string("LIME read error occurred"));}
         gauge_xml.assign(strbuf,nbytes-1);  // don't output end '\0' null char
         //QDPIO::cout << "file_xml = "<<gauge_xml<<endl;
         delete [] strbuf;
         flag=true;
         break;
         }
      }}
  catch(const string& message){
     cerr << "Error: "<<message<<endl;
     flag=false;}
 } // end of primary node branch

 QDPInternal::broadcast(flag);

 if (Layout::primaryNode()){
    if (fp!=(FILE*)NULL){    
      limeDestroyReader(reader);
      DCAP(fclose)(fp); }}

 if (!flag){
    throw(string("Unable to extract file-xml from "+config_file_name));}

 QDPInternal::broadcast_str(gauge_xml);
}


// ******************************************************************
  }
}
