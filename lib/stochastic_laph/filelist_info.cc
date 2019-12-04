//Arjun Singh Gambhir
//Real authors: Ben Hoerz, John Bulava, Colin Morningstar, ??
//Code directly ported from chroma_laph, more acknowledgements coming later.
#include "filelist_info.h"
#include "xml_help.h"

using namespace std;

namespace Chroma {
  namespace LaphEnv {

 // *************************************************************


FileListInfo::FileListInfo(XmlReader& xml_in)
{ 
 xml_tag_assert(xml_in,"FileListInfo");
 XmlReader xmlr(xml_in,"./descendant-or-self::FileListInfo");
 set_info(xmlr);
}


FileListInfo::FileListInfo(XmlReader& xml_in, const string& inpath)
{
 string path(tidyString(inpath));
 if (xml_tag_count(xml_in,path+"/FileListInfo")!=1){
    xml_cerr(xml_in,"Could not find unique path for FileListInfo construction");
    xml_cerr(xml_in,"Path was "+path+"/FileListInfo");
    xmlreadfail(xml_in,"FileListInfo");}
 XmlReader xmlr(xml_in, "./descendant-or-self::"+path+"/FileListInfo");
 set_info(xmlr);
}


FileListInfo::FileListInfo(const std::string& stub, int min_suffix, 
                           int max_suffix, bool over_write)
{
 try{set_info(stub,min_suffix,max_suffix,over_write);}
 catch(string& msg){
    std::cerr << "invalid FileListInfo construction"<<std::endl;
    std::cerr << msg<<std::endl;
    QDP_abort(1);}
}


void FileListInfo::set_info(XmlReader& xmlr)
{
 string stub;
 xmlread(xmlr,"FileNameStub",stub,"FileListInfo");

 int min_suffix=0;
 if (xml_tag_count(xmlr,"MinFileNumber")==1)
    xmlread(xmlr,"MinFileNumber",min_suffix,"FileListInfo");

 int max_suffix;
 xmlread(xmlr,"MaxFileNumber",max_suffix,"FileListInfo");

 bool overwrite = false;  // protect mode
 if (xml_tag_count(xmlr,"FileMode")==1){
    string fmode;
    xmlread(xmlr,"FileMode",fmode,"FileListInfo");
    fmode=tidyString(fmode);
    if (fmode=="overwrite") overwrite=true;}

 try{set_info(stub,min_suffix,max_suffix,overwrite);}
 catch(string& msg){
    xml_cerr(xmlr,msg); xmlreadfail(xmlr,"FileListInfo");}
}


void FileListInfo::set_info(const std::string& stub, int min_suffix, 
                            int max_suffix, bool over_write)
{
 m_file_stub=tidyString(stub);
 if (m_file_stub.empty()){
    throw(string("Blank file name in FileListInfo"));}
 m_max_file_number=max_suffix;
 m_overwrite_mode = over_write;  // protect mode
 m_min_file_number=min_suffix;
 if ((m_min_file_number>m_max_file_number)||(m_min_file_number<0)){
    throw(string("minimum file number > maximum file number in FileListInfo"));}
} 


FileListInfo::FileListInfo(const FileListInfo& fin)
   : m_file_stub(fin.m_file_stub), 
     m_max_file_number(fin.m_max_file_number),
     m_min_file_number(fin.m_min_file_number),
     m_overwrite_mode(fin.m_overwrite_mode) {}


FileListInfo& FileListInfo::operator=(const FileListInfo& fin)
{
 m_file_stub=fin.m_file_stub;
 m_max_file_number=fin.m_max_file_number;
 m_min_file_number=fin.m_min_file_number;
 m_overwrite_mode=fin.m_overwrite_mode;
 return *this;
}


std::string FileListInfo::getFileName(int suffix) const
{
 stringstream fs;
 fs << m_file_stub << "." << suffix;
 return fs.str();
}

int FileListInfo::getFirstAvailableSuffix(bool global_mode) const
{
 for (int suffix=m_min_file_number;suffix<=m_max_file_number;suffix++){
    string filename=getFileName(suffix);
    if (!fileExists(filename,global_mode)) return suffix;}
 std::cerr << "no suffix numbers are available for writing"<<std::endl;
 std::cerr << " ... increase maxFilenumber"<<std::endl;
 throw(string("error"));
}


bool FileListInfo::operator==(const FileListInfo& in) const
{
 return  ((m_file_stub==in.m_file_stub)           
        &&(m_max_file_number==in.m_max_file_number)
        &&(m_min_file_number==in.m_min_file_number));
}

void FileListInfo::output(XmlWriter& xmlout) const
{
 push(xmlout,"FileListInfo");
 write(xmlout,"FileNameStub",m_file_stub);
 write(xmlout,"MinFileNumber",m_min_file_number);
 write(xmlout,"MaxFileNumber",m_max_file_number);
 if (m_overwrite_mode) write(xmlout,"FileMode","overwrite");
 else write(xmlout,"FileMode","protect");
 pop(xmlout);
}

// ***************************************************************
  }
}

