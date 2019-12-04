//Arjun Singh Gambhir
//Real authors: Ben Hoerz, John Bulava, Colin Morningstar, ??
//Code directly ported from chroma_laph, more acknowledgements coming later.
#include "xml_help.h"
#include <string>
#include <unistd.h>
using namespace std;

namespace Chroma {
  namespace LaphEnv {

// *********************************************************


XmlReader::XmlReader(std::istream& is, bool global_mode)
          : derived(false), global(global_mode)
{
 active = (global) ? Layout::primaryNode() : true;
 if (active) BasicXPathReader::open(is);
}

XmlReader::XmlReader(std::string& s, bool global_mode)
          : derived(false), global(global_mode)
{
 stringstream oss;
 oss << s;
 active = (global) ? Layout::primaryNode() : true;
 if (active) BasicXPathReader::open(oss);
}

XmlReader::XmlReader(const XmlBufferWriter& mw, bool global_mode)
          : derived(false), global(global_mode)
{
 if (mw.global!=global) 
    throw("global mode mismatch in XmlBufferWriter conversion to XmlReader");
 active = (global) ? Layout::primaryNode() : true;
 if (active){  
    std::istringstream is(mw.output_stream.str()+"\n");
    BasicXPathReader::open(is);}
}

XmlReader::XmlReader(XMLReader& old, bool global_mode) 
          : derived(false), global(global_mode)
{
 if (!global) 
    throw("must be in global mode for XMLReader conversion to XmlReader");
 active = (global) ? Layout::primaryNode() : true;
 std::stringstream xmlbuf;
 old.print(xmlbuf);
 if (active){  
    BasicXPathReader::open(xmlbuf);}
}

XmlReader::XmlReader(const XMLBufferWriter& mw, bool global_mode)
          : derived(false), global(global_mode)
{
 if (!global) 
    throw("global mode required in XMLBufferWriter conversion to XmlReader");
 active = (global) ? Layout::primaryNode() : true;
 std::istringstream is(mw.str()+"\n");
 if (active){  
    BasicXPathReader::open(is);}
}

XmlReader::XmlReader(XmlReader& start, const string& xpath)
          : BasicXPathReader(), derived(true), active(start.active),
            global(start.global)
{
 if (active) BasicXPathReader::open((BasicXPathReader&)start, xpath);
}


XmlReader::~XmlReader() 
{
 if (active) BasicXPathReader::close();
}


int XmlReader::count(const string& xpath)
{
 int n=0;
 if (active) n = BasicXPathReader::count(xpath);
 if (global) QDPInternal::broadcast(n);
 return n;
}

void XmlReader::read(const std::string& xpath, string& input)
{
 if (active) BasicXPathReader::get(xpath, input);
 if (global) QDPInternal::broadcast_str(input);
}

void XmlReader::read(const std::string& xpath, int& input)
 { readPrimitive<int>(xpath, input); }

void XmlReader::read(const std::string& xpath, unsigned int& input)
 { readPrimitive<unsigned int>(xpath, input); }

void XmlReader::read(const std::string& xpath, short int& input)
 { readPrimitive<short int>(xpath, input); }

void XmlReader::read(const std::string& xpath, unsigned short int& input)
 { readPrimitive<unsigned short int>(xpath, input); }

void XmlReader::read(const std::string& xpath, long int& input)
 { readPrimitive<long int>(xpath, input); }

void XmlReader::read(const std::string& xpath, unsigned long int& input)
 { readPrimitive<unsigned long int>(xpath, input); }

void XmlReader::read(const std::string& xpath, float& input)
 { readPrimitive<float>(xpath, input); }

void XmlReader::read(const std::string& xpath, double& input)
 { readPrimitive<double>(xpath, input); }

void XmlReader::read(const std::string& xpath, bool& input)
 { readPrimitive<bool>(xpath, input); }
   

template<typename T>
void XmlReader::multiReadPrimitive(const std::string& xpath, multi1d<T>& input)
{
 std::string list_string;
 read(xpath, list_string);
 std::istringstream list_stream(list_string);
 int array_size = 0;
 T dummy;
 while (list_stream >> dummy) ++array_size;
 if ((!list_stream.eof()) && list_stream.fail()){
    throw("Error in reading multi1d "+xpath);}
 input.resize(array_size);
 std::istringstream list_stream2(list_string);
 for (int i=0;i<input.size();i++){
    list_stream2 >> input[i];}
}


void XmlReader::read(const std::string& xpath, multi1d<int>& result)
 { multiReadPrimitive(xpath, result); }

void XmlReader::read(const std::string& xpath, multi1d<unsigned int>& result)
 { multiReadPrimitive(xpath, result); }

void XmlReader::read(const std::string& xpath, multi1d<short int>& result)
 { multiReadPrimitive(xpath, result); }

void XmlReader::read(const std::string& xpath, multi1d<unsigned short int>& result)
 { multiReadPrimitive(xpath, result); }

void XmlReader::read(const std::string& xpath, multi1d<long int>& result)
 { multiReadPrimitive(xpath, result); }

void XmlReader::read(const std::string& xpath, multi1d<unsigned long int>& result)
 { multiReadPrimitive(xpath, result); }

void XmlReader::read(const std::string& xpath, multi1d<float>& result)
 { multiReadPrimitive(xpath, result); }

void XmlReader::read(const std::string& xpath, multi1d<double>& result)
 { multiReadPrimitive(xpath, result); }

void XmlReader::read(const std::string& xpath, multi1d<bool>& result)
 { multiReadPrimitive(xpath, result); }


void XmlReader::print(ostream& os)
{
 ostringstream newos;
 std::string s;
 if (active){
    BasicXPathReader::print(newos);
    s = newos.str();}
 if (global) QDPInternal::broadcast_str(s);
 os << s;
}
   
void XmlReader::printCurrentContext(ostream& os)
{
 ostringstream newos;
 std::string s;
 if (active){
    if (derived) BasicXPathReader::printChildren(newos);
    else BasicXPathReader::printRoot(newos);
    s = newos.str();}
 if (global) QDPInternal::broadcast_str(s);
 os << s;
}
   
std::string XmlReader::str()
{
 ostringstream newos;
 std::string s;
 if (active){
    BasicXPathReader::print(newos);
    s = newos.str();}
 if (global) QDPInternal::broadcast_str(s);
 return s;
}

// *********************************************************


XmlWriter::XmlWriter(bool global_mode) : global(global_mode)
 { active=(global_mode) ? Layout::primaryNode() : true; }

void XmlWriter::openTag(const string& tagname)
 { if (active) XMLSimpleWriter::openTag(tagname); }

void XmlWriter::push(const string& tagname)
 { if (active) XMLSimpleWriter::openTag(tagname); }

void XmlWriter::pop()
 { if (active) XMLSimpleWriter::closeTag(); }

void XmlWriter::closeTag()
 { if (active) XMLSimpleWriter::closeTag(); }

void XmlWriter::emptyTag(const string& tagname)
 { if (active) XMLSimpleWriter::emptyTag(tagname); }

void XmlWriter::write(const string& tagname)
 { if (active) XMLSimpleWriter::emptyTag(tagname); }


void XmlWriter::write(const string& tagname, const string& value)
 { writePrimitive<string>(tagname, value); }

void XmlWriter::write(const string& tagname, const char* value)
 { writePrimitive<string>(tagname, string(value)); }

void XmlWriter::write(const string& tagname, const int& value)
 { writePrimitive<int>(tagname, value); }

void XmlWriter::write(const string& tagname, const unsigned int& value)
 { writePrimitive<unsigned int>(tagname, value); }

void XmlWriter::write(const string& tagname, const short int& value)
 { writePrimitive<short int>(tagname, value); }

void XmlWriter::write(const string& tagname, const unsigned short int& value)
 { writePrimitive<unsigned short int>(tagname, value); }

void XmlWriter::write(const string& tagname, const long int& value)
 { writePrimitive<long int>(tagname, value); }

void XmlWriter::write(const string& tagname, const unsigned long int& value)
 { writePrimitive<unsigned long int>(tagname, value); }

void XmlWriter::write(const string& tagname, const float& value)
 { writePrimitive<float>(tagname, value); }

void XmlWriter::write(const string& tagname, const double& value)
 { writePrimitive<double>(tagname, value); }

void XmlWriter::write(const string& tagname, const bool& value)
 { writePrimitive<bool>(tagname, value); }



void XmlWriter::write(const string& tagname, const multi1d<int>& values)
 { multiWritePrimitive<int>(tagname, values); }

void XmlWriter::write(const string& tagname, const multi1d<unsigned int>& values)
 { multiWritePrimitive<unsigned int>(tagname, values); }

void XmlWriter::write(const string& tagname, const multi1d<short int>& values)
 { multiWritePrimitive<short int>(tagname, values); }

void XmlWriter::write(const string& tagname, const multi1d<unsigned short int>& values)
 { multiWritePrimitive<unsigned short int>(tagname, values); }

void XmlWriter::write(const string& tagname, const multi1d<long int>& values)
 { multiWritePrimitive<long int>(tagname, values); }

void XmlWriter::write(const string& tagname, const multi1d<unsigned long int>& values)
 { multiWritePrimitive<unsigned long int>(tagname, values); }

void XmlWriter::write(const string& tagname, const multi1d<float>& values)
 { multiWritePrimitive<float>(tagname, values); }

void XmlWriter::write(const string& tagname, const multi1d<double>& values)
 { multiWritePrimitive<double>(tagname, values); }

void XmlWriter::write(const string& tagname, const multi1d<bool>& values)
 { multiWritePrimitive<bool>(tagname, values); }


void XmlWriter::writeXML(const string& xmlstr)
 { if (active) XMLSimpleWriter::writeXML(xmlstr); }


void push(XmlWriter& xml, const string& tagname) 
 { xml.openTag(tagname); }

void pop(XmlWriter& xml) 
 { xml.closeTag(); }

void write(XmlWriter& xml, const string& tagname) 
 { xml.write(tagname); }


 // *******************************************************************

XmlBufferWriter::XmlBufferWriter(bool global_mode) : XmlWriter(global_mode)
 {indent_level=0;}

XmlBufferWriter::XmlBufferWriter(const std::string& s, bool global_mode) 
                : XmlWriter(global_mode)
 {open(s);}

void XmlBufferWriter::open(const std::string& s) 
 {if (active){ output_stream.clear(); output_stream.str(s);} }

string XmlBufferWriter::str() const
{
 string s;
 if (active){
    s=output_stream.str();}
 if (global) QDPInternal::broadcast_str(s);
 return s;
}

XmlBufferWriter::~XmlBufferWriter() {}


// **************************************************************


    // This returns the number of times that the tag "tagname"
    // is found in the descendents of the current context or the current context.
    // A **tag name** should be input, not an Xpath, since
    // a "./descendant-or-self::" is prepended to form an Xpath.

int xml_tag_count(XmlReader& xmlr, const string& tagname)
{
 return xmlr.count("./descendant-or-self::"+tagname);
}

int xml_tag_count(XMLReader& xmlr, const string& tagname)
{
 return xmlr.count("./descendant-or-self::"+tagname);
}

void xml_tag_assert(XmlReader& xmlr, const std::string& tagname)
{
 if (xml_tag_count(xmlr,tagname)!=1){
    if (xmlr.isActive()){
       std::cerr << "Bad XML input to "<<tagname<<std::endl;
       std::cerr << "Expected one <"<<tagname<<"> tag"<<std::endl;}
    xmlreadfail(xmlr,tagname);}
}

void xml_tag_assert(XmlReader& xmlr, const std::string& tagname,
                    const std::string& infoname)
{
 if (xml_tag_count(xmlr,tagname)!=1){
    if (xmlr.isActive()){
       std::cerr << "Bad XML input to "<<infoname<<std::endl;
       std::cerr << "Expected one <"<<tagname<<"> tag"<<std::endl;}
    xmlreadfail(xmlr,infoname);}
}


void xml_cerr(XmlReader& xmlr, const std::string& msg)
{
 if (xmlr.isActive()) std::cerr << msg << std::endl;
}


 // ************************************************************

   // outputs the current context for a failed xml read
   // (typically in an Info constructor)

void xmlreadfail(XmlReader& xmlr, const std::string& infoname)
{
 ostringstream tempout;
 xmlr.printCurrentContext(tempout);
 if (xmlr.isActive()){
    std::cerr << "*****ERROR**** "<<infoname<<" constructor failed"
              << " on XML content:"<<std::endl;
    std::cerr << tempout.str()<<std::endl;}
 QDP_abort(1);
}

void xmlreadfail(XMLReader& xmlr, const std::string& infoname)
{
 QDPIO::cerr << "*****ERROR**** "<<infoname<<" constructor failed"
             << " on XML content:"<<endl;
 ostringstream tempout;
 xmlr.printCurrentContext(tempout);
 QDPIO::cerr << tempout.str()<<endl;
 QDP_abort(1);
}


 // ************************************************************

  // removes tabs, newline, linefeed characters, then trims
  // leading and trailing blanks.

string tidyString(const string& str)   
{
 if (str.empty()) return "";
 string tmp;
 for (uint i=0;i<str.length();i++)
    if ((str[i]!='\n')&&(str[i]!='\t')&&(str[i]!='\r'))
       tmp.push_back(str[i]);
 uint start=tmp.find_first_not_of(" ");
 if (start==string::npos) return "";
 uint len=tmp.find_last_not_of(" ")-start+1;
 return tmp.substr(start,len);
}

// *************************************************************

  // first tidies the string in "str", then removes any
  // leading path/subdirectory information in the file name

string tidyFileName(const string& str)   
{
 string tmp(tidyString(str));
 uint pos=tmp.rfind('/');
 if (pos==string::npos) return tmp;
 else return tmp.substr(pos+1,tmp.length()-pos);
}

// *************************************************************

  // converts an integer to a string

string int_to_string(int intval)
{ 
 ostringstream oss;
 oss << intval;
 return oss.str();
}


// *************************************************************


bool xmlContentIsEqual(const string& doc1, const string& doc2, 
                       float float_rel_tol)
{
 bool result;
 XMLContentComparer temp;
 result=temp.IsEqual(doc1,doc2,float_rel_tol);
 return result;  
}

bool xmlContentIsEqual(XmlReader& xmlr1, XmlReader& xmlr2, 
                       float float_rel_tol)
{
 stringstream oss1,oss2;
 xmlr1.print(oss1);
 xmlr2.print(oss2);
 return xmlContentIsEqual(oss1.str(),oss2.str(),float_rel_tol);
}


// ************************************************************

    //  tests if file having name "file_name" exists
    //  on the primary node or not; broadcasts results
    //  to all nodes

bool fileExists(const std::string& file_name, bool global_mode)
{
 bool result=false;
 if ((Layout::primaryNode())||(!global_mode)){
    result = (access(file_name.c_str(),F_OK) == 0) ? true : false;}
 if (global_mode) QDPInternal::broadcast(result);
 return result;
}


// ************************************************************



bool XML_NodePtr::operator<(const XML_NodePtr& rhs)
{
 return (m_ctl->node_cmp(m_ptr,rhs.m_ptr)<0);
}

bool XML_AttrPtr::operator<(const XML_AttrPtr& rhs)
{
 return (m_ctl->attr_cmp(m_ptr,rhs.m_ptr)<0);
}


bool XMLContentComparer::IsEqual(const string& doc1, const string& doc2, 
                                 float float_rel_tol)
{
 m_float_rel_tol=abs(float_rel_tol);
 if (m_float_rel_tol<1e-12) m_float_rel_tol=1e-12;

    // Parse XML document from a character string

 xmlDocPtr xml1 = xmlParseMemory(doc1.c_str(),doc1.size());
 if (xml1 == 0){
    QDPIO::cerr << "Error in XMLContentComparer::IsEqual"<<endl;
    QDPIO::cerr << "   ...could not parse first XML document"<<endl;
    xmlFreeDoc(xml1);
    return false;}

 xmlDocPtr xml2 = xmlParseMemory(doc2.c_str(),doc2.size());
 if (xml2 == 0){
    QDPIO::cerr << "Error in XMLContentComparer::IsEqual"<<endl;
    QDPIO::cerr << "   ...could not parse second XML document"<<endl;
    xmlFreeDoc(xml1);
    xmlFreeDoc(xml2);
    return false;}

    //Get the root element nodes
 xmlNodePtr root1 = xmlDocGetRootElement(xml1);
 xmlNodePtr root2 = xmlDocGetRootElement(xml2);

 int result=sibling_node_cmp(root1,root2);

 xmlFreeDoc(xml1);
 xmlFreeDoc(xml2);

 return (result==0);
}


  // Beginning with character at "start", moves the "character"
  // pointer to the first non-ignorable character (or NULL if
  // no such characters found).  "nchar" will contain the 
  // number of non-ignorable characters in the next token. 

void XMLContentComparer::get_next_token(const char*& start, int& nchar)
{
 char delimit[] = " \t\n\r";
 start += strspn(start, delimit); 
 if (*start=='\0'){ start=0; nchar=0; return;}
 const char *stop = strpbrk(start, delimit);
 if (stop==0) nchar=strlen(start);
 else{
    nchar=0; const char *tmp=start;
    while (tmp!=stop){ tmp++; nchar++;}}
}

   //  Generalization of strcmp but not using null characters
   //  as string ending; uses numbers of characters in na,nb
   //  Return 0 if na-character string in "a" equals nb-character
   //  string in "b".  Negative return if a-string precedes b-string,
   //  and >0 return if a-string comes after b-string.

int XMLContentComparer::token_strcmp(const char *a, int na, 
                                     const char *b, int nb)
{
 if (na==nb) return strncmp(a,b,na);
 else if (na>nb){
    int k=strncmp(a,b,nb);
    if (k!=0) return k;
    else return 1;}
 else{
    int k=strncmp(a,b,na);
    if (k!=0) return k;
    else return -1;}
}

   // Check if the "nchar"-character token starting at "a"
   // is a valid integer; if so, return integer value in "ivalue".
   // We assume the token has no blank characters.

bool XMLContentComparer::is_integer_token(const char *const a, int nchar, 
                                          int& ivalue)
{
 if (nchar==0) return false;
 const char* ap=a;
 int k=0;
 if ((*ap=='+')||(*ap=='-')){ap++; k++;}
 if (k==nchar) return false;
 while (k<nchar){
    if (!isdigit(*ap)) return false;
    ap++; k++;}
 stringstream oss; oss << "%" << nchar << "d";
 return sscanf(a,oss.str().c_str(),&ivalue);
}

   // Check if the "nchar"-character token starting at "a"
   // is a valid floating-point number; if so, return value in "fvalue".
   // We assume the token has no blank characters.

bool XMLContentComparer::is_float_token(const char *const a, int nchar,
                                        float& fvalue)
{
 if (nchar==0) return false;
 const char* ap=a;
 int k=0;
 if ((*ap=='+')||(*ap=='-')){ap++; k++;}
 if (k==nchar) return false;
 bool dotflag=false, eflag=false;
 while (k<(nchar-1)){
    if (*ap=='.'){
       if (dotflag) return false;
       dotflag=true;}
    else if ((*ap=='e')||(*ap=='E')){
       if (eflag) return false;
       dotflag=true;
       eflag=true;  ap++; k++;
       if ((*ap!='+')&&(*ap!='-')){ ap--; k--;}
       else if (k==nchar) return false;}
    else if (!isdigit(*ap)) return false;
    ap++; k++;
    }
 if (!((isdigit(*ap))||((*ap=='.')&&(!dotflag)))) return false; 
     // last character must be digit or . (if not dots yet)
 stringstream oss; oss << "%" << nchar << "g";
 return sscanf(a,oss.str().c_str(),&fvalue);
}

    // This routine compares the content in the token strings "a" and "b"
    // having "na" and "nb" number of characters.  First the token are
    // compared as integers, then as floats (to within "float_rel_tol"),
    // then finally as strings.  A zero value is returned if they are
    // equal, a number > 0 if a>b and a number <0 is a<b.

int XMLContentComparer::token_content_cmp(const char *a, int na, 
                                          const char *b, int nb)
{
 int ia,ib;
 float fa,fb;
 if (is_integer_token(a,na,ia)&&is_integer_token(b,nb,ib)){
    if (ia>ib) return 1;
    else if (ia==ib) return 0;
    else return -1;}
 else if (is_float_token(a,na,fa) && is_float_token(b,nb,fb)){
    float favg=0.5*abs(fa+fb);
    if ((fa-fb)>(favg*m_float_rel_tol)) return 1;
    else if ((fa-fb)<-(favg*m_float_rel_tol)) return -1;
    else return 0;}
 return token_strcmp(a,na,b,nb);
}


int XMLContentComparer::content_cmp(const char *a, const char *b)
{
 int result;
 const char *ap=a; int na;
 const char *bp=b; int nb;
 get_next_token(ap,na);
 get_next_token(bp,nb);
 while ((ap!=0)&&(bp!=0)){
    result=token_content_cmp(ap,na,bp,nb);
    if (result!=0) return result;
    ap+=na; get_next_token(ap,na);
    bp+=nb; get_next_token(bp,nb);
    }
 if ((ap==0)&&(bp==0)) return 0;
 if (ap!=0) return 1;
 else return -1;
}


int XMLContentComparer::attr_cmp(xmlAttrPtr a_attr, 
                                 xmlAttrPtr b_attr)
{
    // first compare by attribute name
 int result = strcmp((char*) a_attr->name, (char*) b_attr->name);
 if (result!=0) return result;

    // compare by value
 if  ((!xmlNodeIsText(a_attr->children))
    ||(!xmlNodeIsText(b_attr->children))) QDP_abort(1);  // quick check
 return strcmp((char*) a_attr->children->content,
               (char*) b_attr->children->content);
}


int XMLContentComparer::attrlist_cmp(xmlAttrPtr a_attr, 
                                     xmlAttrPtr b_attr)
{
 xmlAttrPtr cur_attr;

 list<XML_AttrPtr> alist;
 for (cur_attr=a_attr; cur_attr; cur_attr=cur_attr->next)
    alist.push_back(XML_AttrPtr(cur_attr,this));
 alist.sort();

 list<XML_AttrPtr> blist;
 for (cur_attr=b_attr; cur_attr; cur_attr=cur_attr->next)
    blist.push_back(XML_AttrPtr(cur_attr,this));
 blist.sort();

 int na=alist.size();
 int nb=blist.size();
 int n = (na>=nb)?nb:na;

 list<XML_AttrPtr>::const_iterator at=alist.begin();
 list<XML_AttrPtr>::const_iterator bt=blist.begin();

 for (int k=0;k<n;k++){
    int result=attr_cmp(at->m_ptr,bt->m_ptr);
    if (result!=0) return result;
    at++; bt++;}
 if ((at==alist.end())&&(bt==blist.end())) return 0;
 if (at==alist.end()) return -1;
 else return 1;
}

int XMLContentComparer::node_cmp(xmlNodePtr a_node,
                                 xmlNodePtr b_node)
{
    // compare node types first
 if (a_node->type < b_node->type) return -1;
 else if (a_node->type > b_node->type) return 1;

    // compare node names next
 int result=strcmp((char*) a_node->name, (char*) b_node->name);
 if (result!=0) return result;

    // compare by attributes thirdly
 result=attrlist_cmp(a_node->properties,b_node->properties);
 if (result!=0) return result;

    // if text node, compare by textual content
 if (xmlNodeIsText(a_node)){
    if (!xmlNodeIsText(b_node)) QDP_abort(1); // something went wrong
    return content_cmp((char*) a_node->content, (char*) b_node->content);}

    // now compare by children
 return sibling_node_cmp(a_node->children,b_node->children);
}


int XMLContentComparer::sibling_node_cmp(xmlNodePtr a_node,
                                         xmlNodePtr b_node)
{
 xmlNodePtr cur_node;

 list<XML_NodePtr> alist;
 for (cur_node=a_node; cur_node; cur_node=cur_node->next)
    alist.push_back(XML_NodePtr(cur_node,this));
 alist.sort();

 list<XML_NodePtr> blist;
 for (cur_node=b_node; cur_node; cur_node=cur_node->next)
    blist.push_back(XML_NodePtr(cur_node,this));
 blist.sort();

 int na=alist.size();
 int nb=blist.size();
 int n = (na>=nb)?nb:na;

 list<XML_NodePtr>::const_iterator at=alist.begin();
 list<XML_NodePtr>::const_iterator bt=blist.begin();

 for (int k=0;k<n;k++){
    int result=node_cmp(at->m_ptr,bt->m_ptr);
    if (result!=0) return result;
    at++; bt++;}
 if ((at==alist.end())&&(bt==blist.end())) return 0;
 if (at==alist.end()) return -1;
 else return 1;
}

// **********************************************************
  }
}
