//Arjun Singh Gambhir
//Real authors: Ben Hoerz, John Bulava, Colin Morningstar, ??
//Code directly ported from chroma_laph, more acknowledgements coming later.
#ifndef XML_HELP_H
#define XML_HELP_H

#include "qdp.h"
#include "chromabase.h"
#include <libxml/tree.h>
#include <string>

namespace Chroma {
  namespace LaphEnv {

    //  Objects of classes XmlReader, XmlWriter, XmlBufferWriter can
    //  operate in global mode or local mode, unlike objects of the
    //  QDP++ classes XMLReader, XMLWriter, XMLBufferWriter.
    //  In global mode, input is done by the primary node and broadcast
    //  to all nodes, and output is just done by the primary node.
    //  In local mode, each node can separately do input and output,
    //  with no communication with other nodes.  The reader throws an
    //  exception if a failure occurs.


class XmlBufferWriter;

// *********************************************************

class XmlReader : protected XMLXPathReader::BasicXPathReader
{

    bool derived; // is this reader derived from another reader?
    bool active;  // does this node participate?
    bool global;  // broadcast results globally

         // disallow copies
    XmlReader& operator=(const XmlReader&);
    XmlReader(const XmlReader&);

  public:

    XmlReader(std::istream& is, bool global_mode=true);
    XmlReader(std::string& s, bool global_mode=true);
    XmlReader(const XmlBufferWriter& mw, bool global_mode=true);
    XmlReader(XMLReader& old, bool global_mode=true);
    XmlReader(const XMLBufferWriter& mw, bool global_mode=true);
    XmlReader(XmlReader& start, const std::string& xpath);
    ~XmlReader();
    
    bool isGlobal() const {return global;}
    bool isLocal() const {return !global;}
    bool isActive() const {return active;}

    int count(const std::string& xpath);

        // throw exception if failure
    void read(const std::string& xpath, std::string& input);
    void read(const std::string& xpath, int& input);
    void read(const std::string& xpath, unsigned int& input);
    void read(const std::string& xpath, short int& input);
    void read(const std::string& xpath, unsigned short int& input);
    void read(const std::string& xpath, long int& input);
    void read(const std::string& xpath, unsigned long int& input);
    void read(const std::string& xpath, float& input);
    void read(const std::string& xpath, double& input);
    void read(const std::string& xpath, bool& input);

    void read(const std::string& xpath, multi1d<int>& input);
    void read(const std::string& xpath, multi1d<unsigned int>& input);
    void read(const std::string& xpath, multi1d<short int>& input);
    void read(const std::string& xpath, multi1d<unsigned short int>& input);
    void read(const std::string& xpath, multi1d<long int>& input);
    void read(const std::string& xpath, multi1d<unsigned long int>& input);
    void read(const std::string& xpath, multi1d<float>& input);
    void read(const std::string& xpath, multi1d<double>& input);
    void read(const std::string& xpath, multi1d<bool>& input);
    
    void print(std::ostream& is);
    void printCurrentContext(std::ostream& is);
    std::string str();

  protected:

    template <typename T>
    void readPrimitive(const std::string& xpath, T& input)
    { if (active) BasicXPathReader::get(xpath, input);
      if (global) QDPInternal::broadcast(input); }

    template <typename T>
    void multiReadPrimitive(const std::string& xpath, multi1d<T>& input);

};


template <typename T>
void read(XmlReader& xml, const std::string& xpath, T& value)
{
 xml.read(xpath,value);
}

template <typename T>
void read(XmlReader& xml, const std::string& xpath, multi1d<T>& values)
{
 xml.read(xpath,values);
}

// *********************************************************

   //  In global mode, results are written only by the primary
   //  node.  In local mode, each node writes the results.

class XmlWriter : protected XMLWriterAPI::XMLSimpleWriter
{

  public:

    XmlWriter(bool global_mode=true);
    virtual ~XmlWriter() {}

    bool isGlobal() const {return global;}
    bool isLocal() const {return !global;}

    void push(const std::string& tagname);
    void openTag(const std::string& tagname);
    void pop();
    void closeTag();

    void emptyTag(const std::string& tagname);
    void write(const std::string& tagname);

    void write(const std::string& tagname, const std::string& value);
    void write(const std::string& tagname, const char* value);
    void write(const std::string& tagname, const int& value);
    void write(const std::string& tagname, const unsigned int& value);
    void write(const std::string& tagname, const short int& value);
    void write(const std::string& tagname, const unsigned short int& value);
    void write(const std::string& tagname, const long int& value);
    void write(const std::string& tagname, const unsigned long int& value);
    void write(const std::string& tagname, const float& value);
    void write(const std::string& tagname, const double& value);
    void write(const std::string& tagname, const bool& value);

    void write(const std::string& tagname, const multi1d<int>& values);
    void write(const std::string& tagname, const multi1d<unsigned int>& values);
    void write(const std::string& tagname, const multi1d<short int>& values);
    void write(const std::string& tagname, const multi1d<unsigned short int>& values);
    void write(const std::string& tagname, const multi1d<long int>& values);
    void write(const std::string& tagname, const multi1d<unsigned long int>& values);
    void write(const std::string& tagname, const multi1d<float>& values);
    void write(const std::string& tagname, const multi1d<double>& values);
    void write(const std::string& tagname, const multi1d<bool>& values);
    
    void writeXML(const std::string& xml);

  protected:

    bool active;
    bool global;

    template<typename T>
    void writePrimitive(const std::string& tagname, const T& value)
     { if (active){
         XMLSimpleWriter::openTag(tagname);
         XMLSimpleWriter::write(value);
         XMLSimpleWriter::closeTag();} }


    template<typename T>
    void multiWritePrimitive(const std::string& tagname, const multi1d<T>& values)
     { if (active){
         XMLSimpleWriter::openTag(tagname);
         if (values.size()>0){
            XMLSimpleWriter::write(values[0]);
            std::string sp(" ");
            for (int k=1;k<values.size();k++){
               XMLSimpleWriter::write(sp);
               XMLSimpleWriter::write(values[k]);}}
         XMLSimpleWriter::closeTag();} }

};





void push(XmlWriter& xml, const std::string& tagname);

void pop(XmlWriter& xml);

void write(XmlWriter& xml, const std::string& tagname);


template <typename T>
void write(XmlWriter& xml, const std::string& tagname, const T& value)
{
 xml.write(tagname,value);
}

template <typename T>
void write(XmlWriter& xml, const std::string& tagname, const multi1d<T>& values)
{
 xml.write(tagname,values);
}



 // ****************************************************************

class XmlBufferWriter : public XmlWriter
{
  public:

    XmlBufferWriter(bool global_mode=true);
    explicit XmlBufferWriter(const std::string& s, bool global_mode=true);
    ~XmlBufferWriter();
    void open(const std::string& s);
    std::string str() const;
    void flush() {output_stream.flush();}
    bool fail() const {return output_stream.fail();}

  private:

  std::ostringstream output_stream;
  std::ostream& getOstream() {return output_stream;}

  friend class XmlReader;

};


// **************************************************************

    // This routine calls an xml read inside a try block
    // and outputs an informative message if the read fails.
    // The added parameter is a string which should be the
    // name of the class which called this function.
    // Note that a **tag name** should be used. A "./descendant-or-self::"
    // is prepended to form an appropriate Xpath.

template <typename T>
void xmlread(XmlReader& xmlr, const std::string& tagname, T& val,
             const std::string& callingClass)
{
 try{
    read(xmlr,"./descendant-or-self::"+tagname,val);}
 catch(const std::string& err_msg){
    if (xmlr.isActive())
       std::cerr <<"Invalid read of "<<tagname<<" in "<<callingClass<<std::endl;
    throw(std::string("error"));}
} 

template <typename T>
void xmlread(XMLReader& xmlr, const std::string& tagname, T& val,
             const std::string& callingClass)
{
 try{
    read(xmlr,"./descendant-or-self::"+tagname,val);}
 catch(const std::string& err_msg){
    QDPIO::cerr <<"Invalid read of "<<tagname<<" in "<<callingClass<<std::endl;
    throw(std::string("error"));}
}

// *********************************************************

    // This returns the number of times that the tag "tagname"
    // is found in the descendents of the current context.
    // A **tag name** should be input, not an Xpath, since
    // a "./descendant-or-self::" is prepended to form an Xpath.

int xml_tag_count(XmlReader& xmlr, const std::string& tagname);
int xml_tag_count(XMLReader& xmlr, const std::string& tagname);

void xml_tag_assert(XmlReader& xmlr, const std::string& tagname);
void xml_tag_assert(XmlReader& xmlr, const std::string& tagname,
                    const std::string& infoname);

void xml_cerr(XmlReader& xmlr, const std::string& msg);


// *********************************************************

template <typename T>
bool xmlreadif(XmlReader& xmlr, const std::string& tagname, T& val,
               const std::string& callingClass)
{
 if (xml_tag_count(xmlr,tagname)==1){
    xmlread(xmlr,tagname,val,callingClass);
    return true;}
 return false;
}

template <typename T>
bool xmlreadif(XMLReader& xmlr, const std::string& tagname, T& val,
               const std::string& callingClass)
{
 if (xml_tag_count(xmlr,tagname)==1){
    xmlread(xmlr,tagname,val,callingClass);
    return true;}
 return false;
}

// *********************************************************

template <typename T>
void assertEqual(const T& obj1, const T& obj2, const std::string& callingClass)
{
 try{
    obj1.checkEqual(obj2);}
 catch(const std::string& err_msg){
    std::cerr << err_msg <<" in "<<callingClass<<std::endl;
    QDP_abort(1);}
}
 
// *********************************************************

   // outputs the current context for a failed xml read
   // (typically in an Info constructor)

void xmlreadfail(XmlReader& xmlr, const std::string& infoname);
void xmlreadfail(XMLReader& xmlr, const std::string& infoname);


// *********************************************************

  // Removes tabs and newline characters, then trims
  // leading and trailing blanks.

std::string tidyString(const std::string& str);    


  // first tidies the string in "str", then removes any
  // leading path/subdirectory information in the file name

std::string tidyFileName(const std::string& str);


  // converts an integer to a string

std::string int_to_string(int ival);


// *********************************************************


  // Compares the XML content in "doc1" and "doc2" and returns
  // "true" if they are the same.  In doing the comparison, the
  // order of sibling nodes is irrelevant, and textual content
  // is compared token by token.  Integer tokens are compared as
  // integers, and floating-point tokens are compared as 
  // floats.  If the difference between floats is less than
  // "float_rel_tol", they are considered the same.


bool xmlContentIsEqual(const std::string& doc1, const std::string& doc2, 
                       float float_rel_tol = 1e-6);


  // Same as above, but applied to XmlReaders "xmlr1" and "xmlr2".

bool xmlContentIsEqual(XmlReader& xmlr1, XmlReader& xmlr2, 
                       float float_rel_tol = 1e-6);


// *********************************************************

    //  tests if file having name "file_name" exists
    //  on the primary node or not; broadcasts results
    //  to all nodes

bool fileExists(const std::string& file_name, bool global_mode=true);



// *********************************************************


  // declarations below are NOT available to the end user;
  // these classes are used by xmlContentIsEqual

class XMLContentComparer;  // forward declaration


class XML_NodePtr {

    xmlNodePtr m_ptr;
    XMLContentComparer *m_ctl;

    XML_NodePtr() {}
    XML_NodePtr(xmlNodePtr ptr, XMLContentComparer *ctl) 
            : m_ptr(ptr), m_ctl(ctl) {}

 public:

    XML_NodePtr(const XML_NodePtr& in) 
            : m_ptr(in.m_ptr), m_ctl(in.m_ctl)  {}

    ~XML_NodePtr() {}

    friend class XMLContentComparer;

    bool operator<(const XML_NodePtr& rhs);

};

class XML_AttrPtr {

    xmlAttrPtr m_ptr;
    XMLContentComparer *m_ctl;

    XML_AttrPtr() {}
    XML_AttrPtr(xmlAttrPtr ptr, XMLContentComparer *ctl) 
            : m_ptr(ptr), m_ctl(ctl) {}

 public:

    XML_AttrPtr(const XML_AttrPtr& in) 
            : m_ptr(in.m_ptr), m_ctl(in.m_ctl)  {}

    ~XML_AttrPtr() {}

    friend class XMLContentComparer;

    bool operator<(const XML_AttrPtr& rhs);

};

class XMLContentComparer {

    float m_float_rel_tol;

       // for use only by xmlContentIsEqual
    XMLContentComparer() {}
    XMLContentComparer(const XMLContentComparer& in);
    XMLContentComparer& operator=(const XMLContentComparer& in);


 public:

    ~XMLContentComparer() {}

 private:

    friend bool xmlContentIsEqual(const std::string& doc1, const std::string& doc2, 
                                  float float_rel_tol);
    friend class XML_NodePtr;
    friend class XML_AttrPtr;

    bool IsEqual(const std::string& doc1, const std::string& doc2, 
                 float float_rel_tol);

    void get_next_token(const char*& start, int& nchar);

    int token_strcmp(const char *a, int na, const char *b, int nb);

    bool is_integer_token(const char *const a, int nchar, int& ivalue);

    bool is_float_token(const char *const a, int nchar, float& fvalue);

    int token_content_cmp(const char *a, int na, const char *b, int nb);

    int content_cmp(const char *a, const char *b);

    int attr_cmp(xmlAttrPtr a_attr, xmlAttrPtr b_attr);

    int attrlist_cmp(xmlAttrPtr a_attr, xmlAttrPtr b_attr);

    int sibling_node_cmp(xmlNodePtr a_node, xmlNodePtr b_node);

    int node_cmp(xmlNodePtr a_node, xmlNodePtr b_node);


};



// *********************************************************
  }
}
#endif
