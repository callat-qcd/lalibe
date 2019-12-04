//Arjun Singh Gambhir
//Real authors: Ben Hoerz, John Bulava, Colin Morningstar, ??
//Code directly ported from chroma_laph, more acknowledgements coming later.
#ifndef LAPH_DATA_IO_HANDLER_H
#define LAPH_DATA_IO_HANDLER_H

#include "qdp.h"
#include "chromabase.h"
#include "io_map.h"
#include "filelist_info.h"
#include "xml_help.h"
#include <set>

namespace Chroma {
  namespace LaphEnv {


 // *********************************************************************************
 // *                                                                               *
 // *   The classes defined in this file that are important are                     *
 // *                                                                               *
 // *         DataPutHandlerMF   (multi-file)                                       *
 // *         DataGetHandlerMF   (multi-file)                                       *
 // *         DataPutHandlerSF   (single file)                                      *
 // *         DataGetHandlerSF   (single file)                                      *
 // *                                                                               *
 // *   and                                                                         *
 // *                                                                               *
 // *         DataPutHandlerMemMF                                                   *
 // *         DataGetHandlerMemMF                                                   *
 // *                                                                               *
 // *                                                                               *
 // *   The "DataPutHandlerMF" and "DataGetHandlerMF" (multi-file) classes handle   *
 // *   inserting data into IOMap records in files specified in a FileListInfo and  *
 // *   subsequently reading the data.  The file key determines which file to       *
 // *   access, and the record key determines which IOMap record in a given file    *
 // *   to access.  The file key and record key are combined into a storage key     *
 // *   which is used in an internal map to store the data in memory.  Do not use   *
 // *   the "put" and "get" classes simultaneously;  use "put" to build up the data *
 // *   in the files, then destroy the insert object; after this, a "get" object    *
 // *   can be created to access the  data in the files.  Objects of these classes  *
 // *   can be created globally (all nodes) or locally (one node).                  *
 // *                                                                               *
 // *   Objects of these "put" handlers always assume an "updating" mode. Existing  *
 // *   files are never erased, and new files are created as needed.  New records   *
 // *   are added to the files.  If the key of a record to be put already exists    *
 // *   in a file, the put will only occur if "overwrite" is specified AND the      *
 // *   size of the data to be put does not exceed the size of the data already in  *
 // *   the file for that key.  For the multi-file handler, the "overwrite" bool    *
 // *   must be given in the FileListInfo; for the single-file handler, the         *
 // *   constructor takes an explicit bool parameter.                               *
 // *                                                                               *
 // *   Data*HandlerMemMF provide the same interface as do Data*HandlerMF and       *
 // *   emulate all the aforementioned functionality, but all data is stored in and *
 // *   read from memory only, i.e. the data is stored in Chroma's NamedObjectMap.  *
 // *   This data model is different from the standard approach adopted in the rest *
 // *   of this suite in that it allows one to chain several tasks within one run   *
 // *   without the need for intermittent IO.                                       *
 // *   A possible use case is the computation of meson and baryon line ends        *
 // *   where the needed Laplace eigenvectors on each time slice are only computed  *
 // *   on the fly to avoid expensive IO between the tasks.                         *
 // *                                                                               *
 // *   Every handler that wants to support both the file and memory data model can *
 // *   store a pointer to the base class Data(Get/Put)HandlerBaseMF and            *
 // *   instantiate an object of the corresponding child class.                     *
 // *                                                                               *
 // *   To use these classes, one needs the following ingredients:                  *
 // *                                                                               *
 // *    (a) a FileListInfo object and an IOHandler file ID string                  *
 // *                                                                               *
 // *    (b) a handler class "H" that has members                                   *
 // *                                                                               *
 // *            bool checkHeader(XmlReader& xmlr, int suffix)                      *
 // *                                                                               *
 // *          that checks that the header information in the file is good,         *
 // *          returning a boolean value, and                                       *
 // *                                                                               *
 // *            void writeHeader(XmlWriter&,const F&,int)                          *
 // *                                                                               *
 // *        where "F" is the file key type, that writes out the header string      *
 // *        for each file                                                          *
 // *                                                                               *
 // *    (c) both classes need to extract the file key from the header string, so   *
 // *          the file key type must have a constructor that takes only an         *
 // *          XmlReader and an output(XmlWriter&) member; you need to specify the  *
 // *          header tag that the file key will be enclosed in; the file key must  *
 // *          also have a "<" operator and an "==" operator defined                *
 // *                                                                               *
 // *    (d) the record key class must have all of the features of an IOMap key:    *
 // *                                                                               *
 // *       -- since used in a C++ map, a less than operator must be define         *
 // *              const K& K::operator<(const K& rhs);                             *
 // *       -- a numbytes(ioh,K) function, where ioh is an IOHandler, must be       *
 // *           defined to give the number of bytes each key occupies in an         *
 // *           IOHandler file                                                      *
 // *       -- a copy constructor K(const K& in) must be defined                    *
 // *           (a default constructor is not needed)                               *
 // *       -- a multi_read(ioh, vector<K>&,n) must be defined to read n keys       *
 // *       -- a multi_write(ioh, const vector<K>&) must be defined                 *
 // *                                                                               *
 // *    (e) the value type must have the following features:                       *
 // *                                                                               *
 // *       -- a write(ioh, const V&) must be defined (ioh is an IOHandler object)  *
 // *       -- a read(ioh, V&) must be defined                                      *
 // *       -- a numbytes(ioh,V) must be defined giving number of bytes occupied    *
 // *            by V in an IOHandler file                                          *
 // *                                                                               *
 // *                                                                               *
 // *   Usage:  (inserting)                                                         *
 // *                                                                               *
 // *     FileListInfo files;     // specify overwrite mode in here                 *
 // *     string fid("id string");                                                  *
 // *     string HeaderTag("HeaderTag");                                            *
 // *     SomeHandler H;  //  R=record type, F=file type, D=data type               *
 // *     bool globalmode=true;  // false is local mode                             *
 // *                                                                               *
 // *     DataPutHandlerMF<SomeHandler,F,R,D>                                       *
 // *               DIH(H,files,fid,HeaderTag,globalmode);                          *
 // *                                                                               *
 // *        // best to insert data having same file key                            *
 // *        // in sequence to minimize file open/close                             *
 // *                                                                               *
 // *     F fkey; R rkey; D data;                                                   *
 // *     DIH.openFile(fkey);                                                       *
 // *     DIH.putData(rkey,data);                                                   *
 // *     R rkey2; D data2;                                                         *
 // *     DIH.putData(rkey2,data2);                                                 *
 // *     R rkey3; D data3;                                                         *
 // *     DIH.queryData(fkey,rkey3); // check if exists                             *
 // *                                                                               *
 // *        // other members                                                       *
 // *                                                                               *
 // *     DIH.getFileListInfo();                                                    *
 // *     DIH.setOverWrite();                                                       *
 // *     DIH.setNoOverWrite();                                                     *
 // *     DIH.close();    // manual close current file                              *
 // *                                                                               *
 // *                                                                               *
 // *   Usage:  (reading)                                                           *
 // *                                                                               *
 // *     FileListInfo files;                                                       *
 // *     string fid("id string");                                                  *
 // *     string HeaderTag("HeaderTag");                                            *
 // *     SomeHandler H;  //  R=record type, F=file type, D=data type               *
 // *     bool globalmode=true;  // false is local mode                             *
 // *                                                                               *
 // *     DataGetHandlerMF<SomeHandler,F,R,D>                                       *
 // *               DRH(H,files,fid,HeaderTag,globalmode);                          *
 // *                                                                               *
 // *     F fkey; R rkey; D data;                                                   *
 // *     const D& data=DRH.getData(fkey,rkey);                                     *
 // *     bool flag1=DRH.queryData(fkey,rkey);                                      *
 // *     bool flag2=DRH.queryData(fkey);                                           *
 // *                                                                               *
 // *     DRH.remove(fkey,rkey);                                                    *
 // *     DRH.clearData();                                                          *
 // *                                                                               *
 // *        // other members                                                       *
 // *                                                                               *
 // *     DRH.getFileListInfo();                                                    *
 // *     XmlBufferWriter xmlout;                                                   *
 // *     DRH.getFileMap(xmlout);                                                   *
 // *     list<pair<F,list<R> > > keys=DRH.getKeys();                               *
 // *     DRH.outputKeys(xmlout);                                                   *
 // *                                                                               *
 // *                                                                               *
 // *   The "DataPutHandlerSF" and "DataGetHandlerSF" classes are similar to their  *
 // *   "MF" (multi-file) counterpart but only one single file (SF) is handled.     *
 // *                                                                               *
 // *   The corresponding Data*HandlerMemSF using the memory data model have not    *
 // *   been implemented so far.                                                    *
 // *                                                                               *
 // *                                                                               *
 // *********************************************************************************
 

template <typename H, typename F, typename R, typename D>
class DataPutHandlerBaseMF
{

 public:

    virtual ~DataPutHandlerBaseMF() {}

    virtual void setOverWrite() {}

    virtual void setNoOverWrite() {}

    virtual const FileListInfo& getFileListInfo() const {}

    virtual bool isGlobal() const = 0;

    virtual bool isLocal() const = 0;


    virtual void open(const F& fkey) = 0;

    virtual void putData(const R& rkey, const D& data) = 0;  // insert into current file

    virtual void flush() = 0;   // flush current file

    virtual void close() = 0;   // close current file

             // calls openFile(fkey) first, then inserts
             
    virtual void putData(const F& fkey, const R& rkey, const D& data) = 0;

    virtual void flush(const F& fkey) = 0;

    virtual void close(const F& fkey) = 0;
    
    virtual void flushAll() = 0;

    virtual void closeAll() = 0;



    virtual bool queryData(const F& fkey, const R& rkey) = 0;

    virtual bool queryData(const R& rkey) = 0;  // query in current open file

    virtual bool queryFile(const F& fkey) = 0;



    virtual void getFileMap(XmlWriter& xmlout) const = 0;

    virtual std::map<int,F> getSuffixMap() const = 0;

    virtual std::set<F> getFileKeys() const = 0;

};

template <typename H, typename F, typename R, typename D>
class DataGetHandlerBaseMF
{

 public:

    virtual ~DataGetHandlerBaseMF() {};

    virtual const FileListInfo& getFileListInfo() const {};

    virtual bool isGlobal() const {return true;};

    virtual bool isLocal() const {return false;};


    virtual bool queryData(const F& fkey, const R& rkey) = 0;

    virtual bool queryFile(const F& fkey) = 0;


    virtual const D& getData(const F& fkey, const R& rkey) = 0;

    virtual void removeData(const F& fkey, const R& rkey) = 0;

    virtual void removeData(const F& fkey) = 0;

    virtual void clearData() = 0;


    virtual void getFileMap(XmlWriter& xmlout) const = 0;

    virtual std::map<int,F> getSuffixMap() const = 0;

    virtual std::set<F> getFileKeys() const = 0;

    virtual std::set<R> getKeys(const F& fkey) = 0;

    virtual void outputKeys(XmlWriter& xmlout) = 0;

};

 
   // **************************************************************
   // *                                                            *
   // *                      DataGetHandlerMF                      *
   // *                                                            *
   // **************************************************************


   // "get" class handles the internal storage in "m_storage"
   // which is a map of pointers.  The use of pointers is efficient
   // whenever each data structure is fairly large (100 elements or 
   // more) since only pointers get copied.  
   
   // "fileMap" is a map associating a file key to a FileMapValue, 
   // which contains a suffix and an IOMap pointer.  Upon construction, 
   // "fileMap" is assigned by opening each file one by one, reading the 
   // header string, and extracting the file key.  None of the files is
   // left open.  As data is accessed, the files are opened and left
   // open.  While a file is open, the IOMap keeps all of the record
   // keys in memory.
   
   // "queryData" and "getData" open files as needed, and leave them open.
   // "removeData" removes the data from memory, but the files are left open.
   // "clearData" wipes the data from memory and closes all files.


template <typename H, typename F, typename R, typename D>
class DataGetHandlerMF : public DataGetHandlerBaseMF<H,F,R,D>
{

    struct StorageKey
    {
      F fkey;
      R rkey;
      StorageKey(const F& in_fkey, const R& in_rkey) : fkey(in_fkey), rkey(in_rkey) {}
      StorageKey(const StorageKey& in)  : fkey(in.fkey), rkey(in.rkey) {}
      StorageKey& operator=(const StorageKey& in)
        {fkey=in.fkey; rkey=in.rkey; return *this;}
      bool operator<(const StorageKey& rhs) const
        {return ((fkey<rhs.fkey) || ((fkey==rhs.fkey)&&(rkey<rhs.rkey)));}
    };

    struct FileMapValue
    {
      int suffix;
      IOMap<R,D> *fptr;
      FileMapValue(int in_suff) : suffix(in_suff), fptr(0) {}
      ~FileMapValue() {delete fptr;}
    };

    typedef std::map<StorageKey,D*>   StorageMapType;
    typedef std::map<F,FileMapValue>  FileMapType;
 
    FileListInfo finfo;
    StorageMapType  m_storage;
    FileMapType  fileMap;
    H& handler;
    bool gmode;
    bool checksums;
    std::string fid;


 public:

    DataGetHandlerMF(H& in_handler, const FileListInfo& in_filelist,
                     const std::string& filetype_id, 
                     const std::string& header_tag, 
                     bool global_mode=true, bool use_checksums=false);

    ~DataGetHandlerMF() {clearData(); fileMap.clear();}

    const FileListInfo& getFileListInfo() const {return finfo;}

    bool isGlobal() const {return gmode;}

    bool isLocal() const {return !gmode;}


    bool queryData(const F& fkey, const R& rkey);

    bool queryFile(const F& fkey);


    const D& getData(const F& fkey, const R& rkey);

    void removeData(const F& fkey, const R& rkey);

    void removeData(const F& fkey);

    void clearData();


    void getFileMap(XmlWriter& xmlout) const;

    std::map<int,F> getSuffixMap() const;

    std::set<F> getFileKeys() const;

    std::set<R> getKeys(const F& fkey);

    void outputKeys(XmlWriter& xmlout);


 private:

    void fail(const F& fkey, const R& rkey);
    
    void fail(const F& fkey);

    void fail(const std::string& msg);
    
    IOMap<R,D>* get_file_ptr(const F& fkey);
    
    void open(FileMapValue& fmv);
    
    void close(FileMapValue& fmv);

          // disallow copies
    DataGetHandlerMF(const DataGetHandlerMF& in);
    DataGetHandlerMF& operator=(const DataGetHandlerMF& in);

};



   // Constructor checks that the information in the headers
   // of all existing files is consistent, then sets up the
   // file map.  Each files is opened (one by one), the header string
   // is read, and then the file is closed.

template <typename H, typename F, typename R, typename D>
DataGetHandlerMF<H,F,R,D>::DataGetHandlerMF(H& in_handler,
                                            const FileListInfo& in_filelist,
                                            const std::string& filetype_id,
                                            const std::string& header_tag,
                                            bool global_mode, bool use_checksums)
  :  finfo(in_filelist), handler(in_handler), gmode(global_mode), 
     checksums(use_checksums), fid(tidyString(filetype_id)) 
{
 for (int suffix=finfo.getMinFileNumber();
          suffix<=finfo.getMaxFileNumber();suffix++){

    std::string filename=finfo.getFileName(suffix);

         // open all existing files and check consistency of headers
    std::string headerxml;
    bool exists;
    {IOMap<R,D> iom(gmode);
     exists=iom.peekHeader(headerxml,filename,fid);}
    if (!exists) continue;

/*
    if (!fileExists(filename,gmode)) continue;

         // open all existing files and check consistency of headers
    string headerxml;
    {IOMap<R,D> iom(gmode);
     try { iom.openReadOnly(filename,fid,headerxml,checksums);}
     catch(...) {
        fail("could not open file "+filename+" for reading");}} */

    if ((!gmode)||(Layout::primaryNode())){
       XmlReader xmlr(headerxml,false);
       if (!handler.checkHeader(xmlr,suffix)){
          fail("Header string in file is\n"+headerxml+"header info in file "
                     +filename+" does not match info in current Handler\n\n"
                     +"...execution aborted...\n");}}

            // extract the file key from this file
    try{
       XmlReader xmlr(headerxml,gmode);
       XmlReader xmlf(xmlr,"./descendant-or-self::"+header_tag);
       F fkey(xmlf);
       typename FileMapType::iterator it=fileMap.find(fkey);
       if (it!=fileMap.end()){
          fail(std::string("duplicate keys in fileMap in current Handler\n")
               +" ... too confusing to continue\n file suffix "
               +int_to_string(suffix)+" and suffix "+int_to_string((it->second).suffix)
               +" have same file key\n");}
       fileMap.insert(std::make_pair(fkey, FileMapValue(suffix)));}
    catch(...){
       fail("Could not extract FileKey from file "+filename+"\n");}
    }
}


template <typename H, typename F, typename R, typename D>
IOMap<R,D>* DataGetHandlerMF<H,F,R,D>::get_file_ptr(const F& fkey)
{
 typename FileMapType::iterator it=fileMap.find(fkey);
 if (it==fileMap.end()) return 0;
 if (it->second.fptr==0) open(it->second);
 return it->second.fptr;
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMF<H,F,R,D>::open(FileMapValue& fmv)
{
 std::string filename=finfo.getFileName(fmv.suffix);
 try {
    fmv.fptr=new IOMap<R,D>(gmode);
    fmv.fptr->openReadOnly(filename,fid,checksums);}
 catch(...){
    fail("failure opening file "+filename+" in DataGetHandlerMF");}
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMF<H,F,R,D>::close(FileMapValue& fmv)
{
 delete fmv.fptr; 
 fmv.fptr=0;
}


template <typename H, typename F, typename R, typename D>
const D& DataGetHandlerMF<H,F,R,D>::getData(const F& fkey, const R& rkey)
{
 StorageKey skey(fkey,rkey);
 typename StorageMapType::const_iterator dt=m_storage.find(skey);
 if (dt!=m_storage.end()) return *(dt->second);
 IOMap<R,D> *fptr=get_file_ptr(fkey);
 if (fptr==0) fail(fkey);
 D *result(new D);
 m_storage[skey]=result;
 try {fptr->get(rkey,*result);}
 catch(...){fail(fkey,rkey);}
 return *result;
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMF<H,F,R,D>::fail(const F& fkey, const R& rkey)
{
 if ((Layout::primaryNode())||(!gmode)){
  std::cerr << "DataGetHandlerMF could not find requested record:"<<std::endl;
  std::cerr << " File stub: "<< finfo.getFileStub()<<std::endl;
  XmlBufferWriter xmlout(gmode);
  push(xmlout,"FileRecordKey");
  fkey.output(xmlout);
  rkey.output(xmlout);
  pop(xmlout);
  std::cerr << xmlout.str()<<std::endl;
  std::cerr << "...execution aborted..."<<std::endl;}
 clearData(); fileMap.clear();
 QDP_abort(1);
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMF<H,F,R,D>::fail(const F& fkey)
{
 if ((Layout::primaryNode())||(!gmode)){
  std::cerr << "DataGetHandlerMF could not find requested file key:"<<std::endl;
  std::cerr << " File stub: "<< finfo.getFileStub()<<std::endl;
  XmlBufferWriter xmlout(gmode);
  fkey.output(xmlout);
  pop(xmlout);
  std::cerr << xmlout.str()<<std::endl;
  std::cerr << "...execution aborted..."<<std::endl;}
 clearData(); fileMap.clear();
 QDP_abort(1);
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMF<H,F,R,D>::fail(const std::string& msg)
{
 if ((Layout::primaryNode())||(!gmode)){
    std::cerr << "DataGetHandlerMF error: "<<msg<<std::endl;}
 QDP_abort(1);
 clearData(); fileMap.clear();
}


template <typename H, typename F, typename R, typename D>
bool DataGetHandlerMF<H,F,R,D>::queryData(const F& fkey, const R& rkey)
{
 StorageKey skey(fkey,rkey);
 typename StorageMapType::const_iterator dt=m_storage.find(skey);
 if (dt!=m_storage.end()) return true;
 IOMap<R,D> *fptr=get_file_ptr(fkey);
 if (fptr==0) return false;
 return fptr->exist(rkey);
}


template <typename H, typename F, typename R, typename D>
bool DataGetHandlerMF<H,F,R,D>::queryFile(const F& fkey)
{
 typename FileMapType::iterator it=fileMap.find(fkey);
 return (it!=fileMap.end());
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMF<H,F,R,D>::clearData()
{
 for (typename StorageMapType::iterator
      it=m_storage.begin();it!=m_storage.end();it++) delete it->second;
 m_storage.clear();
 for (typename FileMapType::iterator 
      it=fileMap.begin();it!=fileMap.end();it++) close(it->second);
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMF<H,F,R,D>::removeData(const F& fkey, const R& rkey)
{
 StorageKey skey(fkey,rkey);
 typename StorageMapType::iterator dt=m_storage.find(skey);
 if (dt!=m_storage.end()){
    delete dt->second;
    m_storage.erase(dt);}
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMF<H,F,R,D>::removeData(const F& fkey)
{
 for (typename StorageMapType::iterator
      it=m_storage.begin();it!=m_storage.end();){
    if (it->first.fkey==fkey){
        delete it->second;
        typename StorageMapType::iterator dt=it;
        it++; m_storage.erase(dt);}
    else
        it++;
    }
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMF<H,F,R,D>::getFileMap(XmlWriter& xmlout) const
{
 push(xmlout,"FileMap");
 for (typename FileMapType::const_iterator it=fileMap.begin();
      it!=fileMap.end();it++){
    push(xmlout,"Entry");
    it->first.output(xmlout);
    write(xmlout,"Suffix",(it->second).suffix);
    pop(xmlout);}
 pop(xmlout);
} 

template <typename H, typename F, typename R, typename D>
std::map<int,F> DataGetHandlerMF<H,F,R,D>::getSuffixMap() const
{
 std::map<int,F> filekeys;
 for (typename FileMapType::const_iterator it=fileMap.begin();
      it!=fileMap.end();it++)
    filekeys.insert(std::make_pair((it->second).suffix,it->first));
 return filekeys;
} 


template <typename H, typename F, typename R, typename D>
std::set<F> DataGetHandlerMF<H,F,R,D>::getFileKeys() const
{
 std::set<F> filekeys;
 for (typename FileMapType::const_iterator it=fileMap.begin();
      it!=fileMap.end();it++)
    filekeys.insert(it->first);
 return filekeys;
} 


template <typename H, typename F, typename R, typename D>
std::set<R> DataGetHandlerMF<H,F,R,D>::getKeys(const F& fkey)
{
 std::set<R> keys;
 IOMap<R,D>* fptr=get_file_ptr(fkey);
 if (fptr!=0) fptr->getKeys(keys);
 return keys;
}
 


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMF<H,F,R,D>::outputKeys(XmlWriter& xmlout)
{
 push(xmlout,"AvailableKeys");
 for (typename FileMapType::iterator it=fileMap.begin();
      it!=fileMap.end();it++){
    std::set<R> keys;
    if (it->second.fptr!=0) 
       it->second.fptr->getKeys(keys);
    else{
       open(it->second);
       it->second.fptr->getKeys(keys);
       close(it->second);} 
    for (typename std::set<R>::const_iterator 
            kt=keys.begin();kt!=keys.end();kt++){
       push(xmlout,"Key");
       it->first.output(xmlout);
       kt->output(xmlout);
       pop(xmlout);}}
 pop(xmlout);
}




   // **************************************************************
   // *                                                            *
   // *                      DataPutHandlerMF                      *
   // *                                                            *
   // **************************************************************

   // "fileMap" is a map associating a file key to a FileMapValue, 
   // which contains a suffix and an IOMap pointer.  Upon construction, 
   // "fileMap" is assigned by opening each file one by one, reading the 
   // header string, and extracting the file key.  None of the files is
   // left open.  "current" is an iterator that keeps track of the
   // current file available for data insertion.
   
   // To insert data, you must first call "open" with a file key.  If the
   // file key refers to a file that is already open, then all that happens
   // is that the current file pointer is changed. If this file is not opened,
   // then it is opened for writing.  If no file is associated with the
   // file key, then a new file is created and opened.   After the "open"
   // command, you then use "putData" to insert records into that one file.
   // Subsequently calling "open" with another file key changes the current
   // file pointer.  The previous file is left open.
   
   // "queryData" and "queryFile" opens the files, and they are left open.
   
   // "flush" and "flushAll" can be used to flush the data in one or all 
   // files.  Flushing means writing out the current IOMap record map to
   // file.  If a program crashing, an IOMap that has not been flushed is
   // corrupted.  If flushed, it will be okay for subsequent programs.
   // Note that flushing does not close the file.  Use "close" or "closeAll"
   // to actually close files.  Keep in mind that while a file is open,
   // its IOMap must maintain the map of record keys in memory.


template <typename H, typename F, typename R, typename D>
class DataPutHandlerMF : public DataPutHandlerBaseMF<H,F,R,D>
{

    struct FileMapValue
    {
      int suffix;
      IOMap<R,D> *fptr;
      FileMapValue(int in_suff) : suffix(in_suff), fptr(0) {}
      ~FileMapValue() {delete fptr;}
    };

    typedef std::map<F,FileMapValue>  FileMapType;


    FileListInfo finfo;
    FileMapType fileMap;
    H& handler;
    bool gmode;
    bool checksums;
    std::string fid;
    typename FileMapType::iterator current;
    int strpfactor, strpunit;


 public:

    DataPutHandlerMF(H& in_handler, const FileListInfo& in_filelist,
                     const std::string& filetype_id,  
                     const std::string& header_tag, 
                     bool global_mode=true, bool use_checksums=false,
                     int striping_factor=1, int striping_unit=0);

    ~DataPutHandlerMF() {fileMap.clear();}

    void setOverWrite() {finfo.setOverWrite();}

    void setNoOverWrite() {finfo.setNoOverWrite();}

    const FileListInfo& getFileListInfo() const {return finfo;}

    bool isGlobal() const {return gmode;}

    bool isLocal() const {return !gmode;}


    void open(const F& fkey);

    void putData(const R& rkey, const D& data);  // insert into current file

    void flush();   // flush current file

    void close();   // close current file

             // calls openFile(fkey) first, then inserts
             
    void putData(const F& fkey, const R& rkey, const D& data);

    void flush(const F& fkey);

    void close(const F& fkey);
    
    void flushAll();

    void closeAll();



    bool queryData(const F& fkey, const R& rkey);

    bool queryData(const R& rkey);  // query in current open file

    bool queryFile(const F& fkey);



    void merge(const FileListInfo& infiles,        // adds data from the "infiles"
               const std::string& header_tag);     // into this file list
                                                   // use in serial code ONLY
                                              
    void merge(const std::vector<FileListInfo>& infiles,  
               const std::string& header_tag);     

    void getFileMap(XmlWriter& xmlout) const;

    std::map<int,F> getSuffixMap() const;

    std::set<F> getFileKeys() const;


 private:


    void fail(const std::string& msg);
    
    void fail(const F& fkey, const R& rkey);

    typename FileMapType::iterator get_file_ptr(const F& fkey);
    
    void openUpdate(FileMapValue& fmv);
    
    void openNew(const F& fkey, FileMapValue& fmv);

    void close(typename FileMapType::iterator it);

    void flush(typename FileMapType::iterator it);


          // disallow copies
    DataPutHandlerMF(const DataPutHandlerMF& in);
    DataPutHandlerMF& operator=(const DataPutHandlerMF& in);

};



   // constructor checks that the information in the headers
   // of all existing files is consistent, then sets up the
   // file map

template <typename H, typename F, typename R, typename D>
DataPutHandlerMF<H,F,R,D>::DataPutHandlerMF(H& in_handler,
                                            const FileListInfo& in_filelist,
                                            const std::string& filetype_id,
                                            const std::string& header_tag,
                                            bool global_mode, bool use_checksums,
                                            int striping_factor, int striping_unit)
  :  finfo(in_filelist), handler(in_handler), gmode(global_mode), checksums(use_checksums), 
     fid(tidyString(filetype_id)), current(fileMap.end()),
     strpfactor(striping_factor), strpunit(striping_unit)
{
 for (int suffix=finfo.getMinFileNumber();
          suffix<=finfo.getMaxFileNumber();suffix++){

    std::string filename=finfo.getFileName(suffix);

         // open all existing files and check consistency of headers
    std::string headerxml;
    bool exists;
    {IOMap<R,D> iom(gmode);
     exists=iom.peekHeader(headerxml,filename,fid);}
    if (!exists) continue;

/*
    if (!fileExists(filename,gmode)) continue;

         // open all existing files and check consistency of headers
    string headerxml;
    {IOMap<R,D> iom(gmode);
     try { iom.openReadOnly(filename,fid,headerxml,checksums);}
     catch(...) {
        fail("could not open file "+filename+" for reading");}} */

    if ((!gmode)||(Layout::primaryNode())){
       XmlReader xmlr(headerxml,false);
       if (!handler.checkHeader(xmlr,suffix)){
          fail("Header string in file is\n"+headerxml+"header info in file "
                     +filename+" does not match info in current Handler\n\n"
                     +"...execution aborted...\n");}}

            // extract the file key from this file
    try{
       XmlReader xmlr(headerxml,gmode);
       XmlReader xmlf(xmlr,"./descendant-or-self::"+header_tag);
       F fkey(xmlf);
       typename FileMapType::iterator it=fileMap.find(fkey);
       if (it!=fileMap.end()){
          fail(std::string("duplicate keys in fileMap in current Handler\n")
               +" ... too confusing to continue\n file suffix "
               +int_to_string(suffix)+" and suffix "+int_to_string((it->second).suffix)
               +" have same file key\n");}
       fileMap.insert(std::make_pair(fkey, FileMapValue(suffix)));}
    catch(...){
       fail("Could not extract FileKey from file "+filename+"\n");}
    }
}


     //  Get file pointer corresponding to fkey.  If "fkey" is already in the
     //  fileMap, then check if the file is open (fptr assigned).  If not open,
     //  open the file.  If already open, we are done.  If "fkey" is not in the
     //  fileMap, then we find a new suffix and open a new file.

template <typename H, typename F, typename R, typename D>
typename DataPutHandlerMF<H,F,R,D>::FileMapType::iterator 
        DataPutHandlerMF<H,F,R,D>::get_file_ptr(const F& fkey)
{
 typename FileMapType::iterator it=fileMap.find(fkey);
 if (it!=fileMap.end()){
    if (it->second.fptr==0) openUpdate(it->second);}
 else{
    int findex=0;
    try{ findex=finfo.getFirstAvailableSuffix(gmode);}
    catch(...){
       fail(std::string("could not create file...no more suffix indices")
                  +" available...execution aborted...\n");}
    std::pair< typename FileMapType::iterator, bool > pr;
    pr=fileMap.insert(std::make_pair(fkey,FileMapValue(findex)));
    if (pr.second==true) it=pr.first;
    else{ fail("DataPutHandlerMF insertion failed on index "+int_to_string(findex));}
    openNew(fkey,it->second);}
 return it;
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::openUpdate(FileMapValue& fmv)
{
 std::string filename=finfo.getFileName(fmv.suffix);
 try {
    fmv.fptr=new IOMap<R,D>(gmode);
    std::string header;
    fmv.fptr->openUpdate(filename,fid,header,'L',1,0,checksums,finfo.isModeOverwrite());}
 catch(...){
    fail("failure opening file "+filename+" in DataPutHandlerMF");}
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::openNew(const F& fkey, FileMapValue& fmv)
{
 std::string filename=finfo.getFileName(fmv.suffix);
 if (fileExists(filename,gmode)){
    fail(std::string("file collision: another process has created a file")
         +" during execution of this program..best to abort\n");}
 try {
    fmv.fptr=new IOMap<R,D>(gmode);
    XmlBufferWriter headerxml(gmode);
    handler.writeHeader(headerxml,fkey,fmv.suffix);  // write header info 
    fmv.fptr->openNew(filename,fid,headerxml.str(),false,'L',strpfactor,strpunit,checksums,
                      finfo.isModeOverwrite());}
 catch(...){
    fail("failure opening file "+filename+" in DataPutHandlerMF");}
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::close(typename FileMapType::iterator it)
{
 if (it==fileMap.end()) return;
 delete it->second.fptr; 
 it->second.fptr=0;
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::flush(typename FileMapType::iterator it)
{
 if (it==fileMap.end()) return;
 if (it->second.fptr!=0) it->second.fptr->flush();
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::fail(const std::string& msg)
{
 if ((Layout::primaryNode())||(!gmode)){
    std::cerr << "DataPutHandlerMF error: "<<msg<<std::endl;}
 fileMap.clear();
 QDP_abort(1);
}

template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::fail(const F& fkey, const R& rkey)
{
 if ((Layout::primaryNode())||(!gmode)){
  std::cerr << "DataPutHandlerMF could not insert requested record:"<<std::endl;
  std::cerr << " File stub: "<< finfo.getFileStub()<<std::endl;
  XmlBufferWriter xmlout(gmode);
  push(xmlout,"FileRecordKey");
  fkey.output(xmlout);
  rkey.output(xmlout);
  pop(xmlout);
  std::cerr << xmlout.str()<<std::endl;
  std::cerr << "...execution aborted..."<<std::endl;}
 fileMap.clear();
 QDP_abort(1);
}


    // opens file, resets the current file pointer

template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::open(const F& fkey)
{
 current=get_file_ptr(fkey);
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::putData(const R& rkey, const D& data)
{
 if (current==fileMap.end()) { fail("No current file; cannot insert Data");}
 if (current->second.fptr==0) openUpdate(current->second);
 try{current->second.fptr->put(rkey,data);}
 catch(...){
    fail(current->first,rkey);}
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::flush()
{
 flush(current);
}

template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::close()
{
 close(current);
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::flushAll()
{
 for (typename FileMapType::iterator it=fileMap.begin();it!=fileMap.end();it++)
    flush(it);
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::closeAll()
{
 for (typename FileMapType::iterator it=fileMap.begin();it!=fileMap.end();it++)
    close(it);
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::putData(
                 const F& fkey, const R& rkey, const D& data)
{
 open(fkey);
 try{ current->second.fptr->put(rkey,data);}
 catch(...){
    fail(fkey,rkey);}
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::flush(const F& fkey)
{
 typename FileMapType::iterator it=fileMap.find(fkey);
 flush(it);
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::close(const F& fkey)
{
 typename FileMapType::iterator it=fileMap.find(fkey);
 close(it);
}


template <typename H, typename F, typename R, typename D>
bool DataPutHandlerMF<H,F,R,D>::queryData(
                 const F& fkey, const R& rkey)
{
 typename FileMapType::iterator it=fileMap.find(fkey);
 if (it==fileMap.end()) return false;
 if (it->second.fptr==0) openUpdate(it->second);
 return (it->second.fptr)->exist(rkey);
}


template <typename H, typename F, typename R, typename D>
bool DataPutHandlerMF<H,F,R,D>::queryData(const R& rkey)
{
 if (current==fileMap.end()) { fail("No current file; cannot query Data");}
 if (current->second.fptr==0) openUpdate(current->second);
 return current->second.fptr->exist(rkey);
}


template <typename H, typename F, typename R, typename D>
bool DataPutHandlerMF<H,F,R,D>::queryFile(const F& fkey)
{
 typename FileMapType::iterator it=fileMap.find(fkey);
 return (it!=fileMap.end());
}



template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::getFileMap(XmlWriter& xmlout) const
{
 push(xmlout,"FileMap");
 for (typename FileMapType::const_iterator it=fileMap.begin();
      it!=fileMap.end();it++){
    push(xmlout,"Entry");
    it->first.output(xmlout);
    write(xmlout,"Suffix",(it->second).suffix);
    pop(xmlout);}
 pop(xmlout);
} 

template <typename H, typename F, typename R, typename D>
std::map<int,F> DataPutHandlerMF<H,F,R,D>::getSuffixMap() const
{
 std::map<int,F> filekeys;
 for (typename FileMapType::const_iterator it=fileMap.begin();
      it!=fileMap.end();it++)
    filekeys.insert(std::make_pair((it->second).suffix,it->first));
 return filekeys;
} 


template <typename H, typename F, typename R, typename D>
std::set<F> DataPutHandlerMF<H,F,R,D>::getFileKeys() const
{
 std::set<F> filekeys;
 for (typename FileMapType::const_iterator it=fileMap.begin();
      it!=fileMap.end();it++)
    filekeys.insert(it->first);
 return filekeys;
} 



template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::merge(const FileListInfo& infiles,
                                      const std::string& header_tag)
{
 DataGetHandlerMF<H,F,R,D> inmerge(handler,infiles,fid,header_tag,gmode);
 std::set<F> fkeys=inmerge.getFileKeys();
 QDPIO::cout << "  number of file keys = "<<fkeys.size()<<std::endl;
 for (typename std::set<F>::const_iterator ft=fkeys.begin();ft!=fkeys.end();ft++){
    open(*ft);
    std::set<R> rkeys=inmerge.getKeys(*ft);
    QDPIO::cout << "  file suffix "<<current->second.suffix
                <<": number of records to be added = "<<rkeys.size()<<std::endl;
    for (typename std::set<R>::const_iterator rt=rkeys.begin();rt!=rkeys.end();rt++){
       const D& data=inmerge.getData(*ft,*rt);
       putData(*rt,data);
       inmerge.removeData(*ft,*rt);}
//    close(*ft);
    inmerge.clearData();
    }
} 

/*
template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::merge(const vector<FileListInfo>& infiles,
                                      const std::string& header_tag)
{
 for (int k=0;k<infiles.size();k++){
    QDPIO::cout << "merging file list k = "<<k<<endl;
    DataGetHandlerMF<H,F,R,D> inmerge(handler,infiles[k],fid,header_tag,gmode);
    std::set<F> fkeys=inmerge.getFileKeys();
    QDPIO::cout << "  number of file keys = "<<fkeys.size()<<endl;
    for (typename std::set<F>::const_iterator ft=fkeys.begin();ft!=fkeys.end();ft++){
       open(*ft);
       std::set<R> rkeys=inmerge.getKeys(*ft);
       QDPIO::cout << "  file suffix "<<current->second.suffix
                   <<": number of records to be added = "<<rkeys.size()<<endl;
       for (typename std::set<R>::const_iterator rt=rkeys.begin();rt!=rkeys.end();rt++){
          const D& data=inmerge.getData(*ft,*rt);
          putData(*rt,data);
          inmerge.removeData(*ft,*rt);}
//       close(*ft);
       inmerge.clearData();
       }
    }
} */

   //  The merge subroutine above causes lustre lock up on ranger when several
   //  jobs are running simultaneously.  Hence, it has been revised below in such
   //  a way to slow down directory searches.

template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::merge(const std::vector<FileListInfo>& infiles,
                                      const std::string& header_tag)
{
 for (unsigned int k=0;k<infiles.size();k++){
    QDPIO::cout << "merging file list k = "<<k<<std::endl;
    std::string instub=infiles[k].getFileStub();
    int nfkeys=0;
    for (int suffix=infiles[k].getMinFileNumber();suffix<=infiles[k].getMaxFileNumber();suffix++){
       FileListInfo inf(instub,suffix,suffix);
       DataGetHandlerMF<H,F,R,D> inmerge(handler,inf,fid,header_tag,gmode);
       std::set<F> fkeys=inmerge.getFileKeys();
       nfkeys+=fkeys.size();
       for (typename std::set<F>::const_iterator ft=fkeys.begin();ft!=fkeys.end();ft++){
          open(*ft);
          std::set<R> rkeys=inmerge.getKeys(*ft);
          QDPIO::cout << " file list "<<k<<":  from suffix "<<suffix
                      << " into merged suffix "<<current->second.suffix
                      <<": number of records added = "<<rkeys.size()<<std::endl;
          for (typename std::set<R>::const_iterator rt=rkeys.begin();rt!=rkeys.end();rt++){
             const D& data=inmerge.getData(*ft,*rt);
             putData(*rt,data);
             inmerge.removeData(*ft,*rt);}
//        close(*ft);
          }}
    QDPIO::cout << "  number of file keys merged in list "<<k<<" was "<<nfkeys<<std::endl; 
    }
}

   // **************************************************************
   // *                                                            *
   // *                      DataGetHandlerSF                      *
   // *                                                            *
   // **************************************************************


template <typename H, typename R, typename D>
class DataGetHandlerSF
{

    std::map<R,D*> m_storage;
    IOMap<R,D> *iomptr;
    H& handler;

 public:

    DataGetHandlerSF(H& in_handler, const std::string& file_name, 
                     const std::string& filetype_id, bool global_mode=true,
                     bool use_checksums=false);

    ~DataGetHandlerSF() {clearData(); delete iomptr;}

    const std::string& getFileName() const {return iomptr->GetFileName();}

    bool isGlobal() const {return iomptr->isGlobal();}

    bool isLocal() const {return iomptr->isLocal();}


    bool queryData(const R& rkey);

    const D& getData(const R& rkey);

    void removeData(const R& rkey);

    void clearData();


    std::set<R> getKeys();

    void outputKeys(XmlWriter& xmlout);


 private:

    void fail(const R& rkey);

    void fail(const std::string& msg);

          // disallow copies
    DataGetHandlerSF(const DataGetHandlerSF& in);
    DataGetHandlerSF& operator=(const DataGetHandlerSF& in);

};



   // constructor checks that the information in the header
   // is consistent

template <typename H, typename R, typename D>
DataGetHandlerSF<H,R,D>::DataGetHandlerSF(H& in_handler, const std::string& filename,
                                          const std::string& filetype_id, 
                                          bool global_mode, bool use_checksums)
                       :  handler(in_handler)
{
 std::string headerxml;
 try{
    iomptr=new IOMap<R,D>(global_mode);
    iomptr->openReadOnly(filename,filetype_id,headerxml,use_checksums);}
 catch(...) {
    fail("could not open file "+filename+" for reading");}

 if ((!global_mode)||(Layout::primaryNode())){
    XmlReader xmlr(headerxml,false);
    if (!handler.checkHeader(xmlr)){
       fail("Header string in file is \n"+headerxml+"\n header info in file "+filename
            +" does not match info in current Handler\n ...execution aborted...\n");}}
}


template <typename H, typename R, typename D>
const D& DataGetHandlerSF<H,R,D>::getData(const R& rkey)
{
 typename std::map<R,D*>::const_iterator dt=m_storage.find(rkey);
 if (dt!=m_storage.end()) return *(dt->second);
 D *result(new D);
 m_storage[rkey]=result;
 try {iomptr->get(rkey,*result);}
 catch(...){delete result; fail(rkey);}
 return *result;
}


template <typename H, typename R, typename D>
bool DataGetHandlerSF<H,R,D>::queryData(const R& rkey)
{
 return iomptr->exist(rkey);
}


template <typename H, typename R, typename D>
void DataGetHandlerSF<H,R,D>::clearData()
{
 for (typename std::map<R,D*>::iterator
      it=m_storage.begin();it!=m_storage.end();it++) delete it->second;
 m_storage.clear();
}


template <typename H, typename R, typename D>
void DataGetHandlerSF<H,R,D>::removeData(const R& rkey)
{
 typename std::map<R,D*>::iterator dt=m_storage.find(rkey);
 if (dt!=m_storage.end()){
    delete dt->second;
    m_storage.erase(dt);}
}



template <typename H, typename R, typename D>
std::set<R> DataGetHandlerSF<H,R,D>::getKeys()
{
 std::set<R> keys;
 iomptr->getKeys(keys);
 return keys;
} 


template <typename H, typename R, typename D>
void DataGetHandlerSF<H,R,D>::outputKeys(XmlWriter& xmlout)
{
 push(xmlout,"AvailableKeys");
 std::set<R> keys;
 iomptr->getKeys(keys);
 for (typename std::set<R>::const_iterator 
            kt=keys.begin();kt!=keys.end();kt++){
    push(xmlout,"Key");
    kt->output(xmlout);
    pop(xmlout);}
 pop(xmlout);
}



template <typename H, typename R, typename D>
void DataGetHandlerSF<H,R,D>::fail(const R& rkey)
{
 if ((Layout::primaryNode())||(iomptr->isLocal())){
  std::cerr << "DataGetHandlerSF could not find requested record:"<<std::endl;
  std::cerr << " File name: "<< iomptr->getFileName() <<std::endl;
  XmlBufferWriter xmlout(iomptr->isGlobal());
  push(xmlout,"RecordKey");
  rkey.output(xmlout);
  pop(xmlout);
  std::cerr << xmlout.str()<<std::endl;
  std::cerr << "...execution aborted..."<<std::endl;}
 clearData(); delete iomptr;
 QDP_abort(1);
}



template <typename H, typename R, typename D>
void DataGetHandlerSF<H,R,D>::fail(const std::string& msg)
{
 if ((Layout::primaryNode())||(iomptr->isLocal())){
    std::cerr << "DataGetHandlerSF error: "<<msg<<std::endl;}
 clearData(); delete iomptr;
 QDP_abort(1);
}



   // **************************************************************
   // *                                                            *
   // *                      DataPutHandlerSF                      *
   // *                                                            *
   // **************************************************************


template <typename H, typename R, typename D>
class DataPutHandlerSF
{

    H& handler;
    IOMap<R,D> *iomptr;


 public:

    DataPutHandlerSF(H& inptr, const std::string& file_name, 
                     const std::string& filetype_id,
                     bool overwrite=false, bool global_mode=true,
                     bool use_checksums=false, int striping_factor=1,
                     int striping_unit=0);

    ~DataPutHandlerSF() {delete iomptr;}

    const std::string& getFileName() const {return iomptr->getFileName();}

    bool isGlobal() const {return iomptr->isGlobal();}

    bool isLocal() const {return iomptr->isLocal();}


    void putData(const R& rkey, const D& data);

    void flush();

    bool queryData(const R& rkey);


 private:
 
    void fail(const std::string& msg);

    void fail(const R& rkey);


          // disallow copies
    DataPutHandlerSF(const DataPutHandlerSF& in);
    DataPutHandlerSF& operator=(const DataPutHandlerSF& in);

};



   // constructor checks that the information in the header
   // is consistent

template <typename H, typename R, typename D>
DataPutHandlerSF<H,R,D>::DataPutHandlerSF(H& in_handler,
                                          const std::string& file_name,
                                          const std::string& filetype_id,
                                          bool overwrite, bool global_mode,
                                          bool use_checksums, int striping_factor,
                                          int striping_unit)
                     :  handler(in_handler)
{
 XmlBufferWriter headerxml(global_mode);
 handler.writeHeader(headerxml);  // write header info 
 try{
   iomptr=new IOMap<R,D>(global_mode);
   std::string header(headerxml.str());
   iomptr->openUpdate(file_name,filetype_id,header,'L',
                      striping_factor,striping_unit,use_checksums,overwrite);
   if (!(iomptr->isNewFile())){
      if ((!global_mode)||(Layout::primaryNode())){
         XmlReader xmlr(header,false);
         if (!handler.checkHeader(xmlr)){
            fail("Header string in file is \n"+header+"\n header info in file "
                 +file_name+" does not match info in current Handler\n "
                 +"...execution aborted...\n");}}}}
 catch(...) {
    fail("could not open file "+file_name+" for writing"); }
}


template <typename H, typename R, typename D>
void DataPutHandlerSF<H,R,D>::putData(const R& rkey, const D& data)
{
 try{
    iomptr->put(rkey,data);}
 catch(...){
    fail(rkey);}
}


template <typename H, typename R, typename D>
void DataPutHandlerSF<H,R,D>::flush()
{
 iomptr->flush();
}


template <typename H, typename R, typename D>
bool DataPutHandlerSF<H,R,D>::queryData(const R& rkey)
{
 return iomptr->exist(rkey);
}



template <typename H, typename R, typename D>
void DataPutHandlerSF<H,R,D>::fail(const std::string& msg)
{
 if ((Layout::primaryNode())||(iomptr->isLocal())){
    std::cerr << "DataPutHandlerSF error: "<<msg<<std::endl;}
 delete iomptr;
 QDP_abort(1);
}


template <typename H, typename R, typename D>
void DataPutHandlerSF<H,R,D>::fail(const R& rkey)
{
 if ((Layout::primaryNode())||(iomptr->isLocal())){
  std::cerr << "DataPutHandlerSF could not insert requested record:"<<std::endl;
  std::cerr << " File name: "<< iomptr->getFileName() <<std::endl;
  XmlBufferWriter xmlout(iomptr->isGlobal());
  push(xmlout,"RecordKey");
  rkey.output(xmlout);
  pop(xmlout);
  std::cerr << xmlout.str()<<std::endl;
  std::cerr << "...execution aborted..."<<std::endl;}
 delete iomptr;
 QDP_abort(1);
}


   // **************************************************************
   // *                                                            *
   // *                      DataGetHandlerMemMF                   *
   // *                                                            *
   // **************************************************************


   // "get" class handles the internal storage in an std::map
   // stored in the NamedObjectMap. Ultimately, only pointers are
   // stored. The use of pointers is efficient
   // whenever each data structure is fairly large (100 elements or 
   // more) since only pointers get copied.  
   
   // "fileMap" is a map associating a file key to a FileMapValue, 
   // which is the above-mentioned std::map mapping RecordKeys to
   // the DataType. Upon construction, all header strings are read
   // from the FileHeaderMap and checked for compatibility just as 
   // the usual GetHandlerMF would do.
   // Data access is trivial because all data is in the memory
   // already.
   
   // "removeData" does remove the data from memory, so make sure
   // that this is actually intended.
   // "clearData" wipes all the data from memory.


template <typename H, typename F, typename R, typename D>
class DataGetHandlerMemMF : public DataGetHandlerBaseMF<H,F,R,D>
{

    typedef std::map<R,D>  FileMapValue;
    typedef std::map<F,FileMapValue>  FileMapType;
    typedef std::map<F,std::string>   FileHeaderMapType;


    FileMapType*  fileMap;
    FileHeaderMapType* headerMap;
    H& handler;
    std::string fid;


 public:

    DataGetHandlerMemMF(H& in_handler, const std::string& filetype_id);

    ~DataGetHandlerMemMF() {}

    bool isGlobal() const {return true;}

    bool isLocal() const {return false;}


    bool queryData(const F& fkey, const R& rkey);

    bool queryFile(const F& fkey);


    const D& getData(const F& fkey, const R& rkey);

    void removeData(const F& fkey, const R& rkey);

    void removeData(const F& fkey);

    void clearData();


    void getFileMap(XmlWriter& xmlout) const;

    std::map<int,F> getSuffixMap() const;

    std::set<F> getFileKeys() const;

    std::set<R> getKeys(const F& fkey);

    void outputKeys(XmlWriter& xmlout);


 private:

    void fail(const F& fkey, const R& rkey);
    
    void fail(const F& fkey);

    void fail(const std::string& msg);
    

          // disallow copies
    DataGetHandlerMemMF(const DataGetHandlerMemMF& in);
    DataGetHandlerMemMF& operator=(const DataGetHandlerMemMF& in);

};



   // Constructor checks that the information in the headers
   // of all existing files is consistent, then sets up the
   // file map.  Each files is opened (one by one), the header string
   // is read, and then the file is closed.

template <typename H, typename F, typename R, typename D>
DataGetHandlerMemMF<H,F,R,D>::DataGetHandlerMemMF(H& in_handler,
                                            const std::string& filetype_id)
  :  handler(in_handler), fid(tidyString(filetype_id))
{

 if (!TheNamedObjMap::Instance().check(fid)) {
   QDPIO::cerr << "  Warning: did not find data in the NamedObjectMap for " << fid <<std::endl<<std::endl;
   TheNamedObjMap::Instance().create<FileMapType>(fid);
 }

 if (!TheNamedObjMap::Instance().check(fid+"Header")) {
   QDPIO::cerr << "  Warning: did not find header data in the NamedObjectMap for " << fid <<std::endl<<std::endl;
   TheNamedObjMap::Instance().create<FileHeaderMapType>(fid+"Header");
 }

 fileMap = &(TheNamedObjMap::Instance().getData<FileMapType>(fid));
 headerMap = &(TheNamedObjMap::Instance().getData<FileHeaderMapType>(fid+"Header"));

 std::set<F> fkeys = getFileKeys();

 for (typename std::set<F>::const_iterator fIt=fkeys.begin(); fIt!=fkeys.end(); ++fIt) {
   if (Layout::primaryNode()){
      XmlReader xmlr((*headerMap)[*fIt],false);
      if (!handler.checkHeader(xmlr,0)){
         fail("Header string in file is\n"+(*headerMap)[*fIt]+"header info in file "
                    +"does not match info in current Handler\n\n"
                    +"...execution aborted...\n");}}

 }

}


template <typename H, typename F, typename R, typename D>
const D& DataGetHandlerMemMF<H,F,R,D>::getData(const F& fkey, const R& rkey)
{
 try {return (*fileMap)[fkey][rkey];}
 catch(...){fail(fkey,rkey);}
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMemMF<H,F,R,D>::fail(const F& fkey, const R& rkey)
{
 if (Layout::primaryNode()){
  std::cerr << "DataGetHandlerMemMF could not find requested record:"<<std::endl;
  XmlBufferWriter xmlout(true);
  push(xmlout,"FileRecordKey");
  fkey.output(xmlout);
  rkey.output(xmlout);
  pop(xmlout);
  std::cerr << xmlout.str()<<std::endl;
  std::cerr << "...execution aborted..."<<std::endl;}
 clearData(); fileMap->clear();
 QDP_abort(1);
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMemMF<H,F,R,D>::fail(const F& fkey)
{
 if (Layout::primaryNode()){
  std::cerr << "DataGetHandlerMemMF could not find requested file key:"<<std::endl;
  XmlBufferWriter xmlout(true);
  fkey.output(xmlout);
  pop(xmlout);
  std::cerr << xmlout.str()<<std::endl;
  std::cerr << "...execution aborted..."<<std::endl;}
 clearData(); fileMap->clear();
 QDP_abort(1);
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMemMF<H,F,R,D>::fail(const std::string& msg)
{
 if (Layout::primaryNode()){
    std::cerr << "DataGetHandlerMemMF error: "<<msg<<std::endl;}
 QDP_abort(1);
 clearData(); fileMap->clear();
}


template <typename H, typename F, typename R, typename D>
bool DataGetHandlerMemMF<H,F,R,D>::queryData(const F& fkey, const R& rkey)
{
 typename FileMapType::iterator it=fileMap->find(fkey);
 if (it == fileMap->end()) return false;

 return (it->second.find(rkey) != it->second.end());
}


template <typename H, typename F, typename R, typename D>
bool DataGetHandlerMemMF<H,F,R,D>::queryFile(const F& fkey)
{
 typename FileMapType::iterator it=fileMap->find(fkey);
 return (it!=fileMap->end());
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMemMF<H,F,R,D>::clearData()
{
 fileMap->clear();
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMemMF<H,F,R,D>::removeData(const F& fkey, const R& rkey)
{
 (*fileMap)[fkey].erase(rkey);
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMemMF<H,F,R,D>::removeData(const F& fkey)
{
 (*fileMap).erase(fkey);
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMemMF<H,F,R,D>::getFileMap(XmlWriter& xmlout) const
{
 push(xmlout,"FileMap");
 for (typename FileMapType::const_iterator it=fileMap->begin();
      it!=fileMap->end();it++){
    push(xmlout,"Entry");
    it->first.output(xmlout);
    pop(xmlout);}
 pop(xmlout);
} 

template <typename H, typename F, typename R, typename D>
std::map<int,F> DataGetHandlerMemMF<H,F,R,D>::getSuffixMap() const
{
 int suffix = 0;
 std::map<int,F> filekeys;
 for (typename FileMapType::const_iterator it=fileMap->begin();
      it!=fileMap->end();it++)
    filekeys.insert(std::make_pair(suffix,it->first));
 return filekeys;
} 


template <typename H, typename F, typename R, typename D>
std::set<F> DataGetHandlerMemMF<H,F,R,D>::getFileKeys() const
{
 std::set<F> filekeys;
 for (typename FileMapType::const_iterator it=fileMap->begin();
      it!=fileMap->end();it++)
    filekeys.insert(it->first);
 return filekeys;
} 


template <typename H, typename F, typename R, typename D>
std::set<R> DataGetHandlerMemMF<H,F,R,D>::getKeys(const F& fkey)
{
 std::set<R> keys;
 for (typename FileMapValue::const_iterator it=(*fileMap)[fkey].begin();
      it!=(*fileMap)[fkey].end();it++)
    keys.insert(it->first);

 return keys;
}
 


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMemMF<H,F,R,D>::outputKeys(XmlWriter& xmlout)
{
 push(xmlout,"AvailableKeys");
 for (typename FileMapType::iterator it=fileMap->begin();
      it!=fileMap->end();it++){
    std::set<R> keys = getKeys(it->first);

    for (typename std::set<R>::const_iterator 
            kt=keys.begin();kt!=keys.end();kt++){
       push(xmlout,"Key");
       it->first.output(xmlout);
       kt->output(xmlout);
       pop(xmlout);}}
 pop(xmlout);
}



   // **************************************************************
   // *                                                            *
   // *                      DataPutHandlerMemMF                   *
   // *                                                            *
   // **************************************************************

   // "fileMap" is a map associating a file key to a FileMapValue, 
   // which is a std::map from RecordKey to DataType. Upon construction, 
   // "fileMap" is assigned looking for data in the NamedObjectMap. 
   // "current" is an iterator that keeps track of the current file available
   // for data insertion.
   
   // To insert data, you must first call "open" with a file key.  If the
   // file key refers to something already in the map, then all that happens
   // is that the current file pointer is changed. If this file is not there,
   // then an empty FileMapValue map is created.
   // After the "open" command, you then use "putData" to insert records into 
   // that one file. Subsequently calling "open" with another file key changes
   //  the current file pointer.
   
   // "queryData" and "queryFile" only checks if they are in the maps.
   
   // "flush" and "flushAll" are only in the interface for compatibility, as
   // are "close" or "closeAll".


template <typename H, typename F, typename R, typename D>
class DataPutHandlerMemMF : public DataPutHandlerBaseMF<H,F,R,D>
{

    typedef std::map<R,D>  FileMapValue;
    typedef std::map<F,FileMapValue>  FileMapType;
    typedef std::map<F,std::string>   FileHeaderMapType;


    FileMapType* fileMap;
    FileHeaderMapType* headerMap;
    H& handler;
    std::string fid;
    typename FileMapType::iterator current;


 public:

    DataPutHandlerMemMF(H& in_handler, const std::string& filetype_id);

    ~DataPutHandlerMemMF() {}

    bool isGlobal() const {return true;}

    bool isLocal() const {return false;}


    void open(const F& fkey);

    void putData(const R& rkey, const D& data);  // insert into current file

    void flush() {};   // flush current file

    void close() { current = fileMap->end(); };   // close current file

             // calls openFile(fkey) first, then inserts
             
    void putData(const F& fkey, const R& rkey, const D& data);

    void flush(const F& fkey) {}

    void close(const F& fkey) { current = fileMap->end(); }
    
    void flushAll() {}

    void closeAll() { current = fileMap->end(); }



    bool queryData(const F& fkey, const R& rkey);

    bool queryData(const R& rkey);  // query in current open file

    bool queryFile(const F& fkey);



    void getFileMap(XmlWriter& xmlout) const;

    std::map<int,F> getSuffixMap() const;

    std::set<F> getFileKeys() const;


 private:


    void fail(const std::string& msg);
    
    void fail(const F& fkey, const R& rkey);

    typename FileMapType::iterator get_file_ptr(const F& fkey);
    

          // disallow copies
    DataPutHandlerMemMF(const DataPutHandlerMemMF& in);
    DataPutHandlerMemMF& operator=(const DataPutHandlerMemMF& in);

};



   // constructor sets up the file map in the NamedObjectMap

template <typename H, typename F, typename R, typename D>
DataPutHandlerMemMF<H,F,R,D>::DataPutHandlerMemMF(H& in_handler,
                                            const std::string& filetype_id)
  :  handler(in_handler), fid(tidyString(filetype_id))
{
 if (!TheNamedObjMap::Instance().check(fid))
   TheNamedObjMap::Instance().create<FileMapType>(fid);

 if (!TheNamedObjMap::Instance().check(fid+"Header"))
   TheNamedObjMap::Instance().create<FileHeaderMapType>(fid+"Header");

 fileMap = &(TheNamedObjMap::Instance().getData<FileMapType>(fid));
 headerMap = &(TheNamedObjMap::Instance().getData<FileHeaderMapType>(fid+"Header"));

 current = fileMap->end();

}


     //  Get file pointer corresponding to fkey.  If "fkey" is already in the
     //  fileMap, then check if the file is open (fptr assigned).  If not open,
     //  open the file.  If already open, we are done.  If "fkey" is not in the
     //  fileMap, then we find a new suffix and open a new file.

template <typename H, typename F, typename R, typename D>
typename DataPutHandlerMemMF<H,F,R,D>::FileMapType::iterator 
        DataPutHandlerMemMF<H,F,R,D>::get_file_ptr(const F& fkey)
{

 typename FileMapType::iterator it = fileMap->find(fkey);
 if (it != fileMap->end()) return it;

 XmlBufferWriter headerxml(true);
 handler.writeHeader(headerxml,fkey,0);  // write header info 
 headerMap->insert( std::make_pair(fkey, headerxml.str()) );

 it = fileMap->insert( std::make_pair(fkey, FileMapValue()) ).first;
 return it;
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMemMF<H,F,R,D>::fail(const std::string& msg)
{
 if (Layout::primaryNode()){
    std::cerr << "DataPutHandlerMemMF error: "<<msg<<std::endl;}
 fileMap->clear();
 QDP_abort(1);
}

template <typename H, typename F, typename R, typename D>
void DataPutHandlerMemMF<H,F,R,D>::fail(const F& fkey, const R& rkey)
{
 if (Layout::primaryNode()){
  std::cerr << "DataPutHandlerMemMF could not insert requested record:"<<std::endl;
  XmlBufferWriter xmlout(false);
  push(xmlout,"FileRecordKey");
  fkey.output(xmlout);
  rkey.output(xmlout);
  pop(xmlout);
  std::cerr << xmlout.str()<<std::endl;
  std::cerr << "...execution aborted..."<<std::endl;}
 fileMap->clear();
 QDP_abort(1);
}


    // opens file, resets the current file pointer

template <typename H, typename F, typename R, typename D>
void DataPutHandlerMemMF<H,F,R,D>::open(const F& fkey)
{
 current=get_file_ptr(fkey);
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMemMF<H,F,R,D>::putData(const R& rkey, const D& data)
{
 if (current==fileMap->end()) { fail("No current file; cannot insert Data");}
 try{current->second.insert( std::make_pair(rkey,data) );}
 catch(...){
    fail(current->first,rkey);}
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMemMF<H,F,R,D>::putData(
                 const F& fkey, const R& rkey, const D& data)
{
 open(fkey);
 try{current->second.insert( std::make_pair(rkey,data) );}
 catch(...){
    fail(fkey,rkey);}
}


template <typename H, typename F, typename R, typename D>
bool DataPutHandlerMemMF<H,F,R,D>::queryData(
                 const F& fkey, const R& rkey)
{
 typename FileMapType::iterator it=fileMap->find(fkey);
 if (it==fileMap->end()) return false;
 return (it->second.find(rkey)!=it->second.end());
}


template <typename H, typename F, typename R, typename D>
bool DataPutHandlerMemMF<H,F,R,D>::queryData(const R& rkey)
{
 if (current==fileMap->end()) { fail("No current file; cannot query Data");}
 return (current->second.find(rkey)!=current->second.end());
}


template <typename H, typename F, typename R, typename D>
bool DataPutHandlerMemMF<H,F,R,D>::queryFile(const F& fkey)
{
 typename FileMapType::iterator it=fileMap->find(fkey);
 return (it!=fileMap->end());
}



template <typename H, typename F, typename R, typename D>
void DataPutHandlerMemMF<H,F,R,D>::getFileMap(XmlWriter& xmlout) const
{
 push(xmlout,"FileMap");
 for (typename FileMapType::const_iterator it=fileMap->begin();
      it!=fileMap->end();it++){
    push(xmlout,"Entry");
    it->first.output(xmlout);
    pop(xmlout);}
 pop(xmlout);
} 

template <typename H, typename F, typename R, typename D>
std::map<int,F> DataPutHandlerMemMF<H,F,R,D>::getSuffixMap() const
{
 int suffix = 0;
 std::map<int,F> filekeys;
 for (typename FileMapType::const_iterator it=fileMap->begin();
      it!=fileMap->end();it++)
    filekeys.insert(std::make_pair(suffix,it->first));
 return filekeys;
} 


template <typename H, typename F, typename R, typename D>
std::set<F> DataPutHandlerMemMF<H,F,R,D>::getFileKeys() const
{
 std::set<F> filekeys;
 for (typename FileMapType::const_iterator it=fileMap->begin();
      it!=fileMap->end();it++)
    filekeys.insert(it->first);
 return filekeys;
} 


// **************************************************************
  }
}
#endif
