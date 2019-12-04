//Arjun Singh Gambhir
//Real authors: Ben Hoerz, John Bulava, Colin Morningstar, ??
//Code directly ported from chroma_laph, more acknowledgements coming later.
#ifndef IO_MAP_LAPH_H
#define IO_MAP_LAPH_H

#include "io_handler.h"
#include <map>
#include <set>
#include "output_help.h"

#define UINTSIZE     4
//#define SIZETSIZE    4
#define SIZETSIZE    8
#define ULONGLONG    8

namespace Chroma {
  namespace LaphEnv {


 // *********************************************************************************
 // *                                                                               *
 // *       class IOMap:    parallel random access input/output                     *
 // *                                                                               *
 // *   Author: Colin Morningstar (Carnegie Mellon University)                      *
 // *                                                                               *
 // *  Objects of class IOMap handle the input and output of one type of data of    *
 // *  class "V" using one type of key of class "K". Essentially, an IOMap object   *
 // *  is like a C++ map, but instead of dealing with memory, it deals with data    *
 // *  on a disk.  The steps in using the class are as follows:                     *
 // *                                                                               *
 // *   (1) define a class for the keys (see below)                                 *
 // *   (2) define a class for the values (see below)                               *
 // *   (3) create an IOMap object (in either global or local mode)                 *
 // *   (4) open a file by calling  openReadOnly, openNew, or openUpdate            *
 // *   (5) use "put" to insert data into the file                                  *
 // *   (6) use "get" to retrieve data from the file                                *
 // *   (7) close the file (or destroy the IOMap object)                            *
 // *                                                                               *
 // *    class RecordType{....};                                                    *
 // *    class DataType{....};                                                      *
 // *    global_mode=true;                                                          *
 // *    IOMap<RecordType,DataType> iomap(global_mode);                             *
 // *    K key(...); V val;                                                         *
 // *    iomap.get(key,val);                                                        *
 // *    iomap.put(key,val);                                                        *
 // *                                                                               *
 // *  More details on the usage are below.                                         *
 // *                                                                               *
 // *  The record key class should have the following features:                     *
 // *                                                                               *
 // *    (1) since used in a C++ map, a less than operator must be define           *
 // *              const K& K::operator<(const K& rhs);                             *
 // *    (2) a numbytes(ioh,K) function, where ioh is an IOHandler, must be         *
 // *         defined to give the number of bytes each key occupies in an           *
 // *         IOHandler file                                                        *
 // *    (3) a copy constructor K(const K& in) must be defined                      *
 // *          (a default constructor is not needed)                                *
 // *    (4) a multi_read(ioh, vector<K>&,n) must be defined to read n keys         *
 // *    (5) a multi_write(ioh, const vector<K>&) must be defined                   *
 // *                                                                               *
 // *  The number of bytes must be the same for all key values.  The keys are used  *
 // *  in a C++ map so keeping the keys small makes for a more efficient search.    *
 // *  Also, all of the keys stored in the file are read and written in a single    *
 // *  MPI read/write, so keeping the total size of all keys below the MPI eager    *
 // *  threshold (100KB or less) is reasonably important.                           *
 // *                                                                               *
 // *  The value type must have the following features:                             *
 // *                                                                               *
 // *   (1) a write(ioh, const V&) must be defined (ioh is an IOHandler object)     *
 // *   (2) a read(ioh, V&) must be defined                                         *
 // *   (3) a numbytes(ioh,V) must be defined giving number of bytes occupied       *
 // *           by V in an IOHandler file                                           *
 // *                                                                               *
 // *  The size of the value type does not need to be the same for all values. All  *
 // *  of the basic data types, multi1d<T>, multi2d<T>, multi3d<T>, vector<T>,      *
 // *  lattice OLattice<T> objects, and TimeSliceOf<OLattice<T> > already have      *
 // *  the above functions defined.                                                 *
 // *                                                                               *
 // *  An object of class IOHandler must be defined as global or local in the       *
 // *  constructor,                                                                 *
 // *                    IOMap(bool global_mode)                                    *
 // *  and this property cannot be changed.  If global, then it is assumed that     *
 // *  an object is created on every node.  In global mode, lattice and other       *
 // *  distributed data can be read and written, and common data is written by the  *
 // *  primary node and read by the primary node with subsequent broadcast to all   *
 // *  nodes.  In local mode, the object is assumed to exist only on the current    *
 // *  node.  Only non-distributed data can be read and written.  Read/writes are   *
 // *  done through the current node, and no data is broadcast to other nodes.      *
 // *                                                                               *
 // *  IOMap files can be opened using one of three open routines:                  *
 // *                                                                               *
 // *     (1) openReadOnly  -- fails if the file does not exist or read not allowed *
 // *     (2) openNew -- deletes any existing file then creates a new file          *
 // *     (3) openUpdate -- updates existing file or creates new                    *
 // *                                                                               *
 // *  A 32-character ID string is used to identify an IOMap file. You should       *
 // *  choose a string that is based on the key and value types.  During an open    *
 // *  of an existing file, an exact match of the ID string is needed or the open   *
 // *  fails.                                                                       *
 // *                                                                               *
 // *  During an open of a new file, a header string is written.  This can be of    *
 // *  any length.  During an open of an existing file, the header string is        *
 // *  read and returned (there is one read-only open routine that does not read    *
 // *  the header string).  NOTE: in global mode, only the primary node gets the    *
 // *  header string; all other nodes will only receive the empty string.           *
 // *                                                                               *
 // *  During an open of a new file, you can request little endian ('L') format,    *
 // *  big endian ('B') format, or native ('N') format.  You can use check sums     *
 // *  or not.  If check sums are included in a file, they will continue to be      *
 // *  included in any future insertions, even if not used.  During reading, you    *
 // *  can ignore checksums even if they are included in the file.                  *
 // *                                                                               *
 // *  Inserting new data adds the data to the file.  If you attempt to add data    *
 // *  whose key already exists in the file, the data in the file will be           *
 // *  overwritten if overwrites are allowed and if the size of the new data is     *
 // *  not larger than the size of the data in the file. To simplify matters, no    *
 // *  erase member is available.  If you really need to erase records, read a      *
 // *  file and copy the records you wish to keep to a new file.                    *
 // *                                                                               *
 // *  Upon closing a file, the map locations are written to the end of the file.   *
 // *  If the program aborts, the map locations do not get written and the file     *
 // *  is corrupted.  To prevent this, use the "flush" command which immediately    *
 // *  writes the current record locations at the end of the file.                  *
 // *                                                                               *
 // *  All errors that occur during a "get" or a "put" throw a string (so you can   *
 // *  output more meaning information). All other errors are fatal and cause       *
 // *  an abort.                                                                    *
 // *                                                                               *
 // *                                                                               *
 // *  A completed IOMap file contains (in the following order):                    *
 // *   (1) endian character 'L' or 'B'                                             *
 // *   (2) a 32-character ID string                                                *
 // *   (3) location where the map is stored in the file (8 bytes)                  *
 // *   (4) character 'Y' or 'N' describing if check sums stored in file            *
 // *   (5) header string of any length                                             *
 // *   (6) the data records one after another                                      *
 // *   (7) the file map describing locations of the records and keys               *
 // *   (8) an ending character 'E'                                                 *
 // *                                                                               *
 // *  Each record contains (in the following order):                               *
 // *   (1) the size in bytes of the record (excluding checksum)                    *
 // *   (2) the data itself in binary format                                        *
 // *   (3) the checksum of the data (optionally)                                   *
 // *                                                                               *
 // *                                                                               *
 // *********************************************************************************

/*
   //  Below is a sample class for an IOMap key.  Use this as a starting
   //  point for defining your own key class.


class TwoSpin
{

    unsigned int s1,s2;    // each value between 0 and 3 (say)
    
  public:
  
    TwoSpin(int in1, int in2);   // no default constructor by design !!
    
    TwoSpin(const TwoSpin& rhs) : s1(rhs.s1), s2(rhs.s2) {}
    
    TwoSpin& operator=(const TwoSpin& rhs)
     {s1=rhs.s1; s2=rhs.s2; return *this;}
    
    std::string output() const;
    
    bool operator<(const TwoSpin& rhs) const
    { return ((s1<rhs.s1) || ( (s1==rhs.s1)&&(s2<rhs.s2) ) ); }
    
    friend void multi_write(IOHandler& ioh, const vector<TwoSpin>& output);

};

   // no default constructor is needed
   
inline TwoSpin::TwoSpin(int in1, int in2)
{
 if ((in1<0)||(in1>3)||(in2<0)||(in2>3)){
    QDPIO::cerr << "Invalid TwoSpin initialization"<<endl;
    QDP_abort(1);}
 s1=static_cast<unsigned int>(in1);
 s2=static_cast<unsigned int>(in2);
}

inline std::string  TwoSpin::output() const
{
 stringstream oss;
 oss << "("<< s1 <<", "<< s2 <<")";
 return oss.str();
}

    // the size the key occupies in the IOHandler file (in bytes)
 
size_t numbytes(IOHandler& ioh, const TwoSpin& ss)
{
 return 2*sizeof(unsigned int);
}

   // copies TwoSpin members into an integer array, then does an
   // integer IOHandler multi_write (handles byte swapping)

void multi_write(IOHandler& ioh, const vector<TwoSpin>& output)
{
 int n=output.size();
 if (n<1) return;
 unsigned int *buf=new(nothrow) unsigned int[2*n];
 if (!buf){ 
    QDPIO::cerr << "could not allocate buffer for TwoSpin multiwrite"<<endl; 
    QDP_abort(1);}
 int k=0;
 for (vector<TwoSpin>::const_iterator it=output.begin();it!=output.end();it++){
    buf[k++]=it->s1; buf[k++]=it->s2;}
 ioh.multi_write(buf,2*n);   // int write handles byte swapping if needed
 delete [] buf;
}

   // does an integer IOHandler multi_read (which handle byte swapping,
   // then copies data into the TwoSpin members

void multi_read(IOHandler& ioh, vector<TwoSpin>& input, int n)
{
 input.clear();
 if (n<1) return;
 input.reserve(n);  // no default constructor needed here
 unsigned int *buf=new(nothrow) unsigned int[2*n];
 if (!buf){ 
    QDPIO::cerr << "could not allocate buffer for TwoSpin multiread"<<endl; 
    QDP_abort(1);}
 ioh.multi_read(buf,2*n);   // read into ints handles byte swapping if needed
 for (int k=0;k<n;k++)
    input.push_back(TwoSpin(buf[2*k],buf[2*k+1]));
 delete [] buf;
}
*/

// ******************************************************************

         //   Helper routines for simple key classes that are 
         //   a small number of unsigned ints.  Just ensure that
         //   the key class "T" has member functions
         //      int numints() const;  <- numints in each key  
         //      void copyTo(unsigned int *buf) const;
         //      T(const unsigned int *buf);   // constructor
         //      size_t numbytes() const <- number of bytes of each key


template <typename T>
void multi_write(IOHandler& ioh, const std::vector<T>& output)
{
 int n=output.size();
 if (n<1) return;
 int nint=output[0].numints();
 std::vector<unsigned int> buf(n*nint);
 int k=0;
 for (typename std::vector<T>::const_iterator it=output.begin();it!=output.end();it++){
    it->copyTo(&buf[k]); k+=nint;}
 ioh.multi_write(&buf[0],n*nint);   // int write handles byte swapping if needed
}

template <typename T>
void multi_read(IOHandler& ioh, std::vector<T>& input, int n)
{
 input.clear();
 if (n<1) return;
 input.reserve(n);  // no default constructor needed here
 int nint=input[0].numints();
 std::vector<unsigned int> buf(nint*n);
 ioh.multi_read(&buf[0],nint*n);   // read into ints handles byte swapping if needed
 for (int k=0;k<n*nint;k+=nint)
    input.push_back(T(&buf[k]));
}


// ***************************************************************************

     //  A simple integer key class that you can use right out of the box.

class UIntKey
{
   unsigned int value;

 public:

   UIntKey(int inval) : value(inval) {}
   UIntKey(const UIntKey& in) : value(in.value) {}
   UIntKey& operator=(const UIntKey& in) {value=in.value; return *this;}
   ~UIntKey() {}

   bool operator<(const UIntKey& rhs) const {return (value<rhs.value);}
   bool operator==(const UIntKey& rhs) const {return (value==rhs.value);}
   bool operator!=(const UIntKey& rhs) const {return (value!=rhs.value);}

   unsigned int getValue() const {return value;}

   void output(XmlWriter& xmlw) const 
   {push(xmlw,"Key");
    write(xmlw,"Value",getValue());
    pop(xmlw);}

   explicit UIntKey(const unsigned int* buf) {value=*buf;}
   int numints() const {return 1;} 
   size_t numbytes() const {return sizeof(unsigned int);}
   void copyTo(unsigned int* buf) const { *buf=value;}

};




 // ********************************************************
 // *                                                      *
 // *              The main event: the IOMap               *
 // *                                                      *
 // ********************************************************


template<typename K, typename V>
class IOMap
{

     typedef IOHandler::pos_type pos_type; 
     typedef IOHandler::off_type off_type; 
     typedef unsigned long long  w_pos;
    
     mutable IOHandler ioh;
     mutable std::map<K, pos_type> file_map;
     bool checksums_in_file;
     bool use_checksums;
     bool allow_overwrites;
     bool verbose1,verbose2;

     w_pos new_write_pos;
     bool write_map_on_close;

      // disallow copying
     IOMap(const IOMap&);
     IOMap(IOMap&);
     IOMap& operator=(const IOMap&);
     IOMap& operator=(IOMap&);


  public:
 
    typedef IOHandler::OpenMode OpenMode;

    IOMap(bool global_mode=true);

            // read only open, returns header string
            
    void openReadOnly(const std::string& filename, 
                      const std::string& filetype_id,
                      std::string& header, bool turn_on_checksum=false);

            // read only open, ignores header string

    void openReadOnly(const std::string& filename, 
                      const std::string& filetype_id,
                      bool turn_on_checksum=false);

            // open a new file in read/write mode, writes the header string (fails 
            // if the file exists and "fail_if_exists" is true; if "fail_if_exists"
            // is false, deletes the existing file to start a new file)

    void openNew(const std::string& filename, 
                 const std::string& filetype_id, 
                 const std::string& header,  
                 bool fail_if_exists=true, char endianness='N',
                 int striping_factor=1, int striping_unit=0,
                 bool turn_on_checksum=false, bool overwrites_allowed=false);

            // open a file in read/write mode; if file exists, the header
            // string is read and returned in "header" and writes will update
            // the existing file; otherwise, a new file is created (in which 
            // case, the header string is needed as input so it can be written
            // into the new file)

    void openUpdate(const std::string& filename, 
                    const std::string& filetype_id, 
                    std::string& header, 
                    char endianness='N', int striping_factor=1, 
                    int striping_unit=0, bool turn_on_checksum=false, 
                    bool overwrites_allowed=false);

    ~IOMap() {close();}

    void close();



    std::string getHeader();  // file must be open
    
        // Version that assumes file is not open; file closed afterwards.
        // Returns empty string if file cannot be opened.
        
    bool peekHeader(std::string& header, const std::string& filename, 
                    const std::string& filetype_id);
    
    std::string getFileName() const { return ioh.getFileName(); }

    bool isOpen() const { return ioh.isOpen(); }

    bool isNewFile() const { return ioh.isNewFile(); }

    bool isLocal() const { return ioh.isLocal(); }

    bool isGlobal() const { return ioh.isGlobal(); }

    bool isOverwriteOn() const { return allow_overwrites; }



    void setHighVerbosity();

    void setMediumVerbosity();
   
    void setNoVerbosity() { verbose1=verbose2=false; }

    void setDisallowOverwrites() { allow_overwrites=false; }

    void setAllowOverwrites() { allow_overwrites=true; }



    void put(const K& key, const V& val);
    
    void get(const K& key, V& val);

    bool exist(const K& key) const;
    
    void flush();  // puts file in finalized state so no data loss if abort occurs
     


    unsigned int size() { return file_map.size(); }

    void getKeys(std::vector<K>& keys) const;

    void getKeys(std::set<K>& keys) const;


    void outputContents(const std::string& logfile, bool printdata);


  private:

    void initialize(const std::string& filename, OpenMode mode,
                    const std::string& filetype_id, char endianness,
                    int striping_factor, int striping_unit,
                    bool turn_on_checksum, std::string& header, 
                    bool read_header, bool overwrites);

    void readMap(w_pos mapstart);
    
    void writeMap(w_pos mapstart);
    
    void check_for_failure(bool errcond, const std::string& msg, bool abort=true);

         // for static (compile time) assertion (produces compiler error if not satisfied)
   template <bool b>
   void staticassert()
   { typedef char asserter[b?1:-1]; }

};

    
    
template <typename K, typename V>
IOMap<K,V>::IOMap(bool global_mode) : ioh(global_mode), checksums_in_file(false), use_checksums(false), 
                                      allow_overwrites(false), verbose1(false), verbose2(false), 
                                      new_write_pos(0), write_map_on_close(false)
{
          // rely on specific sizes in the file so check these sizes
 staticassert<(sizeof(unsigned)==UINTSIZE)&&(sizeof(w_pos)==ULONGLONG)
               &&(sizeof(size_t)==SIZETSIZE)>();
}


template <typename K, typename V>
void IOMap<K,V>::openReadOnly(const std::string& filename, 
                              const std::string& filetype_id,
                              std::string& header, bool turn_on_checksum)
{
 check_for_failure(ioh.isOpen(),"Please do not open an already open IOMap!");
 OpenMode mode=IOHandler::ReadOnly;
 initialize(filename,mode,filetype_id,'N',1,0,turn_on_checksum, 
            header,true,false);
}



template <typename K, typename V>
void IOMap<K,V>::openReadOnly(const std::string& filename, 
                              const std::string& filetype_id,
                              bool turn_on_checksum)
{
 check_for_failure(ioh.isOpen(),"Please do not open an already open IOMap!");
 OpenMode mode=IOHandler::ReadOnly;
 std::string header;
 initialize(filename,mode,filetype_id,'N',1,0,turn_on_checksum, 
            header,false,false);
}


template <typename K, typename V>
void IOMap<K,V>::openNew(const std::string& filename, 
                         const std::string& filetype_id, 
                         const std::string& header,  
                         bool fail_if_exists, char endianness,
                         int striping_factor, int striping_unit,
                         bool turn_on_checksum, bool overwrites_allowed)
{
 check_for_failure(ioh.isOpen(),"Please do not open an already open IOMap!");
 OpenMode mode=(fail_if_exists)?IOHandler::ReadWriteFailIfExists : IOHandler::ReadWriteEraseIfExists;
 std::string tmp_header(header);
 initialize(filename,mode,filetype_id,endianness,striping_factor,striping_unit,
            turn_on_checksum,tmp_header,false,overwrites_allowed);
}


template <typename K, typename V>
void IOMap<K,V>::openUpdate(const std::string& filename, 
                            const std::string& filetype_id, 
                            std::string& header, 
                            char endianness, int striping_factor, 
                            int striping_unit, bool turn_on_checksum, 
                            bool overwrites_allowed)
{
 check_for_failure(ioh.isOpen(),"Please do not open an already open IOMap!");
 OpenMode mode=IOHandler::ReadWriteUpdateIfExists;
 initialize(filename,mode,filetype_id,endianness,striping_factor,striping_unit,
            turn_on_checksum,header,true,overwrites_allowed);
}



template <typename K, typename V>
void IOMap<K,V>::initialize(const std::string& filename, OpenMode mode,
                            const std::string& filetype_id, char endianness,
                            int striping_factor, int striping_unit,
                            bool turn_on_checksum, std::string& header, 
                            bool read_header, bool overwrites)
{
 ioh.open(filename,mode,filetype_id,endianness,striping_factor,
          striping_unit,turn_on_checksum);
 allow_overwrites=overwrites;

 if (ioh.isNewFile()){                // this is a new file
    if (verbose1) std::cout << "IOMap: opening new file "<<filename<<std::endl;
    char cksum=(turn_on_checksum)?'Y':'N';
    checksums_in_file=use_checksums=turn_on_checksum;
    new_write_pos=0;
    ioh.write(new_write_pos);  // placeholder for the map start
    ioh.write(cksum);
    write(ioh,header);
    new_write_pos=ioh.tell();
    write_map_on_close=true;}
 else{                         // this is an existing file
    if (verbose1) std::cout << "IOMap: opening already existing file "<<filename<<std::endl;
    ioh.rewind();
    ioh.read(new_write_pos);
    char cksum;
    ioh.read(cksum);
    checksums_in_file=(cksum=='Y')?true:false;
    use_checksums=turn_on_checksum && checksums_in_file;
    if (!use_checksums) ioh.turnOffChecksum();
    if (read_header) ioh.read(header,false);
    write_map_on_close=false;
    readMap(new_write_pos);
    }

 if (verbose2){
    if (ioh.isChecksumOn()) std::cout << "   Checksums are ON"<<std::endl;
    else std::cout << "   Checksums are OFF"<<std::endl;
    if (ioh.isFileLittleEndian()) std::cout << "   File is little endian"<<std::endl;
    else std::cout << "   File is big endian"<<std::endl;
    if (ioh.isEndianConversionOn()) std::cout << "   Endian conversion is ON"<<std::endl;
    else std::cout << "   Endian conversion is OFF"<<std::endl;
    if (ioh.isReadOnly()) std::cout << "   File is opened in read-only mode"<<std::endl;
    else std::cout << "   File is opened in read-write mode"<<std::endl;
    if (ioh.isGlobal()) std::cout << "   IOMap is global"<<std::endl;
    else std::cout << "   IOMap is local on node "<<Layout::nodeNumber()<<std::endl;
    std::cout << "   Number of key-value pairs is "<<file_map.size()<<std::endl;}
}


template <typename K, typename V>
void IOMap<K,V>::close()
{
 if (!ioh.isOpen()) return;
 if (write_map_on_close) writeMap(new_write_pos);
 if (verbose1) std::cout << "IOMap: closed file "<<ioh.getFileName()<<std::endl;
 ioh.close();
 file_map.clear();
}



template <typename K, typename V>
void IOMap<K,V>::flush()
{
 if (!ioh.isOpen()) return;
 if (write_map_on_close) writeMap(new_write_pos);
 if (verbose1) std::cout << "IOMap: flushed file "<<ioh.getFileName()<<std::endl;
 write_map_on_close=false;
}


template<typename K, typename V>
void IOMap<K,V>::writeMap(w_pos mapstart)
{
 if (checksums_in_file) ioh.turnOnChecksum();
 ioh.seekFromStart(static_cast<off_type>(mapstart));
 size_t nrecords=file_map.size();
 ioh.write(nrecords);

 size_t keysize;
 if (nrecords==0) keysize=0;
 else keysize=numbytes(ioh,file_map.begin()->first); 
 ioh.write(keysize);
 
 std::vector<K> keybuf; keybuf.reserve(nrecords);
 w_pos* posbuf=new(std::nothrow) w_pos[nrecords];
 check_for_failure((!posbuf),"allocation error in IOMap");
 w_pos* posptr=posbuf;
 for (typename std::map<K,pos_type>::const_iterator it=file_map.begin();
          it!=file_map.end();it++){
    keybuf.push_back(it->first);
   *posptr=it->second; posptr++;}
 multi_write(ioh,keybuf);
 ioh.multi_write(posbuf,nrecords);
 delete [] posbuf;
 check_for_failure((size_t(ioh.tell())-size_t(mapstart))!=(nrecords*(keysize+sizeof(mapstart))
                    +2*sizeof(size_t)),"bad write of file map");

 if (checksums_in_file){
    QDPUtil::n_uint32_t checksum=ioh.getChecksum();
    ioh.write(checksum);}
 if (!use_checksums) ioh.turnOffChecksum();
 char eof='E';
 ioh.write(eof);
 ioh.rewind(); 
 ioh.write(mapstart);
}




template<typename K, typename V>
void IOMap<K,V>::readMap(w_pos mapstart)
{
 check_for_failure(mapstart==0,"IOMap file "+ioh.getFileName()+" corrupted from prior aborted execution");
 ioh.seekFromStart(static_cast<off_type>(mapstart));
 size_t nrecords;
 ioh.read(nrecords);
 check_for_failure(nrecords>16777216,"Too many records during readMap: bad read?");

 size_t keysize;
 ioh.read(keysize);
 check_for_failure(keysize>1024,"Key size too large during readMap: bad read?");

 if (nrecords>0){
    std::vector<K> keybuf; keybuf.reserve(nrecords);
    w_pos* posbuf=new(std::nothrow) w_pos[nrecords];
    check_for_failure((!posbuf),"allocation error in IOMap");
    multi_read(ioh,keybuf,nrecords);
    ioh.multi_read(posbuf,nrecords);
    for (unsigned int k=0;k<nrecords;k++){
       file_map.insert(std::make_pair(keybuf[k],static_cast<pos_type>(posbuf[k])));}
    delete [] posbuf;
    check_for_failure(keysize!=numbytes(ioh,file_map.begin()->first),
                     "Bad keysize in reading file map");}
 check_for_failure((size_t(ioh.tell())-size_t(mapstart))!=(nrecords*(keysize+sizeof(mapstart))
                    +2*sizeof(size_t)),"bad read of file map");

 if (use_checksums){
    QDPUtil::n_uint32_t checksumA=ioh.getChecksum();
    QDPUtil::n_uint32_t checksumB;
    ioh.read(checksumB);
    check_for_failure(checksumA!=checksumB,"checksum mismatch: bad read of file map");}
 else if (checksums_in_file)
    ioh.seekFromCurr(sizeof(QDPUtil::n_uint32_t));

 char eof;
 ioh.read(eof);
 check_for_failure(eof!='E',"bad read of file map");
}


template<typename K, typename V>
void IOMap<K,V>::setHighVerbosity() 
{ 
 verbose1=verbose2= (ioh.isLocal()) || (Layout::primaryNode());
}


template<typename K, typename V>
void IOMap<K,V>::setMediumVerbosity() 
{
 verbose1= (ioh.isLocal()) || (Layout::primaryNode());
 verbose2=false;
}


template<typename K, typename V>
bool IOMap<K,V>::exist(const K& key) const 
{
 return (file_map.find(key) == file_map.end()) ? false : true;
}
  

template<typename K, typename V>
std::string IOMap<K,V>::getHeader()
{
 check_for_failure(!ioh.isOpen(),"IOMap not open: cannot get header");
 std::string tmp;
 ioh.seekFromStart(sizeof(w_pos)+1);
 ioh.read(tmp,false);
 return tmp;
}


template<typename K, typename V>
bool IOMap<K,V>::peekHeader(std::string& header,
                            const std::string& filename, 
                            const std::string& filetype_id)
{
 return ioh.peekString(header,sizeof(w_pos)+1,filename,filetype_id);
}
 
 
template<typename K, typename V>
void IOMap<K,V>::put(const K& key, const V& val) 
{
 check_for_failure((!ioh.isOpen())||(ioh.isReadOnly()),"Write operation disallowed",false);
 w_pos start_write_pos=new_write_pos;
 size_t sz=numbytes(ioh,val);
 bool flag=exist(key);
 if (flag){
    check_for_failure(!allow_overwrites,"key already in map and overwrites disallowed",false);
    start_write_pos=static_cast<w_pos>(file_map[key]);
    ioh.seekFromStart(static_cast<off_type>(start_write_pos));
    size_t cursize; 
    ioh.read(cursize);
    check_for_failure((cursize<sz),"can only overwrite record if new size is not larger",false);
    if (verbose2) std::cout << "IOMap::put overwriting existing record at file location "
                            <<start_write_pos<<std::endl;
    }
 else{
    if (verbose2) std::cout << "IOMap::put of new record at file location "
                            <<start_write_pos<<std::endl;
    write_map_on_close=true;
    file_map.insert(std::make_pair(key,start_write_pos));}
 if (checksums_in_file) ioh.turnOnChecksum();
 ioh.seekFromStart(static_cast<off_type>(start_write_pos));
 ioh.write(sz);
 write(ioh,val);
 size_t recordsize=sizeof(size_t)+sz;
 if (checksums_in_file){
    QDPUtil::n_uint32_t checksum=ioh.getChecksum();
    ioh.write(checksum);
    recordsize+=sizeof(checksum);}
 if (!use_checksums) ioh.turnOffChecksum();
 pos_type curpos=ioh.tell();
 check_for_failure((size_t(curpos)-size_t(start_write_pos))!=recordsize,
                   "bad value write in IOMap",false);
 if (!flag){
    new_write_pos=static_cast<w_pos>(curpos);
    ioh.rewind(); w_pos zero=0;
    ioh.write(zero);}
} 


template<typename K, typename V>
void IOMap<K,V>::get(const K& key, V& val) 
{
 check_for_failure((!ioh.isOpen())||(!exist(key)),
                   "Read failed: file not open or key not in file",false);
 off_type start = static_cast<off_type>(file_map.find(key)->second); 
 if (verbose2) std::cout << "IOMap::get at file location "<<start<<std::endl;
 ioh.seekFromStart(start);
 size_t sz;
 ioh.read(sz);
 read(ioh,val);
 check_for_failure(sz!=numbytes(ioh,val),"value size mismatch",false);
 check_for_failure((size_t(ioh.tell())-size_t(start))!=(sizeof(size_t)+sz),
                   "bad value read in IOMap",false);
 if (use_checksums){
    QDPUtil::n_uint32_t checksumA=ioh.getChecksum();
    QDPUtil::n_uint32_t checksumB;
    ioh.read(checksumB);
    check_for_failure(checksumA!=checksumB,"checksum mismatch",false);}
}


template<typename K, typename V>
void IOMap<K,V>::getKeys(std::vector<K>& keys) const 
{
 keys.clear();
 keys.reserve(file_map.size());
 for (typename std::map<K,pos_type>::const_iterator it  = file_map.begin();
      it != file_map.end(); ++it){
      keys.push_back(it->first);}
}


template<typename K, typename V>
void IOMap<K,V>::getKeys(std::set<K>& keys) const 
{
 keys.clear();
 for (typename std::map<K,pos_type>::const_iterator it  = file_map.begin();
      it != file_map.end(); ++it){
      keys.insert(it->first);}
}


template<typename K, typename V>
void IOMap<K,V>::check_for_failure(bool errcond, const std::string& mesg,
                                   bool abort)
{
 if (!errcond) return;
 if ((Layout::primaryNode())||(ioh.isLocal())){
    std::cerr << "IOMap error: "<<mesg<<std::endl;}
 if (abort) QDP_abort(1);
 else throw mesg;
}



   // **************************************************************
   // *                                                            *
   // *                         Outputter                          *
   // *                                                            *
   // **************************************************************


   // to use this routine, there must exist a routine
   //    void output(const D& data, TextFileWriter& tout)
   // outputs only the header XML and the record keys by default;
   // if "printdata" is true, then the data itself is output


template <typename R, typename D>
void IOMap<R,D>::outputContents(const std::string& logfile, bool printdata)
{
 check_for_failure(!ioh.isOpen(),"IOMap not open: cannot output contents");

 std::string filename=getFileName();
 std::string headerxml=getHeader();
 
 TextFileWriter tout(logfile);
 tout << "IOMap file name: "<<filename<<"\n";
 tout << "header info:\n"<<headerxml<<"\n";

 std::set<R> dbkeys; 
 getKeys(dbkeys); 

 for (typename std::set<R>::const_iterator 
            kt=dbkeys.begin();kt!=dbkeys.end();kt++){
    tout << "Record Key:\n";
    XmlBufferWriter xmlout;
    kt->output(xmlout);
    tout << xmlout.str() <<"\n";
    if (printdata){
       tout << "Record Data:\n";
       D result;
       get(*kt,result);
       output(result,tout); tout <<"\n\n\n";}
    }

 tout.close(); 
}

// **************************************************************
  }
}
#endif
