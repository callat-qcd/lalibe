//Arjun Singh Gambhir
//Real authors: Ben Hoerz, John Bulava, Colin Morningstar, ??
//Code directly ported from chroma_laph, more acknowledgements coming later.
#ifndef IO_HANDLER_H
#define IO_HANDLER_H

#define MPICH_IGNORE_CXX_SEEK

#include "qdp.h"
#include "chromabase.h"
#include <list>
#include <vector>

#include <unistd.h> // 2014-01-22 new compiler issue?

#ifdef ARCH_SCALAR
#include <fstream>
#else
#include <mpi.h>
#endif

namespace Chroma {
  namespace LaphEnv {



 // *********************************************************************************
 // *                                                                               *
 // *       class IOHandler:    parallel random access input/output                 *
 // *                                                                               *
 // *   Author: Colin Morningstar (Carnegie Mellon University)                      *
 // *                                                                               *
 // *   This class is similar to the BinaryFileReaderWriter class.  It differs      *
 // *   from that class in a few ways: (a) in a parallel architecture, it uses      *
 // *   MPI-IO so that input/output is parallel and much more efficient; (b) it     *
 // *   allows the user to choose between big-endian and little-endian format;      *
 // *   (c) it allows the user to turn check sums on or off; (d) MPI-IO striping    *
 // *   can be used in a parallel architecture; (e) objects are constructed and     *
 // *   subsequently operate in either global or local mode.  Currently, only       *
 // *   blocking I/O operations are supported, but it would be simple to implement  *
 // *   some of MPI-IO's non-blocking operations.                                   *
 // *                                                                               *
 // *   An object of class IOHandler must be defined as global or local in the      *
 // *   constructor, and this property cannot be changed.  If global, then it is    *
 // *   assumed that an object is created on every node.  In global mode, lattice   *
 // *   and other distributed data can be read and written, and common data is      *
 // *   written by the primary node and read by the primary node with subsequent    *
 // *   broadcast to all nodes.  In local mode, the object is assumed to exist      *
 // *   only on the current node.  Only non-distributed data can be read and        *
 // *   written.  Read/writes are done through the current node, and no data is     *
 // *   broadcast to other nodes.                                                   *
 // *                                                                               *
 // *   Files written by this class will start with a single character 'B' or 'L'   *
 // *   to indicate endian-ness.  Then an ID string is output, which always has     *
 // *   the length "ID_string_length"=32 (spaces are padded if needed).  Files read *
 // *   by this class expect this starting behavior.  An error occurs if the        *
 // *   ID string given in the open command on an existing file does not match      *
 // *   that in the file.                                                           *
 // *                                                                               *
 // *   All errors are considered fatal in this class, so the end user need not     *
 // *   perform any error checking.  The ID string is useful for ensuring that      *
 // *   the file contains the kind of information that the user is expecting.       *
 // *                                                                               *
 // *   Files can be opened in one of four modes:                                   *
 // *     (1) ReadOnly  -- fails if the file does not exist or read not allowed     *
 // *     (2) ReadWriteFailIfExists --  fails if the file exists                    *
 // *     (3) ReadWriteEraseIfExists  -- erases the file if it exists               *
 // *     (4) ReadWriteUpdateIfExists  -- updates existing file or creates new      *
 // *   Modes 2,3,4 will create a new file if it does not exist.  When creating     *
 // *   a new file, the endian format choices are 'B', 'L', or 'N' (native).        *
 // *   Random access to the file is employed.                                      *
 // *                                                                               *
 // *   There are separate open routines for the different modes, or you can        *
 // *   use the OpenMode enum in the general open routine.  Alternatively, you      *
 // *   can use iostream openmodes as follows:                                      *
 // *      ReadOnly                  =   ios::in                                    *
 // *      ReadWriteFailIfExists     =   ios::in | ios::out                         *
 // *      ReadWriteEraseIfExists    =   ios::in | ios::out | ios::trunc            *
 // *                                                                               *
 // *   There are seek and tell members for moving around in the file.  It is an    *
 // *   error to seek before the start of the data, but you can seek past the end   *
 // *   of the file.  You can seek relative to the start of the data (33 characters *
 // *   past the start of the file), the current location, or the end of file.      *
 // *   Use negative offsets relative to the end of file to go backwards.           *
 // *                                                                               *
 // *   Input is done with read(..) commands, and output is done with write(..)     *
 // *   commands. For lattice objects, collective I/O is used (in global mode). For *
 // *   non-lattice objects, the read() and write() commands are either done by     *
 // *   the primary node (global mode) or the current node (local mode).            *
 // *   The data types currently supported for read/write are                       *
 // *                                                                               *
 // *       - basic data types (int, bool, float, string, etc.)                     *
 // *       - multi1d, multi2d, multi3d, vector of basic data types                 *
 // *       - OScalar<T> quantities                                                 *
 // *       - multi1d, multi2d, multi3d, vector of OScalar<T> quantities            *
 // *       - OLattice<T> objects                                                   *
 // *       - multi1d and vector of OLattice<T> objects                             *
 // *       - TimeSliceOf<OLattice<T> > objects                                     *
 // *       - vector<TimeSliceOf<OLattice<T> > > objects                            *
 // *                                                                               *
 // *   The quantities above are contiguous in memory and regular (each element is  *
 // *   the same size) and these facts have been used to speed up I/O.              *
 // *                                                                               *
 // *   (a) For lattice objects, the read() and write() commands do collective I/O  *
 // *   in MPI-IO. All nodes participate so all nodes must call the read/write      *
 // *   command.  Before the read/write call, ensure that all file pointers are set *
 // *   to the same point in the file (such as by a seek call by all nodes).  After *
 // *   each read/write, all file pointers are advanced by the size of the          *
 // *   lattice object.  Only applies to global mode.                               *
 // *                                                                               *
 // *   (b) For reading/writing non-lattice objects in global mode, all nodes must  *
 // *   call the read/write subroutine, but I/O is only done on the primary node.   *
 // *   For reads, the results are then broadcast to all nodes.  All file pointers  *
 // *   are advanced by the size of the object after the I/O operation.  Before the *
 // *   read/write, it is safest to have all nodes call a seek operation so that    *
 // *   all file pointers are initially pointing to the same location in the file.  *
 // *                                                                               *
 // *   (c) Reading/writing non-lattice objects in "local" mode must be done        *
 // *   carefully.  Typically, a seek is called only on one node, and the local     *
 // *   read/write is only called by the one node.  The file pointer for that       *
 // *   one node will be advanced afterwards by the size of the object. For local   *
 // *   reads, the results are NOT broadcast to all other nodes.                    *
 // *   Simultaneous local writes to the same file can cause major problems with    *
 // *   consistency!!  The "local" mode is more meant to be used during             *
 // *   simultaneous writes/reads to **different** files.                           *
 // *   Use local mode cautiously, especially when writing data.                    *
 // *                                                                               *
 // *   (d) For a time slice of a lattice object, use the readTimeSlice()           *
 // *   and writeTimeSlice() members or read(TimeSliceOf< >), write(TimeSliceOf<>). *
 // *   Again, collective I/O is used so all nodes must participate by calling      *
 // *   these read/write routines. Global mode only. NOTE: reading time slices      *
 // *   does NOT allocate memory; the full lattice objects must already exist, and  *
 // *   the read puts the time slice into the appropriate memory locations.  All    *
 // *   other reads DO allocate new memory if necessary.                            *
 // *                                                                               *
 // *   Examples of writing and reading:                                            *
 // *                                                                               *
 // *      IOHandler io;     <-- default is global mode                             *
 // *      io.open(....);                                                           *
 // *      int k=5;     io.seekFromStart(...); write(io,k);                         *
 // *      float x=5.4; io.seekFromStart(...); write(io,x);                         *
 // *      int j; float y;  io.seekFromStart(...); read(io,j,y); // in sequence     *
 // *      string str("ABC"); io.seekFromStart(...); write(io,str);                 *
 // *      multi1d<float> g(1024); io.seekFromStart(...); write(io,g);              *
 // *      LatticeComplex zn;  io.seekFromStart(...); write(io,z);                  *
 // *                                                                               *
 // *   When reading and writing distributed arrays, a column-major order is        *
 // *   assumed in the file, and column-major order of the processor is used.       *
 // *                                                                               *
 // *   For more complicated smaller structures consisting of the 11 basic          *
 // *   fixed-size data types, reads/writes can be delayed to combine them into     *
 // *   a single read/write.                                                        *
 // *                                                                               *
 // *      int ik=8;                                                                *
 // *      double x=7.342;                                                          *
 // *      float y=-1.2;                                                            *
 // *      bool b=true;                                                             *
 // *      char c='G';                                                              *
 // *      unsigned int p=5;                                                        *
 // *      IOHandler ioh(...);                                                      *
 // *      ioh.startDelayedWrite();    <-- internal list of pointers emptied        *
 // *      ioh.addDelayedWrite(ik);    <-- pointer to "ik" added to internal list   *
 // *      ioh.addDelayedWrite(x,y);   <-- pointers to "x", "y" added to list       *
 // *      ioh.addDelayedWrite(b,c,p);                                              *
 // *      ioh.finishDelayedWrite();   <-- all data copied to internal character    *
 // *                                       buffer, then written to file            *
 // *   Similarly, reads can be delayed.                                            *
 // *   The "add" functions add pointers to an internal list, so do NOT change      *
 // *   the data being pointed to until after the "finish".  All copying, byte      *
 // *   swapping, check sum updating, is done in the "finish" routine.              *
 // *   These delayed reads/writes are useful for handling the IO for keys in       *
 // *   maps.                                                                       *
 // *                                                                               *
 // *   If check sums are turned on, objects of this class maintain a check sum.    *
 // *   Doing an explicit seek resets the checksum; changing from a read to a       *
 // *   write or vice versa resets the checksum; or an explicit reset can also      *
 // *   be done.  The checksum is updated as bytes are read or written.             *
 // *   Successive reads or successive writes update the checksum.                  *
 // *                                                                               *
 // *   For local read/writes, create an IOHandler object on one processor in       *
 // *   local mode:                                                                 *
 // *                                                                               *
 // *      bool global_mode=false;                                                  *
 // *      IOHandler localhandler(global_mode);    <-  local write                  *
 // *                                                                               *
 // *   Use local I/O **very carefully**.                                           *
 // *                                                                               *
 // *********************************************************************************


 //  IOHandler class is defined later in this file, after the 
 //  helper class DistArrayViewInfo.


#ifdef ARCH_SCALAR
   typedef size_t       IOH_int; 
#else
   typedef MPI_Offset   IOH_int;
#endif

 
  /*  DistArrayViewInfo objects contain info about multi-dimensional arrays 
      distributed across the MPI nodes.  Blocking assumed, column major.   
      Assumed distributed in a regular way...same size on each node.  
      Distribution on nodes is column major.  This class is mainly
      used by IOHandler for creating file views for MPI-IO. It supports
      entire global arrays, and single slices of the most major index
      (the right-most index).  The end user does not need to know
      anything about this class.
      
      The main purpose of this class is to create an appropriate "filetype"
      to pass to MPI-IO.  The routine  MPI_Type_create_subarray is used.  
      
      Warning:  This class stores an MPI quantity, so make sure all
      destructors get called before MPI_finalize().   */

#ifndef ARCH_SCALAR
                     //   parallel version

class DistArrayViewInfo
{

    int ndim;
    multi1d<int> global_sizes, local_sizes, start_indices;
    IOH_int nelem_this_node, gvvol, lexico_start, el_offset_this_node;
    IOH_int nbytes_per_element;
    MPI_Datatype ftype;
    MPI_Datatype etype;
    
  public:
  
    DistArrayViewInfo(IOH_int bytes_per_site);   // assumes a full lattice
    
    DistArrayViewInfo(int time_slice, IOH_int bytes_per_site);   // for a lattice time slice

                 // this lets MPI determine how to distributes array on nodes
    DistArrayViewInfo(const multi1d<int>& gsizes, IOH_int bytes_per_element);
    
    DistArrayViewInfo(const multi1d<int>& gsizes, const multi1d<int>& numnodes,
                      IOH_int bytes_per_element);

    DistArrayViewInfo(const multi1d<int>& gsizes, const multi1d<int>& numnodes,
                      int major_index, IOH_int bytes_per_element);

    void resetBytes(IOH_int bytes_per_site);

    DistArrayViewInfo(const DistArrayViewInfo& in);

    DistArrayViewInfo& operator=(const DistArrayViewInfo& in);

    ~DistArrayViewInfo();

   

    int getArrayDimensions() const
     { return ndim; }
     
    const multi1d<int>& getGlobalArraySizes() const
     { return global_sizes; }

    const multi1d<int>& getLocalArraySizes() const
     { return local_sizes; }

    const multi1d<int>& getViewLocalStarts() const
     { return start_indices; }
     
    IOH_int getBytesPerSite() const
     { return nbytes_per_element; }
     
    IOH_int getBytesPerElement() const
     { return nbytes_per_element; }
     
    IOH_int getViewElementsThisNode() const
     { return nelem_this_node; }

    IOH_int getViewBytesThisNode() const
     { return nbytes_per_element*nelem_this_node; }

    IOH_int getElementOffsetThisNode() const
     { return el_offset_this_node; }

    IOH_int getViewTotalBytes() const
     { return gvvol*nbytes_per_element; }
     
    IOH_int getViewLexicoStart() const   // in terms of elements
     { return lexico_start; }        // helps with check sums
    
    IOH_int getViewLexicoSpan() const    // in terms of elements
     { return gvvol; }               // helps with check sums


    const MPI_Datatype& getFileViewType() const
     { return ftype; }

    const MPI_Datatype& getFileElemType() const
     { return etype; }
    

 
 private:
 
    void setup_full(const multi1d<int>& gsizes, const multi1d<int>& numnodes,
                    const multi1d<int>& logical_coord);

    void setup_major_slice(const multi1d<int>& gsizes, const multi1d<int>& numnodes,
                           const multi1d<int>& logical_coord, int major_index);
    
    void set_bytes(IOH_int bytes_per_element);
    
    bool checker(const multi1d<int>& gsizes, const multi1d<int>& numnodes);
    
    void error_return(bool errcond, const std::string& mesg);


    void create_filetype();
     
};

#else
                     //   serial version
 
class DistArrayViewInfo
{

    int ndim;
    multi1d<int> local_sizes;
    IOH_int gvvol, lexico_start;
    IOH_int nbytes_per_element;
    
  public:
  
    DistArrayViewInfo(IOH_int bytes_per_site);   // assumes a full lattice
    
    DistArrayViewInfo(int time_slice, IOH_int bytes_per_site);   // for a lattice time slice

                 // this lets MPI determine how to distributes array on nodes
    DistArrayViewInfo(const multi1d<int>& gsizes, IOH_int bytes_per_element);
    
             // numnodes ignored below...for compatibility with parallel compiles
    DistArrayViewInfo(const multi1d<int>& gsizes, const multi1d<int>& numnodes,
                      IOH_int bytes_per_element);

    DistArrayViewInfo(const multi1d<int>& gsizes, const multi1d<int>& numnodes,
                      int major_index, IOH_int bytes_per_element);


    void resetBytes(IOH_int bytes_per_site);

    DistArrayViewInfo(const DistArrayViewInfo& in);

    DistArrayViewInfo& operator=(const DistArrayViewInfo& in);

    ~DistArrayViewInfo();

   
   
    int getArrayDimensions() const
     { return ndim; }
     
    const multi1d<int>& getGlobalArraySizes() const
     { return local_sizes; }

    const multi1d<int>& getLocalArraySizes() const
     { return local_sizes; }

    const multi1d<int> getViewLocalStarts() const
     { multi1d<int> tmp(ndim); tmp=0; return tmp; }

    IOH_int getBytesPerSite() const
     { return nbytes_per_element; }
     
    IOH_int getBytesPerElement() const
     { return nbytes_per_element; }
     
    IOH_int getViewElementsThisNode() const
     { return gvvol; }

    IOH_int getViewBytesThisNode() const
     { return nbytes_per_element*gvvol; }

    IOH_int getElementOffsetThisNode() const
     { return lexico_start; }

    IOH_int getViewTotalBytes() const
     { return gvvol*nbytes_per_element; }
     
    IOH_int getViewLexicoStart() const   // in terms of elements
     { return lexico_start; }        // helps with check sums
    
    IOH_int getViewLexicoSpan() const    // in terms of elements
     { return gvvol; }               // helps with check sums



 
 private:
 
    void setup_full(const multi1d<int>& gsizes);

    void setup_major_slice(const multi1d<int>& gsizes, int major_index);
    
    void set_bytes(IOH_int bytes_per_element);
    
    bool checker(const multi1d<int>& gsizes);
    
    void error_return(bool errcond, const std::string& mesg);

     
};

#endif


  // **********************************************************************
  
  
  //   A little helper class for handling lattice time slices
  //   Use only with lattice quantities:  e.g., TimeSliceOf<LatticeComplex>.
  //   Contains a pointer to the underlying full lattice quantity.  You 
  //   can make a vector of these objects.
  

template <typename T>
class TimeSliceOf
{
   T* data;
   int current_time;

 public:

   TimeSliceOf(T& lattice, int time_slice = 0);
   TimeSliceOf(const TimeSliceOf& rhs) : data(rhs.data), current_time(rhs.current_time) {}
   TimeSliceOf& operator=(const TimeSliceOf& rhs)
    {data=rhs.data; current_time=rhs.current_time; return *this;} 
   void setCurrentTime(int time_slice);
   int getCurrentTime() const {return current_time;}
   T& getData() {return *data;}
   const T& getData() const {return *data;}
   ~TimeSliceOf() {}

 private:
 
   struct OLatticeType {template<typename K> OLatticeType(const OLattice<K>& in){} };

};


template <typename T>
TimeSliceOf<T>::TimeSliceOf(T& lattice, int time_slice) : data(&lattice)
{
 OLatticeType dummy(lattice);  // compiler error if not lattice type
 setCurrentTime(time_slice);
}

template <typename T>
void TimeSliceOf<T>::setCurrentTime(int time_slice)
{
 current_time=time_slice;
 if ((current_time<0)||(current_time>=Layout::lattSize()[QDP::Nd-1])){
    throw(std::string("invalid time in Lattice TimeSliceOf"));}
}




 // ***********************************************************
 // *                                                         *
 // *          Now for the main event:  IOHandler             *
 // *                                                         *
 // ***********************************************************


class IOHandler
{

#ifdef ARCH_SCALAR
   std::fstream fh;
#else
   MPI_File fh;        // the MPI-IO file handler
   MPI_Info finfo;     // MPI-IO file hints
#endif

   bool read_only;
   bool openflag;
   bool read_mode;
   bool global_mode;
   
   char endian_format;            // 'B' for big-endian, 'L' for little-endian
   bool endian_convert;

   bool checksum_on;
   std::string m_filename;
   bool is_new_file;
 
   QDPUtil::n_uint32_t checksum;
   DistArrayViewInfo lattinfo;

      // disallow copying
   IOHandler(const IOHandler&);
   IOHandler(IOHandler&);
   IOHandler& operator=(const IOHandler&);
   IOHandler& operator=(IOHandler&);


 public:
 
   enum OpenMode { ReadOnly, ReadWriteFailIfExists, ReadWriteEraseIfExists, 
                   ReadWriteUpdateIfExists };

   explicit IOHandler(bool global=true);

   IOHandler(const std::string& filename, OpenMode mode=ReadOnly,
             const std::string& filetype_id="", char endianness='N',
             int striping_factor=1, int striping_unit=0,
             bool turn_on_checksum=false, bool global=true);

   IOHandler(const std::string& filename, std::ios_base::openmode mode,
             const std::string& filetype_id="", char endianness='N',
             int striping_factor=1, int striping_unit=0,
             bool turn_on_checksum=false, bool global=true);

   void open(const std::string& filename, OpenMode mode=ReadOnly,
             const std::string& filetype_id="", char endianness='N',
             int striping_factor=1, int striping_unit=0,
             bool turn_on_checksum=false);

   void open(const std::string& filename, std::ios_base::openmode mode,
             const std::string& filetype_id="", char endianness='N',
             int striping_factor=1, int striping_unit=0,
             bool turn_on_checksum=false);

   void openReadOnly(const std::string& filename, const std::string& filetype_id="",
                     bool turn_on_checksum=false);

   void openNew(const std::string& filename, bool fail_if_exists=true,
                const std::string& filetype_id="", char endianness='N',
                int striping_factor=1, int striping_unit=0,
                bool turn_on_checksum=false);

   void openUpdate(const std::string& filename, const std::string& filetype_id="", 
                   char endianness='N', int striping_factor=1, int striping_unit=0,
                   bool turn_on_checksum=false);

   ~IOHandler();

          // closes current file if open, otherwise no action taken
   void close();


          // informational routines
          
   bool isOpen() const { return openflag; }
   
   bool isNewFile() const { return is_new_file; }

   std::string getFileName() const { return m_filename; }
   
   bool isChecksumOn() const { return checksum_on; }
   
   bool isEndianConversionOn() const { return endian_convert; }
   
   bool isFileLittleEndian() const { return (endian_format=='L'); }
   
   bool isFileBigEndian() const { return (endian_format=='B'); }

   bool isReadOnly() const { return read_only; }

   bool isGlobal() const { return global_mode; }
   
   bool isLocal() const { return !global_mode; }
   

#ifdef ARCH_SCALAR
   typedef std::iostream::pos_type   pos_type;  // position in buffer
   typedef std::iostream::off_type   off_type;  // offset in buffer
   typedef std::iostream::seekdir    whence_type;
   typedef std::iostream::openmode   openmode_type;
#else
   typedef MPI_Offset   pos_type;    // position in buffer
   typedef MPI_Offset   off_type;    // offset in buffer
   typedef int          whence_type;
   typedef int          openmode_type;
#endif

         //  Set the file pointer relative to start of data in file,
         //  end of file (use negative to go backward), or current 
         //  location.  Checksum is reset if in use.

         //  Caution: make sure to convert sizeof(...) quantities
         //  to an IOHandler::off_type(sizeof(...)) if you wish to use 
         //  negative offsets.  sizeof(...) returns an unsigned 
         //  integer type.

   void seekFromStart(off_type offset);  // start means start of data (not file)
   void seekFromCurr(off_type offset);
   void seekFromEnd(off_type offset);
    
   void seek(pos_type offset);          // from start of data
   void seekBegin(off_type offset);     // from start of data
   void seekRelative(off_type offset);
   void seekEnd(off_type offset);
   void rewind();   // puts pointer at start of data in file

         //  Getting the file pointer location in bytes from start of data

   pos_type tell();
   pos_type currentPosition();


         //  Check summing
   
   void turnOnChecksum();
   void turnOffChecksum();
   void resetChecksum();   
   QDPUtil::n_uint32_t getChecksum();


         // print out MPI-IO file hints or ID on primary node

   void printFileHints();
   void printFileID();

     // Open file, read string at a particular position, then close. Returns
     // true if the file exists and can be opened and its file type matches 
     // "filetype_id", returns false otherwise.  The string is returned in 
     // "stringvalue".  This routine is used by objects in data_io_handler.h 
     // when the multi-file handlers build up their maps of file keys.
           
   bool peekString(std::string& stringvalue, size_t position, 
                   const std::string& filename,
                   const std::string& filetype_id="");

         // write routines  

   void write(const std::string& output);

   template<typename T>
   void write(const T& output);    // basic types

   template<typename T>
   void multi_write(const T* output, int n);    // basic types

   template <typename T>
   void write(const std::vector<T>& output);   // T is basic type

   template <typename T>
   void write(const multi1d<T>& output);  // T is basic type

   template <typename T>
   void write(const multi2d<T>& output);   // basic types

   template <typename T>
   void write(const multi3d<T>& output);   // basic types

   template<typename T>
   void write(const OScalar<T>& output);

   template <typename T>
   void write(const std::vector<OScalar<T> >& output);

   template <typename T>
   void write(const multi1d<OScalar<T> >& output);

   template <typename T>
   void write(const multi2d<OScalar<T> >& output);

   template <typename T>
   void write(const multi3d<OScalar<T> >& output);

   template<typename T>
   void write(const OLattice<T>& output);

   template <typename T>
   void write(const std::vector<OLattice<T> >& output);

   template <typename T>
   void write(const multi1d<OLattice<T> >& output);

   template <typename T>
   void writeTimeSlice(const OLattice<T>& output, int time_slice);

   template <typename T>
   void write(const TimeSliceOf<T>& output);

   template <typename T>
   void write(const std::vector<TimeSliceOf<T> >& output);


          // read routines
 
   void read(std::string& input, bool broadcast=true); 

   template<typename T>
   void read(T& input);    // basic types

   template<typename T>
   void multi_read(T* input, int n);    // basic types

   template <typename T>
   void read(std::vector<T>& input);   // basic types

   template <typename T>
   void read(multi1d<T>& input);  // T is basic type

   template <typename T>
   void read(multi2d<T>& input);  // basic types

   template <typename T>
   void read(multi3d<T>& input);   // basic types

   template<typename T>
   void read(OScalar<T>& input);

   template <typename T>
   void read(std::vector<OScalar<T> >& input);

   template <typename T>
   void read(multi1d<OScalar<T> >& input);

   template <typename T>
   void read(multi2d<OScalar<T> >& input);

   template <typename T>
   void read(multi3d<OScalar<T> >& input);

   template<typename T>
   void read(OLattice<T>& input);

   template <typename T>
   void read(std::vector<OLattice<T> >& input);

   template <typename T>
   void read(multi1d<OLattice<T> >& input);

   template <typename T>
   void readTimeSlice(OLattice<T>& input, int time_slice);

   template <typename T>
   void read(TimeSliceOf<T>& input);

   template <typename T>
   void read(std::vector<TimeSliceOf<T> >& input);


         //  number of bytes routines (will be useful by IOMap for
         //  verifying read/writes and determining whether overwrites
         //  should occur)
  
   size_t numbytes(const std::string& output);

   template<typename T>
   size_t numbytes(const T& output);    // basic types

   template <typename T>
   size_t numbytes(const std::vector<T>& output);  // T is basic type

   template <typename T>
   size_t numbytes(const multi1d<T>& output);  // T is basic type

   template <typename T>
   size_t numbytes(const multi2d<T>& output);  // basic types

   template <typename T>
   size_t numbytes(const multi3d<T>& output);   // basic types

   template<typename T>
   size_t numbytes(const OScalar<T>& output);

   template <typename T>
   size_t numbytes(const std::vector<OScalar<T> >& output);

   template <typename T>
   size_t numbytes(const multi1d<OScalar<T> >& output);

   template <typename T>
   size_t numbytes(const multi2d<OScalar<T> >& output);

   template <typename T>
   size_t numbytes(const multi3d<OScalar<T> >& output);

   template<typename T>
   size_t numbytes(const OLattice<T>& output);

   template <typename T>
   size_t numbytes(const std::vector<OLattice<T> >& output);

   template <typename T>
   size_t numbytes(const multi1d<OLattice<T> >& output);

   template <typename T>
   size_t numbytes(const TimeSliceOf<T>& output);

   template <typename T>
   size_t numbytes(const std::vector<TimeSliceOf<T> >& output);



         //  Delayed read/write of multiple common data (of fixed sizes).
         //  Warning: These routines build up structures contains pointers
         //  so do not change the quantities they point to until the read/write
         //  is finished.  You can add just the basic 11 data types of fixed size.
         //  The "finishDelayed..." functions copy the data to/from a vector<char>
         //  buffer, issuing a single read/write command.

   void startDelayedWrite();
   void startDelayedRead();

   template <typename T>
   void addDelayedWrite(const T& d);

   template <typename T1, typename T2>
   void addDelayedWrite(const T1& d1, const T2& d2);

   template <typename T1, typename T2, typename T3>
   void addDelayedWrite(const T1& d1, const T2& d2, const T3& d3);

   template <typename T1, typename T2, typename T3, typename T4>
   void addDelayedWrite(const T1& d1, const T2& d2, const T3& d3, const T4& d4);

   template <typename T>
   void addDelayedRead(T& d);

   template <typename T1, typename T2>
   void addDelayedRead(T1& d1, T2& d2);

   template <typename T1, typename T2, typename T3>
   void addDelayedRead(T1& d1, T2& d2, T3& d3);

   template <typename T1, typename T2, typename T3, typename T4>
   void addDelayedRead(T1& d1, T2& d2, T3& d3, T4& d4);

   size_t finishDelayedWrite(); // returns the number of bytes written /read
   size_t finishDelayedRead();



 private:

      // private utility routines
       
   void open_existing_file(const std::string& filetype_id,
                           openmode_type access_mode);

   void open_new_file(const std::string& filetype_id,
                      char endianness,
                      int striping_factor, int striping_unit);

   void clear();

   bool fileExists();

   void writeCommon(const char *data, size_t nbytes);

   void writeDistributed(const char *data, const DistArrayViewInfo& dinfo,
                         bool datashift=false);

   void readCommon(char *data, size_t nbytes, bool broadcast=true);

   void readDistributed(char *data, const DistArrayViewInfo& dinfo,
                        bool datashift=false);

   void check_for_failure(int errcode, const std::string& mesg="");

   std::string int_to_string(int intval);

   std::string tidyString(const std::string& str);  

   void readIDstring(char &endianness, std::string& ID_string);

   void writeIDstring(const std::string& ID_string);

   void write_common(const char* output, size_t element_size, size_t nelements);

   void write_lattice(const char* output, size_t bytes_per_word, 
                      const DistArrayViewInfo& dinfo);

   template<typename T>
   void latt_write_prim(const T& output);


   void read_common(char* output, size_t element_size, size_t nelements);
  
   void read_lattice(char* input, size_t bytes_per_word,
                     const DistArrayViewInfo& dinfo);


   void lexify(char *out, const char* in, const DistArrayViewInfo& dinfo);
 
   void unlexify(char *out, const char* in, const DistArrayViewInfo& dinfo);

   void compute_lattice_checksum(char* data, size_t bytes_per_site, 
                                 size_t lexicostart, size_t nsites);


   struct outbasicdatainfo
    {
     const void* dataptr;
     size_t elsize;
     outbasicdatainfo(const char& d) : dataptr(&d), elsize(sizeof(char)) {}
     outbasicdatainfo(const int& d) : dataptr(&d), elsize(sizeof(int)) {}
     outbasicdatainfo(const unsigned int& d) : dataptr(&d), elsize(sizeof(unsigned int)) {}
     outbasicdatainfo(const short int& d) : dataptr(&d), elsize(sizeof(short int)) {}
     outbasicdatainfo(const unsigned short int& d) : dataptr(&d), elsize(sizeof(unsigned short int)) {}
     outbasicdatainfo(const long int& d) : dataptr(&d), elsize(sizeof(long int)) {}
     outbasicdatainfo(const unsigned long int& d) : dataptr(&d), elsize(sizeof(unsigned long int)) {}
     outbasicdatainfo(const long long int& d) : dataptr(&d), elsize(sizeof(long long int)) {}
     outbasicdatainfo(const float& d) : dataptr(&d), elsize(sizeof(float)) {}
     outbasicdatainfo(const double& d) : dataptr(&d), elsize(sizeof(double)) {}
     outbasicdatainfo(const bool& d) : dataptr(&d), elsize(sizeof(bool)) {}
    };

   struct inbasicdatainfo
    {
     void* dataptr;
     size_t elsize;
     inbasicdatainfo(char& d) : dataptr(&d), elsize(sizeof(char)) {}
     inbasicdatainfo(int& d) : dataptr(&d), elsize(sizeof(int)) {}
     inbasicdatainfo(unsigned int& d) : dataptr(&d), elsize(sizeof(unsigned int)) {}
     inbasicdatainfo(short int& d) : dataptr(&d), elsize(sizeof(short int)) {}
     inbasicdatainfo(unsigned short int& d) : dataptr(&d), elsize(sizeof(unsigned short int)) {}
     inbasicdatainfo(long int& d) : dataptr(&d), elsize(sizeof(long int)) {}
     inbasicdatainfo(unsigned long int& d) : dataptr(&d), elsize(sizeof(unsigned long int)) {}
     inbasicdatainfo(long long int& d) : dataptr(&d), elsize(sizeof(long long int)) {}
     inbasicdatainfo(float& d) : dataptr(&d), elsize(sizeof(float)) {}
     inbasicdatainfo(double& d) : dataptr(&d), elsize(sizeof(double)) {}
     inbasicdatainfo(bool& d) : dataptr(&d), elsize(sizeof(bool)) {}
    };


   std::list<outbasicdatainfo> outdelay;
   std::list<inbasicdatainfo> indelay;

   size_t multi_write_common(const std::list<outbasicdatainfo>& dinfo);

   size_t multi_read_common(const std::list<inbasicdatainfo>& dinfo);

   void delete_file();
   void file_open(openmode_type access_mode, const std::string& errmsg);
   void file_close();
   void file_seek(off_type offset, whence_type whence);
   pos_type file_tell();
   void set_striping(int striping_factor, int striping_unit);

   bool peeker(std::string& stringvalue, size_t position, 
               const std::string& filename,
               const std::string& filetype_id="");

       // constants
       
   static const int ID_string_length;
   static const pos_type data_start_pos;

   static const int IO_ERR_NO_SUCH_FILE;
   static const int IO_ERR_ACCESS;
   static const int IO_ERR_OTHER;
   static const int IO_SUCCESS;

   static const openmode_type IO_MODE_RDONLY;
   static const openmode_type IO_MODE_RDWR;
   static const openmode_type IO_MODE_CREATE;
   static const openmode_type IO_MODE_EXCL;

   static const whence_type IO_SEEK_BEG;
   static const whence_type IO_SEEK_CUR;
   static const whence_type IO_SEEK_END;

         // for static (compile time) assertion
//   template <bool b>
//   void static_assert()
//   { typedef char asserter[b?1:-1]; }

       // purpose of this class is to cause compiler errors if
       // constructed with a non basic type 

    class BasicType
     {
      public:
        explicit BasicType( char s ) {}
        explicit BasicType( int s ) {}
        explicit BasicType( unsigned int s ) {}
        explicit BasicType( short int s ) {}
        explicit BasicType( unsigned short int s ) {}
        explicit BasicType( long int s ) {}
        explicit BasicType( unsigned long int s ) {}
        explicit BasicType( long long int s ) {}
        explicit BasicType( unsigned long long int s ) {}
        explicit BasicType( float s ) {}
        explicit BasicType( double s ) {}
        explicit BasicType( bool s ) {}
     };

       // purpose of this class is to cause compiler errors if
       // constructed with a non-OScalar type (will get no access to
       // private constructor errors)
       
    class OScalarType 
     {
        template<typename T> OScalarType(const T& in) {}
        template<typename T> OScalarType(const OLattice<T>& in) {}
      public:
        template<typename T> OScalarType(const OScalar<T>& in){}
     };

       // purpose of this class is to cause compiler errors if
       // constructed with a non-lattice type (will get no access to
       // private constructor errors)
       
    class OLatticeType 
     {
        template<typename T> OLatticeType(const T& in) {}
        template<typename T> OLatticeType(const OScalar<T>& in) {}
      public:
        template<typename T> OLatticeType(const OLattice<T>& in){}
     };


};


 // ***************************************************************


template <typename T>
void IOHandler::write(const T& output)
{ 
 IOHandler::BasicType dummy(output);  // cause compiler error if T is not a valid basic type
 write_common((const char*)&output, sizeof(T), 1);
}


template <typename T>
void IOHandler::multi_write(const T* output, int n)
{ 
 IOHandler::BasicType dummy(*output);  // cause compiler error if T is not a valid basic type
 write_common((const char*)output, sizeof(T), n);
}


template <typename T>
void IOHandler::write(const std::vector<T>& output)
{
 IOHandler::BasicType dummy(output[0]);  // cause compiler error if T is not a valid basic type
 int n=output.size();
 write(n);
 if (n==0) return;
 write_common((const char*)&output[0],sizeof(T),size_t(n));
}


template <typename T>
void IOHandler::write(const multi1d<T>& output)
{
 IOHandler::BasicType dummy(output[0]);  // cause compiler error if T is not a valid basic type
 int n=output.size();
 write(n);
 if (n==0) return;
 write_common((const char*)&output[0],sizeof(T),size_t(n));
}


template <typename T>
void IOHandler::write(const multi2d<T>& output)
{
 IOHandler::BasicType dummy(output[0][0]);  // cause compiler error if T is not a valid basic type
 unsigned int sz[2];
 sz[0]=output.size1();
 sz[1]=output.size2();
 multi_write(sz,2);
 int n=sz[0]*sz[1];
 if (n==0) return;
 write_common((const char*)&output[0][0],sizeof(T),size_t(n));
}


template <typename T>
void IOHandler::write(const multi3d<T>& output)
{
 IOHandler::BasicType dummy(output[0][0][0]);  // cause compiler error if T is not a valid basic type
 unsigned int sz[3];
 sz[0]=output.size1();
 sz[1]=output.size2();
 sz[2]=output.size3();
 multi_write(sz,3);
 int n=sz[0]*sz[1]*sz[2];
 if (n==0) return;
 write_common((const char*)&output[0][0][0],sizeof(T),size_t(n));
}


template <typename T>
void IOHandler::write(const OScalar<T>& output)
{ 
 size_t wsize=sizeof(typename WordType<T>::Type_t);
 write_common((const char *)&(output.elem(0)), wsize, sizeof(T)/wsize);
}


template <typename T>
void IOHandler::write(const std::vector<OScalar<T> >& output)
{
 int n=output.size();
 write(n);
 if (n==0) return;
 size_t wsize=sizeof(typename WordType<T>::Type_t);
 write_common((const char *)&(output[0].elem(0)), wsize, size_t(n)*sizeof(T)/wsize);
}


template <typename T>
void IOHandler::write(const multi1d<OScalar<T> >& output)
{
 int n=output.size();
 write(n);
 if (n==0) return;
 size_t wsize=sizeof(typename WordType<T>::Type_t);
 write_common((const char *)&(output[0].elem(0)), wsize, size_t(n)*sizeof(T)/wsize);
}


template <typename T>
void IOHandler::write(const multi2d<OScalar<T> >& output)
{
 unsigned int sz[2];
 sz[0]=output.size1();
 sz[1]=output.size2();
 multi_write(sz,2);
 int n=sz[0]*sz[1];
 if (n==0) return;
 size_t wsize=sizeof(typename WordType<T>::Type_t);
 write_common((const char *)&(output[0][0].elem(0)), wsize, size_t(n)*sizeof(T)/wsize);
}


template <typename T>
void IOHandler::write(const multi3d<OScalar<T> >& output)
{
 unsigned int sz[3];
 sz[0]=output.size1();
 sz[1]=output.size2();
 sz[2]=output.size3();
 multi_write(sz,3);
 int n=sz[0]*sz[1]*sz[2];
 if (n==0) return;
 size_t wsize=sizeof(typename WordType<T>::Type_t);
 write_common((const char *)&(output[0][0][0].elem(0)), wsize, size_t(n)*sizeof(T)/wsize);
}


template <typename T>
void IOHandler::write(const OLattice<T>& output)
{
 int sizeT=sizeof(T);
 lattinfo.resetBytes(sizeT);
 unsigned int sz[QDP::Nd+2];
 sz[0]=QDP::Nd;
 for (int k=0;k<QDP::Nd;k++)
    sz[k+1]=Layout::lattSize()[k];
 sz[QDP::Nd+1]=sizeT;
 multi_write(sz,QDP::Nd+2);
 write_lattice((const char *)&(output.elem(0)),
               sizeof(typename WordType<T>::Type_t), lattinfo);
}


template <typename T>
void IOHandler::write(const std::vector<OLattice<T> >& output)
{
 int n=output.size();
 write(n);
 if (n==0) return;
 for (int k=0;k<n;k++) write(output[k]);
}


template <typename T>
void IOHandler::write(const multi1d<OLattice<T> >& output)
{
 int n=output.size();
 write(n);
 if (n==0) return;
 for (int k=0;k<n;k++) write(output[k]);
}


template <typename T>
void IOHandler::writeTimeSlice(const OLattice<T>& output, int time_slice)
{
 int sizeT=sizeof(T);
 DistArrayViewInfo dinfo(time_slice,sizeT);
 unsigned int sz[QDP::Nd+1];
 sz[0]=QDP::Nd-1;
 for (int k=0;k<QDP::Nd-1;k++)
    sz[k+1]=Layout::lattSize()[k];
 sz[QDP::Nd]=sizeT;
 multi_write(sz,QDP::Nd+1);
 write_lattice((const char *)&(output.elem(0)),
               sizeof(typename WordType<T>::Type_t), dinfo);
}


template <typename T>
void IOHandler::write(const TimeSliceOf<T>& output)
{
 writeTimeSlice(output.getData(), output.getCurrentTime());
}


template <typename T>
void IOHandler::write(const std::vector<TimeSliceOf<T> >& output)
{
 int n=output.size();
 write(n);
 if (n==0) return;
 for (int k=0;k<n;k++)
    writeTimeSlice(output[k].getData(), output[k].getCurrentTime());
}


  // *************************************************

  
template<typename T>
void write(IOHandler& io, const T& output)   // basic types
 { io.write(output); }

template <typename T>
void write(IOHandler& io, const std::vector<T>& output)
 { io.write(output); }

template <typename T>
void write(IOHandler& io, const multi1d<T>& output)  // T is basic type
 { io.write(output); }

template <typename T>
void write(IOHandler& io, const multi2d<T>& output)
 { io.write(output); }

template <typename T>
void write(IOHandler& io, const multi3d<T>& output)
 { io.write(output); }

template<typename T>
void write(IOHandler& io, const OScalar<T>& output)
 { io.write(output); }

template <typename T>
void write(IOHandler& io, const std::vector<OScalar<T> >& output)
 { io.write(output); }

template <typename T>
void write(IOHandler& io, const multi1d<OScalar<T> >& output)
 { io.write(output); }

template <typename T>
void write(IOHandler& io, const multi2d<OScalar<T> >& output)
 { io.write(output); }

template <typename T>
void write(IOHandler& io, const multi3d<OScalar<T> >& output)
 { io.write(output); }

template<typename T>
void write(IOHandler& io, const OLattice<T>& output)
 { io.write(output); }

template <typename T>
void write(IOHandler& io, const std::vector<OLattice<T> >& output)
 { io.write(output); }

template <typename T>
void write(IOHandler& io, const multi1d<OLattice<T> >& output)
 { io.write(output); }

template <typename T>
void write(IOHandler& io, const TimeSliceOf<T>& output)
 { io.write(output); }
  
template <typename T>
void write(IOHandler& io, const std::vector<TimeSliceOf<T> >& output)
 { io.write(output); }


// ************************************************************  
  

template <typename T>
void IOHandler::read(T& input)
{ 
 IOHandler::BasicType dummy(input);  // cause compiler error if T is not a valid basic type
 read_common((char*)&input, sizeof(T), 1);
}


template <typename T>
void IOHandler::multi_read(T* input, int n)
{ 
 IOHandler::BasicType dummy(*input);  // cause compiler error if T is not a valid basic type
 read_common((char*)input, sizeof(T), n);
}


template <typename T>
void IOHandler::read(std::vector<T>& input)
{
 IOHandler::BasicType dummy(input[0]);  // cause compiler error if T is not a valid basic type
 int n;
 read(n);
 if (n==0) return;
    // if reading in wrong location, could get nonsense here,
    // so place a reasonable limit
 check_for_failure((n<0)||(n>16777216), "multi1d too large...bad file location?");
 input.resize(n);
 read_common((char*)&input[0],sizeof(T),size_t(n));
}


template <typename T>
void IOHandler::read(multi1d<T>& input)
{
 IOHandler::BasicType dummy(input[0]);  // cause compiler error if T is not a valid basic type
 int n;
 read(n);
 if (n==0) return;
    // if reading in wrong location, could get nonsense here,
    // so place a reasonable limit
 check_for_failure((n<0)||(n>16777216), "multi1d too large...bad file location?");
 input.resize(n);
 read_common((char*)&input[0],sizeof(T),size_t(n));
}


template <typename T>
void IOHandler::read(multi2d<T>& input)
{
 IOHandler::BasicType dummy(input[0][0]);  // cause compiler error if T is not a valid basic type
 unsigned int sz[2];
 multi_read(sz,2);
 unsigned int n=sz[0]*sz[1];
 if (n==0) return;
    // if reading in wrong location, could get nonsense here,
    // so place a reasonable limit
 check_for_failure(n>67108864, "multi2d too large...bad file location?");
 input.resize(sz[1],sz[0]);
 read_common((char*)&input[0][0],sizeof(T),size_t(n));
}


template <typename T>
void IOHandler::read(multi3d<T>& input)
{
 IOHandler::BasicType dummy(input[0][0][0]);  // cause compiler error if T is not a valid basic type
 unsigned int sz[3];
 multi_read(sz,3);
 unsigned int n=sz[0]*sz[1]*sz[2];
 if (n==0) return;
    // if reading in wrong location, could get nonsense here,
    // so place a reasonable limit
 check_for_failure(n>67108864, "multi3d too large...bad file location?");
 input.resize(sz[2],sz[1],sz[0]);
 read_common((char*)&input[0][0][0],sizeof(T),size_t(n));
}


template <typename T>
void IOHandler::read(OScalar<T>& input)
{ 
 size_t wsize=sizeof(typename WordType<T>::Type_t);
 read_common((char *)&(input.elem(0)), wsize, sizeof(T)/wsize);
}


template <typename T>
void IOHandler::read(std::vector<OScalar<T> >& input)
{
 int n;
 read(n);
 if (n==0) return;
    // if reading in wrong location, could get nonsense here,
    // so place a reasonable limit
 check_for_failure((n<0)||(n>16777216), "multi1d too large...bad file location?");
 input.resize(n);
 size_t wsize=sizeof(typename WordType<T>::Type_t);
 read_common((char *)&(input[0].elem(0)), wsize, size_t(n)*sizeof(T)/wsize);
}


template <typename T>
void IOHandler::read(multi1d<OScalar<T> >& input)
{
 int n;
 read(n);
 if (n==0) return;
    // if reading in wrong location, could get nonsense here,
    // so place a reasonable limit
 check_for_failure((n<0)||(n>16777216), "multi1d too large...bad file location?");
 input.resize(n);
 size_t wsize=sizeof(typename WordType<T>::Type_t);
 read_common((char *)&(input[0].elem(0)), wsize, size_t(n)*sizeof(T)/wsize);
}


template <typename T>
void IOHandler::read(multi2d<OScalar<T> >& input)
{
 unsigned int sz[2];
 multi_read(sz,2);
 unsigned int n=sz[0]*sz[1];
 if (n==0) return;
    // if reading in wrong location, could get nonsense here,
    // so place a reasonable limit
 check_for_failure(n>67108864, "multi2d too large...bad file location?");
 input.resize(sz[1],sz[0]);
 size_t wsize=sizeof(typename WordType<T>::Type_t);
 read_common((char *)&(input[0][0].elem(0)), wsize, size_t(n)*sizeof(T)/wsize);
}


template <typename T>
void IOHandler::read(multi3d<OScalar<T> >& input)
{
 unsigned int sz[3];
 multi_read(sz,3);
 unsigned int n=sz[0]*sz[1]*sz[2];
 if (n==0) return;
    // if reading in wrong location, could get nonsense here,
    // so place a reasonable limit
 check_for_failure(n>67108864, "multi3d too large...bad file location?");
 input.resize(sz[2],sz[1],sz[0]);
 size_t wsize=sizeof(typename WordType<T>::Type_t);
 read_common((char *)&(input[0][0][0].elem(0)), wsize, size_t(n)*sizeof(T)/wsize);
}


template <typename T>
void IOHandler::read(OLattice<T>& input)
{
 unsigned int sz[QDP::Nd+2];
 multi_read(sz,QDP::Nd+2);
 bool errflag=((sz[0]!=(unsigned int)(QDP::Nd))||(sz[QDP::Nd+1]!=sizeof(T)));
 for (int k=0;k<QDP::Nd;k++)
    errflag=errflag || (sz[k+1]!=(unsigned int)(Layout::lattSize()[k]));
 check_for_failure(errflag, "Invalid lattice read");
 lattinfo.resetBytes(sizeof(T));
 read_lattice((char *)&(input.elem(0)),
               sizeof(typename WordType<T>::Type_t), lattinfo);
}


template <typename T>
void IOHandler::read(std::vector<OLattice<T> >& input)
{
 int n;
 read(n);
 if (n==0) return;
    // if reading in wrong location, could get nonsense here,
    // so place a reasonable limit
 check_for_failure((n<0)||(n>1024), "multi1d too large...bad file location?");
 input.resize(n);
 for (int k=0;k<n;k++) read(input[k]);
}


template <typename T>
void IOHandler::read(multi1d<OLattice<T> >& input)
{
 int n;
 read(n);
 if (n==0) return;
    // if reading in wrong location, could get nonsense here,
    // so place a reasonable limit
 check_for_failure((n<0)||(n>1024), "multi1d too large...bad file location?");
 input.resize(n);
 for (int k=0;k<n;k++) read(input[k]);
}

  
template <typename T>
void IOHandler::readTimeSlice(OLattice<T>& input, int time_slice)
{
 DistArrayViewInfo dinfo(time_slice,sizeof(T));
 unsigned int sz[QDP::Nd+1];
 multi_read(sz,QDP::Nd+1);
 bool errflag=((sz[0]!=(QDP::Nd-1))||(sz[QDP::Nd]!=sizeof(T)));
 for (int k=0;k<QDP::Nd-1;k++)
    errflag=errflag || (int(sz[k+1])!=Layout::lattSize()[k]);
 check_for_failure(errflag, "Invalid lattice read");
 read_lattice((char *)&(input.elem(0)),
               sizeof(typename WordType<T>::Type_t), dinfo);
}


template <typename T>
void IOHandler::read(TimeSliceOf<T>& input)
{
 readTimeSlice(input.getData(), input.getCurrentTime());
}


template <typename T>
void IOHandler::read(std::vector<TimeSliceOf<T> >& input)
{
 int n;
 read(n);
 check_for_failure( (n!= input.size()), "bad vector<TimeSliceOf<> > read");
 for (int k=0;k<n;k++)
    readTimeSlice(input[k].getData(), input[k].getCurrentTime());
}


// **********************************************

template<typename T>
void read(IOHandler& io, T& input)   // basic types
 { io.read(input); }

template <typename T>
void read(IOHandler& io, std::vector<T>& input)
 { io.read(input); }

template <typename T>
void read(IOHandler& io, multi1d<T>& input)  // T is basic type
 { io.read(input); }

template <typename T>
void read(IOHandler& io, multi2d<T>& input)
 { io.read(input); }

template <typename T>
void read(IOHandler& io, multi3d<T>& input)
 { io.read(input); }

template<typename T>
void read(IOHandler& io, OScalar<T>& input)
 { io.read(input); }

template <typename T>
void read(IOHandler& io, std::vector<OScalar<T> >& input)
 { io.read(input); }

template <typename T>
void read(IOHandler& io, multi1d<OScalar<T> >& input)
 { io.read(input); }

template <typename T>
void read(IOHandler& io, multi2d<OScalar<T> >& input)
 { io.read(input); }

template <typename T>
void read(IOHandler& io, multi3d<OScalar<T> >& input)
 { io.read(input); }

template<typename T>
void read(IOHandler& io, OLattice<T>& input)
 { io.read(input); }

template <typename T>
void read(IOHandler& io, std::vector<OLattice<T> >& input)
 { io.read(input); }

template <typename T>
void read(IOHandler& io, multi1d<OLattice<T> >& input)
 { io.read(input); }

template <typename T>
void read(IOHandler& io, TimeSliceOf<T>& input)
 { io.read(input); }
  
template <typename T>
void read(IOHandler& io, std::vector<TimeSliceOf<T> >& input)
 { io.read(input); }


// ************************************************************  


template <typename T>
size_t IOHandler::numbytes(const T& data)
{ 
 IOHandler::BasicType dummy(data);  // cause compiler error if T is not a valid basic type
 return sizeof(T);
}


template <typename T>
size_t IOHandler::numbytes(const std::vector<T>& data)
{
 IOHandler::BasicType dummy(data[0]);  // cause compiler error if T is not a valid basic type
 return sizeof(int)+data.size()*sizeof(T);
}


template <typename T>
size_t IOHandler::numbytes(const multi1d<T>& data)
{
 IOHandler::BasicType dummy(data[0]);  // cause compiler error if T is not a valid basic type
 return sizeof(int)+data.size()*sizeof(T);
}


template <typename T>
size_t IOHandler::numbytes(const multi2d<T>& data)
{
 IOHandler::BasicType dummy(data[0][0]);  // cause compiler error if T is not a valid basic type
 return 2*sizeof(int)+data.size1()*data.size2()*sizeof(T);
}


template <typename T>
size_t IOHandler::numbytes(const multi3d<T>& data)
{
 IOHandler::BasicType dummy(data[0][0][0]);  // cause compiler error if T is not a valid basic type
 return 3*sizeof(int)+data.size1()*data.size2()*data.size3()*sizeof(T);
}


template <typename T>
size_t IOHandler::numbytes(const OScalar<T>& data)
{ 
 return sizeof(T);
}


template <typename T>
size_t IOHandler::numbytes(const std::vector<OScalar<T> >& data)
{
 return sizeof(int)+data.size()*sizeof(T);
}


template <typename T>
size_t IOHandler::numbytes(const multi1d<OScalar<T> >& data)
{
 return sizeof(int)+data.size()*sizeof(T);
}


template <typename T>
size_t IOHandler::numbytes(const multi2d<OScalar<T> >& data)
{
 return 2*sizeof(int)+data.size1()*data.size2()*sizeof(T);
}


template <typename T>
size_t IOHandler::numbytes(const multi3d<OScalar<T> >& data)
{
 return 3*sizeof(int)+data.size1()*data.size2()*data.size3()*sizeof(T);
}


template <typename T>
size_t IOHandler::numbytes(const OLattice<T>& data)
{
 return (QDP::Nd+2)*sizeof(int)+Layout::vol()*sizeof(T);
}


template <typename T>
size_t IOHandler::numbytes(const std::vector<OLattice<T> >& data)
{
 return sizeof(int)+data.size()*((QDP::Nd+2)*sizeof(int)+Layout::vol()*sizeof(T));
}


template <typename T>
size_t IOHandler::numbytes(const multi1d<OLattice<T> >& data)
{
 return sizeof(int)+data.size()*((QDP::Nd+2)*sizeof(int)+Layout::vol()*sizeof(T));
}


template <typename T>
size_t IOHandler::numbytes(const TimeSliceOf<T>& data)
{
 IOHandler::off_type fullsize=numbytes(data.getData())-(QDP::Nd+2)*sizeof(int);
 fullsize=fullsize/Layout::lattSize()[QDP::Nd-1];
 return (QDP::Nd+1)*sizeof(int)+fullsize;
}


template <typename T>
size_t IOHandler::numbytes(const std::vector<TimeSliceOf<T> >& data)
{
 size_t res=sizeof(int);
 for (uint k=0;k<data.size();k++) res+=numbytes(data[k]);
 return res;
}

// ************************************************************  


template<typename T>
size_t numbytes(IOHandler& io, const T& data)   // generic routine
 { return data.numbytes(); }

template<> inline size_t numbytes(IOHandler& io, const char& data) 
{ return io.numbytes(data); }

template<> inline size_t numbytes(IOHandler& io, const int& data)  
{ return io.numbytes(data); }

template<> inline size_t numbytes(IOHandler& io, const unsigned int& data) 
{ return io.numbytes(data); }

template<> inline size_t numbytes(IOHandler& io, const short int& data)  
{ return io.numbytes(data); }

template<> inline size_t numbytes(IOHandler& io, const unsigned short int& data) 
{ return io.numbytes(data); }

template<> inline size_t numbytes(IOHandler& io, const long int& data)  
{ return io.numbytes(data); }

template<> inline size_t numbytes(IOHandler& io, const unsigned long int& data) 
{ return io.numbytes(data); }

template<> inline size_t numbytes(IOHandler& io, const long long int& data)  
{ return io.numbytes(data); }

template<> inline size_t numbytes(IOHandler& io, const unsigned long long int& data) 
{ return io.numbytes(data); }

template<> inline size_t numbytes(IOHandler& io, const float& data)  
{ return io.numbytes(data); }

template<> inline size_t numbytes(IOHandler& io, const double& data) 
{ return io.numbytes(data); }

template<> inline size_t numbytes(IOHandler& io, const bool& data)  
{ return io.numbytes(data); }


template <typename T>
size_t numbytes(IOHandler& io, const std::vector<T>& data)
 { return io.numbytes(data); }

template <typename T>
size_t numbytes(IOHandler& io, const multi1d<T>& data)  // T is basic type
 { return io.numbytes(data); }

template <typename T>
size_t numbytes(IOHandler& io, const multi2d<T>& data)
 { return io.numbytes(data); }

template <typename T>
size_t numbytes(IOHandler& io, const multi3d<T>& data)
 { return io.numbytes(data); }

template<typename T>
size_t numbytes(IOHandler& io, const OScalar<T>& data)
 { return io.numbytes(data); }

template <typename T>
size_t numbytes(IOHandler& io, const std::vector<OScalar<T> >& data)
 { return io.numbytes(data); }

template <typename T>
size_t numbytes(IOHandler& io, const multi1d<OScalar<T> >& data)
 { return io.numbytes(data); }

template <typename T>
size_t numbytes(IOHandler& io, const multi2d<OScalar<T> >& data)
 { return io.numbytes(data); }

template <typename T>
size_t numbytes(IOHandler& io, const multi3d<OScalar<T> >& data)
 { return io.numbytes(data); }

template<typename T>
size_t numbytes(IOHandler& io, const OLattice<T>& data)
 { return io.numbytes(data); }

template <typename T>
size_t numbytes(IOHandler& io, const std::vector<OLattice<T> >& data)
 { return io.numbytes(data); }

template <typename T>
size_t numbytes(IOHandler& io, const multi1d<OLattice<T> >& data)
 { return io.numbytes(data); }

template <typename T>
size_t numbytes(IOHandler& io, const TimeSliceOf<T>& data)
 { return io.numbytes(data); }
  
template <typename T>
size_t numbytes(IOHandler& io, const std::vector<TimeSliceOf<T> >& data)
 { return io.numbytes(data); }


  // ************************************************

template <typename T>
void IOHandler::addDelayedWrite(const T& d)
{
 outdelay.push_back(outbasicdatainfo(d));
}

template <typename T1, typename T2>
void IOHandler::addDelayedWrite(const T1& d1, const T2& d2)
{
 outdelay.push_back(outbasicdatainfo(d1));
 outdelay.push_back(outbasicdatainfo(d2));
}

template <typename T1, typename T2, typename T3>
void IOHandler::addDelayedWrite(const T1& d1, const T2& d2, const T3& d3)
{
 outdelay.push_back(outbasicdatainfo(d1));
 outdelay.push_back(outbasicdatainfo(d2));
 outdelay.push_back(outbasicdatainfo(d3));
}

template <typename T1, typename T2, typename T3, typename T4>
void IOHandler::addDelayedWrite(const T1& d1, const T2& d2, const T3& d3, const T4& d4)
{
 outdelay.push_back(outbasicdatainfo(d1));
 outdelay.push_back(outbasicdatainfo(d2));
 outdelay.push_back(outbasicdatainfo(d3));
 outdelay.push_back(outbasicdatainfo(d4));
}

template <typename T>
void IOHandler::addDelayedRead(T& d)
{
 indelay.push_back(inbasicdatainfo(d));
}

template <typename T1, typename T2>
void IOHandler::addDelayedRead(T1& d1, T2& d2)
{
 indelay.push_back(inbasicdatainfo(d1));
 indelay.push_back(inbasicdatainfo(d2));
}

template <typename T1, typename T2, typename T3>
void IOHandler::addDelayedRead(T1& d1, T2& d2, T3& d3)
{
 indelay.push_back(inbasicdatainfo(d1));
 indelay.push_back(inbasicdatainfo(d2));
 indelay.push_back(inbasicdatainfo(d3));
}

template <typename T1, typename T2, typename T3, typename T4>
void IOHandler::addDelayedRead(T1& d1, T2& d2, T3& d3, T4& d4)
{
 indelay.push_back(inbasicdatainfo(d1));
 indelay.push_back(inbasicdatainfo(d2));
 indelay.push_back(inbasicdatainfo(d3));
 indelay.push_back(inbasicdatainfo(d4));
}


// **************************************************************
  }
}
#endif
