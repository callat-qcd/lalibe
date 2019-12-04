//Arjun Singh Gambhir
//Real authors: Ben Hoerz, John Bulava, Colin Morningstar, ??
//Code directly ported from chroma_laph, more acknowledgements coming later.
#include "io_handler.h"
#include <vector>

//#define ROMIO_DS_DISABLE
#define USE_COLLECTIVE_IO

using namespace std;

namespace Chroma {
  namespace LaphEnv {
   

const int IOHandler::ID_string_length = 32;
const IOHandler::pos_type IOHandler::data_start_pos = 1+IOHandler::ID_string_length;



   // A given site on the 4-d lattice is specified by integers (x,y,z,t)
   // where each of the four integers varies from 
   //
   //     x = 0 .. Nx-1   where   Nx = Layout::lattSize()[0],
   //     y = 0 .. Ny-1   where   Ny = Layout::lattSize()[1],
   //     z = 0 .. Nz-1   where   Nz = Layout::lattSize()[2],
   //     t = 0 .. Nt-1   where   Nt = Layout::lattSize()[3]
   //
   // "lexicographic" order is defined such that "x" varies most quickly
   // and "t" varies most slowly.  The lattice is written to disk in
   // lexicographic order.  The lexicographic index for a site is
   // given by
   //           lex_index(x,y,z,t)  =  x + Nx*(y + Ny*(z + Nz*t) )
   //
   // On a parallel machine, the lattice is partitioned up into identical 
   // sublattices.  The "geom" parameter specifies how the lattice is 
   // split up.
   //              -geom Px Py Pz Pt
   //
   // There are Px*Py*Pz*Pt sublattices.  Each sublattice has extents
   //
   //        Nx/Px by Ny/Py by Nz/Pz by Nt/Pt 
   //
   // where each ratio **must** be an integer.  Let "xinc" = Nx/Px.
   // Start at lex index = 0.  Incrementing by "xinc" always moves you to
   // a different sublattice.  Also, the sites corresponding to the
   // next "xinc" lexico indices must all be on the **same** sublattice.
   //
   // In local memory, there are different ways of storing the 
   // sublattice:
   //    (a)  QDP_USE_LEXICO_LAYOUT  ->  lexicographic local lay out
   //    (b)  QDP_USE_CB2_LAYOUT -> red/black checkerboard, each checkerboard
   //                                is lexicographic
   //    (c)  QDP_USE_CB3D_LAYOUT  (not currently supported here)
   //    (d)  QDP_USE_CB32_LAYOUT  (not currently supported here)
   //
   // There are basically two ways to deal with the checker boarding.
   // First, we could write out each checkerboard, using non-trivial
   // MPI etypes to skip the parts in the file.  Secondly, we can
   // use a temporary local array to move the checkerboard into 
   // local lexicographic order, the issue one MPI file read/write.
   // We use the second approach here, assuming that memory is "cheap"
   // and minimizing the number of read/writes, and writing/reading as
   // much contiguous data as possible.

// *************************************************************************

            //  serial code

#ifdef ARCH_SCALAR

void IOHandler::delete_file()
{
 check_for_failure(remove(m_filename.c_str()),"Failure deleting file");
}

void IOHandler::file_open(IOHandler::openmode_type access_mode, 
                          const std::string& errmsg)
{
 fh.open(m_filename.c_str(), std::ios::binary | access_mode);
 check_for_failure(!fh, errmsg);
}
   
void IOHandler::file_close()
{
 fh.close();
 check_for_failure(!fh, "Failure during close");
}

const int IOHandler::IO_ERR_NO_SUCH_FILE= 1;
const int IOHandler::IO_ERR_ACCESS=       1;
const int IOHandler::IO_ERR_OTHER=        1;
const int IOHandler::IO_SUCCESS=          0;

const IOHandler::openmode_type IOHandler::IO_MODE_RDONLY= std::ios::in;
const IOHandler::openmode_type IOHandler::IO_MODE_RDWR=   std::ios::in | std::ios::out;
const IOHandler::openmode_type IOHandler::IO_MODE_CREATE= std::ios::trunc;
const IOHandler::openmode_type IOHandler::IO_MODE_EXCL=   std::ios::trunc;

const IOHandler::whence_type IOHandler::IO_SEEK_BEG  = std::ios::beg;
const IOHandler::whence_type IOHandler::IO_SEEK_CUR  = std::ios::cur;
const IOHandler::whence_type IOHandler::IO_SEEK_END  = std::ios::end;


void IOHandler::check_for_failure(int errcode, const std::string& mesg)
{
 if (errcode==IOHandler::IO_SUCCESS) return;
 QDPIO::cout << "IOHandler error with file "<<m_filename<<":"
             <<std::endl<<"  "<<mesg<<std::endl;
 QDP_abort(errcode);
}

void IOHandler::set_striping(int striping_factor, int striping_unit)
{
}


void IOHandler::file_seek(off_type offset, whence_type whence)
{
 fh.seekg(offset, whence);
 check_for_failure(fh.fail(),"Failure during seek");
}

IOHandler::pos_type IOHandler::file_tell()
{
 pos_type current=fh.tellg();
 if (current<0) check_for_failure(1,"Failure during tell");
 return current;
}


void IOHandler::printFileHints()
{
 QDPIO::cout << "Serial mode: no IOHandler file hints"<<std::endl;
}


void IOHandler::writeCommon(const char *data, size_t nbytes)
{
 pos_type currdisp=fh.tellg();
 fh.seekp(currdisp);  // just to be sure
 fh.write(data,nbytes);
 check_for_failure(fh.fail(),"Failure during common write");
 fh.seekg(currdisp+pos_type(nbytes));
}

 
void IOHandler::writeDistributed(const char *data, const DistArrayViewInfo& dinfo,
                                 bool datashift)
{
 const char* dataptr=data;
 if (datashift) dataptr+=dinfo.getElementOffsetThisNode()*dinfo.getBytesPerElement();
 writeCommon(dataptr,dinfo.getViewTotalBytes());
}


void IOHandler::readCommon(char *data, size_t nbytes, bool broadcast)
{
 pos_type currdisp=fh.tellg();
 fh.read(data,nbytes);
 check_for_failure(fh.fail(),"Failure during common read");
 fh.seekg(currdisp+pos_type(nbytes)); 
}

 
void IOHandler::readDistributed(char *data, const DistArrayViewInfo& dinfo,
                                bool datashift)
{
 char* dataptr=data;
 if (datashift) dataptr+=dinfo.getElementOffsetThisNode()*dinfo.getBytesPerElement();
 readCommon(dataptr,dinfo.getViewTotalBytes());
}

   //  check sum in lexico order in memory, not yet checker boarded
   
void IOHandler::compute_lattice_checksum(char* data, size_t bytes_per_site,
                                         size_t lexicostart, size_t nsites)
{
 checksum = QDPUtil::crc32(checksum, data, bytes_per_site*nsites);
}

// *************************************************************************

              //  parallel MPI-IO  code
              
#else


void IOHandler::delete_file()
{
 int errcode;
 if ((Layout::primaryNode())||(!global_mode)){
    errcode=MPI_File_delete((char*)m_filename.c_str(),MPI_INFO_NULL);}
 if (global_mode) QDPInternal::broadcast(errcode);
 check_for_failure(errcode,"Failure while deleting file");
}

void IOHandler::file_open(IOHandler::openmode_type access_mode, 
                          const std::string& errmsg)
{
 if (global_mode)
    check_for_failure(MPI_File_open(MPI_COMM_WORLD, (char*)m_filename.c_str(), 
                      access_mode, finfo, &fh), errmsg);
 else
    check_for_failure(MPI_File_open(MPI_COMM_SELF, (char*)m_filename.c_str(), 
                      access_mode, finfo, &fh), errmsg);
}

void IOHandler::file_close()
{
 check_for_failure(MPI_File_close(&fh),"Failure during close");
}

const int IOHandler::IO_ERR_NO_SUCH_FILE= MPI_ERR_NO_SUCH_FILE;
const int IOHandler::IO_ERR_ACCESS=       MPI_ERR_ACCESS;
const int IOHandler::IO_ERR_OTHER=        MPI_ERR_OTHER;
const int IOHandler::IO_SUCCESS=          MPI_SUCCESS;

const IOHandler::openmode_type IOHandler::IO_MODE_RDONLY= MPI_MODE_RDONLY;
const IOHandler::openmode_type IOHandler::IO_MODE_RDWR=   MPI_MODE_RDWR;
const IOHandler::openmode_type IOHandler::IO_MODE_CREATE= MPI_MODE_CREATE;
const IOHandler::openmode_type IOHandler::IO_MODE_EXCL=   MPI_MODE_EXCL;

const IOHandler::whence_type IOHandler::IO_SEEK_BEG  = MPI_SEEK_SET;
const IOHandler::whence_type IOHandler::IO_SEEK_CUR  = MPI_SEEK_CUR;
const IOHandler::whence_type IOHandler::IO_SEEK_END  = MPI_SEEK_END;


void IOHandler::check_for_failure(int errcode, const std::string& mesg)
{
 if (errcode==MPI_SUCCESS) return;
 char errstring[MPI_MAX_ERROR_STRING];
 int len;
 MPI_Error_string(errcode, errstring, &len);
#ifdef USE_COLLECTIVE_IO
 if ((Layout::primaryNode())||(!global_mode))
#endif
    {std::cout << "IOHandler error with file "<<m_filename<<":"
               <<std::endl<<"  "<<mesg<<std::endl<<errstring<<std::endl;}
 QDP_abort(errcode);
}

void IOHandler::set_striping(int striping_factor, int striping_unit)
{
 if (striping_factor>1){
    MPI_Info_set(finfo,(char*)"striping_factor",
                (char*) int_to_string(striping_factor).c_str());}
 if (striping_unit>0){
    MPI_Info_set(finfo,(char*)"striping_unit",
                (char*) int_to_string(striping_unit).c_str());}
}

   //  Seek in MPI-IO is relative to the current file view!!
   //  offset 0 refers to the start of the view and is measured
   //  in terms of etypes.  To maintain a similarity with seek
   //  in fstream, it is necessary to return to the default file
   //  view after an file view set in a read/write.

void IOHandler::file_seek(off_type offset, whence_type whence)
{
 check_for_failure(MPI_File_seek(fh,offset,whence),
                   "Failure during seek");
}

IOHandler::pos_type IOHandler::file_tell()
{
 pos_type current;
 check_for_failure(MPI_File_get_position(fh,&current),"Failure during tell");
 return current;
}

         // print out MPI-IO file hints on primary node
 
void IOHandler::printFileHints()
{
 int i,nkeys,flag;
 MPI_Info info_used;
 char key[MPI_MAX_INFO_KEY], value[MPI_MAX_INFO_VAL];

 if ((Layout::primaryNode())||(!global_mode)){
    std::cout << "File = "<<m_filename<<std::endl<<"File hints:"<<std::endl;
    MPI_File_get_info(fh,&info_used);
    MPI_Info_get_nkeys(info_used, &nkeys);
    for (i=0; i<nkeys; i++){
       MPI_Info_get_nthkey(info_used,i,key);
       MPI_Info_get(info_used,key,MPI_MAX_INFO_VAL,value,&flag);
       std::cout << "    Key = "<<key<<", value = "<<value<<std::endl;}
    MPI_Info_free(&info_used);
    }
}


      // Write at current location in file. Write is done on primary 
      // node only, but all file pointers updated.  This routine assumes
      // the default file view (entire file, etype = MPI_BYTE).
      // Should be called by all nodes.

void IOHandler::writeCommon(const char *data, size_t nbytes)
{
 pos_type currdisp=file_tell();
 MPI_Status status;
 int errcode,count;
 if ((Layout::primaryNode())||(!global_mode)){
    errcode=MPI_File_write_at(fh,currdisp,(char*)data,nbytes,MPI_BYTE,&status);
    if (errcode==MPI_SUCCESS){
       errcode=MPI_Get_count(&status,MPI_BYTE,&count);
       if (uint(count)!=nbytes) errcode=MPI_ERR_COUNT;}}
 if (global_mode) QDPInternal::broadcast(errcode);
 check_for_failure(errcode,"Failure during common write");
 check_for_failure(MPI_File_seek(fh,currdisp+nbytes,MPI_SEEK_SET),
                   "Failure during common write");
}


void IOHandler::writeDistributed(const char *data, const DistArrayViewInfo& dinfo,
                                 bool datashift)
{
 if (!openflag)
    check_for_failure(IO_ERR_OTHER,"Write failure--no open file");
 if (read_only)
    check_for_failure(IO_ERR_ACCESS,"Write failure--read only file");
 if (!global_mode)
    check_for_failure(IO_ERR_ACCESS,"Distributed write failure: not global mode");
 if (read_mode){
    read_mode=false; checksum=0;}
 pos_type currdisp=file_tell();
              // set the file view
 check_for_failure(MPI_File_set_view(fh,currdisp,
                   dinfo.getFileElemType(),dinfo.getFileViewType(),
                   (char*)"native",finfo),"Failure during MPI_File_set_view");
 MPI_Status status;
 int errcode,count;
 char *dataptr=const_cast<char *>(data);
 if (datashift) dataptr+=dinfo.getElementOffsetThisNode()*dinfo.getBytesPerElement();
       // do the collective write!
#ifdef USE_COLLECTIVE_IO
 errcode=MPI_File_write_all(fh,dataptr,dinfo.getViewElementsThisNode(),
                            dinfo.getFileElemType(),&status);
#else
 bool sync=true;
 errcode=MPI_File_write(fh,dataptr,dinfo.getViewElementsThisNode(),
                        dinfo.getFileElemType(),&status);
 QDPInternal::broadcast(sync);
#endif
 if (errcode==MPI_SUCCESS){
    errcode=MPI_Get_count(&status,dinfo.getFileElemType(),&count);
    if (count!=dinfo.getViewElementsThisNode()) errcode=MPI_ERR_COUNT;}
 check_for_failure(errcode,"Failure during distributed write");
   // reset view to entire file and etype = bytes so seeks will work like fstream
 check_for_failure(MPI_File_set_view(fh,0,MPI_BYTE,MPI_BYTE,
                   (char*)"native",finfo),"Failure during MPI_File_set_view");
 check_for_failure(MPI_File_seek(fh,currdisp+dinfo.getViewTotalBytes(),
                   MPI_SEEK_SET),"Failure during distributed write");
} 


         // Read from current location in file. Memory for "data" 
         // must already be allocated.  File pointer updated.
         // Read is done on primary node. If "broadcast" is true,
         // data is broadcast to all other nodes. This routine assumes
         // the default file view (entire file, etype = MPI_BYTE).

void IOHandler::readCommon(char *data, size_t nbytes, bool broadcast)
{
 pos_type currdisp=file_tell();
 MPI_Status status;
 int errcode,count;
 if ((Layout::primaryNode())||(!global_mode)){
    errcode=MPI_File_read_at(fh,currdisp,data,nbytes,MPI_BYTE,&status);
    if (errcode==MPI_SUCCESS){
       errcode=MPI_Get_count(&status,MPI_BYTE,&count);
       if (uint(count)!=nbytes) errcode=MPI_ERR_COUNT;}}
 if (global_mode) QDPInternal::broadcast(errcode);
 check_for_failure(errcode,"Failure during primary read");
 if ((global_mode)&&(broadcast)) QDPInternal::broadcast(data,nbytes);
 check_for_failure(MPI_File_seek(fh,currdisp+nbytes,MPI_SEEK_SET),
                   "Failure during primary read");
}


void IOHandler::readDistributed(char *data, const DistArrayViewInfo& dinfo,
                                bool datashift)
{
 if (!openflag)
    check_for_failure(IO_ERR_OTHER,"Read failure--no open file");
 if (!global_mode)
    check_for_failure(IO_ERR_ACCESS,"Distributed read failure: not global mode");
 if (!read_mode){
    read_mode=true; checksum=0;}
 pos_type currdisp=file_tell();
 check_for_failure(MPI_File_set_view(fh,currdisp,
                   dinfo.getFileElemType(),dinfo.getFileViewType(),
                   (char*)"native",finfo),"Failure during MPI_File_set_view");
 MPI_Status status;
 int errcode,count;
 char *dataptr=data;
 if (datashift) dataptr+=dinfo.getElementOffsetThisNode()*dinfo.getBytesPerElement();
      // do the collective read!
#ifdef USE_COLLECTIVE_IO
 errcode=MPI_File_read_all(fh,dataptr,dinfo.getViewElementsThisNode(),
                           dinfo.getFileElemType(),&status);
#else
 bool sync=true;
 errcode=MPI_File_read(fh,dataptr,dinfo.getViewElementsThisNode(),
                       dinfo.getFileElemType(),&status);
 QDPInternal::broadcast(sync);
#endif
 if (errcode==MPI_SUCCESS){
    errcode=MPI_Get_count(&status,dinfo.getFileElemType(),&count);
    if (count!=dinfo.getViewElementsThisNode()) errcode=MPI_ERR_COUNT;}
 check_for_failure(errcode,"Failure during distributed read");
   // reset view to entire file and etype = bytes so seeks will work like fstream
 check_for_failure(MPI_File_set_view(fh,0,MPI_BYTE,MPI_BYTE,
                   (char*)"native",finfo),"Failure during MPI_File_set_view");
 check_for_failure(MPI_File_seek(fh,currdisp+dinfo.getViewTotalBytes(),
                   MPI_SEEK_SET),"Failure during distributed read");
} 

    //  check sum after distribution to nodes, but not yet checker boarded
    
void IOHandler::compute_lattice_checksum(char* data, size_t bytes_per_site,
                                         size_t lexicostart, size_t nsites)
{
 uint lastnode=0;
 uint thisnode=Layout::nodeNumber();
 uint node=0;
 const uint xinc = Layout::subgridLattSize()[0];
 size_t chunksize=bytes_per_site*xinc;
 uint count=0;
 size_t lexicostop=lexicostart+nsites;
      // must loop in order of sites
 for (uint site=lexicostart;site<lexicostop;site += xinc){
    multi1d<int> coord = crtesn(site, Layout::lattSize());
    node   = Layout::nodeNumber(coord);
    if (lastnode!=node){
            // send checksum from last node to new node
       QDPInternal::route(&checksum,lastnode,node,sizeof(QDPUtil::n_uint32_t));
       lastnode=node;}
    if (thisnode==node){
       checksum = QDPUtil::crc32(checksum, data+count*bytes_per_site, chunksize);
       count+=xinc;}
    QMP_barrier();
    }
 if (node!=0)  // send to primary node if did not end on primary node
    QDPInternal::route(&checksum,lastnode,0,sizeof(QDPUtil::n_uint32_t));
 QDPInternal::broadcast(checksum);
}


#endif


// *************************************************************************





IOHandler::IOHandler(bool global) : read_only(true), openflag(false), read_mode(true),
                                    global_mode(global), endian_format('U'), 
                                    endian_convert(false), checksum_on(false), 
                                    is_new_file(false), checksum(0), lattinfo(1)
{
#ifndef ARCH_SCALAR
 check_for_failure(MPI_Info_create(&finfo),"Failure during MPI_Info_create");
#ifdef ROMIO_DS_DISABLE
  MPI_Info_set(finfo,(char*)"romio_ds_read",(char*)"disable");
  MPI_Info_set(finfo,(char*)"romio_ds_write",(char*)"disable");
  MPI_Info_set(finfo,(char*)"romio_cb_read",(char*)"enable");
  MPI_Info_set(finfo,(char*)"romio_cb_write",(char*)"enable");
  MPI_Info_set(finfo,(char*)"romio_lustre_ds_in_coll",(char*)"disable");
#endif
#endif
// static_assert<sizeof(int)==4>();
 static_assert(sizeof(int)==4,"sizeof(int) must be 4");
}


IOHandler::IOHandler(const std::string& filename, OpenMode mode,
                     const std::string& filetype_id, char endianness,
                     int striping_factor, int striping_unit,
                     bool turn_on_checksum, bool global) : global_mode(global), lattinfo(1)
{
// static_assert<sizeof(int)==4>();
 static_assert(sizeof(int)==4,"sizeof(int) must be 4");
#ifndef ARCH_SCALAR
 check_for_failure(MPI_Info_create(&finfo),"Failure during MPI_Info_create");
#ifdef ROMIO_DS_DISABLE
  MPI_Info_set(finfo,(char*)"romio_ds_read",(char*)"disable");
  MPI_Info_set(finfo,(char*)"romio_ds_write",(char*)"disable");
  MPI_Info_set(finfo,(char*)"romio_lustre_ds_in_coll",(char*)"disable");
  MPI_Info_set(finfo,(char*)"romio_cb_read",(char*)"enable");
  MPI_Info_set(finfo,(char*)"romio_cb_write",(char*)"enable");
#endif
#endif
 openflag=false;
 open(filename,mode,filetype_id,endianness,striping_factor,striping_unit,
      turn_on_checksum);
}

IOHandler::IOHandler(const std::string& filename, std::ios_base::openmode mode,
                     const std::string& filetype_id, char endianness,
                     int striping_factor, int striping_unit,
                     bool turn_on_checksum, bool global) : global_mode(global), lattinfo(1)
{
// static_assert<sizeof(int)==4>();
 static_assert(sizeof(int)==4,"sizeof(int) must be 4");
#ifndef ARCH_SCALAR
 check_for_failure(MPI_Info_create(&finfo),"Failure during MPI_Info_create");
#ifdef ROMIO_DS_DISABLE
  MPI_Info_set(finfo,(char*)"romio_ds_read",(char*)"disable");
  MPI_Info_set(finfo,(char*)"romio_ds_write",(char*)"disable");
#endif
#endif
 openflag=false;
 open(filename,mode,filetype_id,endianness,striping_factor,striping_unit,
      turn_on_checksum);
}


void IOHandler::open(const std::string& filename, IOHandler::OpenMode mode,
                     const std::string& filetype_id, char endianness,
                     int striping_factor, int striping_unit,
                     bool turn_on_checksum)
{
 close();
 read_only=(mode==ReadOnly)?true:false;

 m_filename=tidyString(filename);
 if (m_filename.empty())
     check_for_failure(IO_ERR_NO_SUCH_FILE,"Empty file name");

 bool exists=fileExists();
 if ((mode==ReadOnly)&&(!exists))
    check_for_failure(IO_ERR_NO_SUCH_FILE,
           "Failure during ReadOnly open: file does not exist");
 else if ((mode==ReadWriteFailIfExists)&&(exists))
    check_for_failure(IO_ERR_ACCESS,
           "Failure during ReadWrite open: file exists and FailIfExists mode");
 else if ((mode==ReadWriteEraseIfExists)&&(exists)){
    delete_file();
    exists=false;}
 if (exists){
    IOHandler::openmode_type access=(mode==ReadOnly) ? 
                IO_MODE_RDONLY : IO_MODE_RDWR;
    is_new_file=false;
    open_existing_file(filetype_id,access);}
 else{
    is_new_file=true;
    open_new_file(filetype_id,endianness,striping_factor,striping_unit);}

 checksum_on=turn_on_checksum;
 checksum=0;
 read_mode=true;
}


void IOHandler::open(const std::string& filename, std::ios_base::openmode mode,
                     const std::string& filetype_id, char endianness,
                     int striping_factor, int striping_unit,
                     bool turn_on_checksum)
{
 IOHandler::OpenMode iomode;
 if (!(mode & std::ios_base::out)) iomode=ReadOnly;
 else{
    if (mode & std::ios_base::trunc) iomode=ReadWriteEraseIfExists;
    else iomode=ReadWriteFailIfExists;}
 open(filename,iomode,filetype_id,endianness,striping_factor,striping_unit,
      turn_on_checksum);
}


void IOHandler::openReadOnly(const std::string& filename, 
                             const std::string& filetype_id,
                             bool turn_on_checksum)
{
 open(filename,ReadOnly,filetype_id,'N',1,0,turn_on_checksum);
}


void IOHandler::openNew(const std::string& filename, bool fail_if_exists,
                        const std::string& filetype_id, char endianness,
                        int striping_factor, int striping_unit,
                        bool turn_on_checksum)
{
 OpenMode iomode=(fail_if_exists) ? ReadWriteFailIfExists : ReadWriteEraseIfExists;
 open(filename,iomode,filetype_id,endianness,striping_factor,striping_unit,
      turn_on_checksum);
}


void IOHandler::openUpdate(const std::string& filename,
                           const std::string& filetype_id, char endianness,
                           int striping_factor, int striping_unit,
                           bool turn_on_checksum)
{
 open(filename,ReadWriteUpdateIfExists,filetype_id,endianness,
      striping_factor,striping_unit,turn_on_checksum);
}



void IOHandler::open_existing_file(const std::string& filetype_id,
                                   IOHandler::openmode_type access_mode)
{
 file_open(access_mode,"Could not open existing file: "+m_filename);
     // get endian info, and check file id
 openflag=true;
 std::string ID_string;
 readIDstring(endian_format,ID_string);  // on primary node only
 bool flag=false;
 if ((Layout::primaryNode())||(!global_mode)){
    flag=true;
    if ((endian_format!='B')&&(endian_format!='L')){
       flag=false; std::cerr << "Invalid endian format"<<std::endl;}
    if (tidyString(filetype_id)!=ID_string){
       flag=false; 
       std::cerr << "File = "<<m_filename<<std::endl;
       std::cerr << "File ID mismatch:"<<std::endl;
       std::cerr << "File contains ID: <"<<ID_string<<">"<<std::endl;
       std::cerr << "ID requested was: <"<<tidyString(filetype_id)<<">"<<std::endl;}}
 if (global_mode) QDPInternal::broadcast(flag);
 if (!flag)
    check_for_failure(IO_ERR_OTHER,"Error during open");
 if (global_mode) QDPInternal::broadcast(endian_format);
 if (QDPUtil::big_endian())
    endian_convert=(endian_format=='L')?true:false;
 else 
    endian_convert=(endian_format=='B')?true:false;
}


void IOHandler::open_new_file(const std::string& filetype_id,
                              char endianness,
                              int striping_factor, int striping_unit)
{
 if (endianness=='N'){  // native 
    endian_format=QDPUtil::big_endian() ? 'B':'L';
    endian_convert=false;}
 else if ((endianness!='B')&&(endianness!='L'))
    check_for_failure(IO_ERR_OTHER,"Invalid endian format");
 else{
    endian_format=endianness;
    if (QDPUtil::big_endian()) 
       endian_convert=(endian_format=='L')?true:false;
    else 
       endian_convert=(endian_format=='B')?true:false;}

 set_striping(striping_factor,striping_unit);
 IOHandler::openmode_type amode = IO_MODE_RDWR | IO_MODE_CREATE;
 
 file_open(amode,"Could not open file "+m_filename);

 openflag=true;
 writeIDstring(filetype_id);
}


       // Destructor

IOHandler::~IOHandler() 
{
 if (openflag) clear(); 
#ifndef ARCH_SCALAR
 MPI_Info_free(&finfo);
#endif
}


     // Open file, read string at a particular position, then close. Returns
     // true if the file exists and can be opened and its file type matches 
     // "filetype_id", returns false otherwise.  The string is returned in 
     // "stringvalue".  This routine is used by objects in data_io_handler.h 
     // when the multi-file handlers build up their maps of file keys.
           
bool IOHandler::peekString(std::string& stringvalue, size_t byte_offset,
                           const std::string& filename, const std::string& filetype_id)
{
 stringvalue.clear();
 std::string fname=tidyString(filename);
 if (fname.empty()) return false;
 bool flag=false;
 if ((Layout::primaryNode())||(!global_mode)){ 
    flag=peeker(stringvalue,byte_offset,fname,filetype_id);}
 if (global_mode){
    QDPInternal::broadcast(flag);
    QDPInternal::broadcast_str(stringvalue);}
 return flag;
}



bool IOHandler::peeker(std::string& stringvalue, size_t byte_offset,
                       const std::string& fname, const std::string& filetype_id)
{
 ifstream in(fname.c_str(), std::ios::binary | std::ios::in);
 if (!in) return false;
 std::string ID_string(ID_string_length+1,' ');
 in.read((char*)&ID_string[0],ID_string_length+1);
 if (in.fail()) return false;
 char endian=ID_string[0];
 if ((endian!='B')&&(endian!='L')) return false;
 ID_string.erase(0,1);
 ID_string=tidyString(ID_string);
 if (tidyString(filetype_id)!=ID_string) return false;
 in.seekg(pos_type(byte_offset),std::ios::cur);
 if (in.fail()) return false;
 unsigned int n;
 in.read((char*)&n,sizeof(int));
 if (in.fail()) return false;
 if (QDPUtil::big_endian())
    endian_convert=(endian=='L')?true:false;
 else 
    endian_convert=(endian=='B')?true:false;
 if (endian_convert) QDPUtil::byte_swap(&n,sizeof(int),1);
 if (n>16777216) return false;  // too large for string...must be corrupt data
 stringvalue.resize(n);
 in.read((char*)&stringvalue[0],sizeof(char)*n);
 if (in.fail()){ stringvalue.clear(); return false;}
 return true;       // in destructor will close the file
}



      // Clear, then reset finfo

void IOHandler::close()
{ 
 if (!openflag) return;
 clear();
#ifndef ARCH_SCALAR
 MPI_Info_free(&finfo);
 check_for_failure(MPI_Info_create(&finfo),"Failure during MPI_Info reset");
#endif
}

     // Close the file, reset all data members except finfo

void IOHandler::clear()
{ 
 file_close();
 m_filename.clear();
 openflag=false;
 endian_convert=false;
 endian_format='U';  // undefined
 checksum=0;
 read_only=true;
 read_mode=true;
 checksum_on=false;
 is_new_file=false;
}

         //  Set the file pointer relative to start of file,
         //  end of file (positive is backward), or current location

void IOHandler::seekFromStart(off_type offset)
{
 file_seek(offset+data_start_pos,IOHandler::IO_SEEK_BEG);
 checksum=0;
}

void IOHandler::seekFromCurr(off_type offset)
{ 
 file_seek(offset,IOHandler::IO_SEEK_CUR);
 checksum=0;
}
 
void IOHandler::seekFromEnd(off_type offset)
{ 
 file_seek(offset,IOHandler::IO_SEEK_END);
 checksum=0;
}

void IOHandler::seek(pos_type offset)
{ 
 file_seek(offset+data_start_pos,IOHandler::IO_SEEK_BEG);
 checksum=0;
}

void IOHandler::seekBegin(off_type offset)
{ 
 file_seek(offset+data_start_pos,IOHandler::IO_SEEK_BEG);
 checksum=0;
}

void IOHandler::seekRelative(off_type offset)
{ 
 file_seek(offset,IOHandler::IO_SEEK_CUR);
 checksum=0;
}

void IOHandler::seekEnd(off_type offset)
{ 
 file_seek(offset,IOHandler::IO_SEEK_END);
 checksum=0;
}

void IOHandler::rewind()
{ 
 file_seek(data_start_pos,IOHandler::IO_SEEK_BEG);
 checksum=0;
}


      //  Get the file pointer location in bytes from start of data

IOHandler::pos_type IOHandler::currentPosition()
{
 return file_tell()-data_start_pos;
}

IOHandler::pos_type IOHandler::tell()
{
 return file_tell()-data_start_pos;
}


void IOHandler::turnOnChecksum()
{
 if (checksum_on) return;
 checksum_on=true;
 checksum=0;
}

void IOHandler::turnOffChecksum()
{
 checksum_on=false;
}

void IOHandler::resetChecksum()
{
 checksum=0;
}

QDPUtil::n_uint32_t IOHandler::getChecksum()
{
 if (checksum_on) return checksum;
 check_for_failure(true,"Invalid call to getChecksum since checksums not turned on");
 return checksum;  // to avoid compiler warnings
}



        // Print out file ID on primary node. Does not 
        // change file pointers.

void IOHandler::printFileID()
{
 char endian;
 std::string ID_string;
 pos_type curr=file_tell();
 readIDstring(endian,ID_string);
 file_seek(curr,IOHandler::IO_SEEK_BEG);
 if ((Layout::primaryNode())||(!global_mode)){
    std::cout << "File = "<<m_filename<<std::endl<<"ID string = <"
              <<tidyString(ID_string)<<">"<<std::endl;}
}

   // Write ID string at start of file. Length of string
   // is always "ID_string_length".  First character will
   // be 'B' or 'L' to indicate endian format to use.

void IOHandler::writeIDstring(const std::string& ID_string)
{
 check_for_failure(int(ID_string.length())>ID_string_length,
                  "IOHandler file ID string too long: cannot exceed "
                  +int_to_string(ID_string_length)+" characters");
 std::string buf(1,endian_format);
 buf+=ID_string;
 int nblanks=ID_string_length-ID_string.length();
 if (nblanks>0)
    buf+=string(nblanks,' ');
 file_seek(0,IOHandler::IO_SEEK_BEG);
 writeCommon(&buf[0],ID_string_length+1);
}


   // Reads endian format and ID string from current file.
   // Read is done only on primary node and is NOT broadcast.
   // Results are returned in "endianness" and "ID_string"
   // on primary node.  All file pointers updated.

void IOHandler::readIDstring(char &endianness, std::string& ID_string)
{
 if ((Layout::primaryNode())||(!global_mode)){
    ID_string.resize(ID_string_length+1);}
 file_seek(0,IOHandler::IO_SEEK_BEG);
 readCommon(&ID_string[0],ID_string_length+1,false);
 if ((Layout::primaryNode())||(!global_mode)){
    endianness=ID_string[0];
    ID_string.erase(0,1);
    ID_string=tidyString(ID_string);}
}




      // Checks if a file exists.
      
bool IOHandler::fileExists()
{
 bool result=false;
 if ((Layout::primaryNode())||(!global_mode)){
    result = (access(m_filename.c_str(),F_OK) == 0) ? true : false;}
 if (global_mode) QDPInternal::broadcast(result);
 return result;
}


      // Converts an integer to a string
      
std::string IOHandler::int_to_string(int intval)
{ 
 std::ostringstream oss;
 oss << intval;
 return oss.str();
}

      // Removes leading and trailing blanks in a string

std::string IOHandler::tidyString(const std::string& str)   
{
 string tmp;
 for (unsigned int i=0;i<str.length();i++)
    if ((str[i]!='\n')&&(str[i]!='\t')&&(str[i]!='\r'))
       tmp.push_back(str[i]);
 int start=tmp.find_first_not_of(" ");
 if (start==int(std::string::npos)) return "";
 int len=tmp.find_last_not_of(" ")-start+1;
 return tmp.substr(start,len);
}


// *******************************************************************

      //   Main input/output routines that do the byte-swapping (if needed),
      //   update the check sum, and do the read/write.

void IOHandler::write_common(const char* output, size_t elementbytes, size_t nelements)
{
 if (!openflag)
    check_for_failure(IO_ERR_OTHER,"Write failure--no open file");
 if (read_only)
    check_for_failure(IO_ERR_ACCESS,"Write failure--read only file");
 if (read_mode){
    read_mode=false; checksum=0;}
 if (endian_convert){
    QDPUtil::byte_swap(const_cast<char *>(output), elementbytes, nelements);
    if (checksum_on) checksum = QDPUtil::crc32(checksum, output, elementbytes*nelements);
    writeCommon(output, elementbytes*nelements);
    QDPUtil::byte_swap(const_cast<char *>(output), elementbytes, nelements);}
 else{
    if (checksum_on) checksum = QDPUtil::crc32(checksum, output, elementbytes*nelements);
    writeCommon(output, elementbytes*nelements);}
}


void IOHandler::read_common(char* input, size_t elementbytes, size_t nelements)
{
 if (!openflag)
    check_for_failure(IO_ERR_OTHER,"Read failure--no open file");
 if (!read_mode){
    read_mode=true; checksum=0;}
 readCommon(input, elementbytes*nelements);
 if (checksum_on) checksum = QDPUtil::crc32(checksum, input, elementbytes*nelements);
 if (endian_convert) QDPUtil::byte_swap(input, elementbytes, nelements);
}


      //   Multi versions read/write via a char* buffer.
      //   Addresses of memory, sizes of each element, and numbers
      //   of elements are specified in "dinfo".

size_t IOHandler::multi_write_common(const list<IOHandler::outbasicdatainfo>& dinfo)
{
 if (!openflag)
    check_for_failure(IO_ERR_OTHER,"Write failure--no open file");
 if (read_only)
    check_for_failure(IO_ERR_ACCESS,"Write failure--read only file");
 if (read_mode){
    read_mode=false; checksum=0;}
 size_t totalbytes=0;
 for (list<IOHandler::outbasicdatainfo>::const_iterator it=dinfo.begin();it!=dinfo.end();it++)
    totalbytes+=it->elsize;
 char *buf=new(nothrow) char[totalbytes];   
 check_for_failure(!buf,"IOHandler::write---unable to allocate memory");
 char* pt=buf;
 for (list<IOHandler::outbasicdatainfo>::const_iterator it=dinfo.begin();it!=dinfo.end();it++){
    memcpy(pt,it->dataptr,it->elsize);
    if (endian_convert) QDPUtil::byte_swap(pt, it->elsize, 1);
    pt+=it->elsize;}
 if (checksum_on) checksum = QDPUtil::crc32(checksum, buf, totalbytes);
 writeCommon(buf, totalbytes);
 delete [] buf;
 return totalbytes;
}


size_t IOHandler::multi_read_common(const list<IOHandler::inbasicdatainfo>& dinfo)
{
 if (!openflag)
    check_for_failure(IO_ERR_OTHER,"Read failure--no open file");
 if (!read_mode){
    read_mode=true; checksum=0;}
 size_t totalbytes=0;
 for (list<IOHandler::inbasicdatainfo>::const_iterator it=dinfo.begin();it!=dinfo.end();it++)
    totalbytes+=it->elsize;
 char *buf=new(nothrow) char[totalbytes];
 check_for_failure(!buf,"IOHandler::read---unable to allocate memory");
 readCommon(buf, totalbytes);
 if (checksum_on) checksum = QDPUtil::crc32(checksum, buf, totalbytes);
 char* pt=buf;
 for (list<IOHandler::inbasicdatainfo>::const_iterator it=dinfo.begin();it!=dinfo.end();it++){
    if (endian_convert) QDPUtil::byte_swap(pt, it->elsize, 1);
    memcpy(it->dataptr,pt,it->elsize);
    pt+=it->elsize;}
 delete [] buf;
 return totalbytes;
}


      //   Main input/output routines that do the byte-swapping (if needed),
      //   update the check sum, and do the read/write.  For lattice quantities,
      //   must also deal with checkerboarding.  "bytes_per_site" is the total number
      //   of bytes needed to define the entire quantity on one lattice site (includes
      //   color matrix, spin vectors, etc.).  "bytes_per_word" 
      //   is the number of bytes in the underlying numbers or words (int, double, complex, 
      //   etc) of each component of the quantity at each site.  This information is needed
      //   for byte-swapping.

void IOHandler::write_lattice(const char* output, size_t bytes_per_word, 
                              const DistArrayViewInfo& dinfo)
{
 if (!openflag)
    check_for_failure(IO_ERR_OTHER,"Write failure--no open file");
 if (read_only)
    check_for_failure(IO_ERR_ACCESS,"Write failure--read only file");
 if (read_mode){
    read_mode=false; checksum=0;}
 IOH_int bytes_this_node=dinfo.getViewBytesThisNode();
#ifdef QDP_USE_LEXICO_LAYOUT
 char* buf=const_cast<char *>(output);
 bool datashift=true;
 char* b=buf+dinfo.getElementOffsetThisNode()*dinfo.getBytesPerElement();
#else
 char* buf=0;
 bool datashift=false; 
 if (bytes_this_node>0){
    buf=new(nothrow) char[bytes_this_node];
    check_for_failure(!buf,"Could not allocate buffer in IOHandler::write_lattice");
    lexify(buf,output,dinfo);} 
 char* b=buf;
#endif 
 if (endian_convert){
    size_t words_this_node=bytes_this_node/bytes_per_word;
    QDPUtil::byte_swap(b, bytes_per_word, words_this_node);
    if (checksum_on) 
       compute_lattice_checksum(b,dinfo.getBytesPerSite(),dinfo.getViewLexicoStart(),
                                dinfo.getViewLexicoSpan());
    writeDistributed(buf, dinfo, datashift);
    if (datashift) 
       QDPUtil::byte_swap(b, bytes_per_word, words_this_node);
    }
 else{ 
    if (checksum_on){
       compute_lattice_checksum(b,dinfo.getBytesPerSite(),dinfo.getViewLexicoStart(),
                                dinfo.getViewLexicoSpan());}
    writeDistributed(buf,dinfo,datashift);}
#ifndef QDP_USE_LEXICO_LAYOUT
 if (buf!=0) delete [] buf;
#endif
}



void IOHandler::read_lattice(char* input, size_t bytes_per_word, 
                             const DistArrayViewInfo& dinfo)
{
 if (!openflag)
    check_for_failure(IO_ERR_OTHER,"Read failure--no open file");
 if (!read_mode){
    read_mode=true; checksum=0;}
 IOH_int bytes_this_node=dinfo.getViewBytesThisNode();
#if ( QDP_USE_LEXICO_LAYOUT == 1 )
 char* buf=input;
 bool datashift=true;
 char* b=buf+dinfo.getElementOffsetThisNode()*dinfo.getBytesPerElement();
#else
 char* buf=0;
 bool datashift=false;
 if (bytes_this_node>0){
    buf=new(nothrow) char[bytes_this_node];
    check_for_failure(!buf,"Could not allocate buffer in IOHandler::write_lattice");}
 char* b=buf;
#endif
 readDistributed(buf, dinfo, datashift);
 if (checksum_on)
    compute_lattice_checksum(b,dinfo.getBytesPerSite(),dinfo.getViewLexicoStart(),
                             dinfo.getViewLexicoSpan());
 if (endian_convert){
    size_t words_this_node=bytes_this_node/bytes_per_word;
    QDPUtil::byte_swap(b, bytes_per_word, words_this_node);}
#ifndef QDP_USE_LEXICO_LAYOUT
 unlexify(input,buf,dinfo);
 delete [] buf;
#endif 
}



void IOHandler::lexify(char *out, const char* in, const DistArrayViewInfo& dinfo)
{
 IOH_int nelements=dinfo.getViewElementsThisNode();
 if (nelements==0) return;
 IOH_int bytes_per_element=dinfo.getBytesPerElement();
 IOH_int offset=dinfo.getElementOffsetThisNode();
#pragma omp parallel for
 for (IOH_int k=0;k<nelements;k++){
    IOH_int kk=k+offset;
       // get lattice site coordinate from local linear index (lexico ordering)
    multi1d<int> coord(QDP::Nd);
    for (int i=0; i<QDP::Nd; i++) {
       coord[i]=Layout::nodeCoord()[i]*Layout::subgridLattSize()[i]
               +(kk % Layout::subgridLattSize()[i]);
       kk /= Layout::subgridLattSize()[i];}
        // now get the actual local linear index with checkerboarding
    IOH_int lindex=Layout::linearSiteIndex(coord);
    memcpy(out+k*bytes_per_element, in+lindex*bytes_per_element, bytes_per_element);}
}


void IOHandler::unlexify(char *out, const char* in, const DistArrayViewInfo& dinfo)
{
 IOH_int nelements=dinfo.getViewElementsThisNode();
 if (nelements==0) return;
 IOH_int bytes_per_element=dinfo.getBytesPerElement();
 IOH_int offset=dinfo.getElementOffsetThisNode();
#pragma omp parallel for
 for (IOH_int k=0;k<nelements;k++){
    IOH_int kk=k+offset;
       // get lattice site coordinate from local linear index (lexico ordering)
    multi1d<int> coord(QDP::Nd);
    for (int i=0; i<QDP::Nd; i++) {
       coord[i]=Layout::nodeCoord()[i]*Layout::subgridLattSize()[i]
               +(kk % Layout::subgridLattSize()[i]);
       kk /= Layout::subgridLattSize()[i];}
        // now get the actual local linear index with checkerboarding
    IOH_int lindex=Layout::linearSiteIndex(coord);
    memcpy(out+lindex*bytes_per_element, in+k*bytes_per_element, bytes_per_element);}
}



// *******************************************************************


void IOHandler::write(const std::string& output)
{
 int n=output.length();
 write_common((const char*)&n, sizeof(int), 1);
 write_common(output.data(), sizeof(char), n);
}


void IOHandler::read(std::string& input, bool broadcast)
{
 int n;
 read_common((char*)&n, sizeof(int), 1);
    // if reading in wrong location, could get nonsense here,
    // so limit to a 16MB string
 check_for_failure(n>16777216, "string for read too large...bad file location?");
 char* str = new(nothrow) char[n];
 check_for_failure((str==0),"IOHandler::read---unable to allocate memory");
 readCommon(str, sizeof(char)*n, broadcast);
 if ((Layout::primaryNode())||(!global_mode)||(global_mode && broadcast)){
    if (checksum_on) checksum = QDPUtil::crc32(checksum, str, sizeof(char)*n);
    input.assign(str, n);}
 if (global_mode && (!broadcast) && checksum_on) QDPInternal::broadcast(checksum);
 delete[] str;
}

    // the number of bytes that these quantities occupy in an
    // IOHandler file

size_t numbytes(IOHandler& ioh, const std::string& output)
{
 return sizeof(int)+output.length();
}



// *************************************************

void IOHandler::startDelayedWrite()
{
 outdelay.clear();
}

void IOHandler::startDelayedRead()
{
 indelay.clear();
}

size_t IOHandler::finishDelayedWrite()
{
 size_t nbytes=multi_write_common(outdelay);
 outdelay.clear();
 return nbytes;
}

size_t IOHandler::finishDelayedRead()
{
 size_t nbytes=multi_read_common(indelay);
 indelay.clear();
 return nbytes;
}


// ********************************************************


#ifndef ARCH_SCALAR

 // **************************************************************************

  /*  Objects of this class contain info about multi-dimensional arrays 
      distributed across the MPI nodes.  Blocking assumed, column major.   
      Assumed distributed in a regular way...same size on each node.  
      Distribution on nodes is column major.  This class is mainly
      used by BinaryFileHandler for creating file views for MPI-IO. 
      
      The main purpose of this class is to create an appropriate "filetype"
      to pass to MPI-IO.  The routine  MPI_Type_create_subarray is used.  */


         // copy constructors
        
DistArrayViewInfo::DistArrayViewInfo(const DistArrayViewInfo& in) 
                     :    ndim(in.ndim), global_sizes(in.global_sizes),
                          local_sizes(in.local_sizes),
                          start_indices(in.start_indices),
                          nelem_this_node(in.nelem_this_node),
                          gvvol(in.gvvol), lexico_start(in.lexico_start),
                          el_offset_this_node(in.el_offset_this_node),
                          nbytes_per_element(in.nbytes_per_element)
{
 MPI_Type_dup(in.ftype,&ftype);
 MPI_Type_dup(in.etype,&etype);
}
    

DistArrayViewInfo& DistArrayViewInfo::operator=(const DistArrayViewInfo& in)
{
 ndim=in.ndim; 
 global_sizes=in.global_sizes;
 local_sizes=in.local_sizes;
 start_indices=in.start_indices;
 nelem_this_node=in.nelem_this_node;
 gvvol=in.gvvol; 
 lexico_start=in.lexico_start;
 el_offset_this_node=in.el_offset_this_node;
 nbytes_per_element=in.nbytes_per_element;
 MPI_Type_dup(in.ftype,&ftype);
 MPI_Type_dup(in.etype,&etype);
 return *this;
}

         // destructor

DistArrayViewInfo::~DistArrayViewInfo()
{
 MPI_Type_free(&ftype);
 MPI_Type_free(&etype);
}

 

void DistArrayViewInfo::setup_full(const multi1d<int>& gsizes, const multi1d<int>& numnodes,
                                   const multi1d<int>& logical_coord)
{
 ndim=gsizes.size();
 global_sizes=gsizes;
 local_sizes.resize(ndim);
 start_indices.resize(ndim);
 IOH_int lavol=1;
 gvvol=1;
 for (int k=0;k<ndim;k++){
    local_sizes[k]=global_sizes[k]/numnodes[k];
    start_indices[k]=logical_coord[k]*local_sizes[k];
    gvvol*=global_sizes[k];
    lavol*=local_sizes[k];}
 nelem_this_node=lavol;
 lexico_start=0;
 el_offset_this_node=0;
}                      


   //  The major index is the most slowly varying one.  Let NG be the global size of
   //  the major index, and NL be the local size of the major index.  NG/NL is the logical 
   //  size in this direction.  The major logical coordinate ranges from 0..NG/NL.
   //  Let LI be the logical coordinate that contains the global major index of value
   //  "major_index".  This routine is basically the same as "setup_full", except we 
   //  switch major logical coordinate 0 with the major logical coordinate LI.
   //  (This keeps the file views of different nodes non-overlapping.)

void DistArrayViewInfo::setup_major_slice(const multi1d<int>& gsizes, const multi1d<int>& numnodes,
                                          const multi1d<int>& logical_coord, int major_index)
{
 ndim=gsizes.size();
 global_sizes=gsizes;
 local_sizes.resize(ndim);
 start_indices.resize(ndim);
 for (int k=0;k<ndim;k++){
    local_sizes[k]=global_sizes[k]/numnodes[k];
    start_indices[k]=logical_coord[k]*local_sizes[k];}
 IOH_int lavol=1;
 gvvol=1;
 for (int k=0;k<ndim-1;k++){
    gvvol*=global_sizes[k];
    lavol*=local_sizes[k];}
 int nMajor=global_sizes[ndim-1];
   //  determine if this node contains the major slice
 int majorblocksize=nMajor/numnodes[ndim-1];
 int logicalkeep=major_index/majorblocksize;  // contains "major_index"
 if (logical_coord[ndim-1]==logicalkeep){  // this node participates
    start_indices[ndim-1]=0;
    nelem_this_node=lavol; // /majorblocksize;
    el_offset_this_node=nelem_this_node*( major_index % majorblocksize);}
 else{
    if (logical_coord[ndim-1]==0) 
        start_indices[ndim-1]=logicalkeep*local_sizes[ndim-1];
    nelem_this_node=0;
    el_offset_this_node=0;}
 lexico_start=gvvol*major_index; 

/*
 ndim=gsizes.size();
 global_sizes=gsizes;
 local_sizes.resize(ndim);
 start_indices.resize(ndim);
 for (int k=0;k<ndim;k++){
    local_sizes[k]=global_sizes[k]/numnodes[k];
    start_indices[k]=logical_coord[k]*local_sizes[k];}
 start_indices[ndim-1]=0;
 int lavol=1;
 gvvol=1;
 for (int k=0;k<ndim-1;k++){
    gvvol*=global_sizes[k];
    lavol*=local_sizes[k];}
 int nMajor=global_sizes[ndim-1];
   //  determine if this node contains the major slice
 int majorblocksize=nMajor/numnodes[ndim-1];
 if ((major_index/majorblocksize)==logical_coord[ndim-1]){  // this node participates
    nelem_this_node=lavol; // /majorblocksize;
    el_offset_this_node=nelem_this_node*( major_index % majorblocksize);}
 else{
    nelem_this_node=0;
    el_offset_this_node=0;}
 lexico_start=gvvol*major_index; */
}



void DistArrayViewInfo::set_bytes(IOH_int bytes_per_element)
{
 error_return(bytes_per_element<=0,"invalid input set_bytes");
 nbytes_per_element=bytes_per_element;
}



bool DistArrayViewInfo::checker(const multi1d<int>& gsizes, 
                                const multi1d<int>& numnodes)
{
 if ((gsizes.size()==0)||(gsizes.size()!=numnodes.size()))
    return false;
 for (int k=0;k<gsizes.size();k++)
    if ((gsizes[k]<=0)||(numnodes[k]<=0))
       return false;
 for (int k=0;k<gsizes.size();k++){
    int ll=gsizes[k]/numnodes[k];
    if (numnodes[k]*ll != gsizes[k]) return false;}
 int check=1;
 for (int k=0;k<gsizes.size();k++) check*=numnodes[k];
 if (check!=Layout::numNodes()) return false;
 return true;
}


void DistArrayViewInfo::error_return(bool cond, const std::string& msg)
{
 if (cond){
    std::cerr << "Error in DistArrayViewInfo: "<<msg<<std::endl;
    QDP_abort(1);}
}

         // constructor for full lattice quantity
         
DistArrayViewInfo::DistArrayViewInfo(IOH_int bytes_per_site)
{
 setup_full(Layout::lattSize(),Layout::logicalSize(),
            Layout::nodeCoord());
 set_bytes(bytes_per_site);
    // make and store the MPI filetype
 create_filetype();
} 

         // constructor for lattice time-slice
        
DistArrayViewInfo::DistArrayViewInfo(int time_slice, IOH_int bytes_per_site)
{
 setup_major_slice(Layout::lattSize(),Layout::logicalSize(),
                   Layout::nodeCoord(),time_slice);
 set_bytes(bytes_per_site);
    // make and store the MPI filetype
 create_filetype();
} 


         // constructor for distributed array (non-lattice), let MPI
         // distribute over the nodes
         

DistArrayViewInfo::DistArrayViewInfo(const multi1d<int>& gsizes, 
                                     IOH_int bytes_per_element)
{
 int nnodes=Layout::numNodes();
 int nd=gsizes.size();
 error_return(nd==0,"bad global sizes");
 multi1d<int> numnodes(nd); numnodes=0;
 error_return(MPI_Dims_create(nnodes,nd,&numnodes[0]),
          "Could not distribute array among processors");
 error_return(!checker(gsizes,numnodes),"invalid input");
 int rank=Layout::nodeNumber();
 multi1d<int> logical_coord(nd);   
   //  given the rank of this node, determine the logical coordinates 
   //  (assume column major; indices to the left vary faster)
 for (int i=0; i<nd; i++) {
    logical_coord[i] = rank % numnodes[i];
    rank /= numnodes[i];}
 setup_full(gsizes,numnodes,logical_coord);
 set_bytes(bytes_per_element);

    // make and store the MPI filetype
 create_filetype();
}


DistArrayViewInfo::DistArrayViewInfo(const multi1d<int>& gsizes, 
                                     const multi1d<int>& numnodes,
                                     IOH_int bytes_per_element)
{
 error_return(!checker(gsizes,numnodes),"invalid input");
 int rank=Layout::nodeNumber();
 int nd=gsizes.size();
 multi1d<int> logical_coord(nd);   
   //  given the rank of this node, determine the logical coordinates 
   //  (assume column major; indices to the left vary faster)
 for (int i=0; i<nd; i++) {
    logical_coord[i] = rank % numnodes[i];
    rank /= numnodes[i];}
 setup_full(gsizes,numnodes,logical_coord);
 set_bytes(bytes_per_element);

    // make and store the MPI filetype
 create_filetype();
}


DistArrayViewInfo::DistArrayViewInfo(const multi1d<int>& gsizes, 
                                     const multi1d<int>& numnodes,
                                     int major_index, IOH_int bytes_per_element)
{
 error_return(!checker(gsizes,numnodes),"invalid input");
 int rank=Layout::nodeNumber(); 
 int nd=gsizes.size();
 multi1d<int> logical_coord(nd);   
   //  given the rank of this node, determine the logical coordinates 
   //  (assume column major; indices to the left vary faster)
 for (int i=0; i<nd; i++) {
    logical_coord[i] = rank % numnodes[i];
    rank /= numnodes[i];}
 error_return((major_index<0)||(major_index>=gsizes[nd-1]),"bad major index");
 setup_major_slice(gsizes,numnodes,logical_coord,major_index);
 set_bytes(bytes_per_element);

    // make and store the MPI filetype
 create_filetype();
}

void DistArrayViewInfo::create_filetype()
{
 MPI_Type_contiguous(nbytes_per_element,MPI_BYTE,&etype);
 MPI_Type_commit(&etype);
 error_return(MPI_Type_create_subarray(ndim,&global_sizes[0],&local_sizes[0],
                  &start_indices[0],MPI_ORDER_FORTRAN,etype,&ftype),
                  "Could not create filetype");            
 MPI_Type_commit(&ftype);
}


void DistArrayViewInfo::resetBytes(IOH_int bytes_per_element)
{
 if (bytes_per_element==nbytes_per_element) return;
 set_bytes(bytes_per_element);
 MPI_Type_free(&ftype);
 MPI_Type_free(&etype);
 create_filetype();
}


#else
 // ***************************************************************
               // serial code

DistArrayViewInfo::DistArrayViewInfo(const DistArrayViewInfo& in) 
                     :    ndim(in.ndim), local_sizes(in.local_sizes),
                          gvvol(in.gvvol), lexico_start(in.lexico_start),
                          nbytes_per_element(in.nbytes_per_element) {}
    

DistArrayViewInfo& DistArrayViewInfo::operator=(const DistArrayViewInfo& in)
{
 ndim=in.ndim; 
 local_sizes=in.local_sizes;
 gvvol=in.gvvol; 
 lexico_start=in.lexico_start;
 nbytes_per_element=in.nbytes_per_element;
 return *this;
}

         // destructor

DistArrayViewInfo::~DistArrayViewInfo() {}

 

void DistArrayViewInfo::setup_full(const multi1d<int>& gsizes)
{
 ndim=gsizes.size();
 local_sizes=gsizes;
 gvvol=1;
 for (int k=0;k<ndim;k++)
    gvvol*=gsizes[k];
 lexico_start=0;
}                      


void DistArrayViewInfo::setup_major_slice(const multi1d<int>& gsizes,
                                          int major_index)
{
 ndim=gsizes.size();
 local_sizes=gsizes;
 gvvol=1;
 for (int k=0;k<ndim-1;k++)
    gvvol*=gsizes[k];
 lexico_start=gvvol*major_index;
}



void DistArrayViewInfo::set_bytes(IOH_int bytes_per_element)
{
 error_return(bytes_per_element<=0,"invalid input set_bytes");
 nbytes_per_element=bytes_per_element;
}



bool DistArrayViewInfo::checker(const multi1d<int>& gsizes)
{
 if (gsizes.size()==0)
    return false;
 for (int k=0;k<gsizes.size();k++)
    if (gsizes[k]<=0) return false;
 return true;
}


void DistArrayViewInfo::error_return(bool cond, const std::string& msg)
{
 if (cond){
    std::cerr << "Error in DistArrayViewInfo: "<<msg<<std::endl;
    QDP_abort(1);}
}

         // constructor for full lattice quantity
         
DistArrayViewInfo::DistArrayViewInfo(IOH_int bytes_per_site)
{
 setup_full(Layout::lattSize());
 set_bytes(bytes_per_site);
} 

         // constructor for lattice time-slice
        
DistArrayViewInfo::DistArrayViewInfo(int time_slice, IOH_int bytes_per_site)
{
 setup_major_slice(Layout::lattSize(),time_slice);
 set_bytes(bytes_per_site);
} 


         // constructor for distributed array (non-lattice), let MPI
         // distribute over the nodes
         

DistArrayViewInfo::DistArrayViewInfo(const multi1d<int>& gsizes, 
                                     IOH_int bytes_per_element)
{
 error_return(!checker(gsizes),"invalid input");
 setup_full(gsizes);
 set_bytes(bytes_per_element);
}

DistArrayViewInfo::DistArrayViewInfo(const multi1d<int>& gsizes, 
                                     const multi1d<int>& numnodes,
                                     IOH_int bytes_per_element)
{
 error_return(!checker(gsizes),"invalid input");
 setup_full(gsizes);
 set_bytes(bytes_per_element);
}


DistArrayViewInfo::DistArrayViewInfo(const multi1d<int>& gsizes, 
                                     const multi1d<int>& numnodes,
                                     int major_index, IOH_int bytes_per_element)
{
 error_return(!checker(gsizes),"invalid input");
 setup_major_slice(gsizes,major_index);
 set_bytes(bytes_per_element);
} 



void DistArrayViewInfo::resetBytes(IOH_int bytes_per_element)
{
 if (bytes_per_element==nbytes_per_element) return;
 set_bytes(bytes_per_element);
}



#endif


// ***************************************************************
  }
}
