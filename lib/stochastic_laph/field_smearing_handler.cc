//Arjun Singh Gambhir
//Real authors: Ben Hoerz, John Bulava, Colin Morningstar, ??
//Code directly ported from chroma_laph, more acknowledgements coming later.
#include "field_smearing_handler.h"
#include "xml_help.h"
#include <typeinfo>
using namespace std;

#ifdef TESTING
#include "tests.h"
#endif

namespace Chroma {
  namespace LaphEnv {


 // **************************************************************
 // *                                                            *
 // *                                                            *
 // *           GluonSmearingHandler implementation              *
 // *                                                            *
 // *                                                            *
 // **************************************************************

   // constructors

GluonSmearingHandler::GluonSmearingHandler() 
     : smearPtr(0), uPtr(0), dh_ptr(0) {}



GluonSmearingHandler::GluonSmearingHandler(const GluonSmearingInfo& gluon_smearing,
                                           const GaugeConfigurationInfo& gauge,
                                           const string& smearedFieldFileName)
{
 dh_ptr = 0;
 set_info(gluon_smearing,gauge,smearedFieldFileName);
}


void GluonSmearingHandler::setInfo(const GluonSmearingInfo& gluon_smearing,
                                   const GaugeConfigurationInfo& gauge,
                                   const string& smearedFieldFileName)
{
 clear();
 set_info(gluon_smearing,gauge,smearedFieldFileName);
}


void GluonSmearingHandler::set_info(const GluonSmearingInfo& gluon_smearing,
                                    const GaugeConfigurationInfo& gauge,
                                    const string& smearedFieldFileName)
{
 h_filename=tidyString(smearedFieldFileName);
 if (h_filename.empty()){
    QDPIO::cerr << "empty file name in GluonSmearingHandler"<<endl;
    QDP_abort(1);}
 try{
    uPtr=new GaugeConfigurationInfo(gauge);
    smearPtr=new GluonSmearingInfo(gluon_smearing);}
 catch(...){
    QDPIO::cerr << "problem in setInfo in GluonSmearingHandler"<<endl;
    QDP_abort(1);}

#if (QDP_ND == 3)
    // check that the file exists and contains the appropriate data

 QDPIO::cout << "Checking file: "<<h_filename<<endl;
 if (!fileExists(h_filename)){
      QDPIO::cerr << "file "<<h_filename<<" does not exist"<<endl;
      QDPIO::cerr << " reading the smeared gauge time slices"<<endl;
      QDPIO::cerr << "   will not be possible"<<endl;
      clear();
      QDP_abort(1);}

      // Open for reading (reads header info and checks)

 dh_ptr = new DataGetHandlerSF<GluonSmearingHandler,RecordKey,
                               multi1d<LatticeColorMatrix> >(*this,h_filename,
                               "Laph--SmearedGaugeField");

     // check that dimensions of 3-d lattice are the same as the
     // spatial extents of the 4-d lattice

 const multi1d<int>& extents4d = uPtr->getExtents();  // 4-d extents
 int tdir=uPtr->getTimeDir();
 multi1d<int> extents3d = QDP::Layout::lattSize();     // 3-d extents
 int dir3d=0;
 for (int dir4d=0;dir4d<=QDP::Nd;dir4d++){
    if (dir4d!=tdir){
       if (extents4d[dir4d]!=extents3d[dir3d]){
          QDPIO::cerr << "Dimensions of 3-d lattice do not match"
                      <<" the 3 spatial dimensions of the 4-d lattice"<<endl;
          QDPIO::cerr << "3-d lattice dimensions: "<<extents3d[0]<<" x "
                      << extents3d[1]<<" x "<<extents3d[2]<<endl;
          QDPIO::cerr << "4-d lattice dimensions: "<<extents4d[0]<<" x "
                      << extents4d[1]<<" x "<<extents4d[2]<<" x "
                      << extents4d[3]<<" time direction = "<<tdir<<endl;
          QDP_abort(1);}
       dir3d++;}}

 QDPIO::cout << "File "<<h_filename
             <<" successfully opened and header matches"<<endl<<endl;

#endif
}

    // destructor

GluonSmearingHandler::~GluonSmearingHandler()
{
 clear();
}


void GluonSmearingHandler::clear()
{
#if (QDP_ND == 3)
 clearData();
 if (dh_ptr) delete dh_ptr; 
 dh_ptr=0;
#endif

 try{
    delete smearPtr;
    delete uPtr;} 
 catch(...) {QDP_abort(1);}
 h_filename.clear();
 smearPtr=0;
 uPtr=0;
}



 // **********************************************************

bool GluonSmearingHandler::isInfoSet() const
{
 return ((smearPtr!=0)&&(uPtr!=0));
}

void GluonSmearingHandler::check_info_set(const string& name) const
{
 if (!isInfoSet()){
    QDPIO::cerr << "error in GluonSmearingHandler:"<<endl;
    QDPIO::cerr << "  must setInfo before calling "<<name<<endl;
    QDP_abort(1);}
}

const GluonSmearingInfo& GluonSmearingHandler::getGluonSmearingInfo() const
{
 check_info_set("getGluonSmearingInfo");
 return *smearPtr;
}

const GaugeConfigurationInfo& GluonSmearingHandler::getGaugeConfigurationInfo() const
{
 check_info_set("getGaugeConfigurationInfo");
 return *uPtr;
}

void GluonSmearingHandler::filefail(const string& message)
{
 QDPIO::cerr << message << endl;
 clear();
 QDP_abort(1);
}


 // **********************************************************

#if (QDP_ND == 4)

           // Compute the smeared gauge field and put into the file
           // whose name is stored in the GluonSmearingInfo.  This must 
           // be run in four dimensions, but the output file is meant to
           // be read in three dimensions.

void GluonSmearingHandler::computeSmearedGaugeField(const string& gauge_id)
{
 check_info_set("computeSmearedGaugeField");

    // check if output file exists; if so, we do not want to overwrite
    // so quit

 if (fileExists(h_filename)){
      QDPIO::cerr << "file "<<h_filename<<" already exists"<<endl;
      QDPIO::cerr << " will not overwrite; no output written"<<endl<<endl;
      return; }

 START_CODE();
 StopWatch outer;
 outer.start();

    // file does not exist, so it is okay to commence computation and output

 double rho = smearPtr->getLinkStapleWeight();
 int niters = smearPtr->getNumberOfLinkIterations();
 int tdir = uPtr->getTimeDir();

 if (tdir!=QDP::Nd-1){
    QDPIO::cerr << "This program does not support a time index that"
                <<" is not the most slowly varying"<<endl;
    QDP_abort(1);}

 int ndir = uPtr->getNumberOfDirections();

 QDPIO::cout << "Computation of smeared gauge field commencing in GluonSmearingHandler"<<endl;
 QDPIO::cout << smearPtr->output()<<endl;
 
 multi1d<bool> smear_dirs(ndir);
 for (int i=0;i<ndir;i++) smear_dirs[i] = true;
 smear_dirs[tdir]=false;

 multi2d<Real> rhomat(ndir,ndir);
 for(int mu=0; mu < ndir; mu++) 
 for(int nu=0; nu < ndir; nu++){
    if( mu != nu ) rhomat[mu][nu] = Real(rho);
    else rhomat[mu][nu] = 0;}
 for(int mu=0; mu < ndir; mu++)
    if ( ! smear_dirs[mu] ){
       for(int nu=0; nu < ndir; nu++){
          rhomat[mu][nu] = 0;
          rhomat[nu][mu] = 0;}}

 multi1d<LatticeColorMatrix> utemp(ndir), usmear(ndir);
 GaugeConfigurationHandler uHandler(*uPtr,gauge_id);

 if ((niters%2)==0) usmear=uHandler.getData();
 else{
    utemp=uHandler.getData();
    Chroma::Stouting::smear_links(utemp,usmear,smear_dirs,rhomat);}

 for (int k=0;k<(niters/2);k++){
    Chroma::Stouting::smear_links(usmear,utemp,smear_dirs,rhomat);
    Chroma::Stouting::smear_links(utemp,usmear,smear_dirs,rhomat);}

 QDPIO::cout << "  now starting write out to file as time slices"<<endl<<endl;
 StopWatch iotimer;
 iotimer.start();
 int textent = uPtr->getTimeExtent();

    // Prepare header xml data for the output file
 XmlBufferWriter xmlw;
 writeHeader(xmlw);
 string header(xmlw.str());

 IOMap<RecordKey,vector<TimeSliceOf<LatticeColorMatrix> > > dmobj;
 try { dmobj.openNew(h_filename,"Laph--SmearedGaugeField",header);}
 catch(...) { filefail("could not open for writing"); }

 vector<TimeSliceOf<LatticeColorMatrix> > slicer;
 slicer.reserve(3);
 for (unsigned int dir=0;dir<3;dir++){
    slicer.push_back(TimeSliceOf<LatticeColorMatrix>(usmear[dir]));}

 for (unsigned int t=0;t<uint(textent);t++){
    for (unsigned int dir=0;dir<3;dir++){
       slicer[dir].setCurrentTime(t);}
       try {dmobj.put(RecordKey(t),slicer);}
       catch(...){filefail("could not insert");}}

 iotimer.stop();

#ifdef TESTING
 printSmearedLinks(usmear,*uPtr,"smeargauge4d.log");
#endif

 outer.stop();
 END_CODE();

 QDPIO::cout << "computeSmearedGaugeField: total time = "
             << outer.getTimeInSeconds() 
             << " secs" << endl;
 QDPIO::cout << "Time to write out smeared gauge time slices = "
        << iotimer.getTimeInSeconds()<<" sec"<<endl;
 QDPIO::cout << "ran successfully in GluonSmearingHandler" <<endl<< endl;
}



#elif (QDP_ND == 3)

           // Set smeared gauge field from the disk map object if smearing
           // parameters same as current info (error otherwise).  

void GluonSmearingHandler::setSmearedGaugeFieldTimeSlice(int time_slice)
{
 check_info_set("setSmearedGaugeFieldTimeSlice");
 RecordKey rkey(time_slice);
 dh_ptr->getData(rkey);
}


           // Provides access to the smeared gauge field; calls
           // "set" if needed; aborts if cannot find entry

const multi1d<LatticeColorMatrix>& 
    GluonSmearingHandler::getSmearedGaugeFieldTimeSlice(int time_slice)
{
 check_info_set("getSmearedGaugeFieldTimeSlice");
 RecordKey rkey(time_slice);
 return dh_ptr->getData(rkey);
}



void GluonSmearingHandler::removeData(int time_slice)
{
 RecordKey rkey(time_slice);
 return dh_ptr->removeData(rkey);
}


void GluonSmearingHandler::clearData()
{
 dh_ptr->clearData();
}


#endif

void GluonSmearingHandler::writeHeader(XmlWriter& file_xml)
{
 push(file_xml, "SmearedGaugeFieldTimeSlices");
 smearPtr->output(file_xml);
 uPtr->output(file_xml);
 pop(file_xml);
}

bool GluonSmearingHandler::checkHeader(XmlReader& xmlr)
{
 GaugeConfigurationInfo gauge_check(xmlr);
 GluonSmearingInfo gsmear_check(xmlr);
 try { uPtr->checkEqual(gauge_check); 
       smearPtr->checkEqual(gsmear_check);}
 catch(...) { return false;}
 return true;
}




 // **************************************************************
 // *                                                            *
 // *                                                            *
 // *           QuarkSmearingHandler implementation              *
 // *                                                            *
 // *                                                            *
 // **************************************************************


   // constructors

QuarkSmearingHandler::QuarkSmearingHandler()
     : uPtr(0), gSmearPtr(0), qSmearPtr(0),
       m_read_mode(true), dh_ptr(0) {}


QuarkSmearingHandler::QuarkSmearingHandler(const GluonSmearingInfo& gluon_smearing,
                                           const GaugeConfigurationInfo& gauge,
                                           const QuarkSmearingInfo& quark_smearing,
                                           const string smeared_quark_file_stub,
                                           bool read_mode)
{
 set_info(gluon_smearing,gauge,quark_smearing,smeared_quark_file_stub,read_mode);
}

void QuarkSmearingHandler::setInfo(const GluonSmearingInfo& gluon_smearing,
                                   const GaugeConfigurationInfo& gauge,
                                   const QuarkSmearingInfo& quark_smearing,
                                   const string smeared_quark_file_stub,
                                   bool read_mode)
{
 clear();
 set_info(gluon_smearing,gauge,quark_smearing,smeared_quark_file_stub,read_mode);
}



void QuarkSmearingHandler::set_info(const GluonSmearingInfo& gluon_smearing,
                                    const GaugeConfigurationInfo& gauge,
                                    const QuarkSmearingInfo& quark_smearing,
                                    const string smeared_quark_file_stub,
                                    bool read_mode)
{
 smearedQuarkFileStub=smeared_quark_file_stub;
 if (smearedQuarkFileStub.empty()){
    QDPIO::cerr << "empty file name or stub in QuarkSmearingHandler"<<endl;
    QDP_abort(1);}

 try{
    uPtr=new GaugeConfigurationInfo(gauge);
    gSmearPtr=new GluonSmearingInfo(gluon_smearing);
    qSmearPtr=new QuarkSmearingInfo(quark_smearing);}
 catch(...){
    QDPIO::cerr << "problem in setInfo in QuarkSmearingHandler"<<endl;
    QDP_abort(1);}

 m_read_mode=read_mode;

 m_memory_only=false;
 if (smeared_quark_file_stub=="NamedObjectMap")
   m_memory_only=true;

 if (m_read_mode){
#if (QDP_ND == 3)
    FileListInfo files(smearedQuarkFileStub+"_time",
                       0,uPtr->getTimeExtent()-1,false);
    try{ 

      if (m_memory_only)
	dh_ptr=new DataGetHandlerMemMF<QuarkSmearingHandler,TimeKey,LevelKey,
                    LatticeColorVector>(*this,"Laph--SmearedQuarkTimeFile");
      else
	dh_ptr=new DataGetHandlerMF<QuarkSmearingHandler,TimeKey,LevelKey,
                    LatticeColorVector>(*this,files,"Laph--SmearedQuarkTimeFile",
                    "LaphEigenvectors");

       if (dh_ptr->getFileKeys().size()!=0)
          QDPIO::cout << "QuarkSmearingHandler file stub appears to be okay"<<endl;
       else{
          QDPIO::cout << "problem with QuarkSmearingHandler file stub: "<<smeared_quark_file_stub<<endl;
          QDP_abort(1);}}
#elif (QDP_ND == 4)
    FileListInfo files(smearedQuarkFileStub+"_level",
                       0,qSmearPtr->getNumberOfLaplacianEigenvectors()-1,false);
    try{ dh_ptr=new DataGetHandlerMF<QuarkSmearingHandler,LevelKey,LevelKey,
                        LatticeColorVector>(*this,files,"Laph--SmearedQuarkLevelFile",
                        "LaphEigenvectors");
       if (dh_ptr->queryFile(LevelKey(0)))
          QDPIO::cout << "QuarkSmearingHandler file stub appears to be okay"<<endl;
       else{
          QDPIO::cout << "problem with QuarkSmearingHandler file stub: "<<smeared_quark_file_stub<<endl;
          QDP_abort(1);}}
#endif
    catch(...){ 
      QDPIO::cerr << "unable to allocate data handler in QuarkSmearingHandler"<<endl;
      QDP_abort(1);}}
 else
    dh_ptr=0;
}

    // use this to increase the number of Laplacian eigenvectors allowed;
    // need inside QuarkHandler since smearing sub-handler is static.

void QuarkSmearingHandler::updateSmearing(const QuarkSmearingInfo& quark_smearing)
{
 check_info_set("updateQuarkSmearingInfo");
 qSmearPtr->increaseUpdate(quark_smearing);
}


QuarkSmearingHandler::~QuarkSmearingHandler()
{
 clear();
}


void QuarkSmearingHandler::clear()
{
 try{
    delete gSmearPtr;
    delete qSmearPtr;
    delete uPtr;
    delete dh_ptr;} 
 catch(...) {QDP_abort(1);}
 smearedQuarkFileStub.clear();
 gSmearPtr=0;
 qSmearPtr=0;
 uPtr=0;
 dh_ptr=0;
#if (QDP_ND == 3)
 Eigenvalues.resize(0);
#elif (QDP_ND == 4)
 Eigenvalues.resize(0,0);
#endif
}


           // access to the info

bool QuarkSmearingHandler::isInfoSet() const
{
 return ((gSmearPtr!=0)&&(uPtr!=0)&&(qSmearPtr!=0)
        &&(!smearedQuarkFileStub.empty()));
} 


   // check_mode = 0 means no check, 1 means check for read mode,
   // 2 means check for write mode

void QuarkSmearingHandler::check_info_set(const string& name,
                                          int check_mode) const
{
 if (!isInfoSet()){
    QDPIO::cerr << "error in QuarkSmearingHandler:"<<endl;
    QDPIO::cerr << "  must setInfo before calling "<<name<<endl;
    QDP_abort(1);}
 if (check_mode==1){
    if (!m_read_mode){
       QDPIO::cerr << "error in QuarkSmearingHandler:"<<endl;
       QDPIO::cerr << "  must be in read mode when calling "<<name<<endl;
       QDP_abort(1);}}
 else if (check_mode==2){
    if (m_read_mode){
       QDPIO::cerr << "error in QuarkSmearingHandler:"<<endl;
       QDPIO::cerr << "  must not be in read mode when calling "<<name<<endl;
       QDP_abort(1);}}
}

const GluonSmearingInfo& QuarkSmearingHandler::getGluonSmearingInfo() const
{
 check_info_set("getGluonSmearingInfo");
 return *gSmearPtr;
}

const QuarkSmearingInfo& QuarkSmearingHandler::getQuarkSmearingInfo() const
{
 check_info_set("getQuarkSmearingInfo");
 return *qSmearPtr;
}

const GaugeConfigurationInfo& QuarkSmearingHandler::getGaugeConfigurationInfo() const
{
 check_info_set("getGaugeConfigurationInfo");
 return *uPtr;
}

const string& QuarkSmearingHandler::getSmearedQuarkFieldFileStub() const
{
 check_info_set("getSmearedQuarkFieldFileStub");
 return smearedQuarkFileStub;
}

void QuarkSmearingHandler::getFileMap(XmlWriter& xmlout) const
{
 check_info_set("getFileMap",1);
 dh_ptr->getFileMap(xmlout);
}

void QuarkSmearingHandler::outputKeys(XmlWriter& xmlout)
{
 check_info_set("outputKeys",1);
 dh_ptr->outputKeys(xmlout);
}


void QuarkSmearingHandler::failure(const string& message)
{
 QDPIO::cerr << message << " in QuarkSmearingHandler"<<endl;
 QDP_abort(1);
}


 // **********************************************************

#if (QDP_ND == 3)

           // Compute the Laph Eigenvectors and write to file.
           // Will not overwrite.
           // User must do a "clear" first.

void QuarkSmearingHandler::computeLaphEigenvectors(
                         const LaphEigenSolverInfo& solver_info,
                         const string& smeared_gauge_file,
			 set<int>& timeslices,
                         int striping_factor, int striping_unit)
{
 check_info_set("computeLaphEigenvectors",2);

 START_CODE();
 StopWatch rolex,bulava;
 rolex.start();
 double iotime=0.0;

    // check that we won't overwrite any files; abort if so

 int nTime = uPtr->getTimeExtent();
 int nEigvecs = qSmearPtr->getNumberOfLaplacianEigenvectors();

    // fire up the GluonSmearingHandler

 GluonSmearingHandler gHandler(*gSmearPtr,*uPtr,smeared_gauge_file);

 QDPIO::cout << "Computation of Laplacian eigenvectors commencing"
             << " in FieldSmearingHandler"<<endl;
 QDPIO::cout << "  number of requested eigenvectors = "<<nEigvecs<<endl;
 QDPIO::cout << "            time extent of lattice = "<<nTime<<endl;

#ifdef TESTING
 TextFileWriter fout("lapheigvec3d.log");
#endif

    // fire up the data handler

 FileListInfo writeFiles(smearedQuarkFileStub+"_time",0,nTime-1,false);
 DataPutHandlerBaseMF<QuarkSmearingHandler,TimeKey,LevelKey,LatticeColorVector> *DHput;

 if (m_memory_only)
   DHput = new DataPutHandlerMemMF<QuarkSmearingHandler,TimeKey,LevelKey,
       LatticeColorVector>(*this,"Laph--SmearedQuarkTimeFile");
 else
   DHput = new DataPutHandlerMF<QuarkSmearingHandler,TimeKey,LevelKey,
       LatticeColorVector>(*this,writeFiles,"Laph--SmearedQuarkTimeFile",
                                 "LaphEigenvectors",true,false,
                                 striping_factor,striping_unit);
 
     // do the computation for each time slice; put results into
     // files with the stub name given and suffix  "_time.t" where
     // t = 0,1,...  

     // check what time slices are to be done
 if (timeslices.empty()) {
   for (int t=0; t<nTime; t++){

      ostringstream oss;
      oss << smearedQuarkFileStub << "_time." << t;
      if (fileExists(oss.str())){
	 QDPIO::cout << "file "<<oss.str()<<" already exists"<<endl;
	 QDPIO::cout << " will not overwrite"<<endl<<endl;
	 continue;}

      timeslices.insert(t);
   }
 }

 for (set<int>::const_iterator t=timeslices.begin(); t!=timeslices.end(); ++t){

    QDPIO::cout <<endl<<endl<<"Beginning computation for time = "
                <<*t<<" (with "<<*timeslices.rbegin()<<" as last)"<<endl<<endl;
    StopWatch iotimer; iotimer.start();
    const multi1d<LatticeColorMatrix>& usmear
                        = gHandler.getSmearedGaugeFieldTimeSlice(*t);
    iotimer.stop();
    iotime+=iotimer.getTimeInSeconds();

         // use the Krylov-spectral restarted Lanczos method to compute
         // the eigenvalues/eigenvectors of the Laplacian on this time slice

    LaphVectorHandler laphEigvecs(usmear);
    int chebyshevOrder = solver_info.getChebyshevOrder();
    if (chebyshevOrder>=2){
       laphEigvecs.setChebyshevAccelerationOn(
                      solver_info.getMaximumEigenvalue(),
                      solver_info.getCutoffEigenvalue(),chebyshevOrder);}

    bulava.reset();bulava.start();
    KSRLanczosEigenSolver Comp(laphEigvecs, nEigvecs,
                  solver_info.getResidualTolerance(),
                  solver_info.getKrylovDimension(),
                  solver_info.getMaximumIterations(),
                  laphEigvecs.getSpectrumEnd(),
                  solver_info.getOutputVerbosity());

    if (!Comp.wasSuccessful()){
       QDPIO::cout << endl<<endl<<"Convergence to requested tolerance"
                   << " was NOT achieved...no output"<<endl;
       delete DHput;
       QDP_abort(1);}
    bulava.stop();
    Eigenvalues.resize(nEigvecs);

     //  multiply each eigenvector by a
     //  phase so that the zero-th color component at site (0,0,0)
     //  is real and positive; this uniquely specifies the overall
     //  phase for each eigenvector, which is important since a
     //  change in these phases can change the effective Laph noise.
     //  Using this convention makes the results robust against
     //  eigenvector deletion and subsequent reconstruction.

    laphEigvecs.applyPhaseConvention();

    QDPIO::cout << endl<<endl<<" SUCCESS: time = "<<bulava.getTimeInSeconds() 
                << " secs" <<endl;
    if (chebyshevOrder>=2)
       QDPIO::cout << "  Chebyshev acceleration was used"<<endl;
    QDPIO::cout << "  Norm of matrix diagonalized was "
                   << Comp.getMatrixNormEstimate() << endl; 
    for (int i=0; i<nEigvecs; i++){  
       QDPIO::cout << "  Raw eigenvalue("<<i<<") = "
                   << Comp.getEigenvalue(i);             
       QDPIO::cout << "  Raw Residual("<<i<<") = "           
                   << Comp.getResidual(i) << endl;
       Eigenvalues[i]=Comp.getEigenvalue(i);}

         // perform a crucial check of the diagonalization
        // if Chebyshev acceleration was used, we need to get the
        // eigenvalues of -Delta

    double offmaxmag=0.0;
    laphEigvecs.resize(nEigvecs+1);
    laphEigvecs.setChebyshevAccelerationOff();
    for (int level=0;level<nEigvecs;level++){
       laphEigvecs.assignMatrixVectorMultiply(nEigvecs,level);
       double rr=sqrt(toDouble(
               laphEigvecs.InnerProductRealPart(nEigvecs,nEigvecs)));
       QDPIO::cout << "magnitude of diagonal element -Laplacian["
                   <<level<<","<<level<<"] = "<<rr<<endl;
       if (chebyshevOrder>=2) Eigenvalues[level]=rr;
       for (int row=0;row<nEigvecs;row++){
          if (row!=level){
             DComplex z=laphEigvecs.InnerProduct(row,nEigvecs);
             double rr=toDouble(sqrt(QDP::real(z)*QDP::real(z)+QDP::imag(z)*QDP::imag(z)));
             if (rr>offmaxmag) offmaxmag=rr;}}}
    QDPIO::cout << "Maximum magnitude of off-diagonal matrix elements"<<endl
                << "       of -Laplacian = " << offmaxmag <<endl<<endl;

#ifdef TESTING
    printLaphEigenvectors(laphEigvecs,*t,*uPtr,nEigvecs,fout);
#endif

    iotimer.reset();iotimer.start(); 
    bulava.reset();bulava.start();
    DHput->open(TimeKey(*t)); // open file, write header
    for (int k=0;k<nEigvecs;k++)
       DHput->putData(LevelKey(k),laphEigvecs.getVector(k));
    DHput->flush();  // finalize current file

    gHandler.clearData();  // clear smeared gauge time slices
    iotimer.stop(); bulava.stop();
    QDPIO::cout << "Time to write = "<<bulava.getTimeInSeconds()<<endl;
    iotime+=iotimer.getTimeInSeconds();
    QDPIO::cout << "Computation done for time = "<<*t<<endl<<endl;

    }
  
#ifdef TESTING
 fout.close();
#endif

 delete DHput;

 rolex.stop();
 QDPIO::cout << endl<<"computeLaphEigenvectors: total time = "
             << rolex.getTimeInSeconds() 
             << " secs" << endl;
 QDPIO::cout << "Total file I/O time = "<< iotime<<" sec"<<endl;
 QDPIO::cout << "ran successfully" << endl<<endl;

}


// *************************************************************


double QuarkSmearingHandler::estimateLargestLaplacianEigenvalue(
                                   const string& smeared_gauge_file)
{
 check_info_set("estimateLargestLaplacianEigenvalue",2);

 double lambda_max=0.0;
 int nTime = uPtr->getTimeExtent();
 int nEigvecs = 1;

    // fire up the GluonSmearingHandler

 GluonSmearingHandler gHandler(*gSmearPtr,*uPtr,smeared_gauge_file);

 QDPIO::cout << "Estimating largest eigenvalue of -Laplacian commencing"
             << " in FieldSmearingHandler"<<endl;
 QDPIO::cout << "            time extent of lattice = "<<nTime<<endl;

 START_CODE();
 StopWatch rolex;
 rolex.start();

     // do the computation for each time slice; put results into
     // files with the stub name given and suffix  ".time_t" where
     // t = 0,1,...  

 for (int t=0; t<nTime; t++){

    QDPIO::cout <<endl<<endl<<"Beginning computation for time = "<<t<<endl<<endl;
    const multi1d<LatticeColorMatrix>& usmear
                        = gHandler.getSmearedGaugeFieldTimeSlice(t);

         // use the Krylov-spectral restarted Lanczos method to compute
         // the eigenvalues/eigenvectors of the Laplacian on this time slice

    LaphVectorHandler laphEigvecs(usmear);
    laphEigvecs.setInitialToRandom();

    KSRLanczosEigenSolver Comp(laphEigvecs, nEigvecs, 1e-3, 24, 24, 'U');

    if (!Comp.wasSuccessful()){
       QDPIO::cout << endl<<endl<<"Convergence to requested tolerance"
                   << " was NOT achieved...no output"<<endl;}

    double this_max = Comp.getEigenvalue(0);
    QDPIO::cout << "largest eigenvalue on this time slice = "
                << this_max<<endl;

    lambda_max=max(lambda_max,this_max);
    
    }
  

 rolex.stop();
 QDPIO::cout << endl<<"computeLaphEigenvectors: total time = "
             << rolex.getTimeInSeconds() 
             << " secs" << endl;
 QDPIO::cout << "ran successfully" << endl;

 return lambda_max;
}

     // Provides access to the quark smearing eigenvectors
     // while in read mode.

const LatticeColorVector& QuarkSmearingHandler::getLaphEigenvector(
                                   int timeslice, int eigpair_num)
{
 check_info_set("getLaphEigenvector",1);
 return dh_ptr->getData(TimeKey(timeslice),LevelKey(eigpair_num));
}

bool QuarkSmearingHandler::queryLaphEigenvector(
                                   int timeslice, int eigpair_num)
{
 check_info_set("queryLaphEigenvector",1);
 return dh_ptr->queryData(TimeKey(timeslice),LevelKey(eigpair_num));
}

void QuarkSmearingHandler::removeLaphEigenvector(
                                   int timeslice, int eigpair_num)
{
 if (!m_read_mode) return;
 return dh_ptr->removeData(TimeKey(timeslice),LevelKey(eigpair_num));
}

void QuarkSmearingHandler::clearLaphEigenvectors()
{
 if (!m_read_mode || m_memory_only) return;
 dh_ptr->clearData();
 Eigenvalues.resize(0);
}

void QuarkSmearingHandler::readEigenvectorsIntoMemory(int time)
{
 if (!m_read_mode || m_memory_only) return;
 
 TimeKey tKey(time);
 set<LevelKey> recKeys = dh_ptr->getKeys(tKey);

 DataPutHandlerMemMF<QuarkSmearingHandler,TimeKey,LevelKey,
     LatticeColorVector> DHput(*this,"Laph--SmearedQuarkTimeFile");
 DHput.open(tKey);

 for (set<LevelKey>::const_iterator lIt = recKeys.begin();
	lIt != recKeys.end(); ++lIt) {
   DHput.putData(*lIt, dh_ptr->getData(tKey, *lIt));
   dh_ptr->removeData(tKey, *lIt);
 }
}


#elif (QDP_ND == 4)

           // The 3-d computation of the Laph eigenvectors produces
           // one file for each time slice, and each file contains all
           // eigen-levels up to the requested number of eigenvectors.
           // This routine re-organizes these files, producing one file
           // for each level, but each file contains all time slices.

void QuarkSmearingHandler::combineTimeSlices(int striping_factor, int striping_unit)
{
 check_info_set("combineTimeSlices",2);

 double rtimer=0.0,wtimer=0.0;
 START_CODE();
 StopWatch rolex,seiko;
 rolex.start();
 seiko.start();
 
 int nTime = uPtr->getTimeExtent();
 int nEigvecs = qSmearPtr->getNumberOfLaplacianEigenvectors();

    // check that all "time" files are available
 for (int t=0;t<nTime;t++){
    ostringstream oss;
    oss << smearedQuarkFileStub << "_time." << t;
    if (!fileExists(oss.str())){
         QDPIO::cerr << "file "<<oss.str()<<" does not exist"<<endl;
         QDPIO::cerr << " combineTimeSlices fails...."<<endl;
         QDP_abort(1);}}

    // determine the smallest level to start with (based on existing files)
 int levelstart;
 for (levelstart=nEigvecs-1;levelstart>=0;levelstart--){
    ostringstream oss;
    oss << smearedQuarkFileStub << "_level." << levelstart;
    if (fileExists(oss.str())) break;}
 levelstart++;
 if (levelstart==nEigvecs) return;

   // create the data handlers (the "get" is only used for checking
   // the headers and extracting the eigenvalues)

 FileListInfo readFiles(smearedQuarkFileStub+"_time",0,nTime-1);
 FileListInfo writeFiles(smearedQuarkFileStub+"_level",levelstart,nEigvecs-1);
 Eigenvalues.resize(nTime,nEigvecs);  // storage for eigenvalues

 {DataGetHandlerMF<QuarkSmearingHandler,TimeKey,LevelKey,
       LatticeColorVector> DHget(*this,readFiles,"Laph--SmearedQuarkTimeFile",
                                 "LaphEigenvectors");}
 seiko.stop();
 rtimer+=seiko.getTimeInSeconds();
 QDPIO::cout << "Input check done; time="<<rtimer<<" seconds"<<endl;
 seiko.reset();seiko.start();
 DataPutHandlerMF<QuarkSmearingHandler,LevelKey,LevelKey,
       LatticeColorVector> DHput(*this,writeFiles,"Laph--SmearedQuarkLevelFile",
                                 "LaphEigenvectors",true,false,
                                 striping_factor,striping_unit);
 seiko.stop();
 wtimer+=seiko.getTimeInSeconds();
 QDPIO::cout << "Put handler set up; time="<<wtimer<<" seconds"<<endl;
 seiko.reset();seiko.start();


    // use a "Handle" in case of exception
 multi1d< Handle<IOMap<LevelKey,
          TimeSliceOf<LatticeColorVector> > > > rdm_ptr(nTime);
 for (int t=0;t<nTime;t++){
    ostringstream oss;
    oss << smearedQuarkFileStub << "_time." << t;
    rdm_ptr[t]=new IOMap<LevelKey,TimeSliceOf<LatticeColorVector> >;
    try { rdm_ptr[t]->openReadOnly(oss.str(),"Laph--SmearedQuarkTimeFile");}
    catch(...) { QDPIO::cerr << "read failure in combineTimeSlices"; 
                 QDP_abort(1);}}
 seiko.stop();
 rtimer+=seiko.getTimeInSeconds();
 QDPIO::cout << "Ready for read; time="<<rtimer<<" seconds"<<endl;

   //  now collect time slices of each level and output to file

 LatticeColorVector laph_eigvecs;
 TimeSliceOf<LatticeColorVector> buffer(laph_eigvecs);

   // create each level file

 for (int level=levelstart;level<nEigvecs;level++){
    double rtime=0.0, wtime=0.0;
    seiko.reset();seiko.start();
    QDPIO::cout << "creating file for level "<<level<<endl;
    LevelKey key(level);
    DHput.open(key);  // creates files, writes header
    seiko.stop();
    wtime+=seiko.getTimeInSeconds();
    seiko.reset();seiko.start();
    for (int t=0;t<nTime;t++){
       buffer.setCurrentTime(t);
       try{ rdm_ptr[t]->get(key,buffer); }
       catch(...){ QDPIO::cerr << "Lookup error in combineTimeSlices"<<endl;
                   QDP_abort(1);}}
    seiko.stop();
    rtime+=seiko.getTimeInSeconds();
    seiko.reset();seiko.start();
    DHput.putData(key,laph_eigvecs);
    DHput.flush();  // finalize current file
    seiko.stop();
    wtime+=seiko.getTimeInSeconds();
    QDPIO::cout << " level "<<level<<"  total read time = "<< rtime  << " secs" << endl;
    QDPIO::cout << " level "<<level<<" total write time = "<< wtime  << " secs" << endl;
    rtimer+=rtime;
    wtimer+=wtime;
    }

 rolex.stop();
 QDPIO::cout << endl<<"combineTimeSlices: total time = "<< rolex.getTimeInSeconds()  << " secs" << endl;
 QDPIO::cout <<"  total read time = "<< rtimer  << " secs" << endl;
 QDPIO::cout <<" total write time = "<< wtimer  << " secs" << endl;
 QDPIO::cout << "ran successfully" << endl;

}



     // Provides access to the quark smearing eigenvectors
     // while in read mode.

const LatticeColorVector& QuarkSmearingHandler::getLaphEigenvector(
                                   int eigpair_num)
{
 check_info_set("getLaphEigenvector",1);
 LevelKey key(eigpair_num);
 return dh_ptr->getData(key,key);
}

bool QuarkSmearingHandler::queryLaphEigenvector(
                                   int eigpair_num)
{
 check_info_set("queryLaphEigenvector",1);
 LevelKey key(eigpair_num);
 return dh_ptr->queryData(key,key);
}

void QuarkSmearingHandler::removeLaphEigenvector(
                                   int eigpair_num)
{
 if (!m_read_mode) return;
 LevelKey key(eigpair_num);
 return dh_ptr->removeData(key,key);
}

void QuarkSmearingHandler::clearLaphEigenvectors()
{
 if (!m_read_mode) return;
 dh_ptr->clearData();
 Eigenvalues.resize(0,0);
}

bool QuarkSmearingHandler::checkAllLevelFilesExist()
{
 check_info_set("checkAllLevelFilesExist");

 for (int level=0;level< qSmearPtr->getNumberOfLaplacianEigenvectors();level++){
    ostringstream oss;
    oss << smearedQuarkFileStub << "_level." << level;
    if (!fileExists(oss.str())){
       QDPIO::cerr << "needed file "<<oss.str()<<" does not exists"<<endl;
       return false; }}
 return true;
}

#endif


#if (QDP_ND == 3)

void QuarkSmearingHandler::writeHeader(XmlWriter& xmlw, const TimeKey& fkey, 
                                       int suffix)
{
 if (!m_memory_only && suffix!=fkey.value){
    QDPIO::cerr << "time suffix does not match file key"<<endl;
    QDP_abort(1);}
 push(xmlw, "LaphEigenvectors");
 gSmearPtr->output(xmlw);
 uPtr->output(xmlw);
 qSmearPtr->output(xmlw);
 fkey.output(xmlw);
 write(xmlw,"Eigenvalues",Eigenvalues);
 pop(xmlw);
}

bool QuarkSmearingHandler::checkHeader(XmlReader& xml_in, int suffix)
{
 if (xml_tag_count(xml_in,"LaphEigenvectors")!=1) return false;
 XmlReader xmlr(xml_in,"./descendant-or-self::LaphEigenvectors");
 try {
    GaugeConfigurationInfo gauge_check(xmlr);
    uPtr->checkEqual(gauge_check);
    GluonSmearingInfo gsmear_check(xmlr);
    gSmearPtr->checkEqual(gsmear_check);
    QuarkSmearingInfo qsmear_check(xmlr);
    if (m_read_mode){
       qSmearPtr->checkOK(qsmear_check);}
    else{
       qSmearPtr->checkEqual(qsmear_check);}
    TimeKey time_check(xmlr);
    if (!m_memory_only && time_check.value!=suffix){throw("error");}
    if (!m_read_mode) read(xmlr,"Eigenvalues",Eigenvalues);}
 catch(...) { return false;}
 return true;
}

#elif (QDP_ND == 4)

void QuarkSmearingHandler::writeHeader(XmlWriter& xmlw, const LevelKey& fkey, 
                                       int suffix)
{
 if (suffix!=fkey.value){
    QDPIO::cerr << "level suffix does not match file key"<<endl;
    QDP_abort(1);}
 push(xmlw, "LaphEigenvectors");
 gSmearPtr->output(xmlw);
 uPtr->output(xmlw);
 qSmearPtr->output(xmlw);
 fkey.output(xmlw);
 int nTime=uPtr->getTimeExtent();
 multi1d<double> evals(nTime);
 for (int t=0;t<nTime;t++) evals[t]=Eigenvalues(t,suffix);
 write(xmlw,"Eigenvalues",evals);
 pop(xmlw);
}

bool QuarkSmearingHandler::checkHeader(XmlReader& xml_in, int suffix)
{
 if (xml_tag_count(xml_in,"LaphEigenvectors")!=1) return false;
 XmlReader xmlr(xml_in,"./descendant-or-self::LaphEigenvectors");
 try {
    GaugeConfigurationInfo gauge_check(xmlr);
    uPtr->checkEqual(gauge_check);
    GluonSmearingInfo gsmear_check(xmlr);
    gSmearPtr->checkEqual(gsmear_check);
    QuarkSmearingInfo qsmear_check(xmlr);
    if (m_read_mode){
       qSmearPtr->checkOK(qsmear_check);
       LevelKey level_check(xmlr);
       if (level_check.value!=suffix){throw("error");}}
    else{
       qSmearPtr->checkEqual(qsmear_check);
       TimeKey time_check(xmlr);
       if (time_check.value!=suffix){throw("error");}
       multi1d<double> eigvals;
       read(xmlr,"Eigenvalues",eigvals);
       for (int level=0;level<eigvals.size();level++)
          Eigenvalues(suffix,level)=eigvals[level];}}
 catch(...) { return false;}
 return true;
}

#endif


#if (QDP_ND == 3)

 // **************************************************************
 // *                                                            *
 // *                                                            *
 // *              LaphVectorHandler implementation              *
 // *                                                            *
 // *                                                            *
 // **************************************************************


LaphVectorHandler::LaphVectorHandler(
            const multi1d<LatticeColorMatrix>& smeared_gauge_timeslice)
         : usmear_timeslice(smeared_gauge_timeslice), spec_end('L'), 
           start_vector_type('E')
{
 use_acceleration=false;
 matrixPtr=&LaphVectorHandler::apply_minus_laplacian;
}

LaphVectorHandler::LaphVectorHandler(
           const multi1d<LatticeColorMatrix>& smeared_gauge_timeslice,
           int nvec)   : usmear_timeslice(smeared_gauge_timeslice), 
                         spec_end('L'), start_vector_type('E')
{
 use_acceleration=false;
 matrixPtr=&LaphVectorHandler::apply_minus_laplacian;
 resize(nvec);
}


void LaphVectorHandler::setChebyshevAccelerationOn(double largest_eigenvalue,
                                                   double cutoff_eigenvalue,
                                                   int chebyshev_polynomial_order)
{
 if ((chebyshev_polynomial_order < 2)
    ||(cutoff_eigenvalue>=(largest_eigenvalue-0.2))
    ||(largest_eigenvalue<3.0)){
    QDPIO::cerr << "invalid parameters to setChebyshevAccelerationOn"<<endl;
    throw(string("invalid"));}

 chebyshev_order = chebyshev_polynomial_order; 
 ov_coef = 2.0/(largest_eigenvalue-cutoff_eigenvalue);
 id_coef = 0.5*(largest_eigenvalue+cutoff_eigenvalue-12.0);
 two_ov_coef=2.0*ov_coef;
 spec_end = 'U';

 use_acceleration=true;
 matrixPtr=&LaphVectorHandler::apply_transf_chebyshev;
 vbuf.resize(3);
}

void LaphVectorHandler::setChebyshevAccelerationOff()
{
 use_acceleration=false;
 matrixPtr=&LaphVectorHandler::apply_minus_laplacian;
 vbuf.clear();
 spec_end='L';
}

void LaphVectorHandler::setInitialToRandom()
{
 start_vector_type='R';
}

void LaphVectorHandler::setInitialToEqualComponents()
{
 start_vector_type='E';
}



void LaphVectorHandler::clear()
{
 for (vector<LatticeColorVector* >::iterator 
        it=m_vecs.begin();it!=m_vecs.end();it++)
    delete *it;
 m_vecs.clear();
}

void LaphVectorHandler::resize(int ndim)
{
 if (ndim<=0) clear();
 if (ndim==int(m_vecs.size())) return;
 if (ndim>int(m_vecs.size())){
    int oldsize=m_vecs.size();
    m_vecs.resize(ndim);
    for (int k=oldsize;k<int(ndim);k++) m_vecs[k]=new LatticeColorVector;
    return;}
 for (unsigned int k=ndim;k<m_vecs.size();k++) delete m_vecs[k];
 m_vecs.resize(ndim);
}

void LaphVectorHandler::check_index(int index) const
{
 if ((index<0)||(index>=int(m_vecs.size()))){
    QDPIO::cerr << "invalid index in LaphVectorHandler"<<endl;
    throw(string("error"));}
}

void LaphVectorHandler::assignInitialUnitNorm(int index)
{
 if (start_vector_type=='R'){         // random initial vector
    assignRandomUnitNorm(index);
    return;}

 check_index(index);
 LatticeColorVector& v= *m_vecs[index];
 ColorVector cv;
 Complex z=1;                        // equal-component initial vector
 for (int color=0;color<QDP::Nc;color++)
   pokeColor(cv,z,color);
 v=cv;
 double vnorm;
 unitNormalize(index,vnorm);
}

void LaphVectorHandler::assignRandomUnitNorm(int index)
{
 check_index(index);
 LatticeColorVector& v= *m_vecs[index];
 gaussian(v);
 double vnorm;
 unitNormalize(index,vnorm);
}


DComplex LaphVectorHandler::InnerProduct(int lhs_index, int rhs_index)
{
 check_index(rhs_index);
 check_index(lhs_index);
 LatticeColorVector& b= *m_vecs[lhs_index];
 LatticeColorVector& a= *m_vecs[rhs_index];
 return innerProduct(b,a);
}
 
double LaphVectorHandler::InnerProductRealPart(int lhs_index, int rhs_index)
{
 check_index(rhs_index);
 check_index(lhs_index);
 LatticeColorVector& b= *m_vecs[lhs_index];
 LatticeColorVector& a= *m_vecs[rhs_index];
 return toDouble(innerProductReal(b,a));
}


void LaphVectorHandler::addTo(int lhs_index, const DComplex& cf, 
                              int rhs_index)
{
 check_index(rhs_index);
 check_index(lhs_index);
 LatticeColorVector& b= *m_vecs[lhs_index];
 LatticeColorVector& a= *m_vecs[rhs_index];
 b+=cf*a;
} 

void LaphVectorHandler::addTo(int lhs_index, double cf, int rhs_index)
{
 check_index(rhs_index);
 check_index(lhs_index);
 LatticeColorVector& b= *m_vecs[lhs_index];
 LatticeColorVector& a= *m_vecs[rhs_index];
 b+=cf*a;
} 

void LaphVectorHandler::copyTo(int lhs_index, int rhs_index)
{
 check_index(rhs_index);
 check_index(lhs_index);
 LatticeColorVector& b= *m_vecs[lhs_index];
 LatticeColorVector& a= *m_vecs[rhs_index];
 b=a;
} 

void LaphVectorHandler::copyTo(int lhs_index, double coef, int rhs_index)
{
 check_index(rhs_index);
 check_index(lhs_index);
 LatticeColorVector& b= *m_vecs[lhs_index];
 LatticeColorVector& a= *m_vecs[rhs_index];
 b=coef*a;
} 

void LaphVectorHandler::swap(int lhs_index, int rhs_index)
{
 check_index(rhs_index);
 check_index(lhs_index);
 LatticeColorVector* tmp=m_vecs[lhs_index];
 m_vecs[lhs_index]=m_vecs[rhs_index];
 m_vecs[rhs_index]=tmp;
} 

void LaphVectorHandler::unitNormalize(int index, double& vnorm)
{
 check_index(index);
 LatticeColorVector& v= *m_vecs[index];
 vnorm=sqrt(toDouble(norm2(v)));
 double rnorm=1.0/vnorm;
 v*=rnorm;
}
/*
void LaphVectorHandler::applyPhaseConvention()
{
 multi1d<int> origin(3);
 origin[0]=origin[1]=origin[2]=0;
 for (int index=0;index<m_vecs.size();++index){
    ColorVector s=peekSite(*m_vecs[index],origin);
    Complex z0=conj(peekColor(s,0));  double r0=toDouble(norm2(z0));
    Complex z1=conj(peekColor(s,1));  double r1=toDouble(norm2(z1));
    if (r1>r0){ z0=z1; r0=r1;}
    Complex z2=conj(peekColor(s,2));  double r2=toDouble(norm2(z2));
    if (r2>r0){ z0=z2; r0=r2;}
    r0=sqrt(r0);
    if (r0<1e-12){
       QDPIO::cerr << "problem applying phase convention"<<endl;
       QDP_abort(1);}
    z0/=r0;
    *m_vecs[index]*=z0;}
} */

void LaphVectorHandler::applyPhaseConvention()
{
 multi1d<int> origin(3);
 origin[0]=origin[1]=origin[2]=0;
 for (unsigned int index=0;index<m_vecs.size();++index){
    ColorVector s=peekSite(*m_vecs[index],origin);
    Complex z=conj(peekColor(s,0));  double r=sqrt(toDouble(norm2(z)));
    if (r<1e-12){
       QDPIO::cerr <<endl<<endl<< "problem applying phase convention"<<endl;
       QDPIO::cerr <<"0-th color component at site (0,0,0) has very"<<endl;
       QDPIO::cerr <<"magnitude: "<<z<<endl;
       QDP_abort(1);}
    z/=r;
    *m_vecs[index]*=z;}
}

void LaphVectorHandler::apply_minus_laplacian(LatticeColorVector& finish,
                                              const LatticeColorVector& start)
{
 finish = Real(6.0)*start;
 for (int mu=0;mu<QDP::Nd;mu++){
    finish -= usmear_timeslice[mu]*shift(start, FORWARD, mu);
    finish -=shift(adj(usmear_timeslice[mu])*start, BACKWARD, mu);}
}

void LaphVectorHandler::apply_transf(LatticeColorVector& finish,
                                     const LatticeColorVector& start)
{
 finish = Real(id_coef)*start;
 for (int mu=0;mu<QDP::Nd;mu++){
    finish += usmear_timeslice[mu]*shift(start, FORWARD, mu);
    finish +=shift(adj(usmear_timeslice[mu])*start, BACKWARD, mu);}
 finish*=ov_coef;
}

void LaphVectorHandler::apply_transf_chebyshev(LatticeColorVector& finish,
                                               const LatticeColorVector& start)
{
 vbuf[0]=start;
 apply_transf(vbuf[1],vbuf[0]);
 int ncur=2, nprev=1, nprevprev=0;
 for (int Tn=2;Tn<chebyshev_order;Tn++){
    apply_transf_chebyshev_step(vbuf[ncur++],vbuf[nprev++],vbuf[nprevprev++]);
    if (ncur>2) ncur=0;
    else if (nprev>2) nprev=0;
    else nprevprev=0;
    }
 apply_transf_chebyshev_step(finish,vbuf[nprev],vbuf[nprevprev]);
}

      //  applies the Chebyshev recurrence relation:
      //       result = 2*x*Tm1 - Tm2

void LaphVectorHandler::apply_transf_chebyshev_step(LatticeColorVector& result,
                                                    const LatticeColorVector& Tm1,
                                                    const LatticeColorVector& Tm2)
{
 result = Real(id_coef)*Tm1;
 for (int mu=0;mu<QDP::Nd;mu++){
    result += usmear_timeslice[mu]*shift(Tm1, FORWARD, mu);
    result +=shift(adj(usmear_timeslice[mu])*Tm1, BACKWARD, mu);}
 result*=two_ov_coef;
 result-=Tm2;
}


 
void LaphVectorHandler::assignMatrixVectorMultiply(int lhs_index, int rhs_index)
{
 check_index(rhs_index);
 check_index(lhs_index);
 LatticeColorVector& b= *m_vecs[lhs_index];
 const LatticeColorVector& a= *m_vecs[rhs_index];
 (this->*matrixPtr)(b,a);
}


const LatticeColorVector& LaphVectorHandler::getVector(int index) const
{
 check_index(index);
 return *m_vecs[index];
}

#endif

// ******************************************************************

  }
}
