//Arjun Singh Gambhir
//Real authors: Ben Hoerz, John Bulava, Colin Morningstar, ??
//Code directly ported from chroma_laph, more acknowledgements coming later.
#ifndef FIELD_SMEARING_HANDLER_H
#define FIELD_SMEARING_HANDLER_H

#include "qdp.h"
#include "chromabase.h"
#include "meas/inline/io/named_objmap.h"
#include "util/gauge/stout_utils.h"
#include "gauge_configuration_handler.h"
#include "field_smearing_info.h"
#include "laph_eigen_info.h"
#include "ksrlanczos.h"
#include "data_io_handler.h"

namespace Chroma {
  namespace LaphEnv {



// *****************************************************************
// *                                                               *
// *  "GluonSmearingHandler" handles access to the smeared gauge   *
// *  field and "QuarkSmearingHandler" handles access to the       *
// *  eigenvectors of the Laplacian used in the quark field        *
// *  Laplacian Heaviside (Laph) smearing.  See declarations below *
// *  for available member functions.                              *
// *                                                               *
// *  GluonSmearingHandler:                                        *
// *                                                               *
// *     -- computes smeared gauge field in 4-d chroma_laph,       *
// *          saving the results to a single file in a format      *
// *          readable by 3-d chroma_laph; each record in the      *
// *          file is one time slice of the smeared field          *
// *     -- in 3-d, just reads a time slice from file              *
// *                                                               *
// *  QuarkSmearingHandler:                                        *
// *                                                               *
// *     -- must be constructed in read mode or write mode         *
// *     -- in 3-d write mode, computes the Laplacian eigenvectors,*
// *          saving results to "*_time.nn" files, each file       *
// *          containing all eigenvectors for one time slice       *
// *          The time slices for which eigenvectors are to be     *
// *          computed can be passed on to the compute() method,   *
// *          or else EV on all time slices will be computed.      *
// *     -- in 3-d read mode, reads the eigenvectors from the      *
// *          "*_time.nn" files and makes them available for the   *
// *          hadron handler to do quark displacements             *
// *     -- in 4-d write mode, reads the "*_time.nn" files and     *
// *          combines each level on all time slices into one      *
// *          4-d file, writing out "_level.nn" files              *
// *     -- in 4-d read mode, reads the "*_level.nn" files and     *
// *          makes them available for the quark handler to        *
// *          compute the quark sinks                              *
// *                                                               *
// *    Note:                                                      *
// *                                                               *
// *       If "NamedObjectMap" is passed as smeared_quark_filestub *
// *       it is assumed that eigenvectors are to be stored in     *
// *       and read from Chroma's NamedObjectMap instead of        *
// *       the physical disk. This can be used to avoid            *
// *       unnecessary IO.                                         *
// *       This behavior is implemented using Data*HandlerMemMF    *
// *       handlers instead of Data*HandlerMF that write to file.  *
// *       Note that clearLaphEigenvectors() will NOT actually     *
// *       delete the eigenvectors from memory in this case.       *
// *       Instead the EV should be erased from the NamedObjectMap *
// *       using the ERASE_OBJECT task, the NamedObjectId being    *
// *       the file id of EV files: Laph--SmearedQuarkTimeFile     *
// *                                                               *
// *  A note concerning the phases multiplying each eigenvector:   *
// *                                                               *
// *   Due to the way that noise is introduced, changing the       *
// *   overall phase of any given Laph eigenvector changes the     *
// *   value of the quark line for one particular noise.  The      *
// *   effect of the phase change is to change the noise (the      *
// *   noise is effectively U(1)).  This is not a problem, but     *
// *   erroneous results can occur if the original eigenvector     *
// *   files used to determine the quark sinks get deleted and     *
// *   the eigenvectors have to be reconstructed for making the    *
// *   hadrons.  With different run parameters, the eigensolver    *
// *   could produce a different phase.  The introduction of a     *
// *   phase convention eliminates this potential problem.         *
// *   The phase convention adopted here is as follows:            *
// *                                                               *
// *   *** each eigenvector is multiplied by a phase so that       *
// *       the 0-th color component at site (0,0,0)                *
// *       is real and positive (abort if component very small)    *
// *                                                               *
// *                                                               *
// *  All Laph Handlers follow the member naming convention:       *
// *                                                               *
// *    compute....()  to do original computation                  *
// *    set...()       to internally set from file or NamedObjMap  *
// *                                                               *
// *    get...()       provides access to results                  *
// *                                                               *
// *****************************************************************



class GluonSmearingHandler
{

   class RecordKey
   {
      unsigned int timeval;

    public:

      RecordKey(int in_time) {check_for_failure(in_time); timeval=in_time;}
      RecordKey(const RecordKey& in) : timeval(in.timeval) {}
      RecordKey& operator=(const RecordKey& in) {timeval=in.timeval; return *this;}
      ~RecordKey() {}

      bool operator<(const RecordKey& rhs) const {return (timeval<rhs.timeval);}
      bool operator==(const RecordKey& rhs) const {return (timeval==rhs.timeval);}
      bool operator!=(const RecordKey& rhs) const {return (timeval!=rhs.timeval);}

      unsigned int getTime() const {return timeval;}

      void output(XmlWriter& xmlw) const 
      {push(xmlw,"Key");
       write(xmlw,"TimeValue",getTime());
       pop(xmlw);}

      explicit RecordKey(const unsigned int* buf) {check_for_failure(*buf); timeval=*buf;}
      int numints() const {return 1;} 
      size_t numbytes() const {return sizeof(unsigned int);}
      void copyTo(unsigned int* buf) const { *buf=timeval;}

    private:

      void check_for_failure(int in_time)
      {if (in_time<0){
          QDPIO::cerr << "invalid direction/time in GluonSmearingHandler::RecordKey"
                      <<std::endl; QDP_abort(1);}}
 
   };

       // pointers to internal infos (managed by this handler
       // with new and delete)

   const GluonSmearingInfo *smearPtr;
   const GaugeConfigurationInfo *uPtr;
   std::string h_filename;

       // storage and/or references to internal data

   DataGetHandlerSF<GluonSmearingHandler,RecordKey,
                    multi1d<LatticeColorMatrix> > *dh_ptr;

       // prevent copying ... handler might contain large
       // amounts of data

   GluonSmearingHandler(const GluonSmearingHandler&);
   GluonSmearingHandler& operator=(const GluonSmearingHandler&);


 public:


   GluonSmearingHandler();

   GluonSmearingHandler(const GluonSmearingInfo& gluon_smearing,
                        const GaugeConfigurationInfo& gauge,
                        const std::string& smearedFieldFileName);

   void setInfo(const GluonSmearingInfo& gluon_smearing,
                const GaugeConfigurationInfo& gauge,
                const std::string& smearedFieldFileName);

   ~GluonSmearingHandler();

   void clear();  // clears everything in the handler


           // access to the info

   bool isInfoSet() const;

   const GluonSmearingInfo& getGluonSmearingInfo() const;

   const GaugeConfigurationInfo& getGaugeConfigurationInfo() const;

   std::string getFileName() const {return h_filename;}

#if (QDP_ND == 4)

           // Computes the smeared gauge field and puts it into the file
           // whose name is stored in the GluonSmearingInfo.  This must 
           // be run in four dimensions, but the output file is meant to
           // be read in three dimensions.

   void computeSmearedGaugeField(const std::string& gauge_id = "");


#elif (QDP_ND == 3)

           // Set smeared gauge field, reading from the file specified in
           // GluonSmearingInfo, if smearing parameters are the same as 
           // the current info (error otherwise).  Memory is allocated for
           // result and this time slice is added to the data map.

   void setSmearedGaugeFieldTimeSlice(int time_slice);

           // Provides access to the smeared gauge field time slice.  
           // Calls "set" to read from file if needed.

   const multi1d<LatticeColorMatrix>& getSmearedGaugeFieldTimeSlice(
                                       int time_slice);

           // remove one time slice, or clear all of the smeared gauge 
           // fields time slices from memory

   void removeData(int time_slice);

   void clearData();


#endif



 private:

   void set_info(const GluonSmearingInfo& gluon_smearing,
                 const GaugeConfigurationInfo& gauge,
                 const std::string& smearedFieldFileName);

   void filefail(const std::string& message);
   void check_info_set(const std::string& name) const;
 
   bool checkHeader(XmlReader& xmlr);
   void writeHeader(XmlWriter& xmlw);

   friend class DataGetHandlerSF<GluonSmearingHandler,RecordKey,
                                 multi1d<LatticeColorMatrix> >;
   friend void applyGaugeTransformation(XmlReader& xml_in);

};



// *******************************************************************





class QuarkSmearingHandler
{

   struct TimeKey
   {
    int value;

    TimeKey() : value(0) {}
    TimeKey(int in_val) : value(in_val) {}
    TimeKey(XmlReader& xmlr)
    {xmlread(xmlr,"TimeSlice",value,"QuarkSmearingHandler::TimeKey");}
    TimeKey(const TimeKey& in) : value(in.value) {}
    ~TimeKey() {}

    bool operator<(const TimeKey& rhs) const
    {return (value<rhs.value);}

    bool operator==(const TimeKey& rhs) const
    {return (value==rhs.value);}

    void output(XmlWriter& xmlw) const 
    {write(xmlw,"TimeSlice",value);}

    explicit TimeKey(const unsigned int* buf) : value(*buf) {}
    int numints() const {return 1;} 
    size_t numbytes() const {return sizeof(unsigned int);}
    void copyTo(unsigned int* buf) const { *buf=value;}

   };

   struct LevelKey
   {
    int value;

    LevelKey() : value(0) {}
    LevelKey(int in_val) : value(in_val) {}
    LevelKey(XmlReader& xmlr)
    {xmlread(xmlr,"Level",value,"QuarkSmearingHandler::LevelKey");}
    LevelKey(const LevelKey& in) : value(in.value) {}
    ~LevelKey() {}

    bool operator<(const LevelKey& rhs) const
    {return (value<rhs.value);}

    bool operator==(const LevelKey& rhs) const
    {return (value==rhs.value);}

    void output(XmlWriter& xmlw) const 
    {write(xmlw,"Level",value);}

    explicit LevelKey(const unsigned int* buf) : value(*buf) {}
    int numints() const {return 1;} 
    size_t numbytes() const {return sizeof(unsigned int);}
    void copyTo(unsigned int* buf) const { *buf=value;}

   };

       // pointers to internal infos (managed by this handler
       // with new and delete)

   const GaugeConfigurationInfo *uPtr;
   const GluonSmearingInfo *gSmearPtr;
   QuarkSmearingInfo *qSmearPtr;
   std::string smearedQuarkFileStub; 
   bool m_read_mode;
   bool m_memory_only;

       // storage and/or references to internal data, and other handlers

#if (QDP_ND == 3)
   multi1d<double> Eigenvalues;
   DataGetHandlerBaseMF<QuarkSmearingHandler,TimeKey,LevelKey,
                    LatticeColorVector> *dh_ptr;
#elif (QDP_ND == 4)
   multi2d<double> Eigenvalues;
   DataGetHandlerMF<QuarkSmearingHandler,LevelKey,LevelKey,
                    LatticeColorVector> *dh_ptr;
#endif

       // prevent copying ... handler might contain large
       // amounts of data

   QuarkSmearingHandler(const QuarkSmearingHandler&);
   QuarkSmearingHandler& operator=(const QuarkSmearingHandler&);


 public:


   QuarkSmearingHandler();

   QuarkSmearingHandler(const GluonSmearingInfo& gluon_smearing,
                        const GaugeConfigurationInfo& gauge,
                        const QuarkSmearingInfo& quark_smearing,
                        const std::string smeared_quark_file_stub,
                        bool read_mode=true);

   void setInfo(const GluonSmearingInfo& gluon_smearing,
                const GaugeConfigurationInfo& gauge,
                const QuarkSmearingInfo& quark_smearing,
                const std::string smeared_quark_file_stub,
                bool read_mode=true);

          // update if number of eigenvectors needs to be increased
   void updateSmearing(const QuarkSmearingInfo& quark_smearing);

   ~QuarkSmearingHandler();

   void clear();      // clears everything in the handler


           // access to the info

   bool isInfoSet() const;

   const GluonSmearingInfo& getGluonSmearingInfo() const;

   const QuarkSmearingInfo& getQuarkSmearingInfo() const;

   const GaugeConfigurationInfo& getGaugeConfigurationInfo() const;

   const std::string& getSmearedQuarkFieldFileStub() const;

   void getFileMap(XmlWriter& xmlout) const;

   void outputKeys(XmlWriter& xmlout);



#if (QDP_ND == 3)

           // Compute largest eigenvalue of smeared covariant Laplacian.
           // This value is useful for Chebyshev acceleration in
           // computing all Laph eigenvectors.

   double estimateLargestLaplacianEigenvalue(const std::string& smeared_gauge_file);

           // Compute the Laph Eigenvectors and output to file.

   void computeLaphEigenvectors(const LaphEigenSolverInfo& solver_info,
                                const std::string& smeared_gauge_file,
				std::set<int>& timeslices,
                                int striping_factor=1, int striping_unit=0);


           // get and query data when in read mode

   const LatticeColorVector& getLaphEigenvector(int time, int eigpair_num);

   bool queryLaphEigenvector(int time, int eigpair_num);

   void removeLaphEigenvector(int time, int eigpair_num);

   void clearLaphEigenvectors();

	   // read eigenvectors from disk into the NamedObjectMap using
	   // DataPutHandlerMemMF

   void readEigenvectorsIntoMemory(int time);


#elif (QDP_ND == 4)

           // The 3-d computation of the Laph eigenvectors produces
           // one file for each time slice, and each file contains all
           // eigen-levels up to the requested number of eigenvectors.
           // This routine re-organizes these files, producing one file
           // for each level, but each file contains all time slices.

   void combineTimeSlices(int striping_factor=1, int striping_unit=0);


           // get and query data when in read mode

   const LatticeColorVector& getLaphEigenvector(int eigpair_num);

   bool queryLaphEigenvector(int eigpair_num);

   void removeLaphEigenvector(int eigpair_num);

   void clearLaphEigenvectors();

           // checks to see if all _level files exist.

   bool checkAllLevelFilesExist();



#endif

 private:

   void set_info(const GluonSmearingInfo& gluon_smearing,
                 const GaugeConfigurationInfo& gauge,
                 const QuarkSmearingInfo& quark_smearing,
                 const std::string smearedQuarkFileStub,
                 bool read_mode);

   void check_info_set(const std::string& name,
                       int check_mode=0) const;

   void failure(const std::string& message);

   bool checkHeader(XmlReader& xmlr, int suffix);

#if (QDP_ND == 3)
   void writeHeader(XmlWriter& xmlw, const TimeKey& fkey, int suffix);
#elif (QDP_ND == 4)
   void writeHeader(XmlWriter& xmlw, const LevelKey& fkey, int suffix);
#endif

   friend class DataGetHandlerMF<QuarkSmearingHandler,TimeKey,LevelKey,
                                 LatticeColorVector>;
   friend class DataGetHandlerMF<QuarkSmearingHandler,LevelKey,LevelKey,
                                 LatticeColorVector>;
   friend class DataPutHandlerMF<QuarkSmearingHandler,TimeKey,LevelKey,
                                 LatticeColorVector>;
   friend class DataPutHandlerMF<QuarkSmearingHandler,LevelKey,LevelKey,
                                 LatticeColorVector>;

   friend class DataGetHandlerMemMF<QuarkSmearingHandler,TimeKey,LevelKey,
                                 LatticeColorVector>;
   friend class DataPutHandlerMemMF<QuarkSmearingHandler,TimeKey,LevelKey,
                                 LatticeColorVector>;

   friend class GlueballHandler;


};


 // ****************************************************************

#if (QDP_ND == 3)

 
// *****************************************************************
// *                                                               *
// *  "LaphVectorHandler" is used by "KSRLanczosEigenSolver" to    *
// *  handle the eigenvectors and define the matrix (the smeared   *
// *  covariant Laplacian) whose eigenvectors are sought.          *
// *                                                               *
// *  The matrix whose eigenvectors are sought is denoted by "A".  *
// *  "A" is minus one times the smeared covariant Laplacian:      *
// *            A = -Delta                                         *
// *  The eigenvalues of "A" lie between 0 and some maximum value  *
// *  denoted by "largest_eigenvalue" or "L" for short.  All of    *
// *  the eigenvalues of "A" are real and positive.  We wish to    *
// *  determine the eigenvectors corresponding to the lowest-lying *
// *  eigenvalues lying between 0 and "cutoff_eigenvalue" or "C"   *
// *  for short.  So we have                                       *
// *                                                               *
// *    desired part of spectrum:   0 ... "C"                      *
// *   unwanted part of spectrum:  "C" .. "L"                      *
// *                                                               *
// *  The rate of convergence to solution depends on the spacing   *
// *  between the levels.  Convergence is much faster for widely   *
// *  spaced levels.  So convergence can be accelerated by         *
// *  transforming the spectrum so that the desired part of the    *
// *  spectrum is more widely spaced.  The following               *
// *  transformation is applied first:                             *
// *                                                               *
// *         1 - 2*(A-C)/(L-C)                                     *
// *                                                               *
// *  The above transformation maps the unwanted spectrum to the   *
// *  range -1 .. 1, and the desired part lies above 1.            *
// *  Then Chebyshev polynomials can be applied.  Eigenvalues      *
// *  lying between -1 and 1 are suppressed (stay between -1..1),  *
// *  and the desired eigenvalues above 1 get spaced out to        *
// *  large values above 1 to speed up convergence.                *
// *  The lowest-lying eigenvalue becomes the highest-lying        *
// *  transformed eigenvalue.  (Transforming the desired levels    *
// *  to the region above 1 is more convenient since it allows     *
// *  the use of Chebyshev polynomials of any order, both even     *
// *  and odd.)                                                    *
// *                                                               *
// *  The Chebyshev polynomials are applied using the following    *
// *  recurrence relation:                                         *
// *                                                               *
// *      T_0(x) = 1    T_1(x) = x                                 *
// *      T_n(x) = 2*x * T_(n-1)(x) - T_(n-2)(x)                   *
// *                                                               *
// *****************************************************************


class LaphVectorHandler : public KSRLVectorHandler
{

   std::vector<LatticeColorVector* > m_vecs;

   const multi1d<LatticeColorMatrix>& usmear_timeslice;
   
   std::vector<LatticeColorVector> vbuf;
   double ov_coef, id_coef, two_ov_coef;
   int chebyshev_order;   // must be >= 2
   bool use_acceleration;   // if true, Chebyshev acceleration is used,
                            // if false, -Laplacian is used as matrix
   void (LaphVectorHandler::*matrixPtr)(LatticeColorVector&, const LatticeColorVector&);
   char spec_end;
   char start_vector_type;

 public:
  
   LaphVectorHandler(const multi1d<LatticeColorMatrix>& smeared_gauge_timeslice);

   LaphVectorHandler(const multi1d<LatticeColorMatrix>& smeared_gauge_timeslice,
                     int nvec);

   void setChebyshevAccelerationOn(double largest_eigenvalue,
                                   double cutoff_eigenvalue,
                                   int chebyshev_polynomial_order);

   void setChebyshevAccelerationOff();

   void setInitialToRandom();

   void setInitialToEqualComponents();

   ~LaphVectorHandler() {clear();}
    
   void clear();

   void resize(int nvec);

   int getNumberOfVectors() const { return m_vecs.size();}

   char getSpectrumEnd() const { return spec_end;}

   void assignInitialUnitNorm(int index);

   void assignRandomUnitNorm(int index);

   double InnerProductRealPart(int left_index, int right_index);

   DComplex InnerProduct(int left_index, int right_index);

   void addTo(int left_index, double coef, int right_index);

   void addTo(int left_index, const DComplex& coef, int right_index);

   void copyTo(int left_index, int right_index);

   void copyTo(int left_index, double coef, int right_index);

   void swap(int left_index, int right_index);

   void unitNormalize(int index, double& vnorm);

   void assignMatrixVectorMultiply(int lhs_index, int rhs_index);

   const LatticeColorVector& getVector(int index) const;

   void applyPhaseConvention();


 private:
 
   void check_index(int index) const;

   void apply_minus_laplacian(LatticeColorVector& finish,
                              const LatticeColorVector& start);

   void apply_transf(LatticeColorVector& finish,
                     const LatticeColorVector& start);

   void apply_transf_chebyshev(LatticeColorVector& finish,
                               const LatticeColorVector& start);

   void apply_transf_chebyshev_step(LatticeColorVector& result,
                                    const LatticeColorVector& Tm1,
                                    const LatticeColorVector& Tm2);

};


#endif


// ***************************************************************
  }
}
#endif  
