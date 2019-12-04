//Arjun Singh Gambhir
//Real authors: Ben Hoerz, John Bulava, Colin Morningstar, ??
//Code directly ported from chroma_laph, more acknowledgements coming later.
#ifndef KSRLANCZOSEIGENSOLVER_H
#define KSRLANCZOSEIGENSOLVER_H
#include "qdp.h"
#include "chromabase.h"
#include <vector>
#include <string>   // remove when debugged



namespace Chroma {
  namespace LaphEnv {

// *******************************************************************
// *                                                                 *
// *  An object of the class "KSRLanczosEigenSolver" sets up and     *
// *  computes the lowest-lying or highest-lying eigenvalues and     *
// *  eigenvectors of a Hermitian matrix using the Krylov-spectral   *
// *  restarted Lanczos (KSRL) method.  The constructor sets up      *
// *  and performs the computation. Results can then be queried      *
// *  using the member functions of the object.  The method is       *
// *  suitable for very large matrices.  Only the action of the      *
// *  matrix on a vector is needed. The method used is based on      *
// *  the thick-restarted Lanczos method of Wu and Simon (trlan).    *
// *                                                                 *
// *  A "KSRLanczosEigenSolver" object manipulates and stores        *
// *  vectors which are usually quite large (and partitioned onto    *
// *  different compute nodes).  The "KSRLanczosEigenSolver" class   *
// *  uses the base class "KSRLVectorHandler" to actually carry      *
// *  out all of the needed vector manipulations.  Hence, the        *
// *  "KSRLanczosEigenSolver" class itself does not perform any      *
// *  parallel computations.   This assumes that the dimension       *
// *  "krdim" of the Krylov space will be small enough that          *
// *  diagonalizing a "krdim" by "krdim" matrix can be done          *
// *  quickly on the primary node. "KSRLanczosEigenSolver" assumes   *
// *  that all needed parallelizations and optimizations are done    *
// *  in the user-written class derived from "KSRLVectorHandler".    *
// *                                                                 *
// *  Sample usage:                                                  *
// *                                                                 *
// *   int num_eigenpair_desired = 20;                               *
// *   char spectrum_end = 'H';    // highest 'H' or lowest 'L'      *
// *   double tolerance = 1e-6;    // residual < tol * matrix-norm   *
// *   int krylov_dimension = 60;  // must exceed number desired     *
// *   int max_restarts = 80;                                        *
// *   derivedVectorHandler V;                                       *
// *                                                                 *
// *   KSRLanczosEigenSolver Comp(V, num_eigenpair_desired,          *
// *                              tolerance, krylov_dimension,       *
// *                              max_restarts, spectrum_end);       *
// *   if (Comp.wasSuccessful()){                                    *
// *      cout << "matrix norm estimate = "                          *
// *           << Comp.getMatrixNormEstimate() << endl;              *
// *      for (int i=0; i<num_eigenpair_desired; i++){               *
// *         cout << "Eigenvalue("<<i<<") = "                        *
// *              << Comp.getEigenvalue(i);                          *
// *         cout << "  Residual("<<i<<") = "                        *
// *              << Comp.getResidual(i) << endl;                    *
// *         }                                                       *
// *      }                                                          *
// *   Access to the eigenvectors is through "V" and should be       *
// *   integer indexed from 0 to num_eigenpair_desired-1.            *
// *   The user is responsible for this access.                      *
// *                                                                 *
// *                                                                 *
// *  To use "KSRLanczosEigenSolver", the end user must write a      *
// *  class derived from abstract base class "KSRLVectorHandler".    *
// *  This class must have the following members:                    *
// *                                                                 *
// *   void clear();                                                 *
// *                                                                 *
// *     -- deletes all vectors, frees up memory                     *
// *                                                                 *
// *   int getNumberOfVectors() const;                               *
// *                                                                 *
// *     -- returns number of vectors currently in memory            *
// *                                                                 *
// *   void resize(int nvec);                                        *
// *                                                                 *
// *     -- resizes the number of vectors in memory; must work       *
// *        like STL resize; if new size is larger than current      *
// *        size, new vectors are appended without changing the      *
// *        original content; if new size is smaller than current    *
// *        size, vectors are deleted from the end, that is,         *
// *        vectors with the highest indices are deleted, without    *
// *        changing content of vectors with indices 0 to            *
// *        new size-1.                                              *
// *                                                                 *
// *   void assignInitialUnitNorm(int index);                        *
// *                                                                 *
// *     -- assigns an initial value with unit norm to the vector    *
// *        with index "index"; should abort if index < 0  or        *
// *        index >= current size;                                   * 
// *                                                                 *
// *   void assignRandomUnitNorm(int index);                         *
// *                                                                 *
// *     -- assigns an initial value with unit norm to the vector    *
// *        with index "index"; should abort if index < 0  or        *
// *        index >= current size;                                   * 
// *                                                                 *
// *   std::complex<double> InnerProduct(int left_index,             *
// *                                     int right_index);           *
// *                                                                 *
// *     -- returns the inner product of the vector having index     *
// *        "left_index" with the vector of index "right_index".     *
// *        The complex conjugate applies to the left vector.        *
// *                                                                 *
// *   double InnerProductRealPart(int left_index,                   *
// *                               int right_index);                 *
// *                                                                 *
// *     -- same as above, except that only the real part of the     *
// *        inner product is returned                                *
// *                                                                 *
// *   void addTo(int left_index, double coef, int right_index);     *
// *                                                                 *
// *     -- adds "coef" times the vector of index "right_index"      *
// *        to the vector at index "left_index"                      *
// *                                                                 *
// *   void addTo(int left_index, const std::complex<double>& coef,  *
// *              int right_index);                                  *
// *                                                                 *
// *     -- adds "coef" times the vector of index "right_index"      *
// *        to the vector at index "left_index", where "coef" is     *
// *        complex valued                                           *
// *                                                                 *
// *   void copyTo(int left_index, int right_index);                 *
// *                                                                 *
// *     -- copies the vector of index "right_index"                 *
// *        to the vector at index "left_index"                      *
// *                                                                 *
// *   void copyTo(int left_index, double coef, int right_index);    *
// *                                                                 *
// *     -- copies "coef" times the vector of index "right_index"    *
// *        to the vector at index "left_index"                      *
// *                                                                 *
// *   void unitNormalize(int index, double& norm);                  *
// *                                                                 *
// *     -- computes the norm of the vector of index "index",        *
// *        putting the value into "norm", then rescales the         *
// *        vector of index "index" so that it has unit norm         *
// *                                                                 *
// *   void assignMatrixVectorMultiply(int left_index,               *
// *                                   int right_index);             *
// *                                                                 *
// *     -- specification of the matrix whose eigenpairs are         *
// *        sought is done here; this routine multiplies the         *
// *        vector of index "right_index" by the matrix, putting     *
// *        the result into the vector of index "left_index"         *
// *                                                                 *
// *******************************************************************


 //  Base class for manipulating vectors with integer index (zero offset).


class KSRLVectorHandler
{

 public:

   KSRLVectorHandler() {}
   virtual ~KSRLVectorHandler() {}

   virtual void clear() = 0;
   virtual void resize(int nvec) = 0;   
   virtual int getNumberOfVectors() const = 0;

   virtual void assignInitialUnitNorm(int index) = 0;
   virtual void assignRandomUnitNorm(int index) = 0; 

   virtual DComplex InnerProduct(int left_index, int right_index) = 0;
   virtual double InnerProductRealPart(int left_index, int right_index) = 0;

   virtual void addTo(int left_index, double coef, int right_index) = 0;
   virtual void addTo(int left_index, const DComplex& coef, int right_index) = 0;
   virtual void copyTo(int left_index, int right_index) = 0;
   virtual void copyTo(int left_index, double coef, int right_index) = 0;
   virtual void swap(int left_index, int right_index) = 0;
   virtual void unitNormalize(int index, double& vnorm) = 0;

   virtual void assignMatrixVectorMultiply(int left_index, int right_index) = 0;

//   virtual void print(const std::string& name, int index) const = 0;  // remove when debugged

};


 // ******************************************************************


     // Krylov-spectral restarted Lanczos computation of the
     // lowest-lying or highest-lying eigenvalues and eigenvectors.
     // Note: this class assumes that accuracy is more important to
     // the end users than speed, and that obtaining the eigenvectors
     // is important.  Hence, full re-orthogonalization is done.


class KSRLanczosEigenSolver
{

   int ned;                 // number of eigenpairs to find
   bool ascending;          // find lowest (true) or highest (false) eigenvalues
   double tol;              // tolerance in residual norms relative to matrix norm
   int krdim;               // dimension of Krylov space
   int maxits;              // maximum number of iterations (restarts)
   int verbosity;

   bool success;            // status of calculation
   double matrixNorm;       // norm of matrix (highest absolute value of all
                            //  Ritz values encountered)

   KSRLVectorHandler& q;    // reference to the (base) vector handler

   int nlocked, nconverged, nkeep;
   double epsilon;                  // unit roundoff

   std::vector<double> alpha,beta,residual,y;


 public:

           // constructor sets up the computation, then does the
           // computation;  results can then be queried from the object

   KSRLanczosEigenSolver(KSRLVectorHandler& q_handler,
                         int num_eigenpair_desired,
                         double tolerance, int krylov_dimension, 
                         int max_restarts, char spectrum_end = 'L',
                         int output_verbosity = 0 );

   ~KSRLanczosEigenSolver() {}

   int getNumberOfEigenPairsDesired() const {return ned;}

   double getTolerance() const {return tol;}

   int getKrylovDimension() const {return krdim;}

   int getMaximumRestarts() const {return maxits;}

   double getMatrixNormEstimate() const;


   bool wasSuccessful() const {return success;}

   double getEigenvalue(int index) const;

   double getResidual(int index) const;


 private:

         // prohibit copying
   KSRLanczosEigenSolver(const KSRLanczosEigenSolver& lhs);            
   KSRLanczosEigenSolver& operator=(const KSRLanczosEigenSolver& lhs);  


   void computeEigenpairs();
   void extend_Lanczos();
   bool reorthogonalize(int index);
   bool compute_Ritz_values_eigvecs();
   int set_number_to_keep(int itconv);
   void compute_kept_Ritz_vectors(int itkeep);
   void reorder();


}; 


 // **********************************************************************
  }
}

#endif
