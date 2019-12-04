//Arjun Singh Gambhir
//Real authors: Ben Hoerz, John Bulava, Colin Morningstar, ??
//Code directly ported from chroma_laph, more acknowledgements coming later.
#include "ksrlanczos.h"
#include <iostream>
#include <float.h>
#include <cstdlib> // removed when debugged
using namespace std;


  // Prototypes of routines in LAPACK library--to call Fortran
  // routines from a C++ program, use extern "C" to tell the
  // compiler that the external routine is a C routine; then
  // add an underscore to the end of the routine name since
  // the routine is in Fortran.  All parameters must be passed
  // as pointers and two-dimensional arrays must follow
  // Fortran conventions.

extern "C" {


   void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, 
               double *w, double *work, int *lwork,  int *info);

}

namespace Chroma {
  namespace LaphEnv {

// *******************************************************************

KSRLanczosEigenSolver::KSRLanczosEigenSolver(
          KSRLVectorHandler& q_handler, int num_eigenpair_desired, 
          double tolerance, int krylov_dimension, int max_restarts,
          char spectrum_end, int output_verbosity )
        : q(q_handler)
{
 success = true;

 ned=0;
 if (num_eigenpair_desired<1){
    QDPIO::cerr << "Number of desired eigenvalues must exceed zero in "
         << "KSRLanczosEigenSolver" << endl;
    success=false;}
 ned=num_eigenpair_desired;

 if ((spectrum_end=='H')||(spectrum_end=='h')
    ||(spectrum_end=='U')||(spectrum_end=='u')) ascending=false;
 else ascending=true;

 tol=0.0;
 if (tolerance<=0.0){
    QDPIO::cerr << "Invalid tolerance requested in KSRLanczosEigenSolver"<<endl;
    success=false;}
 tol=tolerance;

 krdim=0;
 if (krylov_dimension<(ned+6)){
    QDPIO::cerr << "Krylov dimension must be at least number of desired"
              << " eigenpairs plus 6 in KSRLanczosEigenSolver"<<endl;
    success=false;}
 krdim=krylov_dimension;

 maxits=0;
 if (max_restarts<3){
    QDPIO::cerr << "Maximum restarts must exceed 2 in KSRLanczosEigenSolver"<<endl;
    success=false;}
 maxits=max_restarts;

 verbosity=output_verbosity;
 if (verbosity<0) verbosity=0;
 if (verbosity>2) verbosity=2;

 if (success) computeEigenpairs();
}

// ********************************************************************

double KSRLanczosEigenSolver::getMatrixNormEstimate() const
{
 if (success) return matrixNorm;
 QDPIO::cerr << "Unsuccessful KSRLanczosEigenSolver computation:"<<endl;
 QDPIO::cerr << "Matrix norm query is unreliable"<<endl;
 return 0.0;
}

double KSRLanczosEigenSolver::getEigenvalue(int level) const
{
 if ((level<0)||(level>=ned)){
    QDPIO::cout << "Invalid level in KSRLanczosEigenSolver::getEigenvalue"<<endl;
    throw(string("invalid index"));}
 if (!success){
    QDPIO::cout << "Computation did not succeed: unreliable eigenvalue"<<endl;}
 return alpha[level];
}

double KSRLanczosEigenSolver::getResidual(int level) const
{
 if ((level<0)||(level>=ned)){
    QDPIO::cout << "Invalid level in KSRLanczosEigenSolver::getResidue"<<endl;
    throw(string("invalid index"));}
 if (!success){
    QDPIO::cout << "Computation did not succeed: unreliable residual"<<endl;}
 return residual[level];
}

// ********************************************************************


void KSRLanczosEigenSolver::computeEigenpairs()
{
 QDPIO::cout << endl << "  KSRLanczosEigenSolver starting computation:"<<endl<<endl;
 if (ascending) 
    QDPIO::cout << "   Number of desired lowest-lying eigenpairs = "<< ned <<endl;
 else
    QDPIO::cout << "  Number of desired highest-lying eigenpairs = "<< ned <<endl;
 QDPIO::cout << "  Residual tolerance relative to matrix norm = "<<tol<<endl; 
 QDPIO::cout << "                   Dimension of Krylov space = "<<krdim<<endl;
 QDPIO::cout << "                  Maximum number of restarts = "<<maxits<<endl<<endl;

 matrixNorm = 0.0;  // norm of matrix will be highest absolute value of all
                    //  Ritz values encountered; used for convergence test
 nlocked=nconverged=nkeep=0;
 epsilon=DBL_EPSILON;   // unit roundoff error...defined in float.h

 success = false;

 q.resize(krdim+1);            // storage space for the Lanczos vectors
 q.assignInitialUnitNorm(0);   // starting vector (unit normalized)
 alpha.resize(krdim);
 beta.resize(krdim);
 residual.resize(krdim);
 int iter;

 START_CODE();
 StopWatch seiko;
 seiko.start();

    // the main iteration

 for (iter = 0; iter < maxits; iter++){

      // Extend the Lanczos factorization: construct Lanczos
      // vectors q[nkeep+1]..q[krdim], and assign 
      // alpha[nkeep]..alpha[krdim-1] and beta[nkeep]..beta[krdim-1]

    extend_Lanczos(); 

      // Compute the Ritz values in alpha[nlocked]..alpha[krdim-1]
      // and residuals res[nlocked]..res[krdim-1]
      // Compute eigenvectors of arrowhead+tridiagonal, store in "y"
      // (Does not hurt to compute on every node...as long as
      // boolean returned is the same on all nodes.)

    if (!compute_Ritz_values_eigvecs()){
       QDPIO::cerr << "Could not diagonalize:"<<endl;
       break;}

      // update matrixNorm (keep locked vectors as before, though)

    for (int k=nlocked;k<krdim;k++)
       if (abs(alpha[k])>matrixNorm) matrixNorm=abs(alpha[k]);

      // convergence tests and locking

    int itlock = 0;
    if (Layout::primaryNode()){
       double lockval = epsilon*matrixNorm;
       for (int k=1;k < krdim-nlocked;k++){
          if (residual[nlocked+k] < lockval) itlock=k;
          else break;}}
    QDPInternal::broadcast(itlock);

    int itconv = itlock;
    if (Layout::primaryNode()){
       double convval = tol*matrixNorm;
       for (int k=itlock+1;k < krdim-nlocked;k++){
          if (residual[nlocked+k] < convval) itconv=k;
          else break;}}
    QDPInternal::broadcast(itconv);

      // determine the number of Ritz vectors to keep

    int itkeep=0;
    if (Layout::primaryNode()){
       itkeep=set_number_to_keep(itconv);}
    QDPInternal::broadcast(itkeep);

      // compute the Ritz vectors to be kept and compute
      // beta[nlocked]..beta[nlocked+itkeep-1]

    compute_kept_Ritz_vectors(itkeep);

      // update numbers

    nconverged = nlocked + itconv;
    nkeep = nlocked + itkeep;
    nlocked += itlock;

    if (verbosity>0){
       QDPIO::cout << "iteration number: "<<iter<<endl;
       QDPIO::cout << "          nlocked = "<<nlocked<<endl;
       QDPIO::cout << "            nkeep = "<<nkeep<<endl;
       QDPIO::cout << "       nconverged = "<<nconverged<<endl;}

    if (verbosity>1){
       for (int k=0;k<krdim;k++){
          QDPIO::cout <<" Ritz value["<<k<<"] = "<< alpha[k]
                      << "  residual["<<k<<"] = "<<residual[k]<<endl;}}

    if (nconverged>=ned){
       success=true;
       reorder();   // locking could have changed eigenvalue ordering
       break;}   // done

    }

 q.resize(ned);
 alpha.resize(ned);
 residual.resize(ned);
 beta.clear();
 y.clear();

 seiko.stop();
 QDPIO::cout << endl<<"KSRLanczos computation: total time = "
             << seiko.getTimeInSeconds() << " secs" << endl;

 if (success){
    QDPIO::cout <<endl<< "Successful KSRLanczos computation:"<<endl;
    QDPIO::cout << "       number of iterations = "<<iter+1<<endl<<endl;}
 else{
    QDPIO::cout << endl << endl;
    QDPIO::cout << "Warning: maximum iterations encountered before convergence"<<endl;
    QDPIO::cout << "  to requested tolerance was achieved"<<endl<<endl;}
}

// ********************************************************************

    // When this routine is called,
    //   q[0]..q[nlocked-1] are locked since converged to machine precision,
    //   q[0]..q[nconverged-1] have converged to required tolerance
    //   q[0]..q[nkeep] are to be retained
    //   alpha[0]..alpha[nkeep-1] are the kept Ritz values
    //   beta[0]..beta[nkeep-1] have been set
    //   residual[0]..residual[nkeep-1] are the residuals for the kept vectors
    // This routine computes
    //   q[nkeep+1]..q[krdim]
    //   alpha[nkeep] .. alpha[krdim-1]
    //   beta[nkeep] .. beta[krdim-1]

void KSRLanczosEigenSolver::extend_Lanczos()
{
 for (int i=nkeep;i<krdim;i++){

    q.assignMatrixVectorMultiply(i+1,i);
    alpha[i]=q.InnerProductRealPart(i,i+1);

      // orthogonalize

    q.addTo(i+1,-alpha[i],i);
    int jstart = (i>nkeep) ? i-1 : 0;
    for (int j=jstart;j<i;j++) q.addTo(i+1,-beta[j],j);

      // reorthogonalize (to avoid roundoff error)

    bool no_reset=reorthogonalize(i+1);

      // normalize

    if (no_reset) q.unitNormalize(i+1,beta[i]);
    else{
       QDPIO::cout << "encountered a zero norm"<<endl;
       beta[i]=0.0;
       if (i==krdim-1){
          QDPIO::cout << "Oops, last Lanczos vector was a random set"<<endl;
          QDPIO::cout << "This is bad...better quit until better programmed"<<endl;
          exit(1);}
       }
/*
QDPIO::cout <<endl<< " constructed Lanczos vector "<<i+1<<endl;
//q.print("LV",i+1);
double maxoverlap=0.0;
for (int j=0;j<=i;j++){
    DComplex z=q.InnerProduct(j,i+1);
    double rr=toDouble(sqrt(real(z)*real(z)+imag(z)*imag(z)));
   if (rr>maxoverlap) maxoverlap=rr;}
QDPIO::cout << "Orthogonality: maximum overlap is "<<maxoverlap<<endl;
//for (int j=0;j<=i;j++)
//    QDPIO::cout << "LV["<<j<<"].LV["<<i+1<<"] = "<<q.InnerProduct(j,i+1)<<endl;
//QDPIO::cout << "alpha["<<i<<"] := "<<alpha[i]<<endl;
//QDPIO::cout << "beta["<<i<<"] := "<<beta[i]<<endl; */
    } 
}

// *****************************************************************

    //  This routine re-orthogonalizes Lanzcos vector q[index]
    //  against all previous Lanczos vectors q[0]..q[index-1].
    //  Since accuracy of the eigenvectors is important, we always
    //  globally re-orthogonalize, even though this slows down the
    //  computations.  The decision to re-orthogonalize multiple
    //  times is based on the K-criterion.  If the norm of the
    //  vector decreases 1/K, where K = sqrt(2), then further
    //  re-orthgonalization is done.  A maximum of 4 re-orthogonalizations
    //  is enforced.   alpha[index-1] is adjusted.

bool KSRLanczosEigenSolver::reorthogonalize(int index)
{
 double qq = q.InnerProductRealPart(index,index);
 for (int step=0;step<4;step++){
//    QDPIO::cout << "reorthogonalization step = "<<step<<endl;
//cout.precision(16);
//QDPIO::cout << "alpha[index-1] = "<<alpha[index-1]<<endl;
//    alpha[index-1]+=q.InnerProductRealPart(index-1,index);
//QDPIO::cout << "alpha[index-1] = "<<alpha[index-1]<<"  after "<<endl;

    for (int j=0;j<index;j++) q.addTo(index,-q.InnerProduct(j,index),j);


    double qq2 = q.InnerProductRealPart(index,index);
//QDPIO::cout << "qq = "<<qq<<"  qq2 = "<<qq2<<endl;
    if (qq2 > 0.5*qq) return true;
    qq=qq2;
    }
 QDPIO::cout << "******REORTHOGONALIZED MAXIMUM 4 TIMES******"<<endl;
 QDP_abort(1);

/* double qq = q.InnerProductRealPart(index,index);
 double eta = alpha[index-1]*alpha[index-1];
 int jstart = (index>nkeep+1) ? index-2 : 0;
 for (int j=jstart;j<index-1;j++) eta+=beta[j]*beta[j];

 QDPIO::cout << "qq = "<<qq<<"  eta = "<<eta<<endl;
 if (qq>eta){
    QDPIO::cout << "need to do local reorthogonalization"<<endl;
    alpha[index-1]+=q.InnerProductRealPart(index-1,index);
    for (int j=jstart;j<index;j++) q.addTo(index,-q.InnerProduct(j,index),j);
    }
 else if (qq>epsilon*epsilon*eta){
    QDPIO::cout << "need to do global orthogonalization"<<endl;
    alpha[index-1]+=q.InnerProductRealPart(index-1,index);
    for (int j=0;j<index;j++) q.addTo(index,-q.InnerProduct(j,index),j);
    }
 else QDPIO::cout << "need to replace with random vector"<<endl;
*/
 return false;
}

// ********************************************************************

      // Compute the Ritz values in alpha[nlocked]..alpha[krdim-1]
      // and residuals res[nlocked]..res[krdim-1].  The eigenvectors "y"
      // in the Krylov space are also computed.  If more than a small
      // fraction of eigenvectors are needed, then computing all eigenvectors
      // by the QR algorithm is more efficient.  We assume that this is
      // usually the case.  Besides, manipulation of the large vectors
      // is much more costly.

bool KSRLanczosEigenSolver::compute_Ritz_values_eigvecs()
{
 int nkr=krdim-nlocked;
 int arrowheadsize=max(nkeep-nlocked+1,2);
// QDPIO::cout << "arrowheadsize = "<<arrowheadsize<<endl;

 char jobz = 'V';   // solve for eigenvalue and eigenvectors
 char uplo = 'U';   // must specify 'U' and NOT 'L'

          // put upper right of arrowhead matrix into "y"
          // use Fortran convention:  row + nkr*col  (for zero offset)
          // (if highest-lying eigenvalues wanted, multiply matrix by -1)

 if (!ascending){
    for (int k=nlocked;k<krdim-1;k++){
       alpha[k]=-alpha[k]; beta[k]=-beta[k];}
    alpha[krdim-1]=-alpha[krdim-1];}

 if (int(y.size())<nkr*nkr) y.resize(nkr*nkr);
 y.assign(nkr*nkr,0.0);
 for (int k=0;k<nkr;k++)
    y[k*(nkr+1)]=alpha[k+nlocked];
 for (int k=0;k<arrowheadsize-1;k++)
    y[k+nkr*(arrowheadsize-1)]=beta[k+nlocked];
 for (int k=arrowheadsize-1;k<nkr-1;k++)
    y[k+nkr*(k+1)]=beta[k+nlocked];

// q.resize(krdim+2);
// for (int row=0;row<nkr;row++)
// for (int col=row;col<nkr;col++){
//    q.assignMatrixVectorMultiply(krdim+1,col+nlocked);
//    QDPIO::cout << "T["<<row+1<<","<<col+1<<"] = "
//        <<y[row+nkr*col]<<"  dot product = "
//        <<q.InnerProduct(row+nlocked,krdim+1)<<endl;}
// q.resize(krdim+1);

//QDPIO::cout << "beta[krdim-1] = "<<beta[krdim-1]<<endl;

 int lwork = max(10,10*nkr);
 vector<double> work(lwork);
 int info;

     // diagonalize:  eigenvectors returned in "y", eigenvalues in "alpha"

 dsyev_(&jobz,&uplo,&nkr,&y[0],&nkr,&alpha[nlocked],&work[0],&lwork,&info);
 if (info) return false;

     // compute residuals

 for (int k=0;k<nkr;k++)
    residual[k+nlocked]=abs(beta[krdim-1]*y[nkr*(k+1)-1]);

     // undo negative signs in eigenvalues if highest-lying requested
 if (!ascending){
    for (int k=nlocked;k<krdim;k++)
       alpha[k]=-alpha[k];}

// for (int k=0;k<krdim;k++) QDPIO::cout <<"alpha["<<k<<"] = "<<alpha[k]<<endl;

 return true;
}

// ********************************************************************

  //  The method for choosing the number of Ritz vectors to keep
  //  is Equation 6 in  Kesheng Wu and Horst Simon, 
  //  "Dynamic Restarting Schemes For Eigenvalue Problems", 
  //  LBNL report number LBNL-42982.  However, the number must be
  //  at least as large as itconv, and it cannot be bigger than
  //  krdim-nlocked-6.
 
// Equation 5 is   itconv + (krdim - nconverged)/2
// Equation 7 is   itconv + ned 

int KSRLanczosEigenSolver::set_number_to_keep(int itconv)
{
/*
 int nkr=krdim-nlocked;
 int itned = ned-nlocked;
 double tmp = (0.4+itned/(10.0*nkr));
 tmp*=nkr;
 int itkeep = min(int(tmp), nkr-12);
 return max(itconv,itkeep); */

// double tmp = 0.6*double(ned-nlocked);
// return min(itconv+int(tmp),krdim-nlocked-12);

// return min(itconv+ned-nlocked,krdim-nlocked-12);

 return min(itconv+(krdim-nconverged)/2,krdim-nlocked-12);

}

// ********************************************************************

   // evaluates the Ritz vectors to be kept, then update the beta values

void KSRLanczosEigenSolver::compute_kept_Ritz_vectors(int itkeep)
{
 int tmpoffset=krdim+1;
 int nkr=krdim-nlocked;

   // ensure adequate storage space for new Ritz vectors
 if (q.getNumberOfVectors() < tmpoffset+itkeep) q.resize(tmpoffset+itkeep);

 for (int k=0;k<itkeep;k++){
    int kk=tmpoffset+k;
    q.copyTo(kk,y[nkr*k],nlocked);
    for (int j=1;j<nkr;j++) q.addTo(kk,y[k*nkr+j],j+nlocked);
    }

 for (int k=0;k<itkeep;k++)
    q.copyTo(k+nlocked,tmpoffset+k);
 q.copyTo(nlocked+itkeep,krdim);

 for (int k=0;k<itkeep;k++)
    beta[k+nlocked]=beta[krdim-1]*y[nkr*(k+1)-1];

}

// ****************************************************************************

   //  locking could have caused eigenvalue order to deviate from
   //  ascending or descending ... this routine corrects for this
   //  by doing a gnome-sort

void KSRLanczosEigenSolver::reorder()
{
 if (ascending){
    int i=0;
    while (i<krdim){
       if ((i==0)||(alpha[i-1]<=alpha[i])) i++;
       else{
          double tmp = alpha[i];        
          alpha[i] = alpha[i-1]; 
          alpha[--i] = tmp;
          q.swap(i,i-1);}
          }
    }
 else{
    int i=0;
    while (i<krdim){
       if ((i==0)||(alpha[i-1]>=alpha[i])) i++;
       else{
          double tmp = alpha[i];        
          alpha[i] = alpha[i-1]; 
          alpha[--i] = tmp;
          q.swap(i,i-1);}
          }
    }
}

// **********************************************************
  }
}
 
