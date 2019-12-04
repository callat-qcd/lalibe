//Arjun Singh Gambhir
//Real authors: Ben Hoerz, John Bulava, Colin Morningstar, ??
//Code directly ported from chroma_laph, more acknowledgements coming later.
#ifndef LAPH_EIGEN_INFO_H
#define LAPH_EIGEN_INFO_H

#include "qdp.h"
#include "chromabase.h"
#include "xml_help.h"

namespace Chroma {
  namespace LaphEnv {


// *******************************************************************
// *                                                                 *
// *   Objects of class "LaphEigenSolverInfo" store identifying info *
// *   related to solving for the Laplacian eigenvectors using the   *
// *   restarted Lanczos method.  The XML input must have the format *
// *                                                                 *
// *   <LaphEigenSolverInfo>                                         *
// *      <ResidualTolerance> 1e-8 </ResidualTolerance>              *
// *      <MaxIterations> 200 </MaxIterations>                       *
// *      <KrylovDimension> 90 </KrylovDimension>                    *
// *      <ChebyshevOrder> 8 </ChebyshevOrder> (can be omitted)      *
// *      <MaxEigenvalue> 15.0 </MaxEigenvalue> (can be omitted)     *
// *      <CutoffEigenvalue> 3.4 </CutoffEigenvalue> (optional)      *
// *      <StartingVectorType>equal_components</StartingVectorType>  *
// *      <OutputVerbosity>0</OutputVerbosity>                       *
// *   </LaphEigenSolverInfo>                                        *
// *                                                                 *
// *   Note that the number of eigenvectors to solve for is not      *
// *   given in this class; it is specified in QuarkSmearingInfo.    *
// *                                                                 *
// *   The matrix whose eigenvectors are sought is denoted by "A".   *
// *   "A" is minus one times the smeared covariant Laplacian:       *
// *             A = -Delta                                          *
// *   The eigenvalues of "A" lie between 0 and some maximum value   *
// *   denoted by "largest_eigenvalue" or "L" for short.  All of     *
// *   the eigenvalues of "A" are real and positive.  We wish to     *
// *   determine the eigenvectors corresponding to the lowest-lying  *
// *   eigenvalues lying between 0 and "cutoff_eigenvalue" or "C"    *
// *   for short.  So we have                                        *
// *                                                                 *
// *     desired part of spectrum:   0 ... "C"                       *
// *    unwanted part of spectrum:  "C" .. "L"                       *
// *                                                                 *
// *   The rate of convergence to solution depends on the spacing    *
// *   between the levels.  Convergence is much faster for widely    *
// *   spaced levels.  So convergence can be accelerated by          *
// *   transforming the spectrum so that the desired part of the     *
// *   spectrum is more widely spaced.  The following                *
// *   transformation is applied first:                              *
// *                                                                 *
// *          1 - 2*(A-C)/(L-C)                                      *
// *                                                                 *
// *   The above transformation maps the unwanted spectrum to the    *
// *   range -1 .. 1, and the desired part lies above 1.             *
// *   Then Chebyshev polynomials can be applied.  Eigenvalues       * 
// *   lying between -1 and 1 are suppressed (stay between -1..1),   *
// *   and the desired eigenvalues above 1 get spaced out to         *
// *   large values above 1 to speed up convergence.                 *
// *   The lowest-lying eigenvalue becomes the highest-lying         *
// *   transformed eigenvalue.  (Transforming the desired levels     *
// *   to the region above 1 is more convenient since it allows      *
// *   the use of Chebyshev polynomials of any order, both even      *
// *   and odd.)                                                     *
// *                                                                 *
// *   To apply the Chebyshev acceleration, a fairly reasonable      * 
// *   estimate of the large eigenvalue of -Delta is needed.         *
// *                                                                 *
// *   Setting "ChebyshevOrder" to a value of 1 or less means that   *
// *   Chebyshev polynomial acceleration is not applied. In this     *
// *   case, "MaxEigenvalue" and "CutoffEigenvalue" are ignored.     *
// *                                                                 *
// *   If ChebyshevOrder is used, you MUST specify the               *
// *   "MaxEigenvalue" and "CutoffEigenvalue".  "MaxEigenvalue"      *
// *   should be greater than the maximum eigenvalue of -Delta       *
// *   (which is usually around 12.0) and  "CutoffEigenvalue" should *
// *   be just above the highest eigenvalue of -Delta that you wish  *
// *   to compute the corresponding eigenvectors.                    *
// *                                                                 *
// *   "StartingVectorType" is optional: its default value is        *
// *   "equal_components", although "random" is another valid value. *
// *                                                                 *
// *   "OutputVerbosity" is optional: it must have value 0,1,2.      *
// *   Higher numbers lead to more output from the solver.           *
// *                                                                 *
// *                                                                 *
// *******************************************************************


class LaphEigenSolverInfo
{

  int maxIterations;
  int dimKrylov;
  double tolerance;
  int chebyshevOrder;
  double maxEigenvalue;
  double cutoffEigenvalue;
  std::string startVector;
  int outputVerbosity;


 public:  

  LaphEigenSolverInfo(XmlReader& xml_in);

  LaphEigenSolverInfo(const LaphEigenSolverInfo& in);

  LaphEigenSolverInfo& operator=(const LaphEigenSolverInfo& in);

  ~LaphEigenSolverInfo(){}



    // output functions

  int getMaximumIterations() const { return maxIterations; }

  int getKrylovDimension() const {return dimKrylov;}
  
  double getResidualTolerance() const { return tolerance; }

  int getChebyshevOrder() const { return chebyshevOrder; }

  double getMaximumEigenvalue() const { return maxEigenvalue; }

  double getCutoffEigenvalue() const { return cutoffEigenvalue; }

  std::string getStartingVectorType() const { return startVector; }

  int getOutputVerbosity() const { return outputVerbosity; }

  std::string output(int indent = 0) const;

  void output(XmlWriter& xmlout) const;

 private:

  void extract_info_from_reader(XmlReader& xmlr);


};



// **************************************************
  }
}
#endif
