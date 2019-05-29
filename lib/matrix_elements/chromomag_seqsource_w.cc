// -*- C++ -*-
/*!
 *  \brief Apply the ChromoMagnetic Operator to the left of a propagator as the source for seqsource prop.
 *  If there are questions or bugs, blame David Brantley.
 */

#include "chromabase.h"
#include "tr_less_fields.h"
#include "chromomag_seqsource_w.h"
namespace Chroma
{

  //! Construct the Chromomagnetic operator, and then contract the result onto a prop.
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions!
   *
   * The resulting contracted propagator is used as a sequential prop source.
   *
   *
   * \param source_propagator        The input source prop      ( Read )
   * \param u                        gauge field                ( Read )
   */


  LatticePropagator chromoMagneticSeqSource(const LatticePropagator& source_propagator,
              const multi1d<LatticeColorMatrix>& u)
              {
    START_CODE();


    /* Get the Field Strength Tensor from mesfield */
    // mesfield returns the upper triangular elements of F(mu,nu)
    // totaling Nd*(Nd-1) elements, where Nd = number spacetime dimensions.
    // F(mu,nu) = 0 when mu = nu. The upper Nd*(Nd-1)/2 components is the
    // upper triangular portion, while the lower Nd*(Nd-1)/2 are the lower
    // triangular elements. As F(mu,nu) = -F(nu,mu), these are related.

    multi1d<LatticeColorMatrix> F = trLessFieldST(u);

    /* Iterate over the spacetime indices */
    // The Combination Sigma(mu,nu)*F(mu,nu) is symmetric in (mu,nu), so
    // we need only double the result of the upper triangle to
    // get the full operator.

    /* Get the equivalent Gamma((mu,nu)). Only works for Nd = 4 and Ns = 4*/
    // Manually get the Gamma for the right combos of (mu,nu).
    // ONLY FOR Nd = 4!

    int int_converter[Nd*(Nd-1)] = {3,5,9,6,10,12,3,5,6,9,10,12};

    // Declare summand LatticePropagator
    LatticePropagator contracted_prop=zero;

    // Iterate through the upper triangular components.
    for (int iter = 0; iter < Nd*(Nd - 1)/2; iter++)
      {

       LatticePropagator tmp1 = timesI(F[iter]*source_propagator);

       contracted_prop += Gamma(int_converter[iter])*tmp1;

      }

    Real neg = -1;
    // Iterate through the lower triangular components.
    for (int iter = Nd*(Nd - 1)/2; iter < Nd*(Nd - 1); iter++)
      {

       LatticePropagator tmp1 = timesI(F[iter]*source_propagator);

       LatticePropagator tmp2 = Gamma(int_converter[iter])*tmp1;
       contracted_prop += neg*tmp2;

      }


    return contracted_prop;
    END_CODE();
    }
}  // end namespace Chroma
