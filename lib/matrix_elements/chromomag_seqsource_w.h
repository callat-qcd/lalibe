// -*- C++ -*-
/*!
 *  \brief Apply the ChromoMagnetic Operator to the left of a propagator as the source for seqsource prop.
 *  If there are questions or bugs, blame David Brantley.
 */

#ifndef __chromomag_seqsource_w_h__
#define __chromomag_seqsource_w_h__

#include "chromabase.h"
#include "tr_less_fields.h"

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
              const multi1d<LatticeColorMatrix>& u);
}  // end namespace Chroma

#endif
