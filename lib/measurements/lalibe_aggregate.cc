
// I guess we need fermacts to do physics.
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

// INTESTING PHYSICS
#include "lalibe_aggregate.h"
#include "moments_fh_prop_w.h"
#include "stochastic_fh_prop_w.h"
#include "HP_fh_prop_w.h"
#include "stochastic_four_quark_fh_prop_w.h"
#include "baryon_contractions_w.h"
#include "meson_contractions_w.h"
#include "flavor_conserving_fh_baryon_contractions_w.h"
#include "flavor_changing_fh_baryon_contractions_w.h"
#include "fh_prop_w.h"
#include "ZN_prop_w.h"
#include "HP_prop_w.h"
#include "lalibe_bar3ptfn_w.h"
#include "lalibe_seqsource_w.h"
#include "coherent_seqsource_w.h"
#include "pipi_scattering_w.h"
#include "qedm_product_field.h"

// UTILITIES
#include "multi_prop_add.h"
// Protect everything with a preprocessing directive.
#ifdef BUILD_HDF5
#include "hdf5_read_obj.h"
#include "hdf5_write_obj.h"
#include "hdf5_write_erase_obj.h"
#endif

namespace Chroma
{

  //! Name and registration
  namespace LalibeAggregateEnv
  {
    namespace
    {
      //! Local registration flag
      bool registered = false;
    }

    //! Register all the factories
    bool registerAll()
    {
      bool success = true;
      if (! registered)
      {
	// Fermact stuff from chroma
	success &= WilsonTypeFermActsEnv::registerAll();

  // INTERESTING PHYSICS
	success &= LalibeStochasticFHPropagatorEnv::registerAll() ;
	success &= LalibeHPFHPropagatorEnv::registerAll() ;
	success &= LalibeStochasticFourQuarkFHPropagatorEnv::registerAll() ;
	success &= LalibeMomentsFHPropagatorEnv::registerAll() ;
	success &= LalibeFlavorConservingFHBaryonContractionsEnv::registerAll() ;
	success &= LalibeFlavorChangingFHBaryonContractionsEnv::registerAll() ;
	success &= LalibeBaryonContractionsEnv::registerAll() ;
  success &= LalibeMesonContractionsEnv::registerAll() ;
	success &= LalibeFHPropagatorEnv::registerAll() ;
	success &= LalibeZNPropagatorEnv::registerAll() ;
	success &= LalibeHPPropagatorEnv::registerAll() ;
	success &= LalibeBar3ptfnEnv::registerAll() ;
  success &= LalibeSeqSourceEnv::registerAll() ;
  success &= LalibeCoherentSeqsourceEnv::registerAll() ;
  success &= LalibePipiScatteringEnv::registerAll() ;
  success &= LalibeQEDMProductFieldEnv::registerAll() ;

  // USEFUL UTILITIES
  success &= LalibeMultiPropagatorAddEnv::registerAll() ;

#ifdef BUILD_HDF5
	success &= LalibeHDF5ReadNamedObjEnv::registerAll();
	success &= LalibeHDF5WriteNamedObjEnv::registerAll();
	success &= LalibeHDF5WriteEraseNamedObjEnv::registerAll();
#endif

	registered = true;
      }
      return success;
    }
  }

}
