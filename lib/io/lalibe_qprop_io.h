// -*- C++ -*-
/*! \file
 * Taken from Chroma so we don't have to hack or fake record xmls for new stuff in lalibe.
 * Arjun Singh Gambhir
 */

#ifndef __lalibe_qprop_io_h__
#define __lalibe_qprop_io_h__

#include "io/xml_group_reader.h"
#include "io/enum_io/enum_qdpvolfmt_io.h"
#include "io/enum_io/enum_quarkspintype_io.h"
//We still use this for almost everything.
#include "io/qprop_io.h"

namespace Chroma
{
    //! Sequential source parameters
    struct LalibeSeqSource_t
    {
        LalibeSeqSource_t();        /*!< default constructor */

        GroupXML_t    seqsrc;    /*!< Sequential source xml */

        multi1d<int>  sink_mom;  /*!< sink momentum */
        int           j_decay;   /*!< decay direction */
        bool          t_all;
        int           t_0;       // source time
        int           t_sink;    /*!< time slice of sink */
        int           t_sep;     // sink-src separation time
        //Extra stuff from lalibe we need to read/write
        std::string particle;       //which particle.
        std::string flavor;         //which flavor insertion.
        std::string source_spin;    //Source spin state.
        std::string sink_spin;      //Sink spin state.
    };


  //! LalibeSeqSource header read
  void read(XMLReader& xml, const std::string& path, LalibeSeqSource_t& header);

  //! LalibeSeqSource header writer
  void write(XMLWriter& xml, const std::string& path, const LalibeSeqSource_t& header);

}  // end namespace Chroma

#endif
