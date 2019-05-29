// -*- C++ -*-
/*! \file
 * Taken from Chroma so we don't have to hack or fake record xmls for new stuff in lalibe.
 * Arjun Singh Gambhir
 */

#include "chromabase.h"
#include "io/param_io.h"
#include "io/qprop_io.h"
//Lalibe file below.
#include "lalibe_qprop_io.h"

#include "meas/smear/simple_quark_displacement.h"
#include "meas/smear/ape_link_smearing.h"

namespace Chroma
{
    // Initialize header with default values
    LalibeSeqSource_t::LalibeSeqSource_t()
    {
        j_decay   = -1;
        t_sink    = -1;
        sink_mom.resize(Nd-1);
        sink_mom = 0;
    }


    // Anonymous namespace
    // I don't want an anonymous namespace for lalibe's reads/writes.
    //namespace
    //{

    //! LalibeSeqSource header reader
    void read(XMLReader& xml, const std::string& path, LalibeSeqSource_t& param)
    {
        XMLReader paramtop(xml, path);

        int version;
        read(paramtop, "version", version);

        switch (version)
        {
            case 1:
            {
                XMLReader xml_readback;
                {
                    XMLBufferWriter xml_tmp;
                	push(xml_tmp, "Param");
                	write(xml_tmp, "version", 2);
                	push(xml_tmp, "SeqSource");
                	write(xml_tmp, "version", 1);
                	std::string seq_src;
                	read(paramtop, "seq_src", seq_src);
                	write(xml_tmp, "SeqSourceType", seq_src);

                	read(paramtop, "sink_mom", param.sink_mom);
                	write(xml_tmp, "sink_mom",  param.sink_mom);

                	read(paramtop, "t_sink", param.t_sink);
                	write(xml_tmp, "t_sink",  param.t_sink);
                	//Here's where the new stuff starts.
                	read(paramtop, "SeqSourceType", param.particle);
                	write(xml_tmp, "SeqSourceType", param.particle);
                	read(paramtop, "flavor", param.flavor);
                	write(xml_tmp, "flavor", param.flavor);
                	read(paramtop, "source_spin", param.source_spin);
                	write(xml_tmp, "source_spin", param.source_spin);
                	read(paramtop, "sink_spin", param.sink_spin);
                	write(xml_tmp, "sink_spin", param.sink_spin);

                	write(xml_tmp, "j_decay",  Nd-1);

                	pop(xml_tmp);  // Source
                	pop(xml_tmp);  // Param

	                QDPIO::cout << "seqsrc_xml = XX" << xml_tmp.printCurrentContext() << "XX" << std::endl;

	                xml_readback.open(xml_tmp);
                }

                // Recurse back in to re-read
                //read(xml_readback, "/Param", param);
                //Comment above line due to ambigious overloading.
            }
            break;

            case 2:
            {
                param.seqsrc = readXMLGroup(paramtop, "SeqSource", "SeqSourceType");

                XMLReader xml_tmp(paramtop, "SeqSource");
                read(xml_tmp, "j_decay",       param.j_decay);
                read(xml_tmp, "t_all",         param.t_all);
                read(xml_tmp, "t_0",           param.t_0);
                read(xml_tmp, "t_sink",        param.t_sink);
                read(xml_tmp, "t_sep",         param.t_sep);
                read(xml_tmp, "sink_mom",      param.sink_mom);
                //Here's where the new stuff starts.
                read(xml_tmp, "SeqSourceType", param.particle);
                read(xml_tmp, "flavor",        param.flavor);
                read(xml_tmp, "source_spin",   param.source_spin);
                read(xml_tmp, "sink_spin",     param.sink_spin);
            }
            break;

            default:
            QDPIO::cerr << "SeqSource parameter version " << version
		          << " unsupported." << std::endl;
            QDP_abort(1);
        }

    }

    //! LalibeSeqSource header writer
    void write(XMLWriter& xml, const std::string& path, const LalibeSeqSource_t& param)
    {
        push(xml, path);

        int version = 2;
        write(xml, "version", version);
        xml << param.seqsrc.xml;

        pop(xml);
    }
}  // end namespace Chroma
