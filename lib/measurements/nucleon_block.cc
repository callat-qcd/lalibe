/*
  Authors
  Andre Walker-Loud

  re-factoring of Arjun's port of latscat to lalibe
*/

// Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
#include "fermact.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "util/info/unique_id.h"
#include "io/xml_group_reader.h"
#include "util/ferm/paulitodr.h"

// Lalibe Stuff
#include "../momentum/lalibe_sftmom.h"
#include "nucleon_block.h"
//#include "../NN/blockstuff.h"
#include "../matrix_elements/bilinear_gamma.h"
//TODO: bilinear_gamma above shouldn't be needed once boiler plate is removed.
#include "../contractions/baryon_contractions_func_w.h"

// Latscat Stuff
#ifndef CUFFT
#include "../NN/fourier_cpu.h"
#else
#include "../NN/fourier_cuda.h"
#endif
#include "../NN/NN_LC_w.h"
#include "../NN/spinstuff.h"
#include "../NN/transform.h"
#include "../NN/checkpointstuff.h"
#include "../NN/momstuff.h"

namespace Chroma
{
    namespace LalibeNucleonBlockEnv
    {
        namespace
        {
            AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, const std::string& path)
            {
                return new InlineMeas(NucleonBlockParams(xml_in, path));
            }
            //! Local registration flag
            bool registered = false;
        }
        const std::string name = "NUCLEON_BLOCK";

        //! Register all the factories
        bool registerAll()
        {
            bool success = true;
            if (! registered)
                {
                    success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
                    registered = true;
                }
            return success;
        }

        void read(XMLReader& xml, const std::string& path, NucleonBlockParams::NucleonBlock_t& par)
        {
            XMLReader paramtop(xml, path);
            #if !defined(BUILD_FFTW) && !defined(CUFFT)
                QDPIO::cerr << "\n###############################################################\n" << std::endl;
                QDPIO::cerr << "  You must link against FFTW or CUFFT for this NN measurement" << std::endl;
                QDPIO::cerr << "\n###############################################################\n" << std::endl;
                QDP_abort(1);
            #endif
            // FFT info
            if (paramtop.count("fft_chunksize") != 0)
                read(paramtop, "fft_chunksize", par.fft_chunksize); //originally the only parameter in FFTPar struct
            else
                par.fft_chunksize = 0;
            if (paramtop.count("fft_tune") != 0)
                read(paramtop, "fft_tune", par.fft_tune); //tune fft?
            else
                par.fft_tune = false;
            read(paramtop, "compute_locals", par.compute_locals);
            read(paramtop, "negative_parity", par.negative_parity);
            // list of total momentum boosts
            read(paramtop, "boosts", par.boosts);
            // For now - only support P_tot = (0,0,0)
            if (par.boosts.size() > 1 || par.boosts[0][0] != 0 || par.boosts[0][1] != 0 || par.boosts[0][2] != 0){
                QDPIO::cerr << "\n###############################################################\n" << std::endl;
                QDPIO::cerr << "  For now, only a single boost of (0,0,0) is allowed" << std::endl;
                QDPIO::cerr << "\n###############################################################\n" << std::endl;
                QDP_abort(1);
            }

            // Are propagators already in Dirac basis?
            if (paramtop.count("in_dirac_basis") != 0)
                read(paramtop, "in_dirac_basis", par.in_dirac_basis); //specifies props in dirac basis, this is false by default
            else
                par.in_dirac_basis = false;
        }

        void write(XMLWriter& xml, const std::string& path, NucleonBlockParams::NucleonBlock_t& par)
        {
            push(xml, path);
            write(xml, "fft_chunksize",   par.fft_chunksize);
            write(xml, "fft_tune",        par.fft_tune);
            write(xml, "compute_locals",  par.compute_locals);
            write(xml, "boosts",          par.boosts);
            write(xml, "in_dirac_basis",  par.in_dirac_basis);
            write(xml, "negative_parity", par.negative_parity);
            pop(xml);
        }

        //! NamedObject input
        void read(XMLReader& xml, const std::string& path, NucleonBlockParams::NamedObject_t& input)
        {
            XMLReader inputtop(xml, path);
            read(inputtop, "gauge_id"     , input.gauge_id);
            read(inputtop, "prop0_id"     , input.prop0_id);
            read(inputtop, "prop1_id"     , input.prop1_id);
            read(inputtop, "block_map"    , input.block_map);
        }

        //! NamedObject output
        void write(XMLWriter& xml, const std::string& path, const NucleonBlockParams::NamedObject_t& input)
        {
            push(xml, path);
            write(xml, "gauge_id"     , input.gauge_id );
            write(xml, "prop0_id"     , input.prop0_id);
            write(xml, "prop1_id"     , input.prop1_id);
            write(xml, "block_map"    , input.block_map);
            pop(xml);
        }

        // Create string for total boost
        std::string boost_string(const multi1d<int> &boost){
            std::string boostdir("boost_");
            for(unsigned int d=0; d<(Nd-1); d++){
                if(boost[d]>=0) boostdir+="p";
                else boostdir+="m";
                boostdir+=dirlist[d]+std::to_string(abs(boost[d]));
            }
            if(boost.size() > (Nd-1) && boost[Nd-1] !=0){
                if(boost[Nd-1]>=0) boostdir+="p";
                else boostdir+="m";
                boostdir+=dirlist[Nd-1]+std::to_string(abs(boost[Nd-1]));
            }
            return boostdir;
        }

        // timestamp for simple checkpointing
        const std::string getTimestamp() {
            time_t     now = time(0);
            struct tm  tstruct;
            char       buf[80];
            tstruct = *localtime(&now);
            strftime(buf, sizeof(buf), "%Y-%m-%d-%X", &tstruct);
            std::string result(buf);
            result.erase(std::remove(result.begin(), result.end(), ':'), result.end());

            return result;
        }

        // Param stuff
        NucleonBlockParams::NucleonBlockParams()
        {
            frequency = 0;
        }

        NucleonBlockParams::NucleonBlockParams(XMLReader& xml_in, const std::string& path)
        {
            try
                {
                    XMLReader paramtop(xml_in, path);
                    if (paramtop.count("Frequency") == 1)
                        read(paramtop, "Frequency", frequency);
                    else
                        frequency = 1;

                    // Parameters for source construction
                    read(paramtop, "NucleonBlockParams", nblockparam);

                    // Read in the NamedObject info
                    read(paramtop, "NamedObject", named_obj);
                }
            catch(const std::string& e)
                {
                    QDPIO::cerr << __func__ << ": Caught Exception reading XML: "
                                << e << std::endl;
                    QDP_abort(1);
                }
        }

        void NucleonBlockParams::writeXML(XMLWriter& xml_out, const std::string& path)
        {
            push(xml_out, path);
            write(xml_out, "NucleonBlockParams", nblockparam);
            write(xml_out, "NamedObject"       , named_obj);
            pop(xml_out);
        }

        // Function call
        void InlineMeas::operator()(unsigned long update_no, XMLWriter& xml_out)
        {
            START_CODE();

            StopWatch snoop;
            snoop.reset();
            snoop.start();
            QDPIO::cout << "NUCLEON_BLOCK: start" << std::endl;

            // Test and grab a reference to the gauge field
            XMLBufferWriter gauge_xml;
            try
                {
                    TheNamedObjMap::Instance().getData <multi1d <LatticeColorMatrix> >(params.named_obj.gauge_id);
                    TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
                }
            catch( std::bad_cast )
                {
                    QDPIO::cerr << LalibeNucleonBlockEnv::name
                                << ": caught dynamic cast error" << std::endl;
                    QDP_abort(1);
                }
            catch (const std::string& e)
                {
                    QDPIO::cerr << LalibeNucleonBlockEnv::name
                                << ": map call failed: " << e << std::endl;
                    QDP_abort(1);
                }
            const multi1d<LatticeColorMatrix>& u =
                TheNamedObjMap::Instance().getData <multi1d <LatticeColorMatrix> >(params.named_obj.gauge_id);

            //Code from latscat ported below, with comments for missing or adjusted features.
            //Chroma initialize stuff was here, not needed in lalibe.
            //Latscat has an option for its own singnal handlers that goes here, not being ported.
#ifdef DEBUG
            QDPIO::cout << "Warning, DEBUG mode enabled!" << std::endl;
#endif
            //Stuff about LatticePars is read next, this is already read in our main function.
            //For now, this measurement is only going to crunch stuff if we have built with hdf5.
#ifdef BUILD_HDF5

            //set up fourier stuff:
            QDPIO::cout << "Setting up Communicators for FFT..." << std::flush;

            int j_decay = Nd - 1; // Assume t_dir = j_decay

#ifdef BUILD_FFTW
            Fourier fft(j_decay);
            Fourier fftblock(j_decay);
#ifdef PROFILE
            fft.print_flops(true);
            fftblock.print_flops(true);
#endif
            QDPIO::cout << "done!" << std::endl;

            // Note, the truncation size (2nd arg) is not used in this function anymore.  -1 --> do not truncate
            int truncate = -1;
            // comment out while we only make blocks
            //initTopologies(params.nblockparam.contractions_filename, truncate, j_decay);

            if(params.nblockparam.fft_chunksize !=0 )
                QDPIO::cout << "Using chunksize " << params.nblockparam.fft_chunksize << " for the Baryon-Block FFT!" << std::endl;

            //Do FFT tuning if it's enabled; comment below comes from latscat.
            if(params.nblockparam.fft_tune){
                QDPIO::cout << "Tuning FFT for better performance..." << std::flush;
                StopWatch swatch_fftune;
                swatch_fftune.reset();
                swatch_fftune.start();
                fftblock.tune(sizeof(HalfBaryonblock),true);
                swatch_fftune.stop();
                QDPIO::cout << "done! Time " << swatch_fftune.getTimeInSeconds() << std::endl;
            }

            // Get the source locations for prop_0 and prop_1
            multi1d<int> pos0(Nd), pos1(Nd), disp(Nd);

            //Read/set up prop 0.
            XMLReader prop0_file_xml, prop0_record_xml;
            LatticePropagator prop0;
            QDPIO::cout << "Attempt to read propagator 0" << std::endl;
            try {
                prop0 = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop0_id);
                TheNamedObjMap::Instance().get(params.named_obj.prop0_id).getFileXML(prop0_file_xml);
                TheNamedObjMap::Instance().get(params.named_obj.prop0_id).getRecordXML(prop0_record_xml);
                MakeSourceProp_t  orig_header;
                if (prop0_record_xml.count("/Propagator") != 0){
                    read(prop0_record_xml, "/Propagator", orig_header);
                }
                // Or if we pass a smeared propagator
                else if (prop0_record_xml.count("/SinkSmear") != 0){
                    read(prop0_record_xml, "/SinkSmear", orig_header);
                }
                pos0 = orig_header.source_header.getTSrce();
            }
            catch (std::bad_cast) {
                QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
                QDP_abort(1);
            }
            catch (const std::string& e) {
                QDPIO::cerr << name << ": error reading src prop_header: "
                            << e << std::endl;
                QDP_abort(1);
            }
            //Now read/set up prop 1.
            XMLReader prop1_file_xml, prop1_record_xml;
            LatticePropagator prop1;
            QDPIO::cout << "Attempt to read propagator 1" << std::endl;
            try {
                prop1 = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop1_id);
                TheNamedObjMap::Instance().get(params.named_obj.prop1_id).getFileXML(prop1_file_xml);
                TheNamedObjMap::Instance().get(params.named_obj.prop1_id).getRecordXML(prop1_record_xml);
                MakeSourceProp_t  orig_header;
                if (prop1_record_xml.count("/Propagator") != 0){
                    read(prop1_record_xml, "/Propagator", orig_header);
                }
                else if (prop1_record_xml.count("/SinkSmear") != 0){
                    read(prop1_record_xml, "/SinkSmear", orig_header);
                }
                pos1 = orig_header.source_header.getTSrce();
            }
            catch (std::bad_cast) {
                QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
                QDP_abort(1);
            }
            catch (const std::string& e) {
                QDPIO::cerr << name << ": error reading src prop_header: "
                            << e << std::endl;
                QDP_abort(1);
            }

            //Now back to latscat logic of calculating displacements and stuff.
            QDPIO::cout << "    Propagator 0: " << params.named_obj.prop0_id << std::endl;
            QDPIO::cout << "        Location: "; for(int i = 0; i<Nd; i++){ QDPIO::cout << pos0[i] << " " ;};
            QDPIO::cout << std::endl;

            QDPIO::cout << "    Propagator 1: " << params.named_obj.prop1_id << std::endl;
            QDPIO::cout << "        Location: "; for(int i = 0; i<Nd; i++){ QDPIO::cout << pos1[i] << " " ;};
            QDPIO::cout << std::endl;

            disp = pos1 - pos0;

            for(unsigned int d=0; d<Nd; d++){
                QDPIO::cout << disp[d] << "%" << Layout::lattSize()[d] ;
                disp[d] = (disp[d] + Layout::lattSize()[d]) % Layout::lattSize()[d];
                QDPIO::cout << " --> " << disp[d] ;
                if(disp[d] > Layout::lattSize()[d]/2){ disp[d] = disp[d] - Layout::lattSize()[d]; }
                QDPIO::cout << " --> " << disp[d] << std::endl;
            }
            QDPIO::cout << "         Displacement: "; for(int i = 0; i<Nd; i++){ QDPIO::cout << disp[i] << " " ;};
            QDPIO::cout << std::endl;
            QDPIO::cout << "    Real Displacement: "; for(int i = 0; i<Nd; i++){ QDPIO::cout << disp[i] << " " ;};
            QDPIO::cout << std::endl;

            // displacement string
            std::string displacedir;
            for(unsigned int d=0; d<Nd; d++){
                if(disp[d] >= 0) displacedir+="p";
                else displacedir+="m";
                displacedir+=dirlist[d]+std::to_string(abs(disp[d]));
            }

            multi1d<int> protonpos1(Nd), protonpos2(Nd);
            protonpos1 = pos0;  // I know.
            protonpos2 = pos1;  // I'm sorry.

            //set up map
            Timeshiftmap tshiftmap(protonpos1[j_decay],j_decay,Layout::lattSize()[j_decay]);

            //If the propagators aren't specified to already be in the Dirac basis, we rotate them now.
            if(params.nblockparam.dirac_basis == false)
            {
                QDPIO::cout << "Rotating the propagators to Dirac basis." << std::endl;
                rotate_to_Dirac_Basis(prop0);
                rotate_to_Dirac_Basis(prop1);
            }

            //Necessary Fields:
            multi1d<BaryOp> Nup=get_local_MA_single(0,"DP");
            /* comment out while just doing blocks
            LatticeComplex tmplatcomp_P, tmplatcomp_P_34, token;
            std::map<std::string,LatticeHalfSpinMatrix> tmplatmats, tmplatmats_34;
            for(unsigned int s=0; s<12; s++){
                tmplatmats[contterms[s]]=LatticeHalfSpinMatrix();
                tmplatmats_34[contterms[s]]=LatticeHalfSpinMatrix();
            }
            */

            bool neg_par = params.nblockparam.negative_parity;
            std::string neg_par_str = "POS_PARITY";
            if (neg_par){
                // The positive parity code works on the negative parity props
                // since in Dirac Pauli, g5 flips 1 <--> 3 and 2 <--> 4
                SpinMatrixD g5_Dirac = adj(PauliToDRMat()) * Gamma(15) * PauliToDRMat();
                prop0 = g5_Dirac * prop0 * g5_Dirac;
                prop1 = g5_Dirac * prop1 * g5_Dirac;
                neg_par_str = "NEG_PARITY";
            }

            /*
            Make blocks and FFT them and return a map of mom-space blocks
            so that we can add P and S sinks without needing a 3rd FFT

            make block 000[+1]
            if prop1 = prop0:
                make block 000[-1]
            else:
                make block 111[-1]
                if compute_locals:
                    make block 000[-1]
                    make block 111[+1]
                make blocks
                    001, 010, 100 [+1]
                    101, 110, 011 [-1]
            */
            QDPIO::cout << "Creating timers..." << std::flush;
            // TODO: make timings more meaningful?
            //timings:
            StopWatch swatch_pdir, swatch_ploc, swatch_pdisp, swatch_fft, swatch_bblock, swatch_io_read, swatch_io_write;
            swatch_io_read.reset();
            swatch_io_write.reset();
            swatch_pdir.reset();
            swatch_ploc.reset();
            swatch_pdisp.reset();
            swatch_fft.reset();
            swatch_bblock.reset();
            QDPIO::cout << "done!" << std::endl;

            double block_time = 0.;

            // Declare block map
            typedef std::tuple<std::string, std::string, std::string, int, bool> BlockMapKeyType;
            typedef std::map<BlockMapKeyType, LatticeHalfBaryonblock> BlockMapType;
            // Put it in the NamedObjectMap
            TheNamedObjMap::Instance().create<BlockMapType>(params.named_obj.block_map);
            BlockMapType& blockMap = TheNamedObjMap::Instance().getData<BlockMapType>(params.named_obj.block_map);

            // Create a local instance to store blocks while we make them
            LatticeHalfBaryonblock tmpBlock;

            // string ids
            std::string propId_0 = params.named_obj.prop0_id;
            std::string propId_1 = params.named_obj.prop1_id;

            // Block 000[+1]
            BlockMapKeyType theKey = {propId_0, propId_0, propId_0, 1, neg_par};
            std::string block_str = "000[+1]";
            if (blockMap.count(theKey) == 0){
                block_time += get_barblock(tmpBlock, prop0, prop0, prop0, Nup[0].get_gamma(0));
                swatch_fft.start();
                blockMap[theKey] = fft(tmpBlock, 1);
                swatch_fft.stop();
                QDPIO::cout << "created block " << block_str << " " << neg_par_str << std::endl;
                QDPIO::cout << LalibeNucleonBlockEnv::name << ": 000 FFT time " << swatch_fft.getTimeInSeconds() << std::endl;
            }
            else QDPIO::cout << "exists  block " << block_str << " " << neg_par_str << std::endl;

            // Did the user pass the same prop?
            if (propId_0 == propId_1){
                // 000[-1]
                theKey = {propId_0, propId_0, propId_0, -1, neg_par};
                block_str = "000[-1]";
                if (blockMap.count(theKey) == 0){
                    swatch_fft.start();
                    blockMap[theKey] = fft(tmpBlock, -1);
                    swatch_fft.stop();
                    QDPIO::cout << "created block " << block_str << " " << neg_par_str << std::endl;
                }
                else QDPIO::cout << "exists  block " << block_str << " " << neg_par_str << std::endl;
            }
            else{
                // 111[-1]
                theKey = {propId_1, propId_1, propId_1, -1, neg_par};
                block_str = "111[-1]";
                if (blockMap.count(theKey) == 0){
                    block_time += get_barblock(tmpBlock, prop1, prop1, prop1, Nup[0].get_gamma(0));
                    swatch_fft.start();
                    blockMap[theKey] = fft(tmpBlock, -1);
                    swatch_fft.stop();
                    QDPIO::cout << "created block " << block_str << " " << neg_par_str << std::endl;
                }
                else QDPIO::cout << "exists  block " << block_str << " " << neg_par_str << std::endl;

                // Do we want to compute the local NN?  If so, need opposite sign FFT
                if (params.nblockparam.compute_locals){
                    // 111[+1]
                    theKey = {propId_1, propId_1, propId_1, +1, neg_par};
                    block_str = "111[+1]";
                    if (blockMap.count(theKey) == 0){
                        swatch_fft.start();
                        blockMap[theKey] = fft(tmpBlock, +1);
                        swatch_fft.stop();
                        QDPIO::cout << "created block " << block_str << " " << neg_par_str << std::endl;
                    }
                    else QDPIO::cout << "exists  block " << block_str << " " << neg_par_str << std::endl;
                    // 000[-1]
                    theKey = {propId_0, propId_0, propId_0, -1, neg_par};
                    block_str = "000[-1]";
                    if (blockMap.count(theKey) == 0){
                        block_time += get_barblock(tmpBlock, prop0, prop0, prop0, Nup[0].get_gamma(0));
                        swatch_fft.start();
                        blockMap[theKey] = fft(tmpBlock, -1);
                        swatch_fft.stop();
                        QDPIO::cout << "created block " << block_str << " " << neg_par_str << std::endl;
                    }
                    else QDPIO::cout << "exists  block " << block_str << " " << neg_par_str << std::endl;
                }
                // Now make the mixed blocks
                // 001 [+1]
                theKey = {propId_0, propId_0, propId_1, 1, neg_par};
                block_str = "001[+1]";
                if (blockMap.count(theKey) == 0){
                    block_time += get_barblock(tmpBlock, prop0, prop0, prop1, Nup[0].get_gamma(0));
                    swatch_fft.start();
                    blockMap[theKey] = fft(tmpBlock, 1);
                    swatch_fft.stop();
                    QDPIO::cout << "created block " << block_str << " " << neg_par_str << std::endl;
                }
                else QDPIO::cout << "exists  block " << block_str << " " << neg_par_str << std::endl;
                // 010 [+1]
                theKey = {propId_0, propId_1, propId_0, 1, neg_par};
                block_str = "010[+1]";
                if (blockMap.count(theKey) == 0){
                    block_time += get_barblock(tmpBlock, prop0, prop1, prop0, Nup[0].get_gamma(0));
                    swatch_fft.start();
                    blockMap[theKey] = fft(tmpBlock, 1);
                    swatch_fft.stop();
                    QDPIO::cout << "created block " << block_str << " " << neg_par_str << std::endl;
                }
                else QDPIO::cout << "exists  block " << block_str << " " << neg_par_str << std::endl;
                // 100 [+1]
                theKey = {propId_1, propId_0, propId_0, 1, neg_par};
                block_str = "100[+1]";
                if (blockMap.count(theKey) == 0){
                    block_time += get_barblock(tmpBlock, prop1, prop0, prop0, Nup[0].get_gamma(0));
                    swatch_fft.start();
                    blockMap[theKey] = fft(tmpBlock, 1);
                    swatch_fft.stop();
                    QDPIO::cout << "created block " << block_str << " " << neg_par_str << std::endl;
                }
                else QDPIO::cout << "exists  block " << block_str << " " << neg_par_str << std::endl;
                // And the opposite sign FFT blocks
                // 011 [-1]
                theKey = {propId_0, propId_1, propId_1, -1, neg_par};
                block_str = "011[-1]";
                if (blockMap.count(theKey) == 0){
                    block_time += get_barblock(tmpBlock, prop0, prop1, prop1, Nup[0].get_gamma(0));
                    swatch_fft.start();
                    blockMap[theKey] = fft(tmpBlock, -1);
                    swatch_fft.stop();
                    QDPIO::cout << "created block " << block_str << " " << neg_par_str << std::endl;
                }
                else QDPIO::cout << "exists  block " << block_str << " " << neg_par_str << std::endl;
                // 101 [-1]
                theKey = {propId_1, propId_0, propId_1, -1, neg_par};
                block_str = "101[-1]";
                if (blockMap.count(theKey) == 0){
                    block_time += get_barblock(tmpBlock, prop1, prop0, prop1, Nup[0].get_gamma(0));
                    swatch_fft.start();
                    blockMap[theKey] = fft(tmpBlock, -1);
                    swatch_fft.stop();
                    QDPIO::cout << "created block " << block_str << " " << neg_par_str << std::endl;
                }
                else QDPIO::cout << "exists  block " << block_str << " " << neg_par_str << std::endl;
                // 110 [-1]
                theKey = {propId_1, propId_1, propId_0, -1, neg_par};
                block_str = "110[-1]";
                if (blockMap.count(theKey) == 0){
                    block_time += get_barblock(tmpBlock, prop1, prop1, prop0, Nup[0].get_gamma(0));
                    swatch_fft.start();
                    blockMap[theKey] = fft(tmpBlock, -1);
                    swatch_fft.stop();
                    QDPIO::cout << "created block " << block_str << " " << neg_par_str << std::endl;
                }
                else QDPIO::cout << "exists  block " << block_str << " " << neg_par_str << std::endl;
            }

            QDPIO::cout << LalibeNucleonBlockEnv::name << ": total block time " << block_time << std::endl;
            QDPIO::cout << LalibeNucleonBlockEnv::name << ": total FFT time " << swatch_fft.getTimeInSeconds() << std::endl;



#if 0

            /**********************************************************************
             *                                                                    *
             *   LOOP OVER BOOSTS                                                 *
             *        We don't support boosts now, but we might, so keep the loop *
             **********************************************************************/
            for(unsigned int b=0; b < params.nblockparam.boosts.size(); b++){
                multi1d<int> boost = params.nblockparam.boosts[b];
                std::string boostdir=boost_string(boost);

                momentum mom(params.nblockparam.boosts[b]);

                QDPIO::cout << mom << "    " << boostdir << std::endl;

                LatticeComplex phases=get_phases(mom,j_decay,+1);

                swatch_pdisp.start();
                QDPIO::cout << "Computing contractions..." << std::endl;
                // if displacement between sources of uprop_p{1,2} is nonzero:
                QDPIO::cout << "    Doing displaced only." << std::endl;
                QDPIO::cout << "    For now, if you want local, just pass the same propagator twice." << std::endl;
#if QDP_NC == 3
                contract(tmplatcomp_P, tmplatmats, prop_0, prop_1, Nup[0].get_gamma(0), phases, fftblock, false, weights);
#else
                QDPIO::cout<<"Contractions not yet implemented for NC!=3."<<std::endl;
#endif
                QDPIO::cout << "    Contractions done!" << std::endl;

                QDPIO::cout << "Computing contractions with negative parity blocks ..." << std::endl;
                contract(tmplatcomp_P_34, tmplatmats_34, prop_0_34, prop_1_34, Nup[0].get_gamma(0), phases, fftblock, false, weights);
                QDPIO::cout << "    done!" << std::endl;

                //QDPIO::cout << "Skipping the single proton." << std::endl;
                /*******************************************
                 * SINGLE PROTON
                 *******************************************/
                if(b==0){
                    QDPIO::cout << "Single proton." << std::endl;
                    //positive parity proton
                    token=zero;
                    token+=/* sourcefact * tmpprefact * */ tshiftmap(tmplatcomp_P,true);
                    swatch_io_write.start();
                    chk.set_parameter("proton1",token);
                    swatch_io_write.stop();

                    //negative parity proton
                    token=zero;
                    token+=/* sourcefact * tmpprefact * */ tshiftmap(tmplatcomp_P_34,true);
                    swatch_io_write.start();
                    chk.set_parameter("proton1_34",token);
                    swatch_io_write.stop();
                }

                //TODO: This print comes from latscat, but it also serves as a placeholder to extend this measurement and make it more general in a second or third pass.
                QDPIO::cout << "Skipping the local sources." << std::endl;

                for(std::map<std::string,LatticeHalfSpinMatrix>::iterator it=tmplatmats.begin(); it!=tmplatmats.end(); ++it){
                    std::string idstring=it->first;
                    if(idstring.find("loc")!=std::string::npos) continue;

                    size_t firstpos=idstring.find("_");
                    size_t secondpos=idstring.find(firstpos);
                    std::string conttype=idstring.substr(0,firstpos);
                    std::transform(conttype.begin(),conttype.end(),conttype.begin(),::tolower);

                    std::string srcspin=idstring.substr(firstpos+1,secondpos);
                    std::string srcspinval(&srcspin[srcspin.size()-1]);
                    srcspin=srcspin.substr(0,srcspin.size()-1);

                    for(std::map<std::string,HalfSpinMatrix>::iterator innerit=projectors.begin(); innerit!=projectors.end(); ++innerit){
                        std::string inneridstring=innerit->first;
                        std::string snkspin=inneridstring.substr(0,inneridstring.size()-1);
                        if(snkspin!=srcspin) continue;
                        std::string snkspinval(&inneridstring[inneridstring.size()-1]);

                        std::string corrname=conttype+"corr"+"_"+srcspin+"_"+snkspinval+"_"+srcspinval+"_"+displacedir; //s[0]; // 0 = mu hardcoded for this version.
                        token=zero;
                        // TODO src always is 0 for this version!!
                        // if(src==0) token=zero;
                        // else{
                        //     swatch_io_read.start();
                        //     chk.get_parameter(boostdir+"/"+corrname,token);
                        //     swatch_io_read.stop();
                        // }
                        token+=/* sourcefact * tmpprefact * */ tshiftmap(LatticeComplex(trace((innerit->second)*(it->second))));
                        swatch_io_write.start();
                        chk.set_parameter(boostdir+"/"+corrname,token);
                        swatch_io_write.stop();
                    }
                }

                // other-parity (34 entries)
                for(std::map<std::string,LatticeHalfSpinMatrix>::iterator it=tmplatmats_34.begin(); it!=tmplatmats_34.end(); ++it){
                    std::string idstring=it->first;
                    if(idstring.find("loc")!=std::string::npos) continue;

                    size_t firstpos=idstring.find("_");
                    size_t secondpos=idstring.find(firstpos);
                    std::string conttype=idstring.substr(0,firstpos);
                    std::transform(conttype.begin(),conttype.end(),conttype.begin(),::tolower);

                    std::string srcspin=idstring.substr(firstpos+1,secondpos);
                    std::string srcspinval(&srcspin[srcspin.size()-1]);
                    srcspin=srcspin.substr(0,srcspin.size()-1);

                    for(std::map<std::string,HalfSpinMatrix>::iterator innerit=projectors.begin(); innerit!=projectors.end(); ++innerit){
                        std::string inneridstring=innerit->first;
                        std::string snkspin=inneridstring.substr(0,inneridstring.size()-1);
                        if(snkspin!=srcspin) continue;
                        std::string snkspinval(&inneridstring[inneridstring.size()-1]);

                        std::string corrname=conttype+"corr"+"_"+srcspin+"_"+snkspinval+"_"+srcspinval+"_"+displacedir; //s[0]; // 0 = mu hardcoded for this version.
                        token=zero;
                        // TODO src always is 0 for this version!!
                        // if(src==0) token=zero;
                        // else{
                        //     swatch_io_read.start();
                        //     chk.get_parameter(boostdir+"/"+corrname,token);
                        //     swatch_io_read.stop();
                        // }
                        token+=/* sourcefact * tmpprefact * */ tshiftmap(LatticeComplex(trace((innerit->second)*(it->second))));
                        swatch_io_write.start();
                        chk.set_parameter(boostdir+"/"+corrname+"_34",token);
                        swatch_io_write.stop();
                    }
                }

                //state
                chk.set_parameter("mucurrent",static_cast<unsigned int>(0+1)); // 0 = mu hardcoded for this version.
                // if((mu+1)<sourcepars.displacements.nrows()) chk.set_consistency(true);
                chk.set_consistency(true);  // can hardcode THIS too, because mu is always 0!
                chk.close();
                swatch_pdisp.stop();

                /*#ifndef NO_SIGNAL_HANDLERS
                  sighand::exit_code_on_abort();
                  #endif*/
                //This signal handling stuff is not necessary. I don't want the baggage of Thorsten's specialized XML readexs.

                QDPIO::cout << "NN-corr-disp: time=" << swatch_pdisp.getTimeInSeconds() << std::endl;
                swatch_pdisp.reset();
                // if(src!=0) QDPIO::cout << "NN-corr-io-read: time=" << swatch_io_read.getTimeInSeconds() << std::endl;
                QDPIO::cout << "NN-corr-io-write: time=" << swatch_io_write.getTimeInSeconds() << std::endl;
                swatch_io_read.reset();
                swatch_io_write.reset();

            } //This is manually here to close loop while porting.
            /********************************************************
             *                                                       *
             *   END LOOP OVER BOOSTS                                *
             *                                                       *
             ********************************************************/

            //move the checkpoint file:
            std::string timestamp=getTimestamp();
            rename(std::string(params.nblockparam.output_filename+".NN_w.chk").c_str(),std::string(params.nblockparam.output_filename).c_str());

#endif
            //clear baryon blocks:
            //clearTopologies();

            //This is where latscat's swatch_everything stops. I am going to use the timer that's printed at the end.

#else
            QDPIO::cout << "This measurement only works if we have linked against FFTW. Please rebuild." << std::endl;
#endif

#else
            QDPIO::cout << "This measurement only works if we have enabled HDF5. Please rebuild." << std::endl;
#endif

            snoop.stop();
            QDPIO::cout << LalibeNucleonBlockEnv::name << ": total time = " << snoop.getTimeInSeconds() << " secs" << std::endl;
            QDPIO::cout << LalibeNucleonBlockEnv::name << ": ran successfully" << std::endl;
            END_CODE();

        }// Function Call
    }// LalibeNucleonBlockEnv
};
