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

        void read(XMLReader& xml, const std::string& path, NucleonBlockParams::NucleonBlock_t& par){
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
            // Compute blocks needed for local correlators?
            if (paramtop.count("compute_locals") > 0)
                read(paramtop, "compute_locals", par.compute_locals);
            else
                par.compute_locals = true;
            // Are we doing negative parity blocks?
            if (paramtop.count("negative_parity") != 0)
                read(paramtop, "negative_parity", par.negative_parity);
            else
                par.negative_parity = false;
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

        // Param stuff
        NucleonBlockParams::NucleonBlockParams(){}
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

#ifdef DEBUG
            QDPIO::cout << "Warning, DEBUG mode enabled!" << std::endl;
#endif

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

            // Figure out displacement of props so we know if we need to make
            // multiple sets of blocks or just blocks for local corrs
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
            if (disp[3] != 0)
            {
                QDPIO::cerr << "\n###############################################################\n" << std::endl;
                QDPIO::cerr << LalibeNucleonBlockEnv::name << ": props must be at same time-slice" << std::endl;
                QDPIO::cerr << "                prop0 t0 = " << pos0[3] << std::endl;
                QDPIO::cerr << "                prop1 t0 = " << pos1[3] << std::endl;
                QDPIO::cerr << "\n###############################################################\n" << std::endl;
                QDP_abort(1);
            }

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
            Timeshiftmap tshiftmap(pos0[j_decay],j_decay,Layout::lattSize()[j_decay]);

            //If the propagators aren't specified to already be in the Dirac basis, we rotate them now.
            if(params.nblockparam.in_dirac_basis == false)
            {
                QDPIO::cout << "Rotating the propagators to Dirac basis." << std::endl;
                rotate_to_Dirac_Basis(prop0);
                rotate_to_Dirac_Basis(prop1);
            }

            // Spin projector
            multi1d<BaryOp> Nup=get_local_MA_single(0,"DP");

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
            // FFT timer
            StopWatch swatch_fft;
            swatch_fft.reset();
            double block_time = 0.;

            // Put it in the NamedObjectMap (if it doesn't exist yet)
            if (!TheNamedObjMap::Instance().check(params.named_obj.block_map)){
                TheNamedObjMap::Instance().create<BlockMapType>(params.named_obj.block_map);
                QDPIO::cout << "Creating new nucleon block map " << params.named_obj.block_map << std::endl;
            }
            else QDPIO::cout << "block map exists, adding to it " << params.named_obj.block_map << std::endl;
            BlockMapType& blockMap = TheNamedObjMap::Instance().getData<BlockMapType>(params.named_obj.block_map);

            // Create a local instance to store blocks while we make them
            LatticeHalfBaryonblock tmpBlock;

            // string ids
            std::string propId_0 = params.named_obj.prop0_id;
            std::string propId_1 = params.named_obj.prop1_id;

            // Block 000[+1]
            BlockMapKeyType theKey = {propId_0, propId_0, propId_0, 1, neg_par, pos0, displacedir};
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
                theKey = {propId_0, propId_0, propId_0, -1, neg_par, pos0, displacedir};
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
                theKey = {propId_1, propId_1, propId_1, -1, neg_par, pos0, displacedir};
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
                    theKey = {propId_1, propId_1, propId_1, +1, neg_par, pos0, displacedir};
                    block_str = "111[+1]";
                    if (blockMap.count(theKey) == 0){
                        swatch_fft.start();
                        blockMap[theKey] = fft(tmpBlock, +1);
                        swatch_fft.stop();
                        QDPIO::cout << "created block " << block_str << " " << neg_par_str << std::endl;
                    }
                    else QDPIO::cout << "exists  block " << block_str << " " << neg_par_str << std::endl;
                    // 000[-1]
                    theKey = {propId_0, propId_0, propId_0, -1, neg_par, pos0, displacedir};
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
                theKey = {propId_0, propId_0, propId_1, 1, neg_par, pos0, displacedir};
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
                theKey = {propId_0, propId_1, propId_0, 1, neg_par, pos0, displacedir};
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
                theKey = {propId_1, propId_0, propId_0, 1, neg_par, pos0, displacedir};
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
                theKey = {propId_0, propId_1, propId_1, -1, neg_par, pos0, displacedir};
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
                theKey = {propId_1, propId_0, propId_1, -1, neg_par, pos0, displacedir};
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
                theKey = {propId_1, propId_1, propId_0, -1, neg_par, pos0, displacedir};
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

#else
            QDPIO::cout << "This measurement only works if we have linked against FFTW. Please rebuild." << std::endl;
#endif

            snoop.stop();
            QDPIO::cout << LalibeNucleonBlockEnv::name << ": total time = " << snoop.getTimeInSeconds() << " secs" << std::endl;
            QDPIO::cout << LalibeNucleonBlockEnv::name << ": ran successfully" << std::endl;
            END_CODE();

        }// Function Call
    }// LalibeNucleonBlockEnv
};
