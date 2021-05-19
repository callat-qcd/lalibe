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
#include "utils_named_object.h"
#include "../momentum/lalibe_sftmom.h"
#include "nucleon_block.h"
#include "two_nucleons.h"
#include "../matrix_elements/bilinear_gamma.h"
//TODO: bilinear_gamma above shouldn't be needed once boiler plate is removed.
#include "../contractions/baryon_contractions_func_w.h"

// Latscat Stuff
#include "../NN/NN_LC_w.h"
#include "../NN/spinstuff.h"
#include "../NN/transform.h"
#include "../NN/checkpointstuff.h"
#include "../NN/momstuff.h"

// For checking file existence
#include <unistd.h>

namespace Chroma
{
    namespace LalibeTwoNucleonsEnv
    {
        namespace
        {
            AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, const std::string& path)
            {
                return new InlineMeas(TwoNucleonsParams(xml_in, path));
            }
            //! Local registration flag
            bool registered = false;
        }
        const std::string name = "TWO_NUCLEONS";

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

        bool file_exists( const std::string &filename )
        {
            return access( filename.c_str(), 0 ) == 0;
        }

        void read(XMLReader& xml, const std::string& path, TwoNucleonsParams::TwoNucleons_t& par)
        {
            XMLReader paramtop(xml, path);
            read(paramtop, "output_filename",       par.output_filename); //output file

            // check if file already exists
            if( file_exists(par.output_filename) ){
                QDPIO::cout << "\nyour H5 file already exists - exit before doing work" << std::endl;
                QDPIO::cout << par.output_filename << std::endl;
                QDP_abort(1);
            }

            read(paramtop, "contractions_filename", par.contractions_filename); //hdf5 file containing contractions
            // FFT tuning?
            if (paramtop.count("fft_tune") != 0)
                read(paramtop, "fft_tune", par.fft_tune); //tune fft?
            else
                par.fft_tune = false;
            // LUSTRE striping
            if (paramtop.count("output_stripesize") != 0)
                read(paramtop, "output_stripesize", par.output_stripesize); //output stripesize; default recommended
            else
                par.output_stripesize = 0;
            if (paramtop.count("parities") > 0)
                read(paramtop, "parities", par.parities);
            else
            {
                par.parities.resize(2);
                par.parities[0] = "POS_PAR";
                par.parities[1] = "NEG_PAR";
            }
            if (paramtop.count("origin") > 0)
                read(paramtop, "origin", par.origin);
            else{
                QDPIO::cerr << "You must tell us what the origin is" << std::endl;
                QDP_abort(1);
            }
            // Do local contractions?  111[+1] + 111[1-]
            if (paramtop.count("compute_locals") > 0)
                read(paramtop, "compute_locals", par.compute_locals);
            else{
                QDPIO::cout << "compute_locals not set, default = true" << std::endl;
                par.compute_locals = true;
            }
            // Do local origin contractions? 000[+1] 000[-1]
            if (paramtop.count("compute_loc_o") > 0)
                read(paramtop, "compute_loc_o", par.compute_loc_o);
            else{
                QDPIO::cout << "compute_loc_o not set, default = true" << std::endl;
                par.compute_loc_o = true;
            }
            // Compute the proton correlation function?
            if (paramtop.count("compute_proton") > 0)
                read(paramtop, "compute_proton", par.compute_proton);
            else{
                QDPIO::cout << "compute_proton not set, default = true" << std::endl;
                par.compute_proton = true;
            }
            // write spin-flip correlators?
            if (paramtop.count("spin_flip") > 0)
                read(paramtop, "spin_flip", par.spin_flip);
            else{
                QDPIO::cout << "spin_flip not set, default = false" << std::endl;
                par.spin_flip = false;
            }
            // list of total momentum boosts
            read(paramtop, "boosts", par.boosts);
            // For now - only support P_tot = (0,0,0)
            if (par.boosts.size() > 1 || par.boosts[0][0] != 0 || par.boosts[0][1] != 0 || par.boosts[0][2] != 0){
                QDPIO::cerr << "\n###############################################################\n" << std::endl;
                QDPIO::cerr << "  For now, only a single boost of (0,0,0) is allowed" << std::endl;
                QDPIO::cerr << "\n###############################################################\n" << std::endl;
                QDP_abort(1);
            }

            // Do we want to delete blocks?
            if( paramtop.count("prop_list_delete") > 0){
                read(paramtop, "prop_list_delete", par.prop_list_delete);
            }
        }

        void write(XMLWriter& xml, const std::string& path, TwoNucleonsParams::TwoNucleons_t& par)
        {
            push(xml, path);
            write(xml, "contractions_filename",  par.contractions_filename);
            write(xml, "output_filename",        par.output_filename);
            write(xml, "fft_tune",               par.fft_tune);
            write(xml, "output_stripesize",      par.output_stripesize);
            write(xml, "compute_locals",         par.compute_locals);
            write(xml, "compute_proton",         par.compute_proton);
            write(xml, "spin_flip",              par.spin_flip);
            write(xml, "parities",               par.parities);
            write(xml, "boosts",                 par.boosts);
            pop(xml);
        }

        //! NamedObject input
        void read(XMLReader& xml, const std::string& path, TwoNucleonsParams::NamedObject_t& input)
        {
            XMLReader inputtop(xml, path);
            read(inputtop, "gauge_id"     , input.gauge_id);
            read(inputtop, "prop0_list"   , input.prop0_list);
            read(inputtop, "prop1_list"   , input.prop1_list);
            // read group xml of blocks
            input.nucleon_blocks = readXMLArrayGroup(inputtop, "nucleon_blocks", "block");
        }

        //! NamedObject output
        void write(XMLWriter& xml, const std::string& path, const TwoNucleonsParams::NamedObject_t& input)
        {
            push(xml, path);
            write(xml, "gauge_id"     , input.gauge_id);
            write(xml, "prop0_list"   , input.prop0_list);
            write(xml, "prop1_list"   , input.prop1_list);
            // Group writer
            push(xml,  "nucleon_blocks");
            for( int t(0);t<input.nucleon_blocks.size();t++)
            {
                push(xml,"elem");
                xml << input.nucleon_blocks[t].xml;
                pop(xml);
            }
            pop(xml);
            // End nucleon block group writer
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
        TwoNucleonsParams::TwoNucleonsParams(){}
        TwoNucleonsParams::TwoNucleonsParams(XMLReader& xml_in, const std::string& path)
        {
            try
            {
                XMLReader paramtop(xml_in, path);
                if (paramtop.count("Frequency") == 1)
                    read(paramtop, "Frequency", frequency);
                else
                    frequency = 1;

                // Parameters for source construction
                read(paramtop, "TwoNucleonsParams", twonucleonsparam);

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

        void TwoNucleonsParams::writeXML(XMLWriter& xml_out, const std::string& path)
        {
            push(xml_out, path);
            write(xml_out, "TwoNucleonsParams", twonucleonsparam);
            write(xml_out, "NamedObject"      , named_obj);
            pop(xml_out);
        }

        // Function call
        void InlineMeas::operator()(unsigned long update_no, XMLWriter& xml_out)
        {
            START_CODE();
            StopWatch snoop;
            snoop.reset();
            snoop.start();
            QDPIO::cout << "TWO_NUCLEONS: start" << std::endl;

#ifndef BUILD_HDF5
            QDPIO::cerr << "\n###########################################\n" << std::endl;
            QDPIO::cerr << LalibeTwoNucleonsEnv::name
                        << " only works if we have enabled HDF5. Please rebuild." << std::endl;
            QDPIO::cerr << "\n###########################################\n" << std::endl;
            QDP_abort(1);
#else
#if QDP_NC != 3
            QDPIO::cerr << "\n###########################################\n" << std::endl;
            QDPIO::cerr << "Contractions not yet implemented for NC!=3." << std::endl;
            QDPIO::cerr << "\n###########################################\n" << std::endl;
            QDP_abort(1);
#endif
#ifndef BUILD_FFTW
            QDPIO::cerr << "\n###########################################\n" << std::endl;
            QDPIO::cerr << "You must build with FFT or CUFFT" << std::endl;
            QDPIO::cerr << "\n###########################################\n" << std::endl;
            QDP_abort(1);
#endif

            // Grab a reference to the gauge field
            const multi1d<LatticeColorMatrix>& u =
                TheNamedObjMap::Instance().getData <multi1d <LatticeColorMatrix> >(params.named_obj.gauge_id);

            // The structure of this function is
            // 1. Get list of prop0_Ids and prop1_Ids
            // 2. Get list of nucleon blocks keys
            // 3. Based upon requested set of contractions, check if we have all
            //    required blocks
            // 4. If yes, proceed to do contractions

            // Get propIds and locations

            int                  n_blocks = params.named_obj.prop0_list.size();
            multi1d<std::string> block_names(n_blocks);
            multi1d<std::string> prop0_Ids(n_blocks);
            multi1d<std::string> prop1_Ids(params.named_obj.prop1_list.size());
            multi2d< int >       pos0_list(n_blocks, Nd);
            multi2d< int >       pos1_list(n_blocks, Nd);
            multi1d<int>         pos0(Nd), pos1(Nd), disp(Nd);
            multi1d<std::string> disp_list(n_blocks);
            std::string          displacedir, pos0_str, pos1_str;
            std::string          nodisplace="px0py0pz0";
            multi1d<int>         origin = params.twonucleonsparam.origin;

            if (prop0_Ids.size() != prop1_Ids.size()){
                QDPIO::cerr << "You must pass the same number of prop0 and prop1 files" << std::endl;
                QDP_abort(1);
            }
            else
                QDPIO::cout << "We have " << prop0_Ids.size() << " sets of propagators" << std::endl;
            for (int b=0; b<n_blocks; b++){
                prop0_Ids[b] = params.named_obj.prop0_list[b];
                prop1_Ids[b] = params.named_obj.prop1_list[b];
                pos0         = LalibeUtilsNambedObjEnv::get_prop_position(prop0_Ids[b]);
                pos0_list[b] = pos0;
                pos1         = LalibeUtilsNambedObjEnv::get_prop_position(prop1_Ids[b]);
                //             pos0 = origin so no need to compute disp
                disp         = pos1 - origin;
                displacedir  = "";
                pos0_str     = "";
                pos1_str     = "";
                for(unsigned int d=0; d<Nd; d++){
                    if (pos0[d] - origin[d] != 0){
                        QDPIO::cerr << prop0_Ids[b] << " not located at origin! " << std::endl;
                        QDP_abort(1);
                    }
                    // make disp string
                    // adjust for boundary conditions
                    disp[d] = (disp[d] + Layout::lattSize()[d]) % Layout::lattSize()[d];
                    if(disp[d] > Layout::lattSize()[d]/2){ disp[d] = disp[d] - Layout::lattSize()[d]; }
                    if (d < Nd-1){// only use spatial coords for displacement
                        if(disp[d] >= 0) displacedir+="p";
                        else displacedir+="m";
                        displacedir+=dirlist[d]+std::to_string(abs(disp[d]));
                    }
                    pos0_str += dirlist[d]+std::to_string(pos0[d]);
                    pos1_str += dirlist[d]+std::to_string(pos1[d]);
                }
                disp_list[b] = displacedir;
            }
            // we need an inner map for the props to simplify looping over block choices
            std::map<char, multi1d<std::string>> prop01 { {'0', prop0_Ids}, {'1', prop1_Ids}};

            // Get block keys
            // Start by making list of BlockMapType
            if ( n_blocks != params.named_obj.nucleon_blocks.size()){
                QDPIO::cerr << "You must pass the same number of blocks as props" << std::endl;
                QDP_abort(1);
            }
            multi1d<LalibeNucleonBlockEnv::BlockMapType*> blockMap_list(n_blocks);
            multi1d<Complex> weights(n_blocks);
            for (int b=0; b<n_blocks; b++){
                std::istringstream xml_block(params.named_obj.nucleon_blocks[b].xml);
                XMLReader block(xml_block);
                // blocks
                read(block, "block", block_names[b]);
                blockMap_list[b] = &TheNamedObjMap::Instance().getData<LalibeNucleonBlockEnv::BlockMapType>(block_names[b]);
                // weights
                read(block, "weight", weights[b]);
            }
            QDPIO::cout << "We have " << blockMap_list.size() << " sets of blocks" << std::endl;

            // Now, loop over keys in block_maps to see if all blocks are present
            bool have_all_blocks = true;
            std::string parity;
            LalibeNucleonBlockEnv::BlockMapKeyType key0, key1;

            for (int p=0; p<params.twonucleonsparam.parities.size(); p++){
                parity = params.twonucleonsparam.parities[p];
                QDPIO::cout << "  checking " << parity << std::endl;
                for (int b=0; b<n_blocks; b++){
                    if ( params.twonucleonsparam.compute_locals ){
                        // We need 000[+1] and 000[-1] for all blocks - for these, we have nodisplace = "px0py0pz0"
                        key0 = {prop0_Ids[b], prop0_Ids[b], prop0_Ids[b],  1, parity, origin, nodisplace};
                        key1 = {prop0_Ids[b], prop0_Ids[b], prop0_Ids[b], -1, parity, origin, nodisplace};
                        if (!blockMap_list[b]->count(key0) || !blockMap_list[b]->count(key1)){
                            have_all_blocks = false;
                            QDPIO::cout << "missing 000 " << prop0_Ids[b] << std::endl;
                        }
                        // if prop1 != prop0, we also need 111[+1] and 111[-1]
                        if (prop1_Ids[0] != prop0_Ids[0]) {
                            key0 = {prop1_Ids[b], prop1_Ids[b], prop1_Ids[b],  1, parity, origin, disp_list[b]};
                            key1 = {prop1_Ids[b], prop1_Ids[b], prop1_Ids[b], -1, parity, origin, disp_list[b]};
                            if (!blockMap_list[b]->count(key0) || !blockMap_list[b]->count(key1)){
                                have_all_blocks = false;
                                QDPIO::cout << "missing 111 " << prop1_Ids[b] << std::endl;
                            }
                        }
                    }
                    // if prop1 != prop0, we need 001, 010, 100 [+1] && 011, 101, 110 [-1]
                    if (prop1_Ids[0] != prop0_Ids[0]){
                        // 001 [+1] and 110[1-]
                        key0 = {prop0_Ids[b], prop0_Ids[b], prop1_Ids[b],  1, parity, origin, disp_list[b]};
                        key1 = {prop1_Ids[b], prop1_Ids[b], prop0_Ids[b], -1, parity, origin, disp_list[b]};
                        if (!blockMap_list[b]->count(key0) || !blockMap_list[b]->count(key1)){
                            have_all_blocks = false;
                            QDPIO::cout << "missing 001 or 110 " << prop0_Ids[b] << " " << prop1_Ids[b] << std::endl;
                        }
                        // 010 [+1] and 101[1-]
                        key0 = {prop0_Ids[b], prop1_Ids[b], prop0_Ids[b],  1, parity, origin, disp_list[b]};
                        key1 = {prop1_Ids[b], prop0_Ids[b], prop1_Ids[b], -1, parity, origin, disp_list[b]};
                        if (!blockMap_list[b]->count(key0) || !blockMap_list[b]->count(key1)){
                            have_all_blocks = false;
                            QDPIO::cout << "missing 010 or 101 " << prop0_Ids[b] << " " << prop1_Ids[b] << std::endl;
                        }
                        // 100 [+1] and 0111[1-]
                        key0 = {prop1_Ids[b], prop0_Ids[b], prop0_Ids[b],  1, parity, origin, disp_list[b]};
                        key1 = {prop0_Ids[b], prop1_Ids[b], prop1_Ids[b], -1, parity, origin, disp_list[b]};
                        if (!blockMap_list[b]->count(key0) || !blockMap_list[b]->count(key1)){
                            have_all_blocks = false;
                            QDPIO::cout << "missing 100 or 011 " << prop0_Ids[b] << " " << prop1_Ids[b] << std::endl;
                        }
                    }
                }
            }
            if (have_all_blocks){
                QDPIO::cout << "  we have all blocks" << std::endl;
            }
            else{
                QDPIO::cerr << "  You did not provide all blocks necessary to perform the requested contractions" << std::endl;
                QDP_abort(1);
            }


            /*
                Proceed with NN contractions
            */
            //relevant spin projectors, in DP rep
            std::map<std::string,HalfSpinMatrix> projectors;
            projectors["SING0"]=getHalfProjector("SING0");
            projectors["TRIPP"]=getHalfProjector("TRIPP1");
            projectors["TRIP0"]=getHalfProjector("TRIP0");
            projectors["TRIPM"]=getHalfProjector("TRIPM1");


            // Minimal checkpointing
            //checkpoint chk(params.twonucleonsparam.output_filename+".NN_w.chk",params.twonucleonsparam.output_stripesize);
            /* turn off checkpointing now as we do not use it
               also, we will set the write mode to HDF5Base::ate which will not allow writing to the same file
             */
            HDF5Base::writemode h5mode, mucurr_mode;
            h5mode = HDF5Base::ate;
            mucurr_mode = HDF5Base::trunc;
            checkpoint chk(params.twonucleonsparam.output_filename,params.twonucleonsparam.output_stripesize);
            

            int j_decay = Nd - 1; // Assume t_dir = j_decay
            // We need FFT for the non-local contractions
            Fourier fft(j_decay);
            Fourier fftblock(j_decay);
            // TUNE the FFT?
            if(params.twonucleonsparam.fft_tune){
                QDPIO::cout << "Tuning FFT for better performance..." << std::flush;
                StopWatch swatch_fftune;
                swatch_fftune.reset();
                swatch_fftune.start();
                fftblock.tune(sizeof(HalfBaryonblock),true);
                swatch_fftune.stop();
                QDPIO::cout << "done! Time " << swatch_fftune.getTimeInSeconds() << std::endl;
            }
            /*
                Truncate is used to truncate the mom space correlators
                we do not use this anymore and retain the full mom space
                -1 --> do not truncate
            */
            QDPIO::cout << "Creating timers..." << std::endl;
            StopWatch swatch_topologies, swatch_io_write, swatch_contract_local, swatch_contract_nonlocal, swatch_fft, swatch_blocks;

            swatch_topologies.reset();
            swatch_topologies.start();
            int truncate = -1;
            initTopologies(params.twonucleonsparam.contractions_filename, truncate, j_decay);
            swatch_topologies.stop();
            QDPIO::cout << LalibeTwoNucleonsEnv::name << ": initTopologies time = "
                        << swatch_topologies.getTimeInSeconds() << " secs" << std::endl;

            swatch_io_write.reset();
            swatch_contract_local.reset();
            swatch_contract_nonlocal.reset();
            swatch_fft.reset();
            swatch_blocks.reset();

            // all props have the same source, so just use the first prop0 position
            Timeshiftmap tshiftmap(pos0_list[0][j_decay],j_decay,Layout::lattSize()[j_decay]);

            // Create objects to store results
            LatticeComplex proton0, proton1, proton0_34, proton1_34, token;
            proton0    = zero;
            proton0_34 = zero;
            proton1    = zero;
            proton1_34 = zero;
            std::map<std::string,LatticeHalfSpinMatrix> NN_map, NN_map_34;
            // zero out the memory
            for (int par=0; par < params.twonucleonsparam.parities.size(); par++){
                parity = params.twonucleonsparam.parities[par];
                
                auto& NN_map_iterator = (parity == "POS_PAR") ? NN_map : NN_map_34;
                for(unsigned int s = 0; s < nonlocalterms.size(); s++){
                    NN_map_iterator[nonlocalterms[s]]=LatticeHalfSpinMatrix();
                }
                if( params.twonucleonsparam.compute_locals ){
                    for (unsigned int s = 0; s < localterms1.size(); s++){
                        NN_map_iterator[localterms1[s]]=LatticeHalfSpinMatrix();
                    }
                    if( params.twonucleonsparam.compute_loc_o ){
                        for (unsigned int s = 0; s < localterms0.size(); s++){
                            NN_map_iterator[localterms0[s]]=LatticeHalfSpinMatrix();
                        }
                    }
                }
                // We must zero out the data before adding to it
                for(std::map<std::string,LatticeHalfSpinMatrix>::iterator it=NN_map_iterator.begin(); it!=NN_map_iterator.end(); ++it){
                    it->second=zero;
                }
            }


            /**********************************************************************
             *                                                                    *
             *   LOOP OVER BOOSTS                                                 *
             *        We don't support boosts now, but we might, so keep the loop *
             **********************************************************************/
            for(unsigned int boost=0; boost < params.twonucleonsparam.boosts.size(); boost++){
                multi1d<int> p_tot = params.twonucleonsparam.boosts[boost];
                std::string boostdir=boost_string(p_tot);

                momentum mom(params.twonucleonsparam.boosts[boost]);
                QDPIO::cout << mom << "    " << boostdir << std::endl;

                chk.create_directory(boostdir);
                chk.close();

                LatticeComplex phases=get_phases(mom,j_decay,+1);

                for (int par=0; par < params.twonucleonsparam.parities.size(); par++){
                    parity = params.twonucleonsparam.parities[par];
                    QDPIO::cout << "Starting " << parity << " contractions" << std::endl;
                    // positive parity contraction
                    if( parity == "POS_PAR"){
                        contract(proton0, proton1, NN_map, 
                                 blockMap_list, prop0_Ids, prop1_Ids, 
                                 weights, origin, disp_list, parity,
                                 phases, fft, 
                                 params.twonucleonsparam.compute_locals, params.twonucleonsparam.compute_loc_o);
                    } else if (parity == "NEG_PAR"){
                        contract(proton0_34, proton1_34, NN_map_34, 
                                 blockMap_list, prop0_Ids, prop1_Ids, 
                                 weights, origin, disp_list, parity,
                                 phases, fft, 
                                 params.twonucleonsparam.compute_locals, params.twonucleonsparam.compute_loc_o);
                    }

                    // now write data
                    chk.open();
                    chk.set_consistency(false);
                   
                    //Single proton
                    if (p_tot[0]==0 && p_tot[1] == 0 && p_tot[2] == 0 && params.twonucleonsparam.compute_proton && params.twonucleonsparam.compute_locals){
                        const auto& prot0_obj = (parity == "POS_PAR") ? proton0 : proton0_34;
                        const auto& prot1_obj = (parity == "POS_PAR") ? proton1 : proton1_34;
                        QDPIO::cout << "Single proton1." << std::endl;
                        token=zero;
                        token += tshiftmap(prot1_obj, true);
                        swatch_io_write.start();
                        std::string corrname = "proton1_"+pos1_str+"_"+parity;
                        chk.set_parameter(corrname,token, h5mode);
                        swatch_io_write.stop();

                        if ( params.twonucleonsparam.compute_loc_o){
                            QDPIO::cout << "Single proton0." << std::endl;
                            token=zero;
                            token += tshiftmap(prot0_obj, true);
                            swatch_io_write.start();
                            std::string corrname = "proton0_"+pos0_str+"_"+parity;
                            chk.set_parameter(corrname,token, h5mode);
                            swatch_io_write.stop();
                        }
                    }

                    /*
                        Two nucleons
                    */
                    // make directory for origin of contractions
                    std::string p0_dir = "O0_"+pos0_str;
                    std::string p1_dir = "O1_"+pos1_str;
                    chk.create_directory(boostdir+"/"+p0_dir);
                    if( params.twonucleonsparam.compute_locals ){
                        chk.create_directory(boostdir+"/"+p1_dir);
                    }

                    const auto& NN_map_iterator = (parity == "POS_PAR") ? NN_map : NN_map_34;
                    //for(std::map<std::string,LatticeHalfSpinMatrix>::iterator it=NN_map.begin(); it != NN_map.end(); ++it){
                    for(auto it=NN_map_iterator.cbegin(); it != NN_map_iterator.cend(); ++it){

                        std::string idstring=it->first;
                        size_t firstpos=idstring.find("_");
                        size_t secondpos=idstring.find(firstpos);
                        std::string conttype=idstring.substr(0,firstpos);
                        std::transform(conttype.begin(),conttype.end(),conttype.begin(),::tolower);

                        std::string srcspin, srcspinval;
                        if( idstring.find("SING") != std::string::npos){
                            srcspin    = "SING";
                            srcspinval = "0";
                        } else if (idstring.find("TRIP") != std::string::npos ){
                            srcspin = "TRIP";
                            srcspinval = idstring.substr(idstring.find("TRIP")+4,1);
                        }

                        for(std::map<std::string,HalfSpinMatrix>::iterator innerit=projectors.begin(); innerit!=projectors.end(); ++innerit){

                            std::string inneridstring=innerit->first;
                            std::string snkspin=inneridstring.substr(0,inneridstring.size()-1);
                            if(snkspin!=srcspin) continue;
                            std::string snkspinval(&inneridstring[inneridstring.size()-1]);

                            if( snkspinval == srcspinval || params.twonucleonsparam.spin_flip){
                                std::string corrname;
                                if( idstring.find("loc") != std::string::npos){
                                    corrname = conttype+"_"+srcspin+"_"+snkspinval+"_"+srcspinval+"_"+nodisplace;
                                } else {
                                    corrname = conttype+"_"+srcspin+"_"+snkspinval+"_"+srcspinval+"_"+disp_list[0];
                                }
                                if (parity == "NEG_PAR")
                                    corrname += "_34";

                                token  = zero;
                                token += tshiftmap(LatticeComplex(trace((innerit->second)*(it->second))));

                                swatch_io_write.start();
                                if (idstring.find("loc1") != std::string::npos){
                                    chk.set_parameter(boostdir+"/"+p1_dir+"/"+corrname,token, h5mode);
                                } else {
                                    chk.set_parameter(boostdir+"/"+p0_dir+"/"+corrname,token, h5mode);
                                }
                                swatch_io_write.stop();
                            }
                        }
                    }
                }
            }
            // 0 = mu hardcoded for this version.
            // if((mu+1)<sourcepars.displacements.nrows()) chk.set_consistency(true);
            // we can overwrite the mucurrent - so change this h5mode
            chk.set_parameter("mucurrent",static_cast<unsigned int>(0+1), mucurr_mode);
            chk.set_consistency(true);  // can hardcode THIS too, because mu is always 0!
            chk.close();

            //move the checkpoint file:
            // We are turning off checkpointing for now
            /*
            std::string timestamp=getTimestamp();
            rename(std::string(params.twonucleonsparam.output_filename+".NN_w.chk").c_str(),
                    std::string(params.twonucleonsparam.output_filename).c_str());
            */

            /* these timers are in the contract() function
            QDPIO::cout << LalibeTwoNucleonsEnv::name << " Block add time="
                        << swatch_blocks.getTimeInSeconds() << std::endl;
            QDPIO::cout << LalibeTwoNucleonsEnv::name << " NN local: time="
                        << swatch_contract_local.getTimeInSeconds() << std::endl;
            QDPIO::cout << LalibeTwoNucleonsEnv::name << " NN non-local: time="
                        << swatch_contract_nonlocal.getTimeInSeconds() << std::endl;
            */
            QDPIO::cout << LalibeTwoNucleonsEnv::name << " NN I/O: time="
                        << swatch_io_write.getTimeInSeconds() << std::endl;
            //clear baryon blocks:
            clearTopologies();

#endif //BUILD_HDF5

            // Now - delete any blocks the user requests
            for( int bList=0; bList < blockMap_list.size(); bList++){
                for (int prop_i = 0; prop_i < params.twonucleonsparam.prop_list_delete.size(); prop_i++){
                    for (auto it=blockMap_list[bList]->begin(); it != blockMap_list[bList]->end();) {

                        if (std::get<0>(it->first) == params.twonucleonsparam.prop_list_delete[prop_i] 
                            || std::get<1>(it->first) == params.twonucleonsparam.prop_list_delete[prop_i]
                            || std::get<2>(it->first) == params.twonucleonsparam.prop_list_delete[prop_i]) {
                            QDPIO::cout << "Deleting block with key: " 
                                        << std::get<0>(it->first) << " "
                                        << std::get<1>(it->first) << " "
                                        << std::get<2>(it->first) << " "
                                        << std::get<3>(it->first) << " "
                                        << std::get<4>(it->first) << " "
                                        << "x"<<std::get<5>(it->first)[0] << " "
                                        << "y"<<std::get<5>(it->first)[1] << " "
                                        << "z"<<std::get<5>(it->first)[2] << " "
                                        << "t"<<std::get<5>(it->first)[3] << " "
                                        << std::get<6>(it->first) << std::endl;
                            // Delete the block
                            blockMap_list[bList]->erase(it++);
                            
                        }
                        else {
                            ++it;
                        }
                    }
                }
            }


            snoop.stop();
            QDPIO::cout << LalibeTwoNucleonsEnv::name << ": total time = " << snoop.getTimeInSeconds() << " secs" << std::endl;
            QDPIO::cout << LalibeTwoNucleonsEnv::name << ": ran successfully" << std::endl;
            END_CODE();

        }// Function Call
    }// LalibeTwoNucleonsEnv
};
