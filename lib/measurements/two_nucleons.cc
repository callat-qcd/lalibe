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

        void read(XMLReader& xml, const std::string& path, TwoNucleonsParams::TwoNucleons_t& par)
        {
            XMLReader paramtop(xml, path);
            read(paramtop, "output_filename",       par.output_filename); //output file
            read(paramtop, "contractions_filename", par.contractions_filename); //hdf5 file containing contractions
            if (paramtop.count("output_stripesize") != 0)
                read(paramtop, "output_stripesize", par.output_stripesize); //output stripesize; default recommended
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
            // list of total momentum boosts
            read(paramtop, "boosts", par.boosts);
            // For now - only support P_tot = (0,0,0)
            if (par.boosts.size() > 1 || par.boosts[0][0] != 0 || par.boosts[0][1] != 0 || par.boosts[0][2] != 0){
                QDPIO::cerr << "\n###############################################################\n" << std::endl;
                QDPIO::cerr << "  For now, only a single boost of (0,0,0) is allowed" << std::endl;
                QDPIO::cerr << "\n###############################################################\n" << std::endl;
                QDP_abort(1);
            }
        }

        void write(XMLWriter& xml, const std::string& path, TwoNucleonsParams::TwoNucleons_t& par)
        {
            push(xml, path);
            write(xml, "contractions_filename",  par.contractions_filename);
            write(xml, "output_filename",        par.output_filename);
            write(xml, "output_stripesize",      par.output_stripesize);
            write(xml, "compute_locals",         par.compute_locals);
            write(xml, "compute_proton",         par.compute_proton);
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

        void addWeightedHalfBlock(LatticeHalfBaryonblock& res,
                                  const Complex& weight,
                                  const LatticeHalfBaryonblock& in) {
            int nSites = Layout::sitesOnNode();
            for (int iX=0; iX<nSites; ++iX) {
                for (int iS1=0; iS1<(Ns>>1); ++iS1) {
                    for (int jS1=0; jS1<(Ns>>1); ++jS1) {
                        for (int iC1=0; iC1<Nc; ++iC1) {
                            for (int iS2=0; iS2<(Ns>>1); ++iS2) {
                                for (int jS2=0; jS2<(Ns>>1); ++jS2) {
                                    for (int iC2=0; iC2<Nc; ++iC2) {
                                        for (int jC2=0; jC2<Nc; ++jC2) {
                                            res.elem(iX).elem(iS1, jS1).elem(iC1).elem(iS2, jS2).elem(iC2, jC2)
                                                += weight.elem().elem().elem()
                                                * in.elem(iX).elem(iS1, jS1).elem(iC1).elem(iS2, jS2).elem(iC2, jC2);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
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
            std::string          displacedir;
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
                pos1         = LalibeUtilsNambedObjEnv::get_prop_position(prop1_Ids[b]);
                //             pos0 = origin so no need to compute disp
                disp         = pos1 - origin;
                displacedir  = "";
                for(unsigned int d=0; d<Nd; d++){
                    if (pos0[d] - origin[d] != 0){
                        QDPIO::cerr << prop0_Ids[b] << " not located at origin! " << std::endl;
                        QDP_abort(1);
                    }
                    if (d < Nd-1){// only use spatial coords for displacement
                        // adjust for boundary conditions
                        disp[d] = (disp[d] + Layout::lattSize()[d]) % Layout::lattSize()[d];
                        if(disp[d] > Layout::lattSize()[d]/2){ disp[d] = disp[d] - Layout::lattSize()[d]; }
                        // make disp string
                        if(disp[d] >= 0) displacedir+="p";
                        else displacedir+="m";
                        displacedir+=dirlist[d]+std::to_string(abs(disp[d]));
                    }
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
            multi1d<const LalibeNucleonBlockEnv::BlockMapType*> blockMap_list(n_blocks);
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
            checkpoint chk(params.twonucleonsparam.output_filename+".NN_w.chk",params.twonucleonsparam.output_stripesize);

            int j_decay = Nd - 1; // Assume t_dir = j_decay
            // We need FFT for the non-local contractions
            Fourier fft(j_decay);
            Fourier fftblock(j_decay);
            /*
                Truncate is used to truncate the mom space correlators
                we do not use this anymore and retain the full mom space
                -1 --> do not truncate
            */
            int truncate = -1;
            initTopologies(params.twonucleonsparam.contractions_filename, truncate, j_decay);

            QDPIO::cout << "Creating timers..." << std::endl;
            StopWatch swatch_io_write, swatch_contract_local, swatch_contract_nonlocal, swatch_fft, swatch_blocks;
            swatch_io_write.reset();
            swatch_contract_local.reset();
            swatch_contract_nonlocal.reset();
            swatch_fft.reset();
            swatch_blocks.reset();

            // all props have the same source, so just use the first prop0 position
            Timeshiftmap tshiftmap(pos0_list[0][j_decay],j_decay,Layout::lattSize()[j_decay]);

            // Create two baryon block objects to store blocks for contractions
            LatticeHalfBaryonblock block0, block1;

            // Create objects to store results
            LatticeComplex latt_prot0, latt_prot1, token;
            std::map<std::string,LatticeHalfSpinMatrix> latt_nn_map;
            for(unsigned int s = 0; s < 8; s++){
                latt_nn_map[contterms[s]]=LatticeHalfSpinMatrix();
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
                    contract(latt_prot0, latt_nn_map, 
                             blockMap_list, prop0_Ids, prop1_Ids, 
                             origin, weights, disp_list, parity,
                             phases, fft, 
                             params.twonucleonsparam.compute_locals, params.twonucleonsparam.compute_loc_o);


                    // now write data
                    chk.open();
                    chk.set_consistency(false);
                    /*
                        Single proton
                    if (p_tot[0]==0 && p_tot[1] == 0 && p_tot[2] == 0 && params.twonucleonsparam.compute_proton && params.twonucleonsparam.compute_locals){
                        QDPIO::cout << "Single proton." << std::endl;
                        token=zero;
                        token += tshiftmap(latt_prot0, true);
                        swatch_io_write.start();
                        std::string corrname = "proton_"+parity;
                        chk.set_parameter(corrname,token);
                        swatch_io_write.stop();
                    }
                    */
                    /*
                        Two nucleons
                    */
                    for(std::map<std::string,LatticeHalfSpinMatrix>::iterator it=latt_nn_map.begin();
                        it != latt_nn_map.end(); ++it){

                        std::string idstring=it->first;
                        size_t firstpos=idstring.find("_");
                        size_t secondpos=idstring.find(firstpos);
                        std::string conttype=idstring.substr(0,firstpos);
                        std::transform(conttype.begin(),conttype.end(),conttype.begin(),::tolower);

                        std::string srcspin=idstring.substr(firstpos+1,secondpos);
                        std::string srcspinval(&srcspin[srcspin.size()-1]);
                        srcspin=srcspin.substr(0,srcspin.size()-1);

                        for(std::map<std::string,HalfSpinMatrix>::iterator innerit=projectors.begin();
                            innerit!=projectors.end(); ++innerit){

                            std::string inneridstring=innerit->first;
                            std::string snkspin=inneridstring.substr(0,inneridstring.size()-1);
                            if(snkspin!=srcspin) continue;
                            std::string snkspinval(&inneridstring[inneridstring.size()-1]);

                            // the disp_list[0] needs to get replaced
                            std::string corrname = conttype+"corr"+"_"+srcspin+"_"+snkspinval+"_"+srcspinval+"_"+disp_list[0];
                            token  = zero;
                            token += tshiftmap(LatticeComplex(trace((innerit->second)*(it->second))));
                            // DEBUG PRINT
                            /*
                            for(int nt=0; nt<8; nt++){
                                QDPIO::cout << "nt = " << nt << " " << corrname << " = " << token.elem(nt).elem(0).elem(0).elem(0);
                            }
                            */
                            swatch_io_write.start();
                            if (parity == "POS_PAR")
                                chk.set_parameter(boostdir+"/"+corrname,token);
                            else
                                chk.set_parameter(boostdir+"/"+corrname+"_34",token);
                            swatch_io_write.stop();
                        }
                    }
                }
            }
            // 0 = mu hardcoded for this version.
            // if((mu+1)<sourcepars.displacements.nrows()) chk.set_consistency(true);
            chk.set_parameter("mucurrent",static_cast<unsigned int>(0+1));
            chk.set_consistency(true);  // can hardcode THIS too, because mu is always 0!
            chk.close();

            //move the checkpoint file:
            std::string timestamp=getTimestamp();
            rename(std::string(params.twonucleonsparam.output_filename+".NN_w.chk").c_str(),
                    std::string(params.twonucleonsparam.output_filename).c_str());


            QDPIO::cout << LalibeTwoNucleonsEnv::name << " Block add time="
                        << swatch_blocks.getTimeInSeconds() << std::endl;
            QDPIO::cout << LalibeTwoNucleonsEnv::name << " NN local: time="
                        << swatch_contract_local.getTimeInSeconds() << std::endl;
            QDPIO::cout << LalibeTwoNucleonsEnv::name << " NN non-local: time="
                        << swatch_contract_nonlocal.getTimeInSeconds() << std::endl;
            QDPIO::cout << LalibeTwoNucleonsEnv::name << " NN I/O: time="
                        << swatch_io_write.getTimeInSeconds() << std::endl;
#if 0
            // moved from block code

            /**********************************************************************
             *                                                                    *
             *   LOOP OVER BOOSTS                                                 *
             *        We don't support boosts now, but we might, so keep the loop *
             **********************************************************************/
            for(unsigned int b=0; b < params.twonucleonsparam.boosts.size(); b++){



                swatch_pdisp.start();
                QDPIO::cout << "Computing contractions..." << std::endl;
                // if displacement between sources of uprop_p{1,2} is nonzero:
                QDPIO::cout << "    Doing displaced only." << std::endl;
                QDPIO::cout << "    For now, if you want local, just pass the same propagator twice." << std::endl;

                contract(tmplatcomp_P, tmplatmats, prop_0, prop_1, Nup[0].get_gamma(0), phases, fftblock, false, weights);
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
            rename(std::string(params.twonucleonsparam.output_filename+".NN_w.chk").c_str(),std::string(params.twonucleonsparam.output_filename).c_str());

#endif
            //clear baryon blocks:
            clearTopologies();

            //This is where latscat's swatch_everything stops. I am going to use the timer that's printed at the end.

#endif

            snoop.stop();
            QDPIO::cout << LalibeTwoNucleonsEnv::name << ": total time = " << snoop.getTimeInSeconds() << " secs" << std::endl;
            QDPIO::cout << LalibeTwoNucleonsEnv::name << ": ran successfully" << std::endl;
            END_CODE();

        }// Function Call
    }// LalibeTwoNucleonsEnv
};
