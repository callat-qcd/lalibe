//Arjun Singh Gambhir
//Code ported from latscat, which originally seems to be from Balint.
// $Id: PP_w.cc,v 1.14 2007-02-22 15:58:30 tkurth Exp $
/*! \file
 *  \brief PP scattering
 */

//#include "progs.h"
//Progs has stuff that's needed and also a lot of stuff that's not duwe to how we link to chroma.
//Items that are requires will manually be put here below.
#include "NN_LC_w.h"

#define MAXIMUM(a,b) (a>b ? a : b)
#define MINIMUM(a,b) (a<b ? a : b)

#ifdef MP
#define MPREC 30
#endif

//using namespace Chroma;
//Let's do this properly.
namespace Chroma {

    //global variables:
    sparsearr<Real> single_neutron;
    Topology PP_SING0, PP_TRIP0, PP_TRIPP, PP_TRIPM, PN_SING0, PN_TRIP0, PN_TRIPP, PN_TRIPM;

    //const unsigned int barblocksize=6912;
    //const unsigned int halfbarblocksize=432;
    const unsigned int redblocksize=1728;
    const unsigned int halfredblocksize=216;
    const unsigned int Nshalf=2;
    unsigned int vol3;

    // block function that takes 3 props
    double get_barblock(LatticeHalfBaryonblock& block,
                        const LatticePropagator& prop0,
                        const LatticePropagator& prop1,
                        const LatticePropagator& prop2,
                        const SpinMatrix& diquark_proj){
        StopWatch swatch_block;
        swatch_block.reset();

        swatch_block.start();
        block=getHalfBaryonblock(prop0,prop1,prop2,diquark_proj);
        swatch_block.stop();
        QDPIO::cout << "Build block: time=" << swatch_block.getTimeInSeconds() << std::endl;
        return swatch_block.getTimeInSeconds();
    }

    // block function2 that takes 2 props and mode
    void get_barblock(LatticeHalfBaryonblock& block, const std::string mode, const LatticePropagator& prop0, const LatticePropagator& prop1, const SpinMatrix& diquark_proj){
        StopWatch swatch_block;
        swatch_block.reset();

        //always boost the first propagator:
        swatch_block.start();
        if(mode.compare("000")==0){
            block=getHalfBaryonblock(prop0,prop0,prop0,diquark_proj);
        }
        else if(mode.compare("001")==0){
            block=getHalfBaryonblock(prop0,prop0,prop1,diquark_proj);
        }
        else if(mode.compare("010")==0){
            block=getHalfBaryonblock(prop0,prop1,prop0,diquark_proj);
        }
        else if(mode.compare("011")==0){
            block=getHalfBaryonblock(prop0,prop1,prop1,diquark_proj);
        }
        else if(mode.compare("100")==0){
            block=getHalfBaryonblock(prop1,prop0,prop0,diquark_proj);
        }
        else if(mode.compare("101")==0){
            block=getHalfBaryonblock(prop1,prop0,prop1,diquark_proj);
        }
        else if(mode.compare("110")==0){
            block=getHalfBaryonblock(prop1,prop1,prop0,diquark_proj);
        }
        else if(mode.compare("111")==0){
            block=getHalfBaryonblock(prop1,prop1,prop1,diquark_proj);
        }
        else{
            QDP_error_exit("error, unknown baryon block mode!");
        }
        swatch_block.stop();
        QDPIO::cout << "Build block: time=" << swatch_block.getTimeInSeconds() << std::endl;
    }

    void get_barblock_boost(LatticeHalfBaryonblock& block, const std::string mode, const LatticePropagator& prop0, const LatticePropagator& prop1, const SpinMatrix& diquark_proj, const LatticeComplex& phases){
        StopWatch swatch_block;
        swatch_block.reset();

        //always boost the first propagator:
        swatch_block.start();
        if(mode.compare("000")==0){
            block=getHalfBaryonblock(LatticePropagator(phases*prop0),prop0,prop0,diquark_proj);
        }
        else if(mode.compare("001")==0){
            block=getHalfBaryonblock(LatticePropagator(phases*prop0),prop0,prop1,diquark_proj);
        }
        else if(mode.compare("010")==0){
            block=getHalfBaryonblock(LatticePropagator(phases*prop0),prop1,prop0,diquark_proj);
        }
        else if(mode.compare("011")==0){
            block=getHalfBaryonblock(LatticePropagator(phases*prop0),prop1,prop1,diquark_proj);
        }
        else if(mode.compare("100")==0){
            block=getHalfBaryonblock(LatticePropagator(phases*prop1),prop0,prop0,diquark_proj);
        }
        else if(mode.compare("101")==0){
            block=getHalfBaryonblock(LatticePropagator(phases*prop1),prop0,prop1,diquark_proj);
        }
        else if(mode.compare("110")==0){
            block=getHalfBaryonblock(LatticePropagator(phases*prop1),prop1,prop0,diquark_proj);
        }
        else if(mode.compare("111")==0){
            block=getHalfBaryonblock(LatticePropagator(phases*prop1),prop1,prop1,diquark_proj);
        }
        else{
            QDP_error_exit("error, unknown baryon block mode!");
        }
        swatch_block.stop();
        QDPIO::cout << "Build boosted block: time=" << swatch_block.getTimeInSeconds() << std::endl;
    }

    /* get_barblock but with arbitrary sinks
     *
     *
     */

    void get_barblock(LatticeHalfBaryonblock& block, const std::string mode, const LatticePropagator& prop0, const LatticePropagator& prop1, const SpinMatrix& diquark_proj, const sink& snk){
        StopWatch swatch_block;
        swatch_block.reset();

        //always boost the first propagator:
        swatch_block.start();
        if(mode.compare("000")==0){
            block=getHalfBaryonblock(prop0,prop0,prop0,diquark_proj, snk);
        }
        else if(mode.compare("001")==0){
            block=getHalfBaryonblock(prop0,prop0,prop1,diquark_proj, snk);
        }
        else if(mode.compare("010")==0){
            block=getHalfBaryonblock(prop0,prop1,prop0,diquark_proj, snk);
        }
        else if(mode.compare("011")==0){
            block=getHalfBaryonblock(prop0,prop1,prop1,diquark_proj, snk);
        }
        else if(mode.compare("100")==0){
            block=getHalfBaryonblock(prop1,prop0,prop0,diquark_proj, snk);
        }
        else if(mode.compare("101")==0){
            block=getHalfBaryonblock(prop1,prop0,prop1,diquark_proj, snk);
        }
        else if(mode.compare("110")==0){
            block=getHalfBaryonblock(prop1,prop1,prop0,diquark_proj, snk);
        }
        else if(mode.compare("111")==0){
            block=getHalfBaryonblock(prop1,prop1,prop1,diquark_proj, snk);
        }
        else{
            QDP_error_exit("error, unknown baryon block mode!");
        }
        swatch_block.stop();
        QDPIO::cout << "Build block: time=" << swatch_block.getTimeInSeconds() << std::endl;
    }

    void get_barblock_boost(LatticeHalfBaryonblock& block, const std::string mode, const LatticePropagator& prop0, const LatticePropagator& prop1, const SpinMatrix& diquark_proj, const LatticeComplex& phases, const sink& snk){
        StopWatch swatch_block;
        swatch_block.reset();

        //always boost the first propagator:
        swatch_block.start();
        if(mode.compare("000")==0){
            block=getHalfBaryonblock(LatticePropagator(phases*prop0),prop0,prop0,diquark_proj, snk);
        }
        else if(mode.compare("001")==0){
            block=getHalfBaryonblock(LatticePropagator(phases*prop0),prop0,prop1,diquark_proj, snk);
        }
        else if(mode.compare("010")==0){
            block=getHalfBaryonblock(LatticePropagator(phases*prop0),prop1,prop0,diquark_proj, snk);
        }
        else if(mode.compare("011")==0){
            block=getHalfBaryonblock(LatticePropagator(phases*prop0),prop1,prop1,diquark_proj, snk);
        }
        else if(mode.compare("100")==0){
            block=getHalfBaryonblock(LatticePropagator(phases*prop1),prop0,prop0,diquark_proj, snk);
        }
        else if(mode.compare("101")==0){
            block=getHalfBaryonblock(LatticePropagator(phases*prop1),prop0,prop1,diquark_proj, snk);
        }
        else if(mode.compare("110")==0){
            block=getHalfBaryonblock(LatticePropagator(phases*prop1),prop1,prop0,diquark_proj, snk);
        }
        else if(mode.compare("111")==0){
            block=getHalfBaryonblock(LatticePropagator(phases*prop1),prop1,prop1,diquark_proj, snk);
        }
        else{
            QDP_error_exit("error, unknown baryon block mode!");
        }
        swatch_block.stop();
        QDPIO::cout << "Build boosted block: time=" << swatch_block.getTimeInSeconds() << std::endl;
    }

    /* get_barblock but with linear combination sinks.
     *
     *
     */

    void get_barblock(LatticeHalfBaryonblock& block, const std::string mode, const LatticePropagator& prop0, const LatticePropagator& prop1, const SpinMatrix& diquark_proj, const multi1d<sink*>& snk, const multi1d<Complex>& weights){
        StopWatch swatch_block;
        swatch_block.reset();

        //always boost the first propagator:
        swatch_block.start();
        if(mode.compare("000")==0){
            block=getHalfBaryonblock(prop0,prop0,prop0,diquark_proj, snk, weights);
        }
        else if(mode.compare("001")==0){
            block=getHalfBaryonblock(prop0,prop0,prop1,diquark_proj, snk, weights);
        }
        else if(mode.compare("010")==0){
            block=getHalfBaryonblock(prop0,prop1,prop0,diquark_proj, snk, weights);
        }
        else if(mode.compare("011")==0){
            block=getHalfBaryonblock(prop0,prop1,prop1,diquark_proj, snk, weights);
        }
        else if(mode.compare("100")==0){
            block=getHalfBaryonblock(prop1,prop0,prop0,diquark_proj, snk, weights);
        }
        else if(mode.compare("101")==0){
            block=getHalfBaryonblock(prop1,prop0,prop1,diquark_proj, snk, weights);
        }
        else if(mode.compare("110")==0){
            block=getHalfBaryonblock(prop1,prop1,prop0,diquark_proj, snk, weights);
        }
        else if(mode.compare("111")==0){
            block=getHalfBaryonblock(prop1,prop1,prop1,diquark_proj, snk, weights);
        }
        else{
            QDP_error_exit("error, unknown baryon block mode!");
        }
        swatch_block.stop();
        QDPIO::cout << "Build block: time=" << swatch_block.getTimeInSeconds() << std::endl;
    }

    void get_barblock_boost(LatticeHalfBaryonblock& block, const std::string mode, const LatticePropagator& prop0, const LatticePropagator& prop1, const SpinMatrix& diquark_proj, const LatticeComplex& phases, const multi1d<sink*>& snk, const multi1d<Complex>& weights){
        StopWatch swatch_block;
        swatch_block.reset();

        //always boost the first propagator:
        swatch_block.start();
        if(mode.compare("000")==0){
            block=getHalfBaryonblock(LatticePropagator(phases*prop0),prop0,prop0,diquark_proj, snk, weights);
        }
        else if(mode.compare("001")==0){
            block=getHalfBaryonblock(LatticePropagator(phases*prop0),prop0,prop1,diquark_proj, snk, weights);
        }
        else if(mode.compare("010")==0){
            block=getHalfBaryonblock(LatticePropagator(phases*prop0),prop1,prop0,diquark_proj, snk, weights);
        }
        else if(mode.compare("011")==0){
            block=getHalfBaryonblock(LatticePropagator(phases*prop0),prop1,prop1,diquark_proj, snk, weights);
        }
        else if(mode.compare("100")==0){
            block=getHalfBaryonblock(LatticePropagator(phases*prop1),prop0,prop0,diquark_proj, snk, weights);
        }
        else if(mode.compare("101")==0){
            block=getHalfBaryonblock(LatticePropagator(phases*prop1),prop0,prop1,diquark_proj, snk, weights);
        }
        else if(mode.compare("110")==0){
            block=getHalfBaryonblock(LatticePropagator(phases*prop1),prop1,prop0,diquark_proj, snk, weights);
        }
        else if(mode.compare("111")==0){
            block=getHalfBaryonblock(LatticePropagator(phases*prop1),prop1,prop1,diquark_proj, snk, weights);
        }
        else{
            QDP_error_exit("error, unknown baryon block mode!");
        }
        swatch_block.stop();
        QDPIO::cout << "Build boosted block: time=" << swatch_block.getTimeInSeconds() << std::endl;
    }

    /* We can sink the propagators FIRST
     *
     *
     */

    void get_barblock(LatticeHalfBaryonblock& block, const std::string mode, const multi1d<LatticePropagator>& prop0, const multi1d<LatticePropagator>& prop1, const SpinMatrix& diquark_proj, const multi1d<Complex>& weights){
        StopWatch swatch_block;
        swatch_block.reset();

        //always boost the first propagator:
        swatch_block.start();
        if(mode.compare("000")==0){
            block=getHalfBaryonblock(prop0,prop0,prop0,diquark_proj, weights);
        }
        else if(mode.compare("001")==0){
            block=getHalfBaryonblock(prop0,prop0,prop1,diquark_proj, weights);
        }
        else if(mode.compare("010")==0){
            block=getHalfBaryonblock(prop0,prop1,prop0,diquark_proj, weights);
        }
        else if(mode.compare("011")==0){
            block=getHalfBaryonblock(prop0,prop1,prop1,diquark_proj, weights);
        }
        else if(mode.compare("100")==0){
            block=getHalfBaryonblock(prop1,prop0,prop0,diquark_proj, weights);
        }
        else if(mode.compare("101")==0){
            block=getHalfBaryonblock(prop1,prop0,prop1,diquark_proj, weights);
        }
        else if(mode.compare("110")==0){
            block=getHalfBaryonblock(prop1,prop1,prop0,diquark_proj, weights);
        }
        else if(mode.compare("111")==0){
            block=getHalfBaryonblock(prop1,prop1,prop1,diquark_proj, weights);
        }
        else{
            QDP_error_exit("error, unknown baryon block mode!");
        }
        swatch_block.stop();
        QDPIO::cout << "Build boosted block: time=" << swatch_block.getTimeInSeconds() << std::endl;
    }

    void get_barblock_boost(LatticeHalfBaryonblock& block, const std::string mode, const multi1d<LatticePropagator>& prop0, const multi1d<LatticePropagator>& prop1, const SpinMatrix& diquark_proj, const LatticeComplex& phases, const multi1d<Complex>& weights){
        StopWatch swatch_block;
        swatch_block.reset();

        multi1d<LatticePropagator> boosted(prop0.size());

        //always boost the first propagator:
        swatch_block.start();
        if(mode.compare("000")==0){
            for(int s=0; s<prop0.size(); s++) boosted[s]=LatticePropagator(phases*prop0[s]);
            block=getHalfBaryonblock(boosted,prop0,prop0,diquark_proj, weights);
        }
        else if(mode.compare("001")==0){
            for(int s=0; s<prop0.size(); s++) boosted[s]=LatticePropagator(phases*prop0[s]);
            block=getHalfBaryonblock(boosted,prop0,prop1,diquark_proj, weights);
        }
        else if(mode.compare("010")==0){
            for(int s=0; s<prop0.size(); s++) boosted[s]=LatticePropagator(phases*prop0[s]);
            block=getHalfBaryonblock(boosted,prop1,prop0,diquark_proj, weights);
        }
        else if(mode.compare("011")==0){
            for(int s=0; s<prop0.size(); s++) boosted[s]=LatticePropagator(phases*prop0[s]);
            block=getHalfBaryonblock(boosted,prop1,prop1,diquark_proj, weights);
        }
        else if(mode.compare("100")==0){
            for(int s=0; s<prop0.size(); s++) boosted[s]=LatticePropagator(phases*prop1[s]);
            block=getHalfBaryonblock(boosted,prop0,prop0,diquark_proj, weights);
        }
        else if(mode.compare("101")==0){
            for(int s=0; s<prop0.size(); s++) boosted[s]=LatticePropagator(phases*prop1[s]);
            block=getHalfBaryonblock(boosted,prop0,prop1,diquark_proj, weights);
        }
        else if(mode.compare("110")==0){
            for(int s=0; s<prop0.size(); s++) boosted[s]=LatticePropagator(phases*prop1[s]);
            block=getHalfBaryonblock(boosted,prop1,prop0,diquark_proj, weights);
        }
        else if(mode.compare("111")==0){
            for(int s=0; s<prop0.size(); s++) boosted[s]=LatticePropagator(phases*prop1[s]);
            block=getHalfBaryonblock(boosted,prop1,prop1,diquark_proj, weights);
        }
        else{
            QDP_error_exit("error, unknown baryon block mode!");
        }
        swatch_block.stop();
        QDPIO::cout << "Build boosted block: time=" << swatch_block.getTimeInSeconds() << std::endl;
    }




    /*
     *
     *
     */

    void initTopologies(const std::string& filename, const int& truncsize, const unsigned int& j_decay){
        //compute 3-Volume:
        vol3=1;
        for(int nu = 0; nu < Nd; ++nu){
            //skip time coordinate:
            if (nu == j_decay) continue;
            vol3*=Layout::lattSize()[nu];
        }

        std::vector<std::string> blockmodes_PP, blockmodes_PN;
        //assign blockmodes_PP={"000|111","001|101","001|110","010|110","010|101","100|011"}
        blockmodes_PP.resize(6);
        blockmodes_PP[0]="000|111"; blockmodes_PP[1]="001|101"; blockmodes_PP[2]="001|110"; blockmodes_PP[3]="010|110"; blockmodes_PP[4]="010|101"; blockmodes_PP[5]="100|011";

        blockmodes_PN.resize(9);
        blockmodes_PN[0]="001|110"; blockmodes_PN[1]="001|101"; blockmodes_PN[2]="010|101"; blockmodes_PN[3]="010|110"; blockmodes_PN[4]="011|001"; blockmodes_PN[5]="011|010"; blockmodes_PN[6]="101|100"; blockmodes_PN[7]="110|100"; blockmodes_PN[8]="111|000";


        QDPIO::cout << "Init baryon contraction tensors in double precision..." << std::endl;
        //single neutron:
        single_neutron=readTensor(filename,"/neutron/UDAV/000");

        //set global symmetries for PP topologies:
        PP_SING0.setSymmetry(+1);
        PP_TRIP0.setSymmetry(-1);
        PP_TRIPP.setSymmetry(-1);
        PP_TRIPM.setSymmetry(-1);

        //displaced tensors for PP:
        sparsearr2<Real> tensor;
        for(unsigned int b=0; b<blockmodes_PP.size(); b++){

            //singlet:
            tensor=readTensor(filename,"/neutron2/SING0/"+blockmodes_PP[b]).split_indices(halfredblocksize);
            PP_SING0.addDiagram(blockmodes_PP[b],tensor,+1);

            //triplet:
            tensor=readTensor(filename,"/neutron2/TRIP0/"+blockmodes_PP[b]).split_indices(halfredblocksize);
            PP_TRIP0.addDiagram(blockmodes_PP[b],tensor,+1);

            tensor=readTensor(filename,"/neutron2/TRIPP1/"+blockmodes_PP[b]).split_indices(halfredblocksize);
            PP_TRIPP.addDiagram(blockmodes_PP[b],tensor,+1);

            tensor=readTensor(filename,"/neutron2/TRIPM1/"+blockmodes_PP[b]).split_indices(halfredblocksize);
            PP_TRIPM.addDiagram(blockmodes_PP[b],tensor,+1);
        }
        //local s-wave tensor:
        tensor=readTensor(filename,"/neutron2/SING0/all").split_indices(halfredblocksize);
        PP_SING0.addDiagram("000|000",tensor,+1);


        //set global symmetries for PN topologies (here, the spin-triplet is combined with spatially even wf, i.e. it should be symmetric and the spin-singlet case always anti-symmetric):
        PN_SING0.setSymmetry(-1);
        PN_TRIP0.setSymmetry(+1);
        PN_TRIPP.setSymmetry(+1);
        PN_TRIPM.setSymmetry(+1);

        //displaced tensors for PN:
        for(unsigned int b=0; b<blockmodes_PN.size(); b++){

            //determine the sign of the mode: basically, when
            std::string firstmode=blockmodes_PN[b].substr(0,3);
            std::string secondmode=blockmodes_PN[b].substr(4,3);

            //compare to blockmodes PP and decide whether a minus 1 sign is needed
            int sign=+1;
            for(unsigned int c=0; c<blockmodes_PP.size(); c++){
                if(blockmodes_PP[c].find(firstmode)==4 || blockmodes_PP[c].find(secondmode)==0) sign=-1;
            }

            //singlet:
            tensor=readTensor(filename,"/deuteron/SING0/"+blockmodes_PN[b]).split_indices(halfredblocksize);
            PN_SING0.addDiagram(blockmodes_PN[b],tensor,sign);

            //triplet:
            tensor=readTensor(filename,"/deuteron/TRIP0/"+blockmodes_PN[b]).split_indices(halfredblocksize);
            PN_TRIP0.addDiagram(blockmodes_PN[b],tensor,sign);

            tensor=readTensor(filename,"/deuteron/TRIPP1/"+blockmodes_PN[b]).split_indices(halfredblocksize);
            PN_TRIPP.addDiagram(blockmodes_PN[b],tensor,sign);

            tensor=readTensor(filename,"/deuteron/TRIPM1/"+blockmodes_PN[b]).split_indices(halfredblocksize);
            PN_TRIPM.addDiagram(blockmodes_PN[b],tensor,sign);
        }
        //local s-wave tensors:
        tensor=readTensor(filename,"/deuteron/TRIP0/all").split_indices(halfredblocksize);
        PN_TRIP0.addDiagram("000|000",tensor,+1);
        tensor=readTensor(filename,"/deuteron/TRIPP1/all").split_indices(halfredblocksize);
        PN_TRIPP.addDiagram("000|000",tensor,+1);
        tensor=readTensor(filename,"/deuteron/TRIPM1/all").split_indices(halfredblocksize);
        PN_TRIPM.addDiagram("000|000",tensor,+1);

        QDPIO::cout << "done!" << std::endl;
    }

    void clearTopologies(){
        //PP
        PP_SING0.clear();
        PP_TRIP0.clear();
        PP_TRIPP.clear();
        PP_TRIPM.clear();

        //PN
        PN_SING0.clear();
        PN_TRIPP.clear();
        PN_TRIP0.clear();
        PN_TRIPM.clear();
    }


    //***********************************************************************************************************************************
    //***********************************************************************************************************************************
    //N->N  S = +/- 1/2
    //***********************************************************************************************************************************
    //***********************************************************************************************************************************
    //#ifndef CUDACOMP
    int one_proton(LatticeComplex& result, const LatticeHalfBaryonblock& protonblock){
        START_CODE();

        //protonblock is already in momentum space:
        QDPIO::cout << "Performing one_proton contractions on the CPU..." << std::flush;

        StopWatch swatch_cont;
        swatch_cont.reset();

        //init fields:
        result=zero;

        // Create the time-slice set
        swatch_cont.start();
        unsigned int nsites=Layout::sitesOnNode();
#pragma omp parallel firstprivate(nsites,zero,halfbarblocksize,single_neutron) shared(protonblock,result)
        {
            Complex* bblock;
            unsigned int i;
#pragma omp for schedule(static)
            for(i=0; i<nsites; i++){
                bblock=(Complex*)(&(protonblock.elem(i)));

                Complex cont=zero;
#ifndef MP
                for(unsigned int k=0; k<single_neutron.get_nvals(); k++){
                    cont+=bblock[single_neutron.ind[k]]*single_neutron.val[k];
                }
#else
                ComplexMP cnt(0.,0.,MPREC);
                for(unsigned int k=0; k<single_neutron.get_nvals(); k++){
                    ComplexMP bbl(toFloat(real(bblock[single_neutron.ind[k]])),toFloat(imag(bblock[single_neutron.ind[k]])),MPREC);
                    FloatMP tens(toFloat(single_neutron.val[k]),MPREC);
                    cnt+=bbl*tens;
                }
                cont=cnt.roundDouble();
#endif
                result.elem(i)=cont.elem();
            }
        }
        swatch_cont.stop();

        QDPIO::cout << "done!" << std::endl;
        QDPIO::cout << "P internal-contrac: time=" << swatch_cont.getTimeInSeconds() << std::endl;
        END_CODE();

        return EXIT_SUCCESS;
    }

    //***********************************************************************************************************************************
    //***********************************************************************************************************************************
    //NN->NN  L = 0; S = 0
    //***********************************************************************************************************************************
    //***********************************************************************************************************************************
    //constructs a local source:
    int two_proton_source_local(LatticeHalfSpinMatrix& result, const LatticeHalfBaryonblock& block_plus, const LatticeHalfBaryonblock& block_minus, const sparsearr2<Real>& tensor){
        START_CODE();

        QDPIO::cout << "Performing source-local two-proton contractions on the CPU..." << std::flush;

        //variables:
        unsigned int nsites=Layout::sitesOnNode();
        unsigned int count1, count2, offset1, offset2;
        count1=halfbarblocksize;
        count2=halfbarblocksize;
        offset1=halfredblocksize;
        offset2=halfredblocksize;

        //start contractions
        StopWatch swatch_cont;
        swatch_cont.reset();
        swatch_cont.start();

        //init result to zero:
        result=zero;
#pragma omp parallel firstprivate(Nshalf,nsites,count1,count2,offset1,offset2,zero) shared(tensor,block_plus,block_minus,result)
        {
            //init buffers
            Complex* bblock_p;
            Complex* bblock_m;
            Complex* cont;

            unsigned int i;
#pragma omp parallel for private(i) schedule(static)
            for(i=0; i<nsites; i++){
                bblock_p=(Complex*)(&(block_plus.elem(i)));
                bblock_m=(Complex*)(&(block_minus.elem(i)));
                cont=(Complex*)(&(result.elem(i)));

#ifndef MP
                for(unsigned int l=0; l<Nshalf; l++){
                    for(unsigned int m=0; m<Nshalf; m++){
                        cont[l+Nshalf*m]=zero;
                        for(unsigned int k=0; k<tensor.get_nvals(); k++){
                            unsigned int ind1=tensor.ind1[k]+offset1*l;
                            unsigned int ind2=tensor.ind2[k]+offset2*m;
                            cont[l+Nshalf*m]+=bblock_p[ind1]*bblock_m[ind2]*tensor.val[k];
                        }
                    }
                }
#else
                for(unsigned int l=0; l<Nshalf; l++){
                    for(unsigned int m=0; m<Nshalf; m++){
                        ComplexMP cnt(0.,0.,MPREC);
                        for(unsigned int k=0; k<tensor.get_nvals(); k++){
                            unsigned int ind1=tensor.ind1[k]+offset1*l;
                            unsigned int ind2=tensor.ind2[k]+offset2*m;
                            ComplexMP bblp(toFloat(real(bblock_p[ind1])),toFloat(imag(bblock_p[ind1])),MPREC);
                            ComplexMP bblm(toFloat(real(bblock_m[ind2])),toFloat(imag(bblock_m[ind2])),MPREC);
                            FloatMP tens(toFloat(tensor.val[k]),MPREC);
                            cnt+=bblp*bblm*tens;
                        }
                        cont[l+Nshalf*m]=cnt.roundDouble();
                    }
                }
#endif
            }
        }
        swatch_cont.stop();

        QDPIO::cout << "done!" << std::endl;
        QDPIO::cout << "NN-loc: internal-contrac: time=" << swatch_cont.getTimeInSeconds() << std::endl;

        END_CODE();

        return EXIT_SUCCESS;
    }

    //***********************************************************************************************************************************
    //***********************************************************************************************************************************
    //NN->NN  L = +/- 1; S = -/+ 1
    //***********************************************************************************************************************************
    //***********************************************************************************************************************************
    //constructs displaced sources:
    int two_proton_displaced(LatticeHalfSpinMatrix& result, const LatticeHalfBaryonblock& block0, const LatticeHalfBaryonblock& block1, const sparsearr2<Real>& tensor, const LatticeComplex& phases, Fourier& fft, const int& psign){
        START_CODE();

        if(psign!=0){
            unsigned int nsites=Layout::sitesOnNode();

            //now do contraction:
            StopWatch swatch_cont;
            swatch_cont.reset();

            LatticeHalfSpinMatrix tmpresult=zero;
            swatch_cont.start();
#pragma omp parallel firstprivate(Nshalf,nsites,halfbarblocksize,halfredblocksize,zero) shared(tensor,block0,block1,result)
            {
                //init buffers
                Complex* bblock_p;
                Complex* bblock_m;
                Complex* cont;

                unsigned int i;
#pragma omp for private(i) schedule(static)
                for(i=0; i<nsites; i++){
                    bblock_p=(Complex*)(&(block0.elem(i)));
                    bblock_m=(Complex*)(&(block1.elem(i)));
                    cont=(Complex*)(&(tmpresult.elem(i)));

#ifndef MP
                    for(unsigned int l=0; l<Nshalf; l++){
                        for(unsigned int m=0; m<Nshalf; m++){
                            cont[l+Nshalf*m]=zero;
                            for(unsigned int k=0; k<tensor.get_nvals(); k++){
                                unsigned int ind1=tensor.ind1[k]+halfredblocksize*l;
                                unsigned int ind2=tensor.ind2[k]+halfredblocksize*m;
                                cont[l+Nshalf*m]+=bblock_p[ind1]*bblock_m[ind2]*tensor.val[k];
                            }
                        }
                    }
#else
                    for(unsigned int l=0; l<Nshalf; l++){
                        for(unsigned int m=0; m<Nshalf; m++){
                            ComplexMP cnt(0.,0.,MPREC);
                            for(unsigned int k=0; k<tensor[b].get_nvals(); k++){
                                unsigned int ind1=tensor[b].ind1[k]+halfredblocksize*l;
                                unsigned int ind2=tensor[b].ind2[k]+halfredblocksize*m;
                                ComplexMP bblp(toFloat(real(bblock_p[ind1])),toFloat(imag(bblock_p[ind1])),MPREC);
                                ComplexMP bblm(toFloat(real(bblock_m[ind2])),toFloat(imag(bblock_m[ind2])),MPREC);
                                FloatMP tens(toFloat(tensor[b].val[k]),MPREC);
                                cnt+=bblp*bblm*tens;
                            }
                            cont[l+Nshalf*m]=cnt.roundDouble();
                        }
                    }
#endif
                }
            }
            swatch_cont.stop();

            //fft if -p -> +p:
            StopWatch swatch_fft;
            if(psign==-1){
                swatch_fft.reset();
                LatticeHalfSpinMatrix tmpfft;

                //now invert argument, i.e. go from B1(q+P)*B2(-q) -> B1(x-r)*B2(x)
                swatch_fft.start();
                tmpfft=fft(tmpresult,+1);
                //multiply phase factors and do another FFT, effectively going from B1(x-r)*B2(x) -> B2(q+P)*B1(-q)
                result+=fft(LatticeHalfSpinMatrix(phases*tmpfft),+1);
                swatch_fft.stop();
            }
            else{
                result+=tmpresult;
            }

            QDPIO::cout << "done!" << std::endl;
            QDPIO::cout << "NN-disp: internal-contrac: time=" << swatch_cont.getTimeInSeconds() << std::endl;
            if(psign==-1) QDPIO::cout << "NN-disp: internal-fourier: time=" << swatch_fft.getTimeInSeconds() << std::endl;
        }

        END_CODE();

        return EXIT_SUCCESS;
    }

    //***********************************************************************************************************************************
    //***********************************************************************************************************************************
    //Applies symmetry transformation on matrix
    //***********************************************************************************************************************************
    //***********************************************************************************************************************************
    //for unit phases (zero boost) it does not matter whether (+1,+1) or (-1,-1) is used. However, since phase belongs to +1 and adj(phase) to -1, it makes a difference in a boosted frame
    void symmetrize(LatticeHalfSpinMatrix& result, const LatticeComplex& phases, Fourier& fft, const int& sign){
        START_CODE();

        //symmetrization (sign=+1) or anti-symmetrization (sign=-1):
        LatticeHalfSpinMatrix tmp;
        tmp=fft(result,+1);
        //multiply phase factors accounting for boosted systems
        result+=sign*fft(LatticeHalfSpinMatrix(phases*tmp),+1);

        END_CODE();
    }

    //*********************************************************************************************************
    //*********************************************************************************************************
    //Combined contractions:
    //*********************************************************************************************************
    //*********************************************************************************************************
    int contract(LatticeComplex& result_P, std::map<std::string,LatticeHalfSpinMatrix>& resultmats, const LatticePropagator& prop0, const LatticePropagator& prop1, const SpinMatrix& diquark_proj, const LatticeComplex& phases, Fourier& fft, const bool& compute_locals){

        START_CODE();
        QDPIO::cout << "Performing source-displaced two-proton contractions on the CPU..." << std::flush;

        //measure fft performance:
        StopWatch swatch_fft, swatch_cont;
        swatch_fft.reset();
        swatch_cont.reset();

        //initialize variables and set all input fields to zero:
        result_P=zero;
        for(std::map<std::string,LatticeHalfSpinMatrix>::iterator it=resultmats.begin(); it!=resultmats.end(); ++it){
            it->second=zero;
        }
        LatticeHalfBaryonblock block0, block1;

        if(compute_locals){
            //get baryon block combination: block0 is boosted
            // change barblock to take 3 props
            get_barblock_boost(block0,"000",prop0,prop1,diquark_proj,phases);
            get_barblock(block1,"000",prop0,prop1,diquark_proj);

            //fourier:
            swatch_fft.start();
            block0=fft(block0,+1);
            block1=fft(block1,-1);
            swatch_fft.stop();

            //contract PP
            swatch_cont.start();
            one_proton(result_P,block0);
            two_proton_source_local(resultmats["PP_SING0_loc"],block0,block1,PP_SING0.getTensor("000|000"));

            //contract PN
            two_proton_source_local(resultmats["PN_TRIPP_loc"],block0,block1,PN_TRIPP.getTensor("000|000"));
            two_proton_source_local(resultmats["PN_TRIP0_loc"],block0,block1,PN_TRIP0.getTensor("000|000"));
            two_proton_source_local(resultmats["PN_TRIPM_loc"],block0,block1,PN_TRIPM.getTensor("000|000"));
            swatch_cont.stop();
        }

        //store the blocks 101, 110, 011:
        std::vector<std::string> bmodes(3);
        bmodes[0]="101";
        bmodes[1]="110";
        bmodes[2]="011";
        std::map<std::string, LatticeHalfBaryonblock> bblockmap;
        for(unsigned int b=0; b<bmodes.size(); b++){
            get_barblock(block1,bmodes[b],prop0,prop1,diquark_proj);
            swatch_fft.start();
            block1=fft(block1,-1);
            swatch_fft.stop();
            bblockmap[bmodes[b]]=block1;
        }

        //do 000 with 111 first:
        if(!compute_locals){
            //first block is boosted:
            get_barblock_boost(block0,"000",prop0,prop1,diquark_proj,phases);
            swatch_fft.start();
            block0=fft(block0,+1);
            swatch_fft.stop();
        }
        get_barblock(block1,"111",prop0,prop1,diquark_proj);
        swatch_fft.start();
        block1=fft(block1,-1);
        swatch_fft.stop();

        //contractions
        swatch_cont.start();
        //PP
        //s-wave:
        two_proton_displaced(resultmats["PP_SING0"],block0,block1,PP_SING0.getTensor("000|111"),phases,fft,PP_SING0.getFourierSign("000|111"));

        //P+1
        two_proton_displaced(resultmats["PP_TRIPP"],block0,block1,PP_TRIPP.getTensor("000|111"),phases,fft,PP_TRIPP.getFourierSign("000|111"));

        //P0
        two_proton_displaced(resultmats["PP_TRIP0"],block0,block1,PP_TRIP0.getTensor("000|111"),phases,fft,PP_TRIP0.getFourierSign("000|111"));

        //P-1
        two_proton_displaced(resultmats["PP_TRIPM"],block0,block1,PP_TRIPM.getTensor("000|111"),phases,fft,PP_TRIPM.getFourierSign("000|111"));


        //PN
        //s-wave:
        two_proton_displaced(resultmats["PN_SING0"],block0,block1,PN_SING0.getTensor("111|000"),phases,fft,PN_SING0.getFourierSign("111|000"));

        //P+1
        two_proton_displaced(resultmats["PN_TRIPP"],block0,block1,PN_TRIPP.getTensor("111|000"),phases,fft,PN_TRIPP.getFourierSign("111|000"));

        //P0
        two_proton_displaced(resultmats["PN_TRIP0"],block0,block1,PN_TRIP0.getTensor("111|000"),phases,fft,PN_TRIP0.getFourierSign("111|000"));

        //P-1
        two_proton_displaced(resultmats["PN_TRIPM"],block0,block1,PN_TRIPM.getTensor("111|000"),phases,fft,PN_TRIPM.getFourierSign("111|000"));
        swatch_cont.stop();


        //do the rest with all the others
        std::vector<std::string> bmodes2(3);
        bmodes2[0]="001";
        bmodes2[1]="010";
        bmodes2[2]="100";

        for(unsigned int c=0; c<bmodes2.size(); c++){
            //block0 is boosted:
            get_barblock_boost(block0,bmodes2[c],prop0,prop1,diquark_proj,phases);
            swatch_fft.start();
            block0=fft(block0,+1);
            swatch_fft.stop();

            swatch_cont.start();
            std::string mode, swapmode;
            for(unsigned int b=0; b<bmodes.size(); b++){
                //PP
                //everything is "normalized" to PP, so the string and orders of blocks can stay as it is:
                mode=bmodes2[c]+"|"+bmodes[b];
                swapmode=bmodes[b]+"|"+bmodes2[c];
                //s-wave:
                two_proton_displaced(resultmats["PP_SING0"],block0,bblockmap[bmodes[b]],PP_SING0.getTensor(mode),phases,fft,PP_SING0.getFourierSign(mode));

                //P+1
                two_proton_displaced(resultmats["PP_TRIPP"],block0,bblockmap[bmodes[b]],PP_TRIPP.getTensor(mode),phases,fft,PP_TRIPP.getFourierSign(mode));

                //P0
                two_proton_displaced(resultmats["PP_TRIP0"],block0,bblockmap[bmodes[b]],PP_TRIP0.getTensor(mode),phases,fft,PP_TRIP0.getFourierSign(mode));

                //P-1
                two_proton_displaced(resultmats["PP_TRIPM"],block0,bblockmap[bmodes[b]],PP_TRIPM.getTensor(mode),phases,fft,PP_TRIPM.getFourierSign(mode));

                //PN
                //here we might need to swap the blocks:
                //s-wave:
                if(PN_SING0.getFourierSign(mode)==0){
                    two_proton_displaced(resultmats["PN_SING0"],bblockmap[bmodes[b]],block0,PN_SING0.getTensor(swapmode),phases,fft,PN_SING0.getFourierSign(swapmode));
                }
                else{
                    two_proton_displaced(resultmats["PN_SING0"],block0,bblockmap[bmodes[b]],PN_SING0.getTensor(mode),phases,fft,PN_SING0.getFourierSign(mode));
                }

                //P+1
                if(PN_TRIPP.getFourierSign(mode)==0){
                    two_proton_displaced(resultmats["PN_TRIPP"],bblockmap[bmodes[b]],block0,PN_TRIPP.getTensor(swapmode),phases,fft,PN_TRIPP.getFourierSign(swapmode));
                }
                else{
                    two_proton_displaced(resultmats["PN_TRIPP"],block0,bblockmap[bmodes[b]],PN_TRIPP.getTensor(mode),phases,fft,PN_TRIPP.getFourierSign(mode));
                }

                //P0
                if(PN_TRIP0.getFourierSign(mode)==0){
                    two_proton_displaced(resultmats["PN_TRIP0"],bblockmap[bmodes[b]],block0,PN_TRIP0.getTensor(swapmode),phases,fft,PN_TRIP0.getFourierSign(swapmode));
                }
                else{
                    two_proton_displaced(resultmats["PN_TRIP0"],block0,bblockmap[bmodes[b]],PN_TRIP0.getTensor(mode),phases,fft,PN_TRIP0.getFourierSign(mode));
                }

                //P-1
                if(PN_TRIPM.getFourierSign(mode)==0){
                    two_proton_displaced(resultmats["PN_TRIPM"],bblockmap[bmodes[b]],block0,PN_TRIPM.getTensor(swapmode),phases,fft,PN_TRIPM.getFourierSign(swapmode));
                }
                else{
                    two_proton_displaced(resultmats["PN_TRIPM"],block0,bblockmap[bmodes[b]],PN_TRIPM.getTensor(mode),phases,fft,PN_TRIPM.getFourierSign(mode));
                }
            }
            swatch_cont.stop();
        }


        //do symmetrization or anti-symmetrization of result vectors:
        swatch_fft.start();
        //PP
        symmetrize(resultmats["PP_SING0"],phases,fft,PP_SING0.getSymmetry());
        symmetrize(resultmats["PP_TRIPP"],phases,fft,PP_TRIPP.getSymmetry());
        symmetrize(resultmats["PP_TRIP0"],phases,fft,PP_TRIP0.getSymmetry());
        symmetrize(resultmats["PP_TRIPM"],phases,fft,PP_TRIPM.getSymmetry());

        //PN
        symmetrize(resultmats["PN_SING0"],phases,fft,PN_SING0.getSymmetry());
        symmetrize(resultmats["PN_TRIPP"],phases,fft,PN_TRIPP.getSymmetry());
        symmetrize(resultmats["PN_TRIP0"],phases,fft,PN_TRIP0.getSymmetry());
        symmetrize(resultmats["PN_TRIPM"],phases,fft,PN_TRIPM.getSymmetry());
        swatch_fft.stop();

        QDPIO::cout << "done!" << std::endl;
        QDPIO::cout << "Contract: fourier: time=" << swatch_fft.getTimeInSeconds() << std::endl;
        QDPIO::cout << "Contract: contrac: time=" << swatch_cont.getTimeInSeconds() << std::endl;
        END_CODE();

        return EXIT_SUCCESS;
    }

    // With arbitrary sink
    int contract(LatticeComplex& result_P, std::map<std::string,LatticeHalfSpinMatrix>& resultmats, const LatticePropagator& prop0, const LatticePropagator& prop1, const SpinMatrix& diquark_proj, const LatticeComplex& phases, Fourier& fft, const bool& compute_locals, const sink& snk){

        START_CODE();
        QDPIO::cout << "Performing source-displaced two-proton contractions on the CPU..." << std::flush;

        //measure fft performance:
        StopWatch swatch_fft, swatch_cont;
        swatch_fft.reset();
        swatch_cont.reset();

        //initialize variables and set all input fields to zero:
        result_P=zero;
        for(std::map<std::string,LatticeHalfSpinMatrix>::iterator it=resultmats.begin(); it!=resultmats.end(); ++it){
            it->second=zero;
        }
        LatticeHalfBaryonblock block0, block1;

        if(compute_locals){
            //get baryon block combination: block0 is boosted
            get_barblock_boost(block0,"000",prop0,prop1,diquark_proj,phases,snk);
            get_barblock(block1,"000",prop0,prop1,diquark_proj,snk);

            //fourier:
            swatch_fft.start();
            block0=fft(block0,+1);
            block1=fft(block1,-1);
            swatch_fft.stop();

            //contract PP
            swatch_cont.start();
            one_proton(result_P,block0);
            two_proton_source_local(resultmats["PP_SING0_loc"],block0,block1,PP_SING0.getTensor("000|000"));

            //contract PN
            two_proton_source_local(resultmats["PN_TRIPP_loc"],block0,block1,PN_TRIPP.getTensor("000|000"));
            two_proton_source_local(resultmats["PN_TRIP0_loc"],block0,block1,PN_TRIP0.getTensor("000|000"));
            two_proton_source_local(resultmats["PN_TRIPM_loc"],block0,block1,PN_TRIPM.getTensor("000|000"));
            swatch_cont.stop();
        }

        //store the blocks 101, 110, 011:
        std::vector<std::string> bmodes(3);
        bmodes[0]="101";
        bmodes[1]="110";
        bmodes[2]="011";
        std::map<std::string, LatticeHalfBaryonblock> bblockmap;
        for(unsigned int b=0; b<bmodes.size(); b++){
            get_barblock(block1,bmodes[b],prop0,prop1,diquark_proj,snk);
            swatch_fft.start();
            block1=fft(block1,-1);
            swatch_fft.stop();
            bblockmap[bmodes[b]]=block1;
        }

        //do 000 with 111 first:
        if(!compute_locals){
            //first block is boosted:
            get_barblock_boost(block0,"000",prop0,prop1,diquark_proj,phases,snk);
            swatch_fft.start();
            block0=fft(block0,+1);
            swatch_fft.stop();
        }
        get_barblock(block1,"111",prop0,prop1,diquark_proj,snk);
        swatch_fft.start();
        block1=fft(block1,-1);
        swatch_fft.stop();

        //contractions
        swatch_cont.start();
        //PP
        //s-wave:
        two_proton_displaced(resultmats["PP_SING0"],block0,block1,PP_SING0.getTensor("000|111"),phases,fft,PP_SING0.getFourierSign("000|111"));

        //P+1
        two_proton_displaced(resultmats["PP_TRIPP"],block0,block1,PP_TRIPP.getTensor("000|111"),phases,fft,PP_TRIPP.getFourierSign("000|111"));

        //P0
        two_proton_displaced(resultmats["PP_TRIP0"],block0,block1,PP_TRIP0.getTensor("000|111"),phases,fft,PP_TRIP0.getFourierSign("000|111"));

        //P-1
        two_proton_displaced(resultmats["PP_TRIPM"],block0,block1,PP_TRIPM.getTensor("000|111"),phases,fft,PP_TRIPM.getFourierSign("000|111"));


        //PN
        //s-wave:
        two_proton_displaced(resultmats["PN_SING0"],block0,block1,PN_SING0.getTensor("111|000"),phases,fft,PN_SING0.getFourierSign("111|000"));

        //P+1
        two_proton_displaced(resultmats["PN_TRIPP"],block0,block1,PN_TRIPP.getTensor("111|000"),phases,fft,PN_TRIPP.getFourierSign("111|000"));

        //P0
        two_proton_displaced(resultmats["PN_TRIP0"],block0,block1,PN_TRIP0.getTensor("111|000"),phases,fft,PN_TRIP0.getFourierSign("111|000"));

        //P-1
        two_proton_displaced(resultmats["PN_TRIPM"],block0,block1,PN_TRIPM.getTensor("111|000"),phases,fft,PN_TRIPM.getFourierSign("111|000"));
        swatch_cont.stop();


        //do the rest with all the others
        std::vector<std::string> bmodes2(3);
        bmodes2[0]="001";
        bmodes2[1]="010";
        bmodes2[2]="100";

        for(unsigned int c=0; c<bmodes2.size(); c++){
            //block0 is boosted:
            get_barblock_boost(block0,bmodes2[c],prop0,prop1,diquark_proj,phases,snk);
            swatch_fft.start();
            block0=fft(block0,+1);
            swatch_fft.stop();

            swatch_cont.start();
            std::string mode, swapmode;
            for(unsigned int b=0; b<bmodes.size(); b++){
                //PP
                //everything is "normalized" to PP, so the string and orders of blocks can stay as it is:
                mode=bmodes2[c]+"|"+bmodes[b];
                swapmode=bmodes[b]+"|"+bmodes2[c];
                //s-wave:
                two_proton_displaced(resultmats["PP_SING0"],block0,bblockmap[bmodes[b]],PP_SING0.getTensor(mode),phases,fft,PP_SING0.getFourierSign(mode));

                //P+1
                two_proton_displaced(resultmats["PP_TRIPP"],block0,bblockmap[bmodes[b]],PP_TRIPP.getTensor(mode),phases,fft,PP_TRIPP.getFourierSign(mode));

                //P0
                two_proton_displaced(resultmats["PP_TRIP0"],block0,bblockmap[bmodes[b]],PP_TRIP0.getTensor(mode),phases,fft,PP_TRIP0.getFourierSign(mode));

                //P-1
                two_proton_displaced(resultmats["PP_TRIPM"],block0,bblockmap[bmodes[b]],PP_TRIPM.getTensor(mode),phases,fft,PP_TRIPM.getFourierSign(mode));


                //PN
                //here we might need to swap the blocks:
                //s-wave:
                if(PN_SING0.getFourierSign(mode)==0){
                    two_proton_displaced(resultmats["PN_SING0"],bblockmap[bmodes[b]],block0,PN_SING0.getTensor(swapmode),phases,fft,PN_SING0.getFourierSign(swapmode));
                }
                else{
                    two_proton_displaced(resultmats["PN_SING0"],block0,bblockmap[bmodes[b]],PN_SING0.getTensor(mode),phases,fft,PN_SING0.getFourierSign(mode));
                }

                //P+1
                if(PN_TRIPP.getFourierSign(mode)==0){
                    two_proton_displaced(resultmats["PN_TRIPP"],bblockmap[bmodes[b]],block0,PN_TRIPP.getTensor(swapmode),phases,fft,PN_TRIPP.getFourierSign(swapmode));
                }
                else{
                    two_proton_displaced(resultmats["PN_TRIPP"],block0,bblockmap[bmodes[b]],PN_TRIPP.getTensor(mode),phases,fft,PN_TRIPP.getFourierSign(mode));
                }

                //P0
                if(PN_TRIP0.getFourierSign(mode)==0){
                    two_proton_displaced(resultmats["PN_TRIP0"],bblockmap[bmodes[b]],block0,PN_TRIP0.getTensor(swapmode),phases,fft,PN_TRIP0.getFourierSign(swapmode));
                }
                else{
                    two_proton_displaced(resultmats["PN_TRIP0"],block0,bblockmap[bmodes[b]],PN_TRIP0.getTensor(mode),phases,fft,PN_TRIP0.getFourierSign(mode));
                }

                //P-1
                if(PN_TRIPM.getFourierSign(mode)==0){
                    two_proton_displaced(resultmats["PN_TRIPM"],bblockmap[bmodes[b]],block0,PN_TRIPM.getTensor(swapmode),phases,fft,PN_TRIPM.getFourierSign(swapmode));
                }
                else{
                    two_proton_displaced(resultmats["PN_TRIPM"],block0,bblockmap[bmodes[b]],PN_TRIPM.getTensor(mode),phases,fft,PN_TRIPM.getFourierSign(mode));
                }
            }
            swatch_cont.stop();
        }


        //do symmetrization or anti-symmetrization of result vectors:
        swatch_fft.start();
        //PP
        symmetrize(resultmats["PP_SING0"],phases,fft,PP_SING0.getSymmetry());
        symmetrize(resultmats["PP_TRIPP"],phases,fft,PP_TRIPP.getSymmetry());
        symmetrize(resultmats["PP_TRIP0"],phases,fft,PP_TRIP0.getSymmetry());
        symmetrize(resultmats["PP_TRIPM"],phases,fft,PP_TRIPM.getSymmetry());

        //PN
        symmetrize(resultmats["PN_SING0"],phases,fft,PN_SING0.getSymmetry());
        symmetrize(resultmats["PN_TRIPP"],phases,fft,PN_TRIPP.getSymmetry());
        symmetrize(resultmats["PN_TRIP0"],phases,fft,PN_TRIP0.getSymmetry());
        symmetrize(resultmats["PN_TRIPM"],phases,fft,PN_TRIPM.getSymmetry());
        swatch_fft.stop();

        QDPIO::cout << "done!" << std::endl;
        QDPIO::cout << "Contract: fourier: time=" << swatch_fft.getTimeInSeconds() << std::endl;
        QDPIO::cout << "Contract: contrac: time=" << swatch_cont.getTimeInSeconds() << std::endl;
        END_CODE();

        return EXIT_SUCCESS;
    }


    // With linear combo sink
    int contract(
                 LatticeComplex& result_P,
                 std::map<std::string,
                 LatticeHalfSpinMatrix>& resultmats,
                 const LatticePropagator& prop0,
                 const LatticePropagator& prop1,
                 const SpinMatrix& diquark_proj,
                 const LatticeComplex& phases,
                 Fourier& fft,
                 const bool& compute_locals,
                 const multi1d<sink*>& snk,
                 const multi1d<Complex>& weights){

        START_CODE();
        QDPIO::cout << "Performing source-displaced two-proton contractions on the CPU..." << std::flush;

        //measure fft performance:
        StopWatch swatch_fft, swatch_cont, swatch_proton;
        swatch_fft.reset();
        swatch_cont.reset();
        swatch_proton.reset();

        //initialize variables and set all input fields to zero:
        result_P=zero;
        for(std::map<std::string,LatticeHalfSpinMatrix>::iterator it=resultmats.begin(); it!=resultmats.end(); ++it){
            it->second=zero;
        }
        LatticeHalfBaryonblock block0, block1;

        if(compute_locals){
            //get baryon block combination: block0 is boosted
            get_barblock_boost(block0,"000",prop0,prop1,diquark_proj,phases,snk,weights);
            get_barblock(block1,"000",prop0,prop1,diquark_proj,snk,weights);

            //fourier:
            swatch_fft.start();
            block0=fft(block0,+1);
            block1=fft(block1,-1);
            swatch_fft.stop();

            //contract PP
            swatch_cont.start();
            one_proton(result_P,block0);
            two_proton_source_local(resultmats["PP_SING0_loc"],block0,block1,PP_SING0.getTensor("000|000"));

            //contract PN
            two_proton_source_local(resultmats["PN_TRIPP_loc"],block0,block1,PN_TRIPP.getTensor("000|000"));
            two_proton_source_local(resultmats["PN_TRIP0_loc"],block0,block1,PN_TRIP0.getTensor("000|000"));
            two_proton_source_local(resultmats["PN_TRIPM_loc"],block0,block1,PN_TRIPM.getTensor("000|000"));
            swatch_cont.stop();
        }

        //store the blocks 101, 110, 011:
        std::vector<std::string> bmodes(3);
        bmodes[0]="101";
        bmodes[1]="110";
        bmodes[2]="011";
        std::map<std::string, LatticeHalfBaryonblock> bblockmap;
        for(unsigned int b=0; b<bmodes.size(); b++){
            get_barblock(block1,bmodes[b],prop0,prop1,diquark_proj,snk,weights);
            swatch_fft.start();
            block1=fft(block1,-1);
            swatch_fft.stop();
            bblockmap[bmodes[b]]=block1;
        }

        //do 000 with 111 first:
        if(!compute_locals){
            //first block is boosted:
            get_barblock_boost(block0,"000",prop0,prop1,diquark_proj,phases,snk,weights);
            swatch_fft.start();
            block0=fft(block0,+1);
            swatch_fft.stop();
        }
        get_barblock(block1,"111",prop0,prop1,diquark_proj,snk,weights);
        swatch_fft.start();
        block1=fft(block1,-1);
        swatch_fft.stop();

        //single nucleon from block1 (it has no momentum shift)
        swatch_proton.start();
        one_proton(result_P,block1);
        swatch_proton.stop();

        //contractions
        swatch_cont.start();
        //PP
        //s-wave:
        two_proton_displaced(resultmats["PP_SING0"],block0,block1,PP_SING0.getTensor("000|111"),phases,fft,PP_SING0.getFourierSign("000|111"));

        //P+1
        two_proton_displaced(resultmats["PP_TRIPP"],block0,block1,PP_TRIPP.getTensor("000|111"),phases,fft,PP_TRIPP.getFourierSign("000|111"));

        //P0
        two_proton_displaced(resultmats["PP_TRIP0"],block0,block1,PP_TRIP0.getTensor("000|111"),phases,fft,PP_TRIP0.getFourierSign("000|111"));

        //P-1
        two_proton_displaced(resultmats["PP_TRIPM"],block0,block1,PP_TRIPM.getTensor("000|111"),phases,fft,PP_TRIPM.getFourierSign("000|111"));


        //PN
        //s-wave:
        two_proton_displaced(resultmats["PN_SING0"],block0,block1,PN_SING0.getTensor("111|000"),phases,fft,PN_SING0.getFourierSign("111|000"));

        //P+1
        two_proton_displaced(resultmats["PN_TRIPP"],block0,block1,PN_TRIPP.getTensor("111|000"),phases,fft,PN_TRIPP.getFourierSign("111|000"));

        //P0
        two_proton_displaced(resultmats["PN_TRIP0"],block0,block1,PN_TRIP0.getTensor("111|000"),phases,fft,PN_TRIP0.getFourierSign("111|000"));

        //P-1
        two_proton_displaced(resultmats["PN_TRIPM"],block0,block1,PN_TRIPM.getTensor("111|000"),phases,fft,PN_TRIPM.getFourierSign("111|000"));
        swatch_cont.stop();


        //do the rest with all the others
        std::vector<std::string> bmodes2(3);
        bmodes2[0]="001";
        bmodes2[1]="010";
        bmodes2[2]="100";

        for(unsigned int c=0; c<bmodes2.size(); c++){
            //block0 is boosted:
            get_barblock_boost(block0,bmodes2[c],prop0,prop1,diquark_proj,phases,snk,weights);
            swatch_fft.start();
            block0=fft(block0,+1);
            swatch_fft.stop();

            swatch_cont.start();
            std::string mode, swapmode;
            for(unsigned int b=0; b<bmodes.size(); b++){
                //PP
                //everything is "normalized" to PP, so the string and orders of blocks can stay as it is:
                mode=bmodes2[c]+"|"+bmodes[b];
                swapmode=bmodes[b]+"|"+bmodes2[c];
                //s-wave:
                two_proton_displaced(resultmats["PP_SING0"],block0,bblockmap[bmodes[b]],PP_SING0.getTensor(mode),phases,fft,PP_SING0.getFourierSign(mode));

                //P+1
                two_proton_displaced(resultmats["PP_TRIPP"],block0,bblockmap[bmodes[b]],PP_TRIPP.getTensor(mode),phases,fft,PP_TRIPP.getFourierSign(mode));

                //P0
                two_proton_displaced(resultmats["PP_TRIP0"],block0,bblockmap[bmodes[b]],PP_TRIP0.getTensor(mode),phases,fft,PP_TRIP0.getFourierSign(mode));

                //P-1
                two_proton_displaced(resultmats["PP_TRIPM"],block0,bblockmap[bmodes[b]],PP_TRIPM.getTensor(mode),phases,fft,PP_TRIPM.getFourierSign(mode));


                //PN
                //here we might need to swap the blocks:
                //s-wave:
                if(PN_SING0.getFourierSign(mode)==0){
                    two_proton_displaced(resultmats["PN_SING0"],bblockmap[bmodes[b]],block0,PN_SING0.getTensor(swapmode),phases,fft,PN_SING0.getFourierSign(swapmode));
                }
                else{
                    two_proton_displaced(resultmats["PN_SING0"],block0,bblockmap[bmodes[b]],PN_SING0.getTensor(mode),phases,fft,PN_SING0.getFourierSign(mode));
                }

                //P+1
                if(PN_TRIPP.getFourierSign(mode)==0){
                    two_proton_displaced(resultmats["PN_TRIPP"],bblockmap[bmodes[b]],block0,PN_TRIPP.getTensor(swapmode),phases,fft,PN_TRIPP.getFourierSign(swapmode));
                }
                else{
                    two_proton_displaced(resultmats["PN_TRIPP"],block0,bblockmap[bmodes[b]],PN_TRIPP.getTensor(mode),phases,fft,PN_TRIPP.getFourierSign(mode));
                }

                //P0
                if(PN_TRIP0.getFourierSign(mode)==0){
                    two_proton_displaced(resultmats["PN_TRIP0"],bblockmap[bmodes[b]],block0,PN_TRIP0.getTensor(swapmode),phases,fft,PN_TRIP0.getFourierSign(swapmode));
                }
                else{
                    two_proton_displaced(resultmats["PN_TRIP0"],block0,bblockmap[bmodes[b]],PN_TRIP0.getTensor(mode),phases,fft,PN_TRIP0.getFourierSign(mode));
                }

                //P-1
                if(PN_TRIPM.getFourierSign(mode)==0){
                    two_proton_displaced(resultmats["PN_TRIPM"],bblockmap[bmodes[b]],block0,PN_TRIPM.getTensor(swapmode),phases,fft,PN_TRIPM.getFourierSign(swapmode));
                }
                else{
                    two_proton_displaced(resultmats["PN_TRIPM"],block0,bblockmap[bmodes[b]],PN_TRIPM.getTensor(mode),phases,fft,PN_TRIPM.getFourierSign(mode));
                }
            }
            swatch_cont.stop();
        }


        //do symmetrization or anti-symmetrization of result vectors:
        swatch_fft.start();
        //PP
        symmetrize(resultmats["PP_SING0"],phases,fft,PP_SING0.getSymmetry());
        symmetrize(resultmats["PP_TRIPP"],phases,fft,PP_TRIPP.getSymmetry());
        symmetrize(resultmats["PP_TRIP0"],phases,fft,PP_TRIP0.getSymmetry());
        symmetrize(resultmats["PP_TRIPM"],phases,fft,PP_TRIPM.getSymmetry());

        //PN
        symmetrize(resultmats["PN_SING0"],phases,fft,PN_SING0.getSymmetry());
        symmetrize(resultmats["PN_TRIPP"],phases,fft,PN_TRIPP.getSymmetry());
        symmetrize(resultmats["PN_TRIP0"],phases,fft,PN_TRIP0.getSymmetry());
        symmetrize(resultmats["PN_TRIPM"],phases,fft,PN_TRIPM.getSymmetry());
        swatch_fft.stop();

        QDPIO::cout << "done!" << std::endl;
        QDPIO::cout << "Contract: fourier: time=" << swatch_fft.getTimeInSeconds()    << std::endl;
        QDPIO::cout << "Contract: proton: time="  << swatch_proton.getTimeInSeconds() << std::endl;
        QDPIO::cout << "Contract: contrac: time=" << swatch_cont.getTimeInSeconds()   << std::endl;
        END_CODE();

        return EXIT_SUCCESS;
    }

    // With linear combo sink
    int contract_local(LatticeComplex& result_P, std::map<std::string,LatticeHalfSpinMatrix>& resultmats, const LatticePropagator& prop0, const LatticePropagator& prop1, const SpinMatrix& diquark_proj, const LatticeComplex& phases, Fourier& fft, const multi1d<sink*>& snk, const multi1d<Complex>& weights){

        START_CODE();
        QDPIO::cout << "Performing source-displaced two-proton contractions on the CPU..." << std::flush;

        //measure fft performance:
        StopWatch swatch_fft, swatch_cont;
        swatch_fft.reset();
        swatch_cont.reset();

        //initialize variables and set all input fields to zero:
        result_P=zero;
        for(std::map<std::string,LatticeHalfSpinMatrix>::iterator it=resultmats.begin(); it!=resultmats.end(); ++it){
            it->second=zero;
        }
        LatticeHalfBaryonblock block0, block1;

        //get baryon block combination: block0 is boosted
        get_barblock_boost(block0,"000",prop0,prop1,diquark_proj,phases,snk,weights);
        get_barblock(block1,"000",prop0,prop1,diquark_proj,snk,weights);

        //fourier:
        swatch_fft.start();
        block0=fft(block0,+1);
        block1=fft(block1,-1);
        swatch_fft.stop();

        //contract PP
        swatch_cont.start();
        one_proton(result_P,block0);
        two_proton_source_local(resultmats["PP_SING0_loc"],block0,block1,PP_SING0.getTensor("000|000"));

        //contract PN
        two_proton_source_local(resultmats["PN_TRIPP_loc"],block0,block1,PN_TRIPP.getTensor("000|000"));
        two_proton_source_local(resultmats["PN_TRIP0_loc"],block0,block1,PN_TRIP0.getTensor("000|000"));
        two_proton_source_local(resultmats["PN_TRIPM_loc"],block0,block1,PN_TRIPM.getTensor("000|000"));
        swatch_cont.stop();

        QDPIO::cout << "done!" << std::endl;
        QDPIO::cout << "Contract: fourier: time=" << swatch_fft.getTimeInSeconds() << std::endl;
        QDPIO::cout << "Contract: contrac: time=" << swatch_cont.getTimeInSeconds() << std::endl;
        END_CODE();

        return EXIT_SUCCESS;
    }

    // With pre-sank linear combinations:
    int contract(
                 LatticeComplex& result_P,
                 std::map<std::string,
                 LatticeHalfSpinMatrix>& resultmats,
                 const multi1d<LatticePropagator>& prop0,
                 const multi1d<LatticePropagator>& prop1,
                 const SpinMatrix& diquark_proj,
                 const LatticeComplex& phases,
                 Fourier& fft,
                 const bool& compute_locals,
                 const multi1d<Complex>& weights){

        START_CODE();
        QDPIO::cout << "Performing source-displaced two-proton contractions on the CPU..." << std::flush;

        //measure fft performance:
        StopWatch swatch_fft, swatch_cont, swatch_proton;
        swatch_fft.reset();
        swatch_cont.reset();
        swatch_proton.reset();

        //initialize variables and set all input fields to zero:
        result_P=zero;
        for(std::map<std::string,LatticeHalfSpinMatrix>::iterator it=resultmats.begin(); it!=resultmats.end(); ++it){
            it->second=zero;
        }
        LatticeHalfBaryonblock block0, block1;

        if(compute_locals){
            //get baryon block combination: block0 is boosted
            get_barblock_boost(block0,"000",prop0,prop1,diquark_proj,phases,weights);
            get_barblock(block1,"000",prop0,prop1,diquark_proj,weights);

            //fourier:
            swatch_fft.start();
            block0=fft(block0,+1);
            block1=fft(block1,-1);
            swatch_fft.stop();

            //contract PP
            swatch_cont.start();
            one_proton(result_P,block0);
            two_proton_source_local(resultmats["PP_SING0_loc"],block0,block1,PP_SING0.getTensor("000|000"));

            //contract PN
            two_proton_source_local(resultmats["PN_TRIPP_loc"],block0,block1,PN_TRIPP.getTensor("000|000"));
            two_proton_source_local(resultmats["PN_TRIP0_loc"],block0,block1,PN_TRIP0.getTensor("000|000"));
            two_proton_source_local(resultmats["PN_TRIPM_loc"],block0,block1,PN_TRIPM.getTensor("000|000"));
            swatch_cont.stop();
        }

        //store the blocks 101, 110, 011:
        std::vector<std::string> bmodes(3);
        bmodes[0]="101";
        bmodes[1]="110";
        bmodes[2]="011";
        std::map<std::string, LatticeHalfBaryonblock> bblockmap;
        for(unsigned int b=0; b<bmodes.size(); b++){
            get_barblock(block1,bmodes[b],prop0,prop1,diquark_proj,weights);
            swatch_fft.start();
            block1=fft(block1,-1);
            swatch_fft.stop();
            bblockmap[bmodes[b]]=block1;
        }

        //do 000 with 111 first:
        if(!compute_locals){
            //first block is boosted:
            get_barblock_boost(block0,"000",prop0,prop1,diquark_proj,phases,weights);
            swatch_fft.start();
            block0=fft(block0,+1);
            swatch_fft.stop();
        }
        get_barblock(block1,"111",prop0,prop1,diquark_proj,weights);
        swatch_fft.start();
        block1=fft(block1,-1);
        swatch_fft.stop();

        //single nucleon from block1 (it has no momentum shift)
        swatch_proton.start();
        one_proton(result_P,block1);
        swatch_proton.stop();

        //contractions
        swatch_cont.start();
        //PP
        //s-wave:
        two_proton_displaced(resultmats["PP_SING0"],block0,block1,PP_SING0.getTensor("000|111"),phases,fft,PP_SING0.getFourierSign("000|111"));

        //P+1
        two_proton_displaced(resultmats["PP_TRIPP"],block0,block1,PP_TRIPP.getTensor("000|111"),phases,fft,PP_TRIPP.getFourierSign("000|111"));

        //P0
        two_proton_displaced(resultmats["PP_TRIP0"],block0,block1,PP_TRIP0.getTensor("000|111"),phases,fft,PP_TRIP0.getFourierSign("000|111"));

        //P-1
        two_proton_displaced(resultmats["PP_TRIPM"],block0,block1,PP_TRIPM.getTensor("000|111"),phases,fft,PP_TRIPM.getFourierSign("000|111"));


        //PN
        //s-wave:
        two_proton_displaced(resultmats["PN_SING0"],block0,block1,PN_SING0.getTensor("111|000"),phases,fft,PN_SING0.getFourierSign("111|000"));

        //P+1
        two_proton_displaced(resultmats["PN_TRIPP"],block0,block1,PN_TRIPP.getTensor("111|000"),phases,fft,PN_TRIPP.getFourierSign("111|000"));

        //P0
        two_proton_displaced(resultmats["PN_TRIP0"],block0,block1,PN_TRIP0.getTensor("111|000"),phases,fft,PN_TRIP0.getFourierSign("111|000"));

        //P-1
        two_proton_displaced(resultmats["PN_TRIPM"],block0,block1,PN_TRIPM.getTensor("111|000"),phases,fft,PN_TRIPM.getFourierSign("111|000"));
        swatch_cont.stop();


        //do the rest with all the others
        std::vector<std::string> bmodes2(3);
        bmodes2[0]="001";
        bmodes2[1]="010";
        bmodes2[2]="100";

        for(unsigned int c=0; c<bmodes2.size(); c++){
            //block0 is boosted:
            get_barblock_boost(block0,bmodes2[c],prop0,prop1,diquark_proj,phases,weights);
            swatch_fft.start();
            block0=fft(block0,+1);
            swatch_fft.stop();

            swatch_cont.start();
            std::string mode, swapmode;
            for(unsigned int b=0; b<bmodes.size(); b++){
                //PP
                //everything is "normalized" to PP, so the string and orders of blocks can stay as it is:
                mode=bmodes2[c]+"|"+bmodes[b];
                swapmode=bmodes[b]+"|"+bmodes2[c];
                //s-wave:
                two_proton_displaced(resultmats["PP_SING0"],block0,bblockmap[bmodes[b]],PP_SING0.getTensor(mode),phases,fft,PP_SING0.getFourierSign(mode));

                //P+1
                two_proton_displaced(resultmats["PP_TRIPP"],block0,bblockmap[bmodes[b]],PP_TRIPP.getTensor(mode),phases,fft,PP_TRIPP.getFourierSign(mode));

                //P0
                two_proton_displaced(resultmats["PP_TRIP0"],block0,bblockmap[bmodes[b]],PP_TRIP0.getTensor(mode),phases,fft,PP_TRIP0.getFourierSign(mode));

                //P-1
                two_proton_displaced(resultmats["PP_TRIPM"],block0,bblockmap[bmodes[b]],PP_TRIPM.getTensor(mode),phases,fft,PP_TRIPM.getFourierSign(mode));


                //PN
                //here we might need to swap the blocks:
                //s-wave:
                if(PN_SING0.getFourierSign(mode)==0){
                    two_proton_displaced(resultmats["PN_SING0"],bblockmap[bmodes[b]],block0,PN_SING0.getTensor(swapmode),phases,fft,PN_SING0.getFourierSign(swapmode));
                }
                else{
                    two_proton_displaced(resultmats["PN_SING0"],block0,bblockmap[bmodes[b]],PN_SING0.getTensor(mode),phases,fft,PN_SING0.getFourierSign(mode));
                }

                //P+1
                if(PN_TRIPP.getFourierSign(mode)==0){
                    two_proton_displaced(resultmats["PN_TRIPP"],bblockmap[bmodes[b]],block0,PN_TRIPP.getTensor(swapmode),phases,fft,PN_TRIPP.getFourierSign(swapmode));
                }
                else{
                    two_proton_displaced(resultmats["PN_TRIPP"],block0,bblockmap[bmodes[b]],PN_TRIPP.getTensor(mode),phases,fft,PN_TRIPP.getFourierSign(mode));
                }

                //P0
                if(PN_TRIP0.getFourierSign(mode)==0){
                    two_proton_displaced(resultmats["PN_TRIP0"],bblockmap[bmodes[b]],block0,PN_TRIP0.getTensor(swapmode),phases,fft,PN_TRIP0.getFourierSign(swapmode));
                }
                else{
                    two_proton_displaced(resultmats["PN_TRIP0"],block0,bblockmap[bmodes[b]],PN_TRIP0.getTensor(mode),phases,fft,PN_TRIP0.getFourierSign(mode));
                }

                //P-1
                if(PN_TRIPM.getFourierSign(mode)==0){
                    two_proton_displaced(resultmats["PN_TRIPM"],bblockmap[bmodes[b]],block0,PN_TRIPM.getTensor(swapmode),phases,fft,PN_TRIPM.getFourierSign(swapmode));
                }
                else{
                    two_proton_displaced(resultmats["PN_TRIPM"],block0,bblockmap[bmodes[b]],PN_TRIPM.getTensor(mode),phases,fft,PN_TRIPM.getFourierSign(mode));
                }
            }
            swatch_cont.stop();
        }


        //do symmetrization or anti-symmetrization of result vectors:
        swatch_fft.start();
        //PP
        symmetrize(resultmats["PP_SING0"],phases,fft,PP_SING0.getSymmetry());
        symmetrize(resultmats["PP_TRIPP"],phases,fft,PP_TRIPP.getSymmetry());
        symmetrize(resultmats["PP_TRIP0"],phases,fft,PP_TRIP0.getSymmetry());
        symmetrize(resultmats["PP_TRIPM"],phases,fft,PP_TRIPM.getSymmetry());

        //PN
        symmetrize(resultmats["PN_SING0"],phases,fft,PN_SING0.getSymmetry());
        symmetrize(resultmats["PN_TRIPP"],phases,fft,PN_TRIPP.getSymmetry());
        symmetrize(resultmats["PN_TRIP0"],phases,fft,PN_TRIP0.getSymmetry());
        symmetrize(resultmats["PN_TRIPM"],phases,fft,PN_TRIPM.getSymmetry());
        swatch_fft.stop();

        QDPIO::cout << "done!" << std::endl;
        QDPIO::cout << "Contract: fourier: time=" << swatch_fft.getTimeInSeconds() << std::endl;
        QDPIO::cout << "Contract: proton: time="  << swatch_proton.getTimeInSeconds() << std::endl;
        QDPIO::cout << "Contract: contrac: time=" << swatch_cont.getTimeInSeconds() << std::endl;
        END_CODE();

        return EXIT_SUCCESS;
    }

    // With linear combo sink
    int contract_local(LatticeComplex& result_P, std::map<std::string,LatticeHalfSpinMatrix>& resultmats, const multi1d<LatticePropagator>& prop0, const multi1d<LatticePropagator>& prop1, const SpinMatrix& diquark_proj, const LatticeComplex& phases, Fourier& fft, const multi1d<Complex>& weights){

        START_CODE();
        QDPIO::cout << "Performing source-displaced two-proton contractions on the CPU..." << std::flush;

        //measure fft performance:
        StopWatch swatch_fft, swatch_cont;
        swatch_fft.reset();
        swatch_cont.reset();

        //initialize variables and set all input fields to zero:
        result_P=zero;
        for(std::map<std::string,LatticeHalfSpinMatrix>::iterator it=resultmats.begin(); it!=resultmats.end(); ++it){
            it->second=zero;
        }
        LatticeHalfBaryonblock block0, block1;

        //get baryon block combination: block0 is boosted
        get_barblock_boost(block0,"000",prop0,prop1,diquark_proj,phases,weights);
        get_barblock(block1,"000",prop0,prop1,diquark_proj,weights);

        //fourier:
        swatch_fft.start();
        block0=fft(block0,+1);
        block1=fft(block1,-1);
        swatch_fft.stop();

        //contract PP
        swatch_cont.start();
        one_proton(result_P,block0);
        two_proton_source_local(resultmats["PP_SING0_loc"],block0,block1,PP_SING0.getTensor("000|000"));

        //contract PN
        two_proton_source_local(resultmats["PN_TRIPP_loc"],block0,block1,PN_TRIPP.getTensor("000|000"));
        two_proton_source_local(resultmats["PN_TRIP0_loc"],block0,block1,PN_TRIP0.getTensor("000|000"));
        two_proton_source_local(resultmats["PN_TRIPM_loc"],block0,block1,PN_TRIPM.getTensor("000|000"));
        swatch_cont.stop();

        QDPIO::cout << "done!" << std::endl;
        QDPIO::cout << "Contract: fourier: time=" << swatch_fft.getTimeInSeconds() << std::endl;
        QDPIO::cout << "Contract: contrac: time=" << swatch_cont.getTimeInSeconds() << std::endl;
        END_CODE();

        return EXIT_SUCCESS;
    }

}
