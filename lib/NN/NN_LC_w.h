//Arjun Singh Gambhir
//Code ported from latscat, which originally seems to be from Balint.
// -*- C++ -*-
// $Id: examples.h,v 1.8 2007-02-24 01:00:29 bjoo Exp $
//
// Include file for test suite
#ifndef _NN_W
#define _NN_W
//sparsearr2 needs to be recognized
#include "blockstuff.h"
//GetWithDef is in mapstuff
#include "mapstuff.h"
#include "../measurements/nucleon_block.h"

namespace Chroma {

    class Topology{
    private:
        std::map<std::string, sparsearr2<Real> > tensors;
        std::map<std::string, int> fourierSigns;
        int globalFourierSign;
                
    public:
        //constructors and destructors
        Topology() : globalFourierSign(static_cast<int>(+1)) {}
        ~Topology(){
            globalFourierSign=+1;
            tensors.clear();
            fourierSigns.clear();
        }
                
        //member functions
        void setSymmetry(const int& globalsign){
            globalFourierSign=globalsign;
        }
                
        void addDiagram(const std::string& mode, const sparsearr2<Real>& tensor, const int& fouriersign){
            tensors[mode]=tensor;
            fourierSigns[mode]=fouriersign;
        }
                
        int getSymmetry()const{
            return globalFourierSign;
        }
                
        int getFourierSign(const std::string& mode)const{
            return GetWithDef(fourierSigns,mode,static_cast<int>(0));
        }
                
        sparsearr2<Real> getTensor(const std::string& mode)const{
            return GetWithDef(tensors,mode,sparsearr2<Real>());
        }
                
        void clear(){
            globalFourierSign=+1;
            tensors.clear();
            fourierSigns.clear();
        }
    };

    void initTopologies(const std::string& filename, const int& truncsize, const unsigned int& j_decay);
    void clearTopologies();

    // get_barblock that returns time
    double get_barblock(LatticeHalfBaryonblock& block, 
                        const LatticePropagator& prop0, 
                        const LatticePropagator& prop1, 
                        const LatticePropagator& prop2, 
                        const SpinMatrix& diquark_proj);

    // Thorsten's get_barblocks
    void get_barblock(LatticeHalfBaryonblock& block, const std::string mode, const LatticePropagator& prop0, const LatticePropagator& prop1, const SpinMatrix& diquark_proj);
    void get_barblock(LatticeHalfBaryonblock& block, const std::string mode, const multi1d<LatticePropagator>& prop0, const multi1d<LatticePropagator>& prop1, const SpinMatrix& diquark_proj, const multi1d<Complex>& weights);
    void get_barblock(LatticeHalfBaryonblock& block, const std::string mode, const LatticePropagator& prop0, const LatticePropagator& prop1, const SpinMatrix& diquark_proj, const sink& snk);
    void get_barblock(LatticeHalfBaryonblock& block, const std::string mode, const LatticePropagator& prop0, const LatticePropagator& prop1, const SpinMatrix& diquark_proj, const multi1d<sink*>& snk, const multi1d<Complex>& weights);


    int one_proton(LatticeComplex& result, const LatticeHalfBaryonblock& protonblock);

    int two_proton_source_local(LatticeHalfSpinMatrix& result, const LatticeHalfBaryonblock& block_plus, const LatticeHalfBaryonblock& block_minus, const sparsearr2<Real>& tensor);

    int two_proton_displaced(LatticeHalfSpinMatrix& result, const LatticeHalfBaryonblock& block0, const LatticeHalfBaryonblock& block1, const sparsearr2<Real>& tensor, const LatticeComplex& phases, Fourier& fft, const int& psign);

    void symmetrize(LatticeHalfSpinMatrix& result, const LatticeComplex& phases, Fourier& fft, const int& sign);

    int contract(LatticeComplex& result_P, std::map<std::string,LatticeHalfSpinMatrix>& resultmats, const LatticePropagator& prop0, const LatticePropagator& prop1, const SpinMatrix& diquark_proj, const LatticeComplex& phases, Fourier& fft, const bool& compute_locals);
    int contract(LatticeComplex& result_P, std::map<std::string,LatticeHalfSpinMatrix>& resultmats, const LatticePropagator& prop0, const LatticePropagator& prop1, const SpinMatrix& diquark_proj, const LatticeComplex& phases, Fourier& fft, const bool& compute_locals, const sink& snk);
    int contract(LatticeComplex& result_P, std::map<std::string,LatticeHalfSpinMatrix>& resultmats, const LatticePropagator& prop0, const LatticePropagator& prop1, const SpinMatrix& diquark_proj, const LatticeComplex& phases, Fourier& fft, const bool& compute_locals, const multi1d<sink*>& snk, const multi1d<Complex>& weights);
    int contract(LatticeComplex& result_P, std::map<std::string,LatticeHalfSpinMatrix>& resultmats, const multi1d<LatticePropagator>& prop0, const multi1d<LatticePropagator>& prop1, const SpinMatrix& diquark_proj, const LatticeComplex& phases, Fourier& fft, const bool& compute_locals, const multi1d<Complex>& weights);

    int contract_local(LatticeComplex& result_P, std::map<std::string,LatticeHalfSpinMatrix>& resultmats, const LatticePropagator& prop0, const LatticePropagator& prop1, const SpinMatrix& diquark_proj, const LatticeComplex& phases, Fourier& fft, const multi1d<sink*>& snk, const multi1d<Complex>& weights);

    // contract blocks
    int contract(LatticeComplex& result_N,
                 std::map<std::string, LatticeHalfSpinMatrix>& result_NN,
                 const multi1d<const LalibeNucleonBlockEnv::BlockMapType*>& blockMap_list,
                 const multi1d<std::string>& prop0_Ids,
                 const multi1d<std::string>& prop1_Ids,
                 const multi1d<int>& origin,
                 const multi1d<Complex>& weights,
                 const multi1d<std::string>& disp_list,
                 const std::string& parity,
                 const LatticeComplex& phases,
                 Fourier& fft,
                 const bool& compute_locals,
                 const bool& compute_loc_origin);


}
#endif
