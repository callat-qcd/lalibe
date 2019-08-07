//Arjun Singh Gambhir
//Code ported from latscat, with many changes included as comments.
#ifndef _BLOCKSTUFF
#define _BLOCKSTUFF

#include "chromabase.h"
//We include chroma types this way.
#include "sinks.h"

namespace Chroma {
        const unsigned int barblocksize=6912;
        const unsigned int halfbarblocksize=432;

        //sparse array storage class:
        template<class T>
        class sparsearr2;

        template<class T>
        class sparsearr{
        private:
                unsigned int nvals;
          
        public:
                multi1d<unsigned int> ind;
                multi1d<unsigned int> filter;
                multi1d<T> val;
          
                sparsearr(unsigned int nnvals) : nvals(nnvals), ind(nnvals), val(nnvals){
                        for(unsigned int i=0; i<nvals; i++){
                                ind[i]=0;
                                val[i]=0;
                        }
                };
                sparsearr() : nvals(0), ind(), val() {};

                unsigned int get_nvals()const{
                        return nvals;
                }
          
                void update(){
                        nvals=ind.size();
                        if(nvals!=val.size()) QDP_error_exit("sparsearr: malformed sparse vector!");
                }

                void resize(unsigned int nnvals){
                        nvals=nnvals;
                        ind.resize(nnvals);
                        val.resize(nnvals);
                        for(unsigned int i=0; i<nvals; i++){
                                ind[i]=0;
                                val[i]=0;
                        }
                }

                template<class S>
                void convert(sparsearr<S>& result)const{
                        result.resize(nvals);
                        for(unsigned int i=0; i<nvals; i++){
                                result.ind[i]=ind[i];
                                result.val[i]=static_cast<S>(val[i]);
                        }
                }
                
                //compresses the vector
                sparsearr<T> compress()const{
                  
                        //determine the number of unique-elements:
                        std::vector<unsigned int> indices;
                        for(unsigned int i=0; i<nvals; i++){
                                indices.push_back(ind[i]);
                        }
                
                        //determine unique elements
                        std::sort(indices.begin(),indices.end());
                        indices.erase(std::unique(indices.begin(),indices.end()),indices.end());
                        
                        //create new vector of outputs, but this time the indices are shifted with respect to the left-out-ones
                        sparsearr<T> result(nvals);
                        //add the arrays to the filter: these indices are needed from the blocks
                        result.filter.resize(indices.size());
                        for(unsigned int i=0; i<indices.size(); i++){
                                result.filter[i]=indices[i];
                        }
                        //update the new indices
                        for(unsigned int i=0; i<nvals; i++){
                                result.val[i]=val[i];
                                for(unsigned int j=0; j<indices.size(); j++){
                                        if(ind[i]==indices[j]){
                                                result.ind[i]=j;
                                                continue;
                                        }
                                }
                        }
                        return result;
                }

                sparsearr2<T> split_indices(const unsigned int& blocksize)const{
                        sparsearr2<T> result(nvals);
                        for(unsigned int i=0; i<nvals; i++){
                                result.val[i]=val[i];
                                result.ind1[i]=ind[i]%blocksize;
                                result.ind2[i]=(ind[i]-result.ind1[i])/blocksize;
                        }
                        return result;
                }
        };

        //special sparse array class for two-baryon stuff:
        template<class T>
        class sparsearr2{
        private:
                unsigned int nvals;

        public:
                multi1d<unsigned int> ind1, ind2;
                multi1d<unsigned int> filter1, filter2;
                multi1d<T> val;

                sparsearr2(unsigned int nnvals) : nvals(nnvals), ind1(nnvals), ind2(nnvals), val(nnvals){
                        for(unsigned int i=0; i<nvals; i++){
                                ind1[i]=0;
                                ind2[i]=0;
                                val[i]=0;
                        }
                };

                sparsearr2() : nvals(0), ind1(), ind2(), val() {};

                unsigned int get_nvals()const{
                        return nvals;
                }

                void update(){
                        nvals=val.size();
                        if( (nvals!=ind1.size()) ||  (nvals!=ind2.size())) QDP_error_exit("sparsearr2: malformed sparse vector!");
                }

                void resize(unsigned int nnvals){
                        nvals=nnvals;
                        ind1.resize(nnvals);
                        ind2.resize(nnvals);
                        val.resize(nnvals);
                        for(unsigned int i=0; i<nvals; i++){
                                ind1[i]=0;
                                ind2[i]=0;
                                val[i]=0;
                        }
                }
                
                //compresses the vector
                sparsearr2<T> compress()const{
                  
                        //determine the number of unique-elements:
                        std::vector<unsigned int> indices1, indices2;
                        for(unsigned int i=0; i<nvals; i++){
                                indices1.push_back(ind1[i]);
                                indices2.push_back(ind2[i]);
                        }
                
                        //determine unique elements
                        std::sort(indices1.begin(),indices1.end());
                        indices1.erase(std::unique(indices1.begin(),indices1.end()),indices1.end());
                        std::sort(indices2.begin(),indices2.end());
                        indices2.erase(std::unique(indices2.begin(),indices2.end()),indices2.end());
                        
                        //create new vector of outputs, but this time the indices are shifted with respect to the left-out-ones
                        sparsearr2<T> result(nvals);
                        //add the arrays to the filter: these indices are needed from the blocks
                        result.filter1.resize(indices1.size());
                        for(unsigned int i=0; i<indices1.size(); i++){
                                result.filter1[i]=indices1[i];
                        }
                        result.filter2.resize(indices2.size());
                        for(unsigned int i=0; i<indices2.size(); i++){
                                result.filter2[i]=indices2[i];
                        }
                        //update the new indices
                        for(unsigned int i=0; i<nvals; i++){
                                result.val[i]=val[i];
                                for(unsigned int j=0; j<indices1.size(); j++){
                                        if(ind1[i]==indices1[j]){
                                                result.ind1[i]=j;
                                                continue;
                                        }
                                }
                                for(unsigned int j=0; j<indices2.size(); j++){
                                        if(ind2[i]==indices2[j]){
                                                result.ind2[i]=j;
                                                continue;
                                        }
                                }
                        }
                        return result;
                }
        };

        template<> template<> void sparsearr<RealD>::convert(sparsearr<RealF>&)const;

        //sparse array storage class for CUDA
#ifdef CUDACOMP
        //structure to store cuda device memory pointers, used for sparse array storage:
        class CuSparse{
        private:
                CuReal* val;
                unsigned int* ind;
                unsigned int nvals;
                unsigned int gridsize;

                //hide that guy, so that the user does not allocate more memory than is available:
                template<class T>
                void assgn(sparsearr<T> arr){
                        if( (val!=NULL) && (ind!=NULL) ) free();

                        nvals=arr.get_nvals();

                        //allocate storage and copy data 
                        //indices:
                        CU_ALLOC_DEV(ind,nvals*sizeof(unsigned int));
                        CU_MEMCPY(ind,&arr.ind[0],nvals*sizeof(unsigned int),cudaMemcpyHostToDevice);
                        //values:
                        CU_ALLOC_DEV(val,nvals*sizeof(CuReal));
                        CU_MEMCPY(val,&arr.val[0],nvals*sizeof(CuReal),cudaMemcpyHostToDevice);
                }

        public:
                CuSparse(){
                        val=NULL;
                        ind=NULL;
                        nvals=0;
                        gridsize=0;
                };

                CuSparse(CuReal* const vval, unsigned int* const iind, const unsigned int& nnvals){
                        val=vval;
                        ind=iind;
                        nvals=nnvals;
                }
          
                ~CuSparse(){
                        free();
                }

                CuReal* get_val(){
                        return val;
                }

                unsigned int* get_ind(){
                        return ind;
                }

                unsigned int get_nvals(){
                        return nvals;
                }
          
                unsigned int get_gridsize(){
                        return gridsize;
                }

#ifdef CUDADOUBLE
                CuSparse(sparsearr<RealD> arr){
                        val=NULL;
                        ind=NULL;
                        nvals=0;

                        assgn(arr);
                }

                void assign(sparsearr<RealD> arr){
                        assgn(arr);
                }
#else
                CuSparse(sparsearr<RealF> arr){
                        val=NULL;
                        ind=NULL;
                        nvals=0;

                        assgn(arr);
                }

                void assign(sparsearr<RealF> arr){
                        assgn(arr);
                }
#endif

                //compute optimal gridsize for contractions
                void compute_gridsize(unsigned int blocksize){
                        gridsize=ceil(nvals/(2*blocksize));
                }

                void free(){
                        if(val!=NULL){
                                CU_FREE_DEV(val);
                                val=NULL;
                        }
                        if(ind!=NULL){
                                CU_FREE_DEV(ind);
                                ind=NULL;
                        }
                        nvals=0;
                }

        }; //CuSparse

        class CuSparse2{
        private:
                CuReal* val;
                unsigned int *ind1, *ind2;
                unsigned int nvals;
                unsigned int gridsize;

                //hide that guy, so that the user does not allocate more memory than is available:
                template<class T>
                void assgn(sparsearr2<T> arr){
                        free();

                        nvals=arr.get_nvals();

                        //allocate storage and copy data 
                        //indices1:
                        CU_ALLOC_DEV(ind1,nvals*sizeof(unsigned int));
                        CU_MEMCPY(ind1,&arr.ind1[0],nvals*sizeof(unsigned int),cudaMemcpyHostToDevice);
                        //indices2:
                        CU_ALLOC_DEV(ind2,nvals*sizeof(unsigned int));
                        CU_MEMCPY(ind2,&arr.ind2[0],nvals*sizeof(unsigned int),cudaMemcpyHostToDevice);
                        //values:
                        CU_ALLOC_DEV(val,nvals*sizeof(CuReal));
                        CU_MEMCPY(val,&arr.val[0],nvals*sizeof(CuReal),cudaMemcpyHostToDevice);
                }

        public:
                CuSparse2(){
                        val=NULL;
                        ind1=NULL;
                        ind2=NULL;
                        nvals=0;
                };

                CuSparse2(CuReal* const vval, unsigned int* const iind1, unsigned int* const iind2, const unsigned int& nnvals){
                        val=vval;
                        ind1=iind1;
                        ind2=iind2;
                        nvals=nnvals;
                }

                ~CuSparse2(){
                        free();
                }

                CuReal* get_val(){
                        return val;
                }

                unsigned int* get_ind1(){
                        return ind1;
                }

                unsigned int* get_ind2(){
                        return ind2;
                }

                unsigned int get_nvals(){
                        return nvals;
                }

                unsigned int get_gridsize(){
                        return gridsize;
                }

#ifdef CUDADOUBLE
                CuSparse2(sparsearr2<RealD> arr){
                        val=NULL;
                        ind1=NULL;
                        ind2=NULL;
                        nvals=0;

                        assgn(arr);
                }

                void assign(sparsearr2<RealD> arr){
                        assgn(arr);
                }
#else
                CuSparse2(sparsearr2<RealF> arr){
                        val=NULL;
                        ind1=NULL;
                        ind2=NULL;
                        nvals=0;

                        assgn(arr);
                }

                void assign(sparsearr2<RealF> arr){
                        assgn(arr);
                }
#endif
          
                void compute_gridsize(unsigned int blocksize){
                        gridsize=ceil(nvals/(2*blocksize));
                }

                void free(){
                        if(val!=NULL){
                                CU_FREE_DEV(val);
                                val=NULL;
                        }
                        if(ind1!=NULL){
                                CU_FREE_DEV(ind1);
                                ind1=NULL;
                        }
                        if(ind2!=NULL){
                                CU_FREE_DEV(ind2);
                                ind2=NULL;
                        }
                        nvals=0;
                }
        }; //CuSparse2
#endif

        //creates baryon block from three propagators: SSS_GAMMA^{eps,BDF}(x)=sum_{a,c,e} prop1^{AB}(x,x1) prop2^{CD}(x,x2) prop3^{EF}(x,x3) eps_{a,c,e}
        Baryonblock getBaryonblock(const Propagator& prop1, const Propagator& prop2, const Propagator& prop3, const SpinMatrix& diquark_proj);
        HalfBaryonblock getHalfBaryonblock(const Propagator& prop1, const Propagator& prop2, const Propagator& prop3, const SpinMatrix& diquark_proj, bool upper=true);
        Complex peekElement(Baryonblock block, const unsigned int& sext, const unsigned int& s1, const unsigned int& c1, const unsigned int& s2, const unsigned int&c2, const unsigned int& s3, const unsigned int& c3);

        //lattice blocks
        LatticeBaryonblock getBaryonblock(const LatticePropagator& propp1, const LatticePropagator& propp2, const LatticePropagator& propp3, SpinMatrix diquark_proj);
        LatticeHalfBaryonblock getHalfBaryonblock(const LatticePropagator& propp1, const LatticePropagator& propp2, const LatticePropagator& propp3, SpinMatrix diquark_proj, bool upper=true);

        // With an arbitrary sink
        LatticeBaryonblock getBaryonblock(const LatticePropagator& propp1, const LatticePropagator& propp2, const LatticePropagator& propp3, SpinMatrix diquark_proj, const sink& snk);
        LatticeHalfBaryonblock getHalfBaryonblock(const LatticePropagator& propp1, const LatticePropagator& propp2, const LatticePropagator& propp3, SpinMatrix diquark_proj, const sink& snk, bool upper=true);

        // A linear combination of baryon blocks
        LatticeBaryonblock getBaryonblock(const LatticePropagator& propp1, const LatticePropagator& propp2, const LatticePropagator& propp3, SpinMatrix diquark_proj, const multi1d<sink*>& snk, const multi1d<Complex>& weights);
        LatticeHalfBaryonblock getHalfBaryonblock(const LatticePropagator& propp1, const LatticePropagator& propp2, const LatticePropagator& propp3, SpinMatrix diquark_proj, const multi1d<sink*>& snk, const multi1d<Complex>& weights, bool upper=true);

        // A linear combination of baryon blocks, with already-smeared propagators
        LatticeBaryonblock      getBaryonblock(     const multi1d<LatticePropagator>& propp1, const multi1d<LatticePropagator>& propp2, const multi1d<LatticePropagator>& propp3, SpinMatrix diquark_proj, const multi1d<Complex>& weights);
        LatticeHalfBaryonblock  getHalfBaryonblock( const multi1d<LatticePropagator>& propp1, const multi1d<LatticePropagator>& propp2, const multi1d<LatticePropagator>& propp3, SpinMatrix diquark_proj, const multi1d<Complex>& weights, bool upper=true);


        //contracton tensor generation
        sparsearr<Real> readTensor(const std::string& filename, const std::string& path);


        //some linalg for half-baryonblocks
        //single blocks
        //unary
        //due to the PETE-stuff, we have to be careful that we do not overload QDP++ operators. We could put them in another NS,
        //but in this case, we just brute-forcefully re-define them:
        HalfBaryonblock& operator*=(HalfBaryonblock& lhs, const Real& scalar);
        HalfBaryonblock& operator*=(HalfBaryonblock& lhs, const Complex& scalar);
        HalfBaryonblock& operator+=(HalfBaryonblock& lhs, const HalfBaryonblock& rhs);
        //BLAS-like: those guys can be templated as there is no possible ambiguity with PETE
        //return a*x+b*y
        template<class T>
        HalfBaryonblock& axpby(const HalfBaryonblock& rhs1, const T& a, const HalfBaryonblock& rhs2, const T& b){
            HalfBaryonblock result;
            Complex* bblockres=(Complex*)(&(result.elem()));
            Complex* bblockr1=(Complex*)(&(rhs1.elem()));
            Complex* bblockr2=(Complex*)(&(rhs2.elem()));
        //#pragma omp parallel simd firstprivate(a,b) shared(bblockres,bblockr1,bblockr2) safelen(halfbarblocksize)
            for(unsigned int s=0; s<halfbarblocksize; s++) bblockres[s]=a*bblockr1[s]+b*bblockr2[s];
            return result;
        }


        //lattice structures
        //unary
        LatticeHalfBaryonblock& operator*=(LatticeHalfBaryonblock& lhs, const Real& scalar);
        LatticeHalfBaryonblock& operator*=(LatticeHalfBaryonblock& lhs, const Complex& scalar);
        LatticeHalfBaryonblock& operator+=(LatticeHalfBaryonblock& lhs, const LatticeHalfBaryonblock& rhs);
        //BLAS-like: those guys can be templated as there is no possible ambiguity with PETE
        //return a*x+b*y
        template<class T>
        LatticeHalfBaryonblock& axpby(const LatticeHalfBaryonblock& rhs1, const T& a, const LatticeHalfBaryonblock& rhs2, const T& b){
            LatticeHalfBaryonblock result;
            unsigned int nsites=Layout::sitesOnNode();
#pragma omp parallel for firstprivate(nsites) shared(rhs1,rhs2,a,b)
            for(unsigned int i=0; i<nsites; i++){
                HalfBaryonblock(result.elem(i))=axpby(HalfBaryonblock(rhs1.elem(i)),a,HalfBaryonblock(rhs2.elem(i),b));
            }
            return result;
        }

}

#endif
