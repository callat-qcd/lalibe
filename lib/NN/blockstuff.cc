//Arjun Singh Gambhir
//Code ported from latscat, with many changes included as comments.
//#include "progs.h"
//Trying to avoid using progs.h due to chroma redundancies.
#include "blockstuff.h"

//this file contains spin projectors for baryons and helpful routines, some of them are reimplementations of chroma routines::
//using namespace Chroma;
//We don't want to use chroma namespace in this way.
namespace Chroma {

        //sparse array stuff:
        template<> template<> 
        void sparsearr<RealD>::convert(sparsearr<RealF>& result)const{
                result.resize(nvals);
                for(unsigned int i=0; i<nvals; i++){
                        result.ind[i]=ind[i];
                        result.val[i]=static_cast<RealF>(toFloat(val[i]));
                }
        }

        //the first two quarks are inside the diquark, the third gives the external spin. Thus for a proton, we have duu and udd for neutron.
        Baryonblock getBaryonblock(const Propagator& propp1, const Propagator& propp2, const Propagator& propp3, const SpinMatrix& diquark_proj){
                //first, multiply the qiquark projector with prop2 from the left:
                Propagator prop1(propp1), prop3(propp3);
                Propagator Gprop2(diquark_proj*propp2);

                //this will be the result:
                Baryonblock result;

                //now, fill the results:
                for(unsigned int eps=0; eps<Ns; eps++){
                        for(unsigned int phi=0; phi<Ns; phi++){
                                ColorMatrix tmp3(peekSpin(prop3,eps,phi));
                                for(unsigned int f=0; f<Nc; f++){

                                        //add this stuff up
                                        Propagator res=zero, tmpres;
                                        Complex tmp3comp;

                                        //eps012=+1 components:
                                        tmp3comp=peekColor(tmp3,2,f);
                                        tmpres=zero;
                                        for(unsigned int b=0; b<Nc; b++){
                                                SpinMatrix tmp1(transpose(SpinMatrix(peekColor(prop1,0,b))));
                                                for(unsigned int d=0; d<Nc; d++){
                                                        SpinMatrix tmp2(peekColor(Gprop2,1,d));
                                                        tmp2=tmp1*tmp2*tmp3comp;
                                                        pokeColor(tmpres,tmp2,b,d); //tmpres hat nun indices transpose(tmpres^{b,d}_{beta,delta})
                                                }
                                        }
                                        res+=tmpres;

                                        //eps102=-1 components:
                                        tmpres=zero;
                                        for(unsigned int b=0; b<Nc; b++){
                                                SpinMatrix tmp1(transpose(SpinMatrix(peekColor(prop1,1,b))));
                                                for(unsigned int d=0; d<Nc; d++){
                                                        SpinMatrix tmp2(peekColor(Gprop2,0,d));
                                                        tmp2=tmp1*tmp2*tmp3comp;
                                                        pokeColor(tmpres,tmp2,b,d);
                                                }
                                        }
                                        res-=tmpres;

                                        //eps201=+1 components:
                                        tmp3comp=peekColor(tmp3,1,f);
                                        tmpres=zero;
                                        for(unsigned int b=0; b<Nc; b++){
                                                SpinMatrix tmp1(transpose(SpinMatrix(peekColor(prop1,2,b))));
                                                for(unsigned int d=0; d<Nc; d++){
                                                        SpinMatrix tmp2(peekColor(Gprop2,0,d));
                                                        tmp2=tmp1*tmp2*tmp3comp;
                                                        pokeColor(tmpres,tmp2,b,d); //tmpres hat nun indices transpose(tmpres^{b,d}_{beta,delta})
                                                }
                                        }
                                        res+=tmpres;

                                        //eps021=-1 components:
                                        tmpres=zero;
                                        for(unsigned int b=0; b<Nc; b++){
                                                SpinMatrix tmp1(transpose(SpinMatrix(peekColor(prop1,0,b))));
                                                for(unsigned int d=0; d<Nc; d++){
                                                        SpinMatrix tmp2(peekColor(Gprop2,2,d));
                                                        tmp2=tmp1*tmp2*tmp3comp;
                                                        pokeColor(tmpres,tmp2,b,d); //tmpres hat nun indices transpose(tmpres^{b,d}_{beta,delta})
                                                }
                                        }
                                        res-=tmpres;

                                        //eps120=+1 components:
                                        tmp3comp=peekColor(tmp3,0,f);
                                        tmpres=zero;
                                        for(unsigned int b=0; b<Nc; b++){
                                                SpinMatrix tmp1(transpose(SpinMatrix(peekColor(prop1,1,b))));
                                                for(unsigned int d=0; d<Nc; d++){
                                                        SpinMatrix tmp2(peekColor(Gprop2,2,d));
                                                        tmp2=tmp1*tmp2*tmp3comp;
                                                        pokeColor(tmpres,tmp2,b,d); //tmpres hat nun indices transpose(tmpres^{b,d}_{beta,delta})
                                                }
                                        }
                                        res+=tmpres;

                                        //eps210=-1 components:
                                        tmpres=zero;
                                        for(unsigned int b=0; b<Nc; b++){
                                                SpinMatrix tmp1(transpose(SpinMatrix(peekColor(prop1,2,b))));
                                                for(unsigned int d=0; d<Nc; d++){
                                                        SpinMatrix tmp2(peekColor(Gprop2,1,d));
                                                        tmp2=tmp1*tmp2*tmp3comp;
                                                        pokeColor(tmpres,tmp2,b,d); //tmpres hat nun indices transpose(tmpres^{b,d}_{beta,delta})
                                                }
                                        }
                                        res-=tmpres;

                                        //res is ready, inject the stuff:
                                        result.elem().elem(eps,phi).elem(f)=Propagator(transpose(res)).elem(); //undo the transposition
                                }
                        }
                }

                return result;
        }

        //the first two quarks are inside the diquark, the third gives the external spin. Thus for a proton, we have duu and udd for neutron.
        HalfBaryonblock getHalfBaryonblock(const Propagator& propp1, const Propagator& propp2, const Propagator& propp3, const SpinMatrix& diquark_proj, bool upper){
                //first, multiply the qiquark projector with prop2 from the left:
                Propagator prop1(propp1), prop3(propp3);
                Propagator Gprop2(diquark_proj*propp2);
                HalfPropagator tmpprop;

                //this will be the result:
                HalfBaryonblock result;

                //now, fill the results:
                for(unsigned int eps=0; eps<(Ns>>1); eps++){
                        for(unsigned int phi=0; phi<(Ns>>1); phi++){
                                ColorMatrix tmp3(peekSpin(prop3,eps,phi));
                                for(unsigned int f=0; f<Nc; f++){

                                        //add this stuff up
                                        Propagator res=zero, tmpres;
                                        Complex tmp3comp;

                                        //eps012=+1 components:
                                        tmp3comp=peekColor(tmp3,2,f);
                                        tmpres=zero;
                                        for(unsigned int b=0; b<Nc; b++){
                                                SpinMatrix tmp1(transpose(SpinMatrix(peekColor(prop1,0,b))));
                                                for(unsigned int d=0; d<Nc; d++){
                                                        SpinMatrix tmp2(peekColor(Gprop2,1,d));
                                                        tmp2=tmp1*tmp2*tmp3comp;
                                                        pokeColor(tmpres,tmp2,b,d); //tmpres hat nun indices transpose(tmpres^{b,d}_{beta,delta})
                                                }
                                        }
                                        res+=tmpres;

                                        //eps102=-1 components:
                                        tmpres=zero;
                                        for(unsigned int b=0; b<Nc; b++){
                                                SpinMatrix tmp1(transpose(SpinMatrix(peekColor(prop1,1,b))));
                                                for(unsigned int d=0; d<Nc; d++){
                                                        SpinMatrix tmp2(peekColor(Gprop2,0,d));
                                                        tmp2=tmp1*tmp2*tmp3comp;
                                                        pokeColor(tmpres,tmp2,b,d);
                                                }
                                        }
                                        res-=tmpres;

                                        //eps201=+1 components:
                                        tmp3comp=peekColor(tmp3,1,f);
                                        tmpres=zero;
                                        for(unsigned int b=0; b<Nc; b++){
                                                SpinMatrix tmp1(transpose(SpinMatrix(peekColor(prop1,2,b))));
                                                for(unsigned int d=0; d<Nc; d++){
                                                        SpinMatrix tmp2(peekColor(Gprop2,0,d));
                                                        tmp2=tmp1*tmp2*tmp3comp;
                                                        pokeColor(tmpres,tmp2,b,d); //tmpres hat nun indices transpose(tmpres^{b,d}_{beta,delta})
                                                }
                                        }
                                        res+=tmpres;

                                        //eps021=-1 components:
                                        tmpres=zero;
                                        for(unsigned int b=0; b<Nc; b++){
                                                SpinMatrix tmp1(transpose(SpinMatrix(peekColor(prop1,0,b))));
                                                for(unsigned int d=0; d<Nc; d++){
                                                        SpinMatrix tmp2(peekColor(Gprop2,2,d));
                                                        tmp2=tmp1*tmp2*tmp3comp;
                                                        pokeColor(tmpres,tmp2,b,d); //tmpres hat nun indices transpose(tmpres^{b,d}_{beta,delta})
                                                }
                                        }
                                        res-=tmpres;

                                        //eps120=+1 components:
                                        tmp3comp=peekColor(tmp3,0,f);
                                        tmpres=zero;
                                        for(unsigned int b=0; b<Nc; b++){
                                                SpinMatrix tmp1(transpose(SpinMatrix(peekColor(prop1,1,b))));
                                                for(unsigned int d=0; d<Nc; d++){
                                                        SpinMatrix tmp2(peekColor(Gprop2,2,d));
                                                        tmp2=tmp1*tmp2*tmp3comp;
                                                        pokeColor(tmpres,tmp2,b,d); //tmpres hat nun indices transpose(tmpres^{b,d}_{beta,delta})
                                                }
                                        }
                                        res+=tmpres;

                                        //eps210=-1 components:
                                        tmpres=zero;
                                        for(unsigned int b=0; b<Nc; b++){
                                                SpinMatrix tmp1(transpose(SpinMatrix(peekColor(prop1,2,b))));
                                                for(unsigned int d=0; d<Nc; d++){
                                                        SpinMatrix tmp2(peekColor(Gprop2,1,d));
                                                        tmp2=tmp1*tmp2*tmp3comp;
                                                        pokeColor(tmpres,tmp2,b,d); //tmpres hat nun indices transpose(tmpres^{b,d}_{beta,delta})
                                                }
                                        }
                                        res-=tmpres;

                                        //res is ready, truncate the propagator:
                                        tmpprop=zero;
                                        if(upper){
                                                for(unsigned int s1=0; s1<(Ns>>1); s1++){
                                                        for(unsigned int s2=0; s2<(Ns>>1); s2++){
                                                                ColorMatrix colmat(peekSpin(res,s1,s2));
                                                                pokeSpin(tmpprop,colmat,s1,s2);
                                                        }
                                                }	
                                        }
                                        else{
                                                for(unsigned int s1=2; s1<Ns; s1++){
                                                        for(unsigned int s2=2; s2<Ns; s2++){
                                                                ColorMatrix colmat(peekSpin(res,s1,s2));
                                                                pokeSpin(tmpprop,colmat,(s1-2),(s2-2));
                                                        }
                                                }
                                        }
                                        result.elem().elem(eps,phi).elem(f)=HalfPropagator(transpose(tmpprop)).elem(); //undo the transposition
                                }
                        }
                }

                return result;
        }

        //peek routines for blocks
        Complex peekElement(Baryonblock block, const unsigned int& sext, const unsigned int& s1, const unsigned int& c1, const unsigned int& s2, const unsigned int& c2, const unsigned int& s3, const unsigned int& c3){
                Propagator tmp;
                tmp.elem()=block.elem().elem(sext,s3).elem(c3);
                return peekSpin(SpinMatrix(peekColor(tmp,c2,c1)),s2,s1);
        }

        Complex peekElement(HalfBaryonblock block, const unsigned int& sext, const unsigned int& s1, const unsigned int& c1, const unsigned int& s2, const unsigned int& c2, const unsigned int& s3, const unsigned int& c3){
                HalfPropagator tmp;
                tmp.elem()=block.elem().elem(sext,s3).elem(c3);
                return peekSpin(HalfSpinMatrix(peekColor(tmp,c2,c1)),s2,s1);
        }

        //full block
        LatticeBaryonblock getBaryonblock(const LatticePropagator& propp1, const LatticePropagator& propp2, const LatticePropagator& propp3, SpinMatrix diquark_proj){
                
                LatticeBaryonblock result;
                unsigned int i, nsites=Layout::sitesOnNode();

#pragma omp parallel for firstprivate(nsites,diquark_proj) private(i) shared(propp1,propp2,propp3,result)
                for(unsigned int i=0; i<nsites; i++){
                        Propagator prop1, prop2, prop3;
                        prop1.elem()=propp1.elem(i);
                        prop2.elem()=propp2.elem(i);
                        prop3.elem()=propp3.elem(i);

                        Baryonblock tmp=getBaryonblock(prop1,prop2,prop3,diquark_proj);
                        result.elem(i)=tmp.elem();
                }
                return result;
        }

        //half block
        LatticeHalfBaryonblock getHalfBaryonblock(const LatticePropagator& propp1, const LatticePropagator& propp2, const LatticePropagator& propp3, SpinMatrix diquark_proj, bool upper){
                
                LatticeHalfBaryonblock result;
                unsigned int i, nsites=Layout::sitesOnNode();

#pragma omp parallel for firstprivate(nsites,diquark_proj) private(i) shared(propp1,propp2,propp3,result)
                for(unsigned int i=0; i<nsites; i++){
                        Propagator prop1, prop2, prop3;
                        prop1.elem()=propp1.elem(i);
                        prop2.elem()=propp2.elem(i);
                        prop3.elem()=propp3.elem(i);

                        HalfBaryonblock tmp=getHalfBaryonblock(prop1,prop2,prop3,diquark_proj,upper);
                        result.elem(i)=tmp.elem();
                }
                return result;
        }

        //with arbitrary sink:
        LatticeBaryonblock getBaryonblock(const LatticePropagator& propp1, const LatticePropagator& propp2, const LatticePropagator& propp3, SpinMatrix diquark_proj, const sink& snk){
            LatticePropagator p1, p2, p3;
            p1=snk(propp1);
            p2=snk(propp2);
            p3=snk(propp3);
            
            return getBaryonblock(p1, p2, p3, diquark_proj);
        }

        LatticeHalfBaryonblock getHalfBaryonblock(const LatticePropagator& propp1, const LatticePropagator& propp2, const LatticePropagator& propp3, SpinMatrix diquark_proj, const sink& snk, bool upper){
            LatticePropagator p1, p2, p3;
            p1=snk(propp1);
            p2=snk(propp2);
            p3=snk(propp3);
            
            return getHalfBaryonblock(p1, p2, p3, diquark_proj, upper);
        }

        // Linear combinations
        LatticeBaryonblock getBaryonblock(const LatticePropagator& propp1, const LatticePropagator& propp2, const LatticePropagator& propp3, SpinMatrix diquark_proj, const multi1d<sink*>& snks, const multi1d<Complex>& weights){
            if( snks.size() != weights.size()){
                QDPIO::cout << "The number of sinks and weights must be the same when constructing a linear combination baryon block." << std::endl;
                QDP_error_exit("linear combination baryon block construction failed.");
            }
            
            LatticeBaryonblock result=zero;
            for(int s = 0; s<snks.size(); s++){
                result+=getBaryonblock(LatticePropagator(weights[s]*propp1), propp2, propp3, diquark_proj, *(snks[s]));
            }
            
            return result;
        }

        LatticeHalfBaryonblock getHalfBaryonblock(const LatticePropagator& propp1, const LatticePropagator& propp2, const LatticePropagator& propp3, SpinMatrix diquark_proj, const multi1d<sink*>& snks, const multi1d<Complex>& weights, bool upper){
            if( snks.size() != weights.size()){
                QDPIO::cout << "The number of sinks and weights must be the same when constructing a linear combination baryon block." << std::endl;
                QDP_error_exit("linear combination baryon block construction failed.");
            }
            
            LatticeHalfBaryonblock result=zero;
            for(int s = 0; s<snks.size(); s++){
                result+=getHalfBaryonblock(LatticePropagator(weights[s]*propp1), propp2, propp3, diquark_proj, *(snks[s]));
            }

            return result;
        }

        // Linear combinations of already-smeared propagators
        LatticeBaryonblock      getBaryonblock(     const multi1d<LatticePropagator>& propp1, const multi1d<LatticePropagator>& propp2, const multi1d<LatticePropagator>& propp3, SpinMatrix diquark_proj, const multi1d<Complex>& weights){
            if( ! (propp1.size() == propp2.size() && propp1.size() == propp3.size() && propp1.size()==weights.size()) ){
                QDPIO::cout << "The number of sinks and weights must be the same when constructing a linear combination baryon block." << std::endl;
                QDP_error_exit("linear combination baryon block construction from already-smeared propagators failed.");
            }
            
            LatticeBaryonblock result=zero;
            for( int s=0; s<weights.size(); s++){
                result+=getBaryonblock(LatticePropagator(weights[s]*propp1[s]), propp2[s], propp3[s], diquark_proj);
            }
            return result;
        }

        LatticeHalfBaryonblock  getHalfBaryonblock( const multi1d<LatticePropagator>& propp1, const multi1d<LatticePropagator>& propp2, const multi1d<LatticePropagator>& propp3, SpinMatrix diquark_proj, const multi1d<Complex>& weights, bool upper){
            if( ! (propp1.size() == propp2.size() && propp1.size() == propp3.size() && propp1.size()==weights.size()) ){
                QDPIO::cout << "The number of sinks and weights must be the same when constructing a linear combination baryon block." << std::endl;
                QDP_error_exit("linear combination half-baryon block construction from already-smeared propagators failed.");
            }
            
            LatticeHalfBaryonblock result=zero;
            for( int s=0; s<weights.size(); s++){
                result+=getHalfBaryonblock(LatticePropagator(weights[s]*propp1[s]), propp2[s], propp3[s], diquark_proj);
            }
            return result;
        }


        //contraction tensor handling:
        sparsearr<Real> readTensor(const std::string& filename, const std::string& path){
                sparsearr<Real> result;

                HDF5Reader input(filename);
                input.read(path+"/indices",result.ind);
                input.read(path+"/values",result.val);
                if(result.ind.size()!=result.val.size()){
                        std::cerr << "Error, inconsistent array sizes in " + filename + "/" + path << std::endl;
                        QDP_error_exit("abort!");
                }
                QDPIO::cout << "Tensor " << filename << "/" << path << " read with " << result.ind.size() << " indices!" << std::endl;
                input.close();
                result.update();

                return result;
        }

        //Linalg stuff
        //single site objects:
        //unary
        HalfBaryonblock& operator*=(HalfBaryonblock& lhs, const Real& scalar){
            Complex* bblockl=(Complex*)(&(lhs.elem()));
            //#pragma omp parallel simd safelen(halfbarblocksize)
            for(unsigned int s=0; s<halfbarblocksize; s++) bblockl[s]*=scalar;
            return lhs;
        }

        HalfBaryonblock& operator*=(HalfBaryonblock& lhs, const Complex& scalar){
            Complex* bblockl=(Complex*)(&(lhs.elem()));
            //#pragma omp parallel simd safelen(halfbarblocksize)
            for(unsigned int s=0; s<halfbarblocksize; s++) bblockl[s]*=scalar;
            return lhs;
        }

        HalfBaryonblock& operator+=(HalfBaryonblock& lhs, const HalfBaryonblock& rhs){
            Complex* bblockl=(Complex*)(&(lhs.elem()));
            Complex* bblockr=(Complex*)(&(rhs.elem()));
        //#pragma omp parallel simd safelen(halfbarblocksize)
            for(unsigned int s=0; s<halfbarblocksize; s++) bblockl[s]+=bblockr[s];
            return lhs;
        }

        //lattice fields:
        //unary
        LatticeHalfBaryonblock& operator*=(LatticeHalfBaryonblock& lhs, const Real& scalar){
            unsigned int nsites=Layout::sitesOnNode();
#pragma omp parallel for firstprivate(nsites) shared(lhs,scalar)
            for(unsigned int i=0; i<nsites; i++){
                        Complex* bblockl=(Complex*)(&(lhs.elem(i)));
                for(unsigned int s=0; s<halfbarblocksize; s++) bblockl[s]*=scalar;
            }
            return lhs;
        }

        LatticeHalfBaryonblock& operator*=(LatticeHalfBaryonblock& lhs, const Complex& scalar){
            unsigned int nsites=Layout::sitesOnNode();
#pragma omp parallel for firstprivate(nsites) shared(lhs,scalar)
            for(unsigned int i=0; i<nsites; i++){
                        Complex* bblockl=(Complex*)(&(lhs.elem(i)));
                for(unsigned int s=0; s<halfbarblocksize; s++) bblockl[s]*=scalar;
            }
            return lhs;
        }

        LatticeHalfBaryonblock& operator+=(LatticeHalfBaryonblock& lhs, const LatticeHalfBaryonblock& rhs){
            unsigned int nsites=Layout::sitesOnNode();
#pragma omp parallel for firstprivate(nsites) shared(lhs,rhs)
            for(unsigned int i=0; i<nsites; i++){
                Complex* bblockl=(Complex*)(&(lhs.elem(i)));
                Complex* bblockr=(Complex*)(&(rhs.elem(i)));
                //#pragma omp parallel simd safelen(halfbarblocksize)
                for(unsigned int s=0; s<halfbarblocksize; s++) bblockl[s]+=bblockr[s];
            }
            return lhs;
        }
}
