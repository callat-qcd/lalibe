//Arjun Singh Gambhir
//Code ported from latscat.
//progs.h was included here, but we are trying to avoid that.
#include "spinstuff.h"
//In the original latscat file, spinstuff.h is not included, I do so here so if there are issues, there are caught before linking.
#include <iostream>
//This is the bare minimum stuff in progs that is needed.
#include "util/ferm/paulitodr.h"
//Additional stuff needed from chroma.

//this file contains spin projectors for baryons and helpful routines, some of them are reimplementations of chroma routines::
namespace Chroma {

  std::ostream& operator<<(std::ostream& os, const SpinMatrix& obj){
    for(unsigned int i=0; i<Ns; i++){
      for(unsigned int j=0; j<Ns; j++){
        os << peekSpin(obj,i,j) << "\t";
      }
      os << std::endl;
    }
    return os;
  }

  StandardOutputStream& operator<<(StandardOutputStream& os, const SpinMatrix& obj){
    for(unsigned int i=0; i<Ns; i++){
      for(unsigned int j=0; j<Ns; j++){
        os << peekSpin(obj,i,j) << "\t";
      }
      os << std::endl;
    }
    return os;
  }

  //returns several types of AS gamma matrices:
  SpinMatrix get_AS_gamma(const unsigned int& row, const unsigned int& column){
    SpinMatrix U = zero;
    Real     prefac = static_cast<Real>(1) / sqrt(static_cast<Real>(2));
    Complex  one = cmplx(prefac,static_cast<Real>(0));
    Complex mone = cmplx(-prefac,static_cast<Real>(0));

    pokeSpin(U,one,row,column);
    pokeSpin(U,mone,column,row);

    return U;
  }

  //get matrix which projects out a single component:
  SpinVector get_proj(const unsigned int& id){
    SpinVector V=zero;
    Complex one=cmplx(static_cast<Real>(1),static_cast<Real>(0));
    
    pokeSpin(V,one,id);

    return V;
  }

  SpinMatrix getProjector(const std::string mode, const std::string& rep){
    SpinMatrix result=zero;
    
    //create spin projectors in dirac pauli (DP) basis:
    if(mode.compare("TRIPP1")==0){
      pokeSpin(result,cmplx(static_cast<Real>(1.),static_cast<Real>(0.)),0,0);
    }
    else if(mode.compare("TRIPM1")==0){
      pokeSpin(result,cmplx(static_cast<Real>(1.),static_cast<Real>(0.)),1,1);
    }
    else if(mode.compare("TRIP0")==0){
      pokeSpin(result,cmplx(static_cast<Real>(1./sqrt(2.)),static_cast<Real>(0.)),0,1);
      pokeSpin(result,cmplx(static_cast<Real>(1./sqrt(2.)),static_cast<Real>(0.)),1,0);
    }
    else if(mode.compare("SING0")==0){
      pokeSpin(result,cmplx(static_cast<Real>(1./sqrt(2.)),static_cast<Real>(0.)),0,1);
      pokeSpin(result,cmplx(static_cast<Real>(-1./sqrt(2.)),static_cast<Real>(0.)),1,0);
    }
    
    //convert to other bases if necessary:
    if(rep.compare("DR")==0){
      SpinMatrix P2DR=PauliToDRMat();
      result=P2DR*result*adj(P2DR);
    }
    
    return result;
  }

  //only meaningful in DP
  HalfSpinMatrix getHalfProjector(const std::string mode){
    HalfSpinMatrix result=zero;
    
    //create spin projectors in dirac pauli basis:
    if(mode.compare("TRIPP1")==0){
      pokeSpin(result,cmplx(static_cast<Real>(1.),static_cast<Real>(0.)),0,0);
    }
    else if(mode.compare("TRIPM1")==0){
      pokeSpin(result,cmplx(static_cast<Real>(1.),static_cast<Real>(0.)),1,1);
    }
    else if(mode.compare("TRIP0")==0){
      pokeSpin(result,cmplx(static_cast<Real>(1./sqrt(2.)),static_cast<Real>(0.)),0,1);
      pokeSpin(result,cmplx(static_cast<Real>(1./sqrt(2.)),static_cast<Real>(0.)),1,0);
    }
    else if(mode.compare("SING0")==0){
      pokeSpin(result,cmplx(static_cast<Real>(1./sqrt(2.)),static_cast<Real>(0.)),0,1);
      pokeSpin(result,cmplx(static_cast<Real>(-1./sqrt(2.)),static_cast<Real>(0.)),1,0);
    }
      
    return result;
  }

  HalfSpinMatrix getHalfSpinMatrix(const SpinMatrix& matrix, const std::string& mode){
    HalfSpinMatrix result=zero;
    unsigned int Nshalf=(Ns>>1);
    
    if(mode=="UPPER"){
      for(unsigned int s1=0; s1<Nshalf; s1++){
        for(unsigned int s2=0; s2<Nshalf; s2++){
          Complex entry=peekSpin(matrix,s1,s2);
          pokeSpin(result,entry,s1,s2);
        }
      }
    }
    else if(mode=="LOWER"){
      for(unsigned int s1=Nshalf; s1<Ns; s1++){
        for(unsigned int s2=Nshalf; s2<Ns; s2++){
          Complex entry=peekSpin(matrix,s1,s2);
          pokeSpin(result,entry,s1-Nshalf,s2-Nshalf);
        }
      }
    }
    
    return result;
  }

  LatticeHalfSpinMatrix getHalfSpinMatrix(const LatticeSpinMatrix& matrix, const std::string& mode){
    LatticeHalfSpinMatrix result=zero;
    unsigned int Nshalf=(Ns>>1);
    
    if(mode=="UPPER"){
      for(unsigned int s1=0; s1<Nshalf; s1++){
        for(unsigned int s2=0; s2<Nshalf; s2++){
          LatticeComplex entry=peekSpin(matrix,s1,s2);
          pokeSpin(result,entry,s1,s2);
        }
      }
    }
    else if(mode=="LOWER"){
      for(unsigned int s1=Nshalf; s1<Ns; s1++){
        for(unsigned int s2=Nshalf; s2<Ns; s2++){
          LatticeComplex entry=peekSpin(matrix,s1,s2);
          pokeSpin(result,entry,s1-Nshalf,s2-Nshalf);
        }
      }
    }
    
    return result;
  }

  //returns the kronecker product between lhs and rhs spin vectors:
  SpinMatrix kron(const SpinVector& lhs, const SpinVector& rhs){
    SpinMatrix result=zero;

    for(unsigned int i=0; i<4; i++){
      for(unsigned int j=0; j<4; j++){
        pokeSpin(result,peekSpin(lhs,i)*peekSpin(rhs,j),i,j);
      }
    }
    
    return result;
  }

  //baryon operator class stuff:
  BaryOp::BaryOp(const std::string& namee, const multi1d<SpinMatrix>& diquarkk, const multi1d<SpinVector>& projj, const multi1d<Complex>& coeffss, const std::string& rep) : name(namee), coeffs(coeffss), proj(projj), diquark(diquarkk), representation(rep){};

  BaryOp::BaryOp(const std::string& namee, const SpinMatrix& diquarkk, const SpinVector& projj, const Complex& coeffss, const std::string& rep) : name(namee), representation(rep){
    coeffs.resize(1);
    proj.resize(1);
    diquark.resize(1);

    diquark[0]=diquarkk;
    proj[0]=projj;
    coeffs[0]=coeffss;
  }

  BaryOp::~BaryOp(){
    diquark.resize(0);
    proj.resize(0);
    coeffs.resize(0);
  }

  //operators:
  std::ostream& operator<<(std::ostream &os,const BaryOp &obj){
    os << "Name: " << obj.name << std::endl;
    for(unsigned int i=0; i<obj.get_count(); i++){
      os << "representation: " << obj.representation << std::endl; 
      os << "coeff[" << i << "]: " << obj.coeffs[i] << std::endl;
      os << "gamma[" << i << "]: " << std::endl;
      os << obj.diquark[i] << std::endl;
      os << "proj[" << i << "]: " << std::endl; 
      os << obj.proj[i] << std::endl << std::endl; 
    }
    return os;
  }

  StandardOutputStream& operator<<(StandardOutputStream &os,const BaryOp &obj){
    os <<"Name: " << obj.name << std::endl;
    for(unsigned int i=0; i<obj.get_count(); i++){
      os << "representation: " << obj.representation << std::endl;
      os << "coeff[" << i << "]: " << obj.coeffs[i] << std::endl;
      os << "gamma[" << i << "]: " << std::endl;
      os << obj.diquark[i] << std::endl;
      os << "proj[" << i << "]: " << std::endl;
      os << obj.proj[i] << std::endl << std::endl;
    }
    return os;
  }


  //convert between DP (Dirac-Pauli) and DR (DeGrand-Rossi)
  void BaryOp::convert(){
    SpinMatrix P2DR=PauliToDRMat();
    if(representation.compare("DR")==0) P2DR=adj(P2DR);
    
    for(unsigned int i=0; i<diquark.size(); i++){
      proj[i]=P2DR*proj[i];
      diquark[i]=P2DR*diquark[i]*adj(P2DR);
    }

    representation=(representation.compare("DR")==0 ? "DP" : "DR");
  }

  void BaryOp::push_back(const SpinMatrix& diquarkk, const SpinVector& projj, const Complex& coeffss){
    //push back into multi1d-arrays:
    multi1d<SpinMatrix> tmp(diquark.size()+1);
    for(unsigned int i=0; i<diquark.size(); i++) tmp[i]=diquark[i];
    tmp[diquark.size()]=diquarkk;
    diquark=tmp;

    multi1d<SpinVector> tmp2(proj.size()+1);
    for(unsigned int i=0; i<proj.size(); i++) tmp2[i]=proj[i];
    tmp2[proj.size()]=projj;
    proj=tmp2;

    multi1d<Complex> tmp3(coeffs.size()+1);
    for(unsigned int i=0; i<coeffs.size(); i++) tmp3[i]=coeffs[i];
    tmp3[coeffs.size()]=coeffss;
    coeffs=tmp3;
  }

  unsigned int BaryOp::get_count()const{
    return diquark.size();
  }

  SpinVector BaryOp::get_proj(const unsigned int& id)const{
    return proj[id];
  }

  SpinMatrix BaryOp::get_gamma(const unsigned int& id)const{
    return diquark[id];
  }

  Complex BaryOp::get_coeff(const unsigned int& id)const{
    return coeffs[id];
  }

  std::string BaryOp::get_name()const{
    return name;
  }

  void spin_basis_transform(LatticePropagator& prop, const std::string new_basis, const std::string old_basis){
      if( new_basis == old_basis){ 
          QDPIO::cout << "Already in the " << new_basis << " basis!" << std::endl;
          return; 
      }
      
      SpinMatrix rotator;
      
      if      ( new_basis == "DIRAC_PAULI" ){
          if  ( old_basis == "DEGRAND_ROSSI")     rotator = adj(PauliToDRMat());
          else{                                   QDPIO::cout << "Don't know how to reach the " << new_basis << " basis from the " << old_basis << "basis." << std::endl;
                                                  QDP_abort(EXIT_FAILURE); }
      }
      else if ( new_basis == "DEGRAND_ROSSI"){
          if  ( old_basis == "DIRAC_PAULI"  )     rotator = PauliToDRMat();
          else{                                   QDPIO::cout << "Don't know how to reach the " << new_basis << " basis from the " << old_basis << "basis." << std::endl;
                                                  QDP_abort(EXIT_FAILURE); }
      }
      else{                                       QDPIO::cout << "Don't know how to reach the " << new_basis << "." << std::endl;
                                                  QDP_abort(EXIT_FAILURE);
      }
          
      QDPIO::cout << "Transforming from " << old_basis << " basis to " << new_basis << " basis...";
      LatticePropagator temp;
      temp = rotator * prop * adj(rotator);
      prop = temp;
      QDPIO::cout << "done!" << std::endl;
  }


  //Class functions:
  //return baryon Operator for local Mixed-antisymmetric baryons and spin up or down (0,1):
  multi1d<BaryOp> get_local_MA(const unsigned int& spinid, const std::string& rep){
    multi1d<BaryOp> result(3);

    if(spinid>1){
      QDPIO::cerr << "get_local_MA: error, the spin should be 0 (up) or 1 (down)!" << std::endl;
      return result;
    }

    std::string spinsign="+";
    if(spinid==1) spinsign="-";

    //spin in DP basis:
    //121:
    result[0]=BaryOp("G1g(1,"+spinsign+"1/2)",get_AS_gamma(0,1),get_proj(0+spinid),cmplx(static_cast<Real>(1),static_cast<Real>(0)),"DP");

    //(143+323+341)/sqrt(3):
    Real     norm = static_cast<Real>(1) / sqrt(static_cast<Real>(3));
    Complex  cnorm = cmplx(norm,static_cast<Real>(0));
    result[1]=BaryOp("G1g(2,"+spinsign+"1/2)",get_AS_gamma(0,3),get_proj(2+spinid),cnorm,"DP");
    result[1].push_back(get_AS_gamma(2,1),get_proj(2+spinid),cnorm);
    result[1].push_back(get_AS_gamma(2,3),get_proj(0+spinid),cnorm);

    //(134+323-341)/sqrt(3):
    result[2]=BaryOp("G1g(3,"+spinsign+"1/2)",get_AS_gamma(0,2+spinid),get_proj(3),cnorm,"DP");
    result[2].push_back(get_AS_gamma(2+spinid,1),get_proj(2),cnorm);
    result[2].push_back(get_AS_gamma(2,3),get_proj(0+spinid),-cnorm);

    //convert to DR if necessary:
    if(rep.compare("DR")==0) for(unsigned int i=0; i<3; i++) result[i].convert();

    return result;
  }

  multi1d<BaryOp> get_local_MA_single(const unsigned int& spinid, const std::string& rep){
    multi1d<BaryOp> result;
    if(spinid>1){
      QDPIO::cerr << "get_local_MA_single: error, the spin should be 0 (up) or 1 (down)!" << std::endl;
      return result;
    }
    result.resize(5);
    
    Complex one=cmplx(static_cast<Real>(1),static_cast<Real>(0));
    
    //spin in DP basis:
    //121/122:
    result[0]=BaryOp("N_12"+std::to_string(1+spinid),get_AS_gamma(0,1),get_proj(0+spinid),one,"DP");
    
    //143/144:
    result[1]=BaryOp("N_14"+std::to_string(3+spinid),get_AS_gamma(0,3),get_proj(2+spinid),one,"DP");

    //323/324:
    result[2]=BaryOp("N_32"+std::to_string(3+spinid),get_AS_gamma(2,1),get_proj(2+spinid),one,"DP");

    //341/342:
    result[3]=BaryOp("N_34"+std::to_string(1+spinid),get_AS_gamma(2,3),get_proj(0+spinid),one,"DP");

    if(spinid==0){
      //134:
      result[4]=BaryOp("N_134",get_AS_gamma(0,2),get_proj(3),one,"DP");
    }
    else{
      //423:
      result[4]=BaryOp("N_423",get_AS_gamma(3,1),get_proj(2),one,"DP");
    }

    //convert to DR if necessary:
    if(rep.compare("DR")==0) for(unsigned int i=0; i<result.size(); i++) result[i].convert();

    return result;
  }

  multi1d<BaryOp> get_local_MA_single_ungerade(const unsigned int& spinid, const std::string& rep){
    multi1d<BaryOp> result;
    if(spinid>1){
      QDPIO::cerr << "get_local_MA_single_ungerade: error, the spin should be 0 (up) or 1 (down)!" << std::endl;
      return result;
    }
    result.resize(5);
    
    Complex one=cmplx(static_cast<Real>(1),static_cast<Real>(0));

    //spin in DP basis:
    //343/344:
    result[0]=BaryOp("N_34"+std::to_string(3+spinid),get_AS_gamma(2,3),get_proj(2+spinid),one,"DP");

    //123/124:
    result[1]=BaryOp("N_12"+std::to_string(3+spinid),get_AS_gamma(0,1),get_proj(2+spinid),one,"DP");

    //141/142:
    result[2]=BaryOp("N_14"+std::to_string(1+spinid),get_AS_gamma(0,3),get_proj(0+spinid),one,"DP");

    //321/322:
    result[3]=BaryOp("N_32"+std::to_string(1+spinid),get_AS_gamma(2,1),get_proj(0+spinid),one,"DP");

    if(spinid==0){
      //312:
      result[4]=BaryOp("N_312",get_AS_gamma(2,0),get_proj(1),one,"DP");
    }
    else{
      //241:
      result[4]=BaryOp("N_241",get_AS_gamma(1,3),get_proj(0),one,"DP");
    }

    //convert to DR if necessary:
    if(rep.compare("DR")==0) for(unsigned int i=0; i<result.size(); i++) result[i].convert();

    return result;
  }
}
