//Arjun Singh Gambhir
//Code ported from latscat.
#ifndef _SPINSTUFF
#define _SPINSTUFF

#include <string>

#include "chromabase.h"
//We include chroma types this way.

//To match cc file, chroma namespace also going here.
namespace Chroma {

  class BaryOp;
  std::ostream& operator<<(std::ostream&,const BaryOp&);
  StandardOutputStream& operator<<(StandardOutputStream &os,const BaryOp &obj);

  class BaryOp{
   private:
    //class members:
    multi1d<ComplexD> coeffs;
    multi1d<SpinVectorD> proj;
    multi1d<SpinMatrixD> diquark;
    std::string representation;
    std::string name;

   public:
    //constructors:
    BaryOp(){};
    BaryOp(const std::string& namee, const multi1d<SpinMatrixD>& diquarkk, const multi1d<SpinVectorD>& projj, const multi1d<Complex>& coeffss, const std::string& rep="DP");
    BaryOp(const std::string& namee, const SpinMatrixD& diquarkk, const SpinVectorD& projj, const Complex& coeffss=cmplx(static_cast<Real>(1),static_cast<Real>(0)), const std::string& rep="DP");
    
    //member functions:
    void convert();
    void push_back(const SpinMatrixD&, const SpinVectorD& projj, const Complex& coeffss);
    SpinVectorD get_proj(const unsigned int&)const;
    SpinMatrixD get_gamma(const unsigned int&)const;
    Complex get_coeff(const unsigned int&)const;
    unsigned int get_count()const;
    std::string get_name()const;

    //destructors:
    ~BaryOp();

    //friends:
    friend std::ostream& operator<<(std::ostream&,const BaryOp&);
    friend StandardOutputStream& operator<<(StandardOutputStream &os,const BaryOp &obj);
  };

  //spin matrix functions:
  SpinMatrixD get_AS_gamma(const unsigned int& row, const unsigned int& column);
  SpinVectorD get_proj(const unsigned int& id);
  SpinMatrix getProjector(const std::string mode, const std::string& rep="DP");
  HalfSpinMatrix getHalfProjector(const std::string mode);
  HalfSpinMatrix getHalfSpinMatrix(const SpinMatrix& matrix, const std::string& mode="UPPER");
  LatticeHalfSpinMatrix getHalfSpinMatrix(const LatticeSpinMatrix& matrix, const std::string& mode="UPPER");
  SpinMatrixD kron(const SpinVectorD& lhs, const SpinVectorD& rhs);
  std::ostream& operator<<(std::ostream&,const SpinMatrixD&);
  StandardOutputStream& operator<<(StandardOutputStream&,const SpinMatrixD&);

  const std::string DESIRED_SPIN_BASIS="DIRAC_PAULI";
  void spin_basis_transform(LatticePropagator& prop, const std::string new_basis, const std::string old_basis);

  //bary op functions:
  multi1d<BaryOp> get_local_MA(const unsigned int& spinid, const std::string& rep="DP");
  multi1d<BaryOp> get_local_MA_single(const unsigned int& spinid, const std::string& rep="DP");
  multi1d<BaryOp> get_local_MA_single_ungerade(const unsigned int& spinid, const std::string& rep="DP");

}
#endif
