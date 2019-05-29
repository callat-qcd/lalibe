// -*- C++ -*-
/*! 
 *  \brief Construct the traceless Field Strength Tensor.
 *  If there are questions or bugs, blame David Brantley.
 */

#include "chromabase.h"
#include "tr_less_fields.h"
namespace Chroma
{
  /* Construct the field strength tensor using the definition in Zweig et. al. (1983)*/
  // This program reads in a gauge field, and constructs the resulting traceless
  // field strength tensor. The output is then returned as a 2d array of lattice color
  // matrices.
  //         **** COMPUTES BOTH THE UPPER AND LOWER TRIANGULAR PORTIONS ****


  multi1d<LatticeColorMatrix> trLessFieldST(const multi1d<LatticeColorMatrix>& u)
  {
    START_CODE();
    
    multi1d<LatticeColorMatrix> F(Nd*(Nd-1));
    
    LatticeColorMatrix t1;
    LatticeColorMatrix t2;
    LatticeColorMatrix t3;
    LatticeColorMatrix t4;

    LatticeComplex tr;

    Real eye = 1;
    
    Real one_eighth = 0.125;
    Real one_over_Nc = 1/Nc;
    
    int ctr = 0;

    for(int mu=0; mu<Nd-1;++mu)
    {
        for(int nu=mu+1; nu<Nd;++nu)
        {

//  U_1 = u(x,mu)*u(x+mu,nu)*u^dag(x+nu,mu)*u^dag(x,nu)
            t1 = shift(u[nu],FORWARD,mu);
            t2 = u[mu]*t1;
            t3 = shift(u[mu],FORWARD,nu);
            t4 = u[nu]*t3;
            F[ctr] = t2*adj(t4);
            
//  U_2 = u(x,nu)*u^dag(x-mu+nu,mu)*u^dag(x-mu,nu)*u(x-mu,mu)
            t1 = shift(u[mu],BACKWARD,mu);
            t2 = shift(shift(u[mu],BACKWARD,mu),FORWARD,nu);
            t3 = shift(u[nu],BACKWARD,mu);
            t4 = t3*t2;
            t2 = adj(t4);
            t3 = t2*t1;
            F[ctr] += u[nu]*t3;
//  U_3 = u^dag(x-mu,mu)*u^dag(x-mu-nu,nu)*u(x-mu-nu,mu)*u(x-nu,nu)
            t1 = shift(u[mu],BACKWARD,mu);
            t2 = shift(shift(u[nu],BACKWARD,mu),BACKWARD,nu);
            t3 = t2*t1;
            t4 = adj(t3);
            t1 = shift(shift(u[mu],BACKWARD,mu),BACKWARD,nu);
            t2 = shift(u[nu],BACKWARD,nu);
            t3 = t1*t2;
            F[ctr] += t4*t3;
//  U_4 = u^dag(x-nu,nu)*u(x-nu,mu)*u(x-nu+mu,nu)*u^dag(x,mu)
            t1 = shift(u[mu],BACKWARD,nu);
            t2 = shift(shift(u[nu],BACKWARD,nu),FORWARD,mu);
            t3 = t1*t2;
            t4 = adj(u[mu]);
            t1 = t3*t4;
            t2 = shift(u[nu],BACKWARD,nu);
            F[ctr] += adj(t2)*t1;

            t1 = adj(F[ctr]);
            F[ctr] -= t1;
            F[ctr] *= one_eighth;

// Enforce tracelessness.
            LatticeColorMatrix I = eye;

            tr = traceColor(F[ctr]);

            t2 = I*tr;
            t2 *= one_over_Nc;

            F[ctr] -= t2;

            ctr++;
        }
   }
 
    for(int mu=1; mu<Nd;++mu)
    {
        for(int nu=0; nu<mu;++nu)
        {


//  U_1 = u(x,mu)*u(x+mu,nu)*u^dag(x+nu,mu)*u^dag(x,nu)
            t1 = shift(u[nu],FORWARD,mu);
            t2 = u[mu]*t1;
            t3 = shift(u[mu],FORWARD,nu);
            t4 = u[nu]*t3;
            F[ctr] = t2*adj(t4);
            
//  U_2 = u(x,nu)*u^dag(x-mu+nu,mu)*u^dag(x-mu,nu)*u(x-mu,mu)
            t1 = shift(u[mu],BACKWARD,mu);
            t2 = shift(shift(u[mu],BACKWARD,mu),FORWARD,nu);
            t3 = shift(u[nu],BACKWARD,mu);
            t4 = t3*t2;
            t2 = adj(t4);
            t3 = t2*t1;
            F[ctr] += u[nu]*t3;
//  U_3 = u^dag(x-mu,mu)*u^dag(x-mu-nu,nu)*u(x-mu-nu,mu)*u(x-nu,nu)
            t1 = shift(u[mu],BACKWARD,mu);
            t2 = shift(shift(u[nu],BACKWARD,mu),BACKWARD,nu);
            t3 = t2*t1;
            t4 = adj(t3);
            t1 = shift(shift(u[mu],BACKWARD,mu),BACKWARD,nu);
            t2 = shift(u[nu],BACKWARD,nu);
            t3 = t1*t2;
            F[ctr] += t4*t3;
//  U_4 = u^dag(x-nu,nu)*u(x-nu,mu)*u(x-nu+mu,nu)*u^dag(x,mu)
            t1 = shift(u[mu],BACKWARD,nu);
            t2 = shift(shift(u[nu],BACKWARD,nu),FORWARD,mu);
            t3 = t1*t2;
            t4 = adj(u[mu]);
            t1 = t3*t4;
            t2 = shift(u[nu],BACKWARD,nu);
            F[ctr] += adj(t2)*t1;

            t1 = adj(F[ctr]);
            F[ctr] -= t1;
            F[ctr] *= one_eighth;

// Enforce tracelessness.
            LatticeColorMatrix I = eye;

            tr = traceColor(F[ctr]);

            t2 = I*tr;
            t2 *= one_over_Nc;

            F[ctr] -= t2;

            ctr++;
        }
   }

  return F;

  END_CODE();

  }
}

