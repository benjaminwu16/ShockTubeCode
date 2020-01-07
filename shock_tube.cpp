//========================================================================================
// athena++ astrophysical mhd code
// copyright(c) 2014 james m. stone <jmstone@princeton.edu> and other code contrIButors
// licensed under the 3-clause bsd license, see license FILE for details
//========================================================================================
//! \FILE shock_tube.cpp
//  \brief problem generator for shock tube problems.
//
// problem generator for shock tube (1-d riemann) problems. initializes plane-parallel
// shock along x1 (in 1d, 2d, 3d), along x2 (in 2d, 3d), and along x3 (in 3d).
//========================================================================================

// c headers
#include <stdio.h>

// c++ headers
#include <cmath>      // sqrt()
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>
#include <random>

// athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../cless/cless.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../bvals/bvals.hpp"

const Real cospiover8 = 0.923879532511286;
const Real pi = 3.141592653589793238;

Real wl[NHYDRO+NFIELD];
Real wr[NHYDRO+NFIELD];
Real phase_angles[1000]; //phase angles for the wave summation to generate turbulence

// fixes bcs on l-x1 (left edge) of grid to postshock flow.
void stinner_ix1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void reflect_ox1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
int is, int ie, int js, int je, int ks, int ke, int ngh);

Real pos(MeshBlock *pmb, int iout);
Real firstjump_onedimension(MeshBlock *pmb, int iout);
Real firstjump_twodimensions(MeshBlock *pmb, int iout);
Real secjump_onedimension(MeshBlock *pmb, int iout);
Real secjump_twodimensions(MeshBlock *pmb, int iout);
Real divergenceb(MeshBlock *pmb, int iout);
Real mhd_jump(MeshBlock *pmb, int iout);
Real bperp_jc(MeshBlock *pmb, int iout);
Real bpara_jc(MeshBlock *pmb, int iout);
//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief function to initialize problem-specific data in Mesh class.  can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this FILE.  called in Mesh constructor.
//========================================================================================
//give me an error

void Mesh::InitUserMeshData(ParameterInput *pin) {
// set iib value function pointer
//	give me an error  
  EnrollUserBoundaryFunction(INNER_X1, stinner_ix1);
  EnrollUserBoundaryFunction(OUTER_X1, reflect_ox1);
  AllocateUserHistoryOutput(5);
  EnrollUserHistoryOutput(0, pos, "shockpos");
//  EnrollUserHistoryOutput(1, firstjump_twodimensions, "momflux");
//  EnrollUserHistoryOutput(2, secjump_twodimensions, "pressurediff");
  EnrollUserHistoryOutput(1, bpara_jc, "bpara_jc");
  EnrollUserHistoryOutput(2, bperp_jc, "bperp_jc"); 
  EnrollUserHistoryOutput(3, mhd_jump, "mhd_jump");
  EnrollUserHistoryOutput(4, divergenceb, "divb");
  
  return;
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(2);
  AllocateRealUserMeshBlockDataField(3);
  ruser_meshblock_data[0].NewAthenaArray(260, 260);
  ruser_meshblock_data[1].NewAthenaArray(4, 260);
  ruser_meshblock_data[2].NewAthenaArray(1000);
 

  for(int i=0; i<260; ++i) {
    for(int j=0; j<260; ++j) {
      ruser_meshblock_data[0](i, j) = 0.0;
    }
  }


  for(int i=0; i<4; ++i) {
    for(int j=0; j<260; ++j) {
      ruser_meshblock_data[1](i, j) = 0.0;
    }
  }

  return;
}


//probability power law distrIBution
Real power_law(Real l, Real a) {
  return 1/(1+pow(l*a, 2.6666667));
}
 
void genAngles() {
  std::random_device rd;
  std::mt19937 gen(rd());  
  std::uniform_real_distribution<> dis(0.0, 1.0);

  for(int i=0; i<1000; i++) {
      phase_angles[i] = dis(gen)*2*pi;
  }
}  

//returns density value for given x, y
Real return_d(Real x, Real y, int steps, Real l, Real avg) {
  std::random_device rd;
  std::mt19937 gen(rd());  
  std::uniform_real_distribution<> dis(0.0, 1.0);

  Real k = (2*pi)/(10*l);
  Real res = 0.0;
  Real ratio = pow(16.0, 1.0/(steps));
  Real norm_constant = 0.0;
  Real theta, xwave, ywave, amplitude;
  for(int i=0; i<steps; i++) {
      amplitude = 1/(1+pow(k, 2.66666667));
      norm_constant += amplitude;
      k *= ratio;
  }
  k = (2*pi)/(10*l);
  //sum over wave modes
  for(int i=0; i<steps; i++) {
     //random angle
      theta = dis(gen)*2*pi;
      xwave = cos(k*cos(theta)*x + phase_angles[i]);
      ywave = cos(k*sin(theta)*y + phase_angles[i]);
      amplitude = 1/(1+pow(k, 2.66666667));
      res += xwave*ywave*amplitude;
      k *= ratio;
  }
  return exp(res/norm_constant);  
}

//returns magnitude potential function at each (x, y)
Real b_potential(Real x, Real y, int steps, Real l) {
  std::random_device rd;
  std::mt19937 gen(rd());  
  std::uniform_real_distribution<> dis(0.0, 1.0);

  Real k = (2*pi)/(10*l);
  Real res = 0.0;
  Real ratio = pow(16.0, 1.0/(steps));
  Real theta, xwave, ywave, amplitude;
  //sum over wave modes
  for(int i=0; i<steps; i++) {
     //random angle
      theta = dis(gen)*2*pi;
      xwave = cos(k*cos(theta)*x + phase_angles[i]);
      ywave = cos(k*sin(theta)*y + phase_angles[i]);
      amplitude = 1/(1+pow(k, 2.66666667));
      res += xwave*ywave*amplitude;
      k *= ratio;
  }
  return res;
}


//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief problem generator for the shock tube tests
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  wl[IDN] = pin->GetReal("problem","dl");
  wl[IVX] = pin->GetReal("problem","ul");
  wl[IVY] = pin->GetReal("problem","vl");
  wl[IVZ] = pin->GetReal("problem","wl");
  if (NON_BAROTROPIC_EOS) wl[IPR] = pin->GetReal("problem","pl");
  if (MAGNETIC_FIELDS_ENABLED) {
    wl[NHYDRO  ] = pin->GetReal("problem","bxl");
    wl[NHYDRO+1] = pin->GetReal("problem","byl");
    wl[NHYDRO+2] = pin->GetReal("problem","bzl");
  }
  genAngles();
  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie; ++i) {
	Real xcoord = pcoord->x1v(i);
	Real ycoord = pcoord->x2v(j);
	phydro->u(IDN,k,j,i) = return_d(xcoord,ycoord,1000,1.0,wl[IDN]);
	//phydro->u(IDN,k,j,i) = wl[IDN];
        phydro->u(IM1,k,j,i) = wl[IVX]*phydro->u(IDN,k,j,i);
        phydro->u(IM2,k,j,i) = wl[IVY]*phydro->u(IDN,k,j,i);
        phydro->u(IM3,k,j,i) = wl[IVZ]*phydro->u(IDN,k,j,i);
        if (NON_BAROTROPIC_EOS) phydro->u(IEN,k,j,i) =
          wl[IPR]/(peos->GetGamma() - 1.0)
          + 0.5*phydro->u(IDN,k,j,i)*(wl[IVX]*wl[IVX] + wl[IVY]*wl[IVY] + wl[IVZ]*wl[IVZ]);
      }
    }
   }

// now set face-centered (interface) magnetic fields -----------------------------------

  if (MAGNETIC_FIELDS_ENABLED) {
      
      for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
    	  Real xcoord = pcoord->x1v(i);
	  Real ycoord = pcoord->x2v(j);
	  ruser_meshblock_data[0](i,j) = b_potential(xcoord,ycoord,1000,1.0);
        
	}
    }

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
     	  
	  int xdiff = (i < ie) ? 1 : -1;
          int ydiff = (j < je) ? 1 : -1;	
          Real xcoord = pcoord->x1v(i);
          Real ycoord = pcoord->x2v(j);
          Real xcoord2 = pcoord->x1v(i+xdiff);
      
          Real ycoord2 = pcoord->x2v(j+ydiff);
          Real orig = ruser_meshblock_data[0](i, j);
          Real ychange = ruser_meshblock_data[0](i, j+ydiff);
          Real xchange = ruser_meshblock_data[0](i+xdiff, j);
	
          //discrete approximation of partial derivative 
          pfield->b.x1f(k,j,i) = (ychange - orig) / (ycoord2 - ycoord); 
          pfield->b.x2f(k,j,i) = (orig - xchange) / (xcoord2 - xcoord);
          pfield->b.x3f(k,j,i) = wl[NHYDRO+2];

//	  pfield->b.x1f(k,j,i) = wl[NHYDRO];
//        pfield->b.x2f(k,j,i) = wl[NHYDRO+1];
//        pfield->b.x3f(k,j,i) = wl[NHYDRO+2];
        
       if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) += 0.5*(SQR(pfield->b.x1f(k,j,i))
            + SQR(pfield->b.x2f(k,j,i)) + SQR(pfield->b.x3f(k,j,i)));
        }
      }
    }}
       
   
    // end by adding bi.x1 at ie+1, bi.x2 at je+1, and bi.x3 at ke+1

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      pfield->b.x1f(k,j,ie+1) = pfield->b.x1f(k,j,ie);
    }}
    for (int k=ks; k<=ke; ++k) {
    for (int i=is; i<=ie; ++i) {
      pfield->b.x2f(k,je+1,i) = pfield->b.x2f(k,je,i);
    }}
    for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie; ++i) {
      pfield->b.x3f(ke+1,j,i) = pfield->b.x3f(ke,j,i);
    }}
  }
  return;
}

//output total magnetic field strength
void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  for(int k=ks; k<=ke; ++k) {
    for(int j=js; j<=je; ++j) {
      for(int i=is; i<=ie; ++i) {
	user_out_var(0,k,j,i) = std::sqrt(SQR(pfield->b.x1f(k,j,i))
				 + SQR(pfield->b.x2f(k,j,i))
				 + SQR(pfield->b.x3f(k,j,i)));
      }
    }
  }

}

//----------------------------------------------------------------------------------------
//! \fn void stinner_ix1()
//  \brief sets boundary condition on left x boundary (iib)
void stinner_ix1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                       FaceField &b, Real time, Real dt,
                       int is, int ie, int js, int je, int ks, int ke, int ngh){
  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for (int i=1; i<=ngh; ++i) {
      Real xcoord = pmb->pcoord->x1v(is-i);
      Real ycoord = pmb->pcoord->x2v(j);
	  
      prim(IDN,k,j,is-i) = return_d(xcoord,ycoord,1000,1.0,1.0);
      prim(IVX,k,j,is-i) = 1.0;
      prim(IVY,k,j,is-i) = 0.0;
      prim(IVZ,k,j,is-i) = 0.0;
      prim(IPR,k,j,is-i) = 1.0;

      if(MAGNETIC_FIELDS_ENABLED) {
	
        pmb->ruser_meshblock_data[1](i-1,j) = b_potential(xcoord,ycoord,1000,1.0);
      }
     
    }
  }}
  if(MAGNETIC_FIELDS_ENABLED) {
    for(int k=ks; k<=ke; ++k) {
    for(int j=js; j<=je; ++j) {
    for(int i=1; i<=ngh; ++i) {

      Real xcoord = pmb->pcoord->x1v(is-i);
      Real ycoord = pmb->pcoord->x2v(j);
      int xdiff = (i < ngh) ? 1 : -1;
      int ydiff = (j < je) ? 1 : -1;	
      Real ycoord2 = pmb->pcoord->x2v(j+ydiff);
      Real xcoord2 = pmb->pcoord->x1v(is-i+xdiff);
      Real orig = pmb->ruser_meshblock_data[1](i-1, j);
      Real ychange = pmb->ruser_meshblock_data[1](i-1, j+ydiff); 
      Real xchange = pmb->ruser_meshblock_data[1](i+xdiff-1, j);
      //taking curl of magnetic potential function   	
      //discrete approxIMation of partial derivative 
      b.x1f(k,j,is-i) = (ychange - orig) / (ycoord2 - ycoord); 
      b.x2f(k,j,is-i) = (orig - xchange) / (xcoord2 - xcoord);
      b.x3f(k,j,is-i) = wl[NHYDRO+2];

//      b.x1f(k,j,i) = wl[NHYDRO];
//      b.x2f(k,j,i) = wl[NHYDRO+1];
//      b.x3f(k,j,i) = wl[NHYDRO+2];
        
      if (NON_BAROTROPIC_EOS) {
         pmb->phydro->u(IEN,k,j,is-i) += 0.5*(SQR(b.x1f(k,j,is-i))
           + SQR(b.x2f(k,j,is-i)) + SQR(b.x3f(k,j,is-i)));
      }

      }
    }}
  }
}


void reflect_ox1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int is, int ie, int js, int je, int ks, int ke, int ngh) {
 // copy hydro variables into ghost zones, reflecting v1
  for (int n=0; n<(NHYDRO); ++n) {
    if (n==(IVX)) {
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          prim(IVX,k,j,ie+i) = -prim(IVX,k,j,(ie-i+1));  // reflect 1-velocity
        }
      }}
    } else {
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          prim(n,k,j,ie+i) = prim(n,k,j,(ie-i+1));
        }
      }}
    }
  }


  // copy face-centered magnetic fields into ghost zones, reflecting b1
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        b.x1f(k,j,(ie+i+1)) = -b.x1f(k,j,(ie-i+1  ));  // reflect 1-field
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        b.x2f(k,j,(ie+i)) =  b.x2f(k,j,(ie-i+1));
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        b.x3f(k,j,(ie+i)) =  b.x3f(k,j,(ie-i+1));
      }
    }}
  }

 
}

//1d/2d
//get position of shock
Real pos(MeshBlock *pmb, int iout) { 
  Real m1 = 1.0, m2 = 1.0, lsum = 0, rsum = 0;
  int pos = 0;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  for (int k=ks; k<=ke; ++k) {
  for (int i=ie; i>=is; --i) {
    rsum = 0;
    for (int j=js; j<=je; ++j) {
	rsum += pmb->phydro->w(IVX,k,j,i);     
    }
    if(i == ie) {
	lsum = rsum;
    }
    if (abs(lsum - rsum) > 5.0) {
      return pmb->pcoord->x1v(i+1);
    }
    lsum = rsum;
  }}
  return 5.0;
}

//1d
//checking first hydrodynamic jump condition
Real firstjump_onedimension(MeshBlock *pmb, int iout) {
  Real m1 = 1.0, m2 = 1.0, rd, ld;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie-4; ++i) {
      ld = pmb->phydro->u(IDN,k,j,i);
      rd = pmb->phydro->u(IDN,k,j,i+4);
      if(rd >= ld + 0.3) {
 	Real lu = pmb->phydro->w(IVX,k,j,i)+1.16;
	Real ru = pmb->phydro->w(IVX,k,j,i+4)+1.16;
	m1 = lu*ld;
	m2 = ru*pmb->phydro->u(IDN,k,j,i+4);
	break;
      }
    }
  }}

  return m1/m2;

}

//1d
//checking second hydrodynamic jump condition
Real secjump_onedimension(MeshBlock *pmb, int iout) {
  Real m1 = 1.0, m2 = 1.0, rd, ld;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie-4; ++i) {
      ld = pmb->phydro->u(IDN,k,j,i);
      rd = pmb->phydro->u(IDN,k,j,i+4);
      if(rd >= ld + 0.3) {
 	Real lu = pmb->phydro->w(IVX,k,j,i)+1.16;
	Real ru = pmb->phydro->w(IVX,k,j,i+4)+1.16;
	Real lp = pmb->phydro->u(IPR,k,j,i);
	Real rp = pmb->phydro->u(IPR,k,j,i+4);
	m1 = lu*lu*ld + lp;
	m2 = ru*ru*rd + rp;
//	m1 = pmb->phydro->u(IDN,k,j,i);
//	m2 = pmb->phydro->w(IDN,k,j,i);
	break;
      }
    }
  }}

  return m1/m2;


}

//2d
//1st jump condition
Real firstjump_twodimensions(MeshBlock *pmb, int iout) {
  Real m1 = 1.0, m2 = 1.0, lsum = 0, rsum = 0;
  int pos = 0;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  for (int k=ks; k<=ke; ++k) {
  for (int i=ie; i>=is; --i) {
    rsum = 0;
    for (int j=js; j<=je; ++j) {
	rsum += pmb->phydro->w(IVX,k,j,i);     
    }
    if(i == ie) {
	lsum = rsum;
    }
    if (abs(lsum - rsum) > 5.0) {
      pos = i+1;
      break;
    }
    lsum = rsum;
  }}
  lsum = 0, rsum = 0;
  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for(int i=std::max(is, pos-2); i<=pos; ++i) {
	lsum += (pmb->phydro->w(IVX,k,j,i)+1.16)*pmb->phydro->w(IDN,k,j,i);
	}
  }}
  Real diff1 = pos - std::max(is, pos-2) + 1;
  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for(int i=std::min(pos+3, ie); i<=std::min(pos+5, ie); ++i) {
	rsum += (pmb->phydro->w(IVX,k,j,i)+1.16)*pmb->phydro->w(IDN,k,j,i);
	}
  }}
  Real diff2 = std::min(pos+5, ie) - std::min(pos+3, ie) + 1;
  lsum /= (js-je+1)*diff1;
  rsum /= (js-je+1)*diff2;
  return (rsum == 0) ? 0.0 : lsum/rsum;

}

Real secjump_twodimensions(MeshBlock *pmb, int iout) {
  Real m1 = 1.0, m2 = 1.0, lsum = 0, rsum = 0;
  int pos = 0;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  for (int k=ks; k<=ke; ++k) {
  for (int i=ie; i>=is; --i) {
    rsum = 0;
    for (int j=js; j<=je; ++j) {
	rsum += pmb->phydro->w(IVX,k,j,i);     
    }
    if(i == ie) {
	lsum = rsum;
    }
    if (abs(lsum - rsum) > 5.0) {
      pos = i+1;
      break;
    }
    lsum = rsum;
  }}
  lsum = 0, rsum = 0;
  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for(int i=std::max(is, pos-2); i<=pos; ++i) {
	lsum += (pmb->phydro->w(IVX,k,j,i)+1.16)*pmb->phydro->w(IDN,k,j,i)*(pmb->phydro->w(IVX,k,j,i)+1.16) + 
		pmb->phydro->w(IPR,k,j,i);
	}
  }}
  Real diff1 = pos - std::max(is, pos-2) + 1;
  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for(int i=std::min(pos+3, ie); i<=std::min(pos+5, ie); ++i) {		       rsum += (pmb->phydro->w(IVX,k,j,i)+1.16)*pmb->phydro->w(IDN,k,j,i)*(pmb->phydro->w(IVX,k,j,i)+1.16) + 
		pmb->phydro->w(IPR,k,j,i);
    }
  }}
  Real diff2 = std::min(pos+5, ie) - std::min(pos+3, ie) + 1;
  lsum /= (js-je+1)*diff1;
  rsum /= (js-je+1)*diff2;
  return (rsum == 0) ? 0.0 : lsum/rsum;

}

//check divb
Real divergenceb(MeshBlock *pmb, int iout) {
  Real divb=0;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  AthenaArray<Real> face1, face2p, face2m, face3p, face3m;
  FaceField &b = pmb->pfield->b;

  face1.NewAthenaArray((ie-is)+2*NGHOST+2);
  face2p.NewAthenaArray((ie-is)+2*NGHOST+1);
  face2m.NewAthenaArray((ie-is)+2*NGHOST+1);
  face3p.NewAthenaArray((ie-is)+2*NGHOST+1);
  face3m.NewAthenaArray((ie-is)+2*NGHOST+1);

  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      pmb->pcoord->Face1Area(k,   j,   is, ie+1, face1);
      pmb->pcoord->Face2Area(k,   j+1, is, ie,   face2p);
      pmb->pcoord->Face2Area(k,   j,   is, ie,   face2m);
      pmb->pcoord->Face3Area(k+1, j,   is, ie,   face3p);
      pmb->pcoord->Face3Area(k,   j,   is, ie,   face3m);
      for(int i=is; i<=ie; i++) {
        divb+=(face1(i+1)*b.x1f(k,j,i+1)-face1(i)*b.x1f(k,j,i)
              +face2p(i)*b.x2f(k,j+1,i)-face2m(i)*b.x2f(k,j,i)
              +face3p(i)*b.x3f(k+1,j,i)-face3m(i)*b.x3f(k,j,i));
      }
    }
  }

  return divb;
}

//checking equation 25.20 in shu "the physics of astrophysics" vol. ii
Real bpara_jc(MeshBlock *pmb, int iout) {
  Real m1 = 1.0, m2 = 1.0, lsum = 0, rsum = 0;
  int pos = 0;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  FaceField &b = pmb->pfield->b;
  for (int k=ks; k<=ke; ++k) {
  for (int i=ie; i>=is; --i) {
    rsum = 0;
    for (int j=js; j<=je; ++j) {
	rsum += pmb->phydro->w(IVX,k,j,i);     
    }
    if(i == ie) {
	lsum = rsum;
    }
    if (abs(lsum - rsum) > 5.0) {
      pos = i+1;
      break;
    }
    lsum = rsum;
  }}
  lsum = 0, rsum = 0;
  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for(int i=std::max(is, pos-2); i<=pos; ++i) {
	lsum += b.x1f(k,j,i)*(pmb->phydro->w(IVY,k,j,i)) 
		- b.x2f(k,j,i)*(pmb->phydro->w(IVX,k,j,i)+1.16);
    }
  }}
  Real diff1 = pos - std::max(is, pos-2) + 1;
  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for(int i=std::min(pos+3, ie); i<=std::min(pos+5, ie); ++i) {		       rsum += b.x1f(k,j,i)*(pmb->phydro->w(IVY,k,j,i)) 
		- b.x2f(k,j,i)*(pmb->phydro->w(IVX,k,j,i)+1.16);
    }
  }}
  Real diff2 = std::min(pos+5, ie) - std::min(pos+3, ie) + 1;
  lsum /= (js-je+1)*diff1;
  rsum /= (js-je+1)*diff2;
  return lsum - rsum;
}

//checking equation 25.21
Real bperp_jc(MeshBlock *pmb, int iout) {
  Real m1 = 1.0, m2 = 1.0, lsum = 0, rsum = 0;
  int pos = 0;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  FaceField &b = pmb->pfield->b;
  for (int k=ks; k<=ke; ++k) {
  for (int i=ie; i>=is; --i) {
    rsum = 0;
    for (int j=js; j<=je; ++j) {
	rsum += pmb->phydro->w(IVX,k,j,i);     
    }
    if(i == ie) {
	lsum = rsum;
    }
    if (abs(lsum - rsum) > 5.0) {
      pos = i+1;
      break;
    }
    lsum = rsum;
  }}
  lsum = 0, rsum = 0;
  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for(int i=std::max(is, pos-2); i<=pos; ++i) {
	lsum += b.x1f(k,j,i);
    }
  }}
  Real diff1 = pos - std::max(is, pos-2) + 1;
  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for(int i=std::min(pos+3, ie); i<=std::min(pos+5, ie); ++i) {		       rsum += b.x1f(k,j,i);
    }
  }}
  Real diff2 = std::min(pos+5, ie) - std::min(pos+3, ie) + 1;
  lsum /= (js-je+1)*diff1;
  rsum /= (js-je+1)*diff2;
  return lsum - rsum;


}

//checking equation 25.22
Real mhd_jump(MeshBlock *pmb, int iout) {
  Real m1 = 1.0, m2 = 1.0, lsum = 0, rsum = 0;
  int pos = 0;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  FaceField &b = pmb->pfield->b;
  for (int k=ks; k<=ke; ++k) {
  for (int i=ie; i>=is; --i) {
    rsum = 0;
    for (int j=js; j<=je; ++j) {
	rsum += pmb->phydro->w(IVX,k,j,i);     
    }
    if(i == ie) {
	lsum = rsum;
    }
    if (abs(lsum - rsum) > 5.0) {
      pos = i+1;
      break;
    }
    lsum = rsum;
  }}
  lsum = 0, rsum = 0;
  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for(int i=std::max(is, pos-2); i<=pos; ++i) {
	lsum += 4*pi*pmb->phydro->w(IDN,k,j,i)*(pmb->phydro->w(IVX,k,j,i)+1.16)*(pmb->phydro->w(IVY,k,j,i));
	lsum -= b.x1f(k,j,i)*b.x2f(k,j,i);	
    }
  }}
  Real diff1 = pos - std::max(is, pos-2) + 1;
  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for(int i=std::min(pos+3, ie); i<=std::min(pos+5, ie); ++i) {		       rsum += 4*pi*pmb->phydro->w(IDN,k,j,i)*(pmb->phydro->w(IVX,k,j,i)+1.16)*(pmb->phydro->w(IVY,k,j,i));
	rsum -= b.x1f(k,j,i)*b.x2f(k,j,i);	
    }
  }}
  Real diff2 = std::min(pos+5, ie) - std::min(pos+3, ie) + 1;
  lsum /= (js-je+1)*diff1;
  rsum /= (js-je+1)*diff2;
  return lsum - rsum;

}


