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

static Real d, u, p;
Real wl[NHYDRO+NFIELD];
Real wr[NHYDRO+NFIELD];

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

//probability power law distrIBution
Real power_law(Real l, Real a) {
  return 1/(1+pow(l*a, 8.0/3.0));
}

//returns density value for given x, y
Real return_d(Real x, Real y, int steps, Real variance, Real l, Real avg) {
  std::random_device rd;
  std::mt19937 gen(rd());  
  std::uniform_real_distribution<> dis(0.0, 1.0);

  Real k = 0.005*l;
  Real res = 0;
  Real ratio = pow(1000, 1.0/steps);
  Real delta_k = k;
  Real prev_k = k;
  //return mult;
  Real c, temp = 0;

  //find normalization constant
  for(int i=0; i<steps; i++) {
    temp += (2 * pi * k * delta_k * power_law(l, k));
    prev_k = k;
    k *= ratio;
    delta_k = k - prev_k;
  }

  c = variance/temp;
  //return c;  
  
  k = 0.005*l;
  delta_k = k;
  prev_k = k;

  //sum over wave modes
  for(int i=0; i<steps; i++) {
    //generate random angles between 0 and 2pi	  
    Real theta = dis(gen);
    theta *= (2*pi);
    Real phi = dis(gen);
    phi *= (2*pi);
    
    //calculate
    Real f1 = sqrt(c * 4 * pi * k * delta_k * power_law(l, k));
    Real f2 = cos(k * x * cos(theta) + k * y * sin(theta) + phi);
    res += (f1 * f2);
    //f_values.push_back(f1*f2);
    prev_k = k;
    k *= ratio;
    delta_k = k - prev_k;
  }
  return avg*exp(-0.11957+res);  
}

//returns magnitude potential function at each (x, y)
Real b_potential(Real x, Real y, int steps, Real variance, Real l, Real avg) {
  std::random_device rd;
  std::mt19937 gen(rd());  
  std::uniform_real_distribution<> dis(0.0, 1.0);

  Real k = 0.005*l;
  Real res = 0;
  Real ratio = pow(1000, 1.0/steps);
  Real delta_k = k;
  Real prev_k = k;
  //return mult;
  Real c, temp = 0;

  //find normalization constant
  for(int i=0; i<steps; i++) {
    temp += (2 * pi * k * delta_k * power_law(l, k));
    prev_k = k;
    k *= ratio;
    delta_k = k - prev_k;
  }

  c = variance/temp;
  //return c;  
  
  k = 0.005*l;
  delta_k = k;
  prev_k = k;

  //sum over wave modes
  for(int i=0; i<steps; i++) {
    //generate random angles between 0 and 2pi	  
    Real theta = dis(gen);
    theta *= (2*pi);
    Real phi = dis(gen);
    phi *= (2*pi);
    
    //calculate
    Real f1 = sqrt(c * 4 * pi * k * delta_k * power_law(l, k));
    Real f2 = cos(k * x * cos(theta) + k * y * sin(theta) + phi);
    res += (f1 * f2);
    //f_values.push_back(f1*f2);
    prev_k = k;
    k *= ratio;
    delta_k = k - prev_k;
  }
  return res;
}


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
  AllocateRealUserMeshBlockDataField(2);
  ruser_meshblock_data[0].NewAthenaArray(260, 260);
  ruser_meshblock_data[1].NewAthenaArray(4, 260);
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
//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief calculate l1 errors in sod (hydro) and rj2a (mhd) tests
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  MeshBlock *pmb = pblock;

  if (!pin->GetOrAddBoolean("problem","compute_error",false)) return;

  // read shock direction and set array indices
  int shk_dir = pin->GetInteger("problem","shock_dir");
  int IM1,IM2,IM3,IB1,IB2,IB3;
  if (shk_dir == 1) {
    IM1 = IM1; IM2 = IM2; IM3 = IM3;
    IB1 = IB1; IB2 = IB2; IB3 = IB3;
  } else if (shk_dir == 2) {
    IM1 = IM2; IM2 = IM3; IM3 = IM1;
    IB1 = IB2; IB2 = IB3; IB3 = IB1;
  } else {
    IM1 = IM3; IM2 = IM1; IM3 = IM2;
    IB1 = IB3; IB2 = IB1; IB3 = IB2;
  }

  // initialize errors to zero
  Real err[NHYDRO+NFIELD];
  for (int i=0; i<(NHYDRO+NFIELD); ++i) err[i]=0.0;

  // errors in rj2a test (dai & woodward 1994 tables ia and IB)
  if (MAGNETIC_FIELDS_ENABLED) {
    Real xfp = 2.2638*tlim;
    Real xrp = (0.53432 + 1.0/std::sqrt(pi*1.309))*tlim;
    Real xsp = (0.53432 + 0.48144/1.309)*tlim;
    Real xc = 0.57538*tlim;
    Real xsm = (0.60588 - 0.51594/1.4903)*tlim;
    Real xrm = (0.60588 - 1.0/std::sqrt(pi*1.4903))*tlim;
    Real xfm = (1.2 - 2.3305/1.08)*tlim;
    Real gm1 = pmb->peos->GetGamma() - 1.0;
    for (int k=pmb->ks; k<=pmb->ke; k++) {
    for (int j=pmb->js; j<=pmb->je; j++) {
      for (int i=pmb->is; i<=pmb->ie; i++) {
        Real r, d0, mx, my, mz, e0, bx, by, bz;
        if (shk_dir == 1) r = pmb->pcoord->x1v(i);
        if (shk_dir == 2) r = pmb->pcoord->x2v(j);
        if (shk_dir == 3) r = pmb->pcoord->x3v(k);

        bx = 2.0/std::sqrt(4.0*pi);
        if (r > xfp) {
          d0 = 1.0;
          mx = 0.0;
          my = 0.0;
          mz = 0.0;
          by = 4.0/std::sqrt(4.0*pi);
          bz = 2.0/std::sqrt(4.0*pi);
          e0 = 1.0/gm1 + 0.5*((mx*mx+my*my+mz*mz)/d0 + (bx*bx+by*by+bz*bz));
        } else if (r > xrp) {
          d0 = 1.3090;
          mx = 0.53432*d0;
          my = -0.094572*d0;
          mz = -0.047286*d0;
          by = 5.3452/std::sqrt(4.0*pi);
          bz = 2.6726/std::sqrt(4.0*pi);
          e0 = 1.5844/gm1 + 0.5*((mx*mx+my*my+mz*mz)/d0 + (bx*bx+by*by+bz*bz));
        } else if (r > xsp) {
          d0 = 1.3090;
          mx = 0.53432*d0;
          my = -0.18411*d0;
          mz = 0.17554*d0;
          by = 5.7083/std::sqrt(4.0*pi);
          bz = 1.7689/std::sqrt(4.0*pi);
          e0 = 1.5844/gm1 + 0.5*((mx*mx+my*my+mz*mz)/d0 + (bx*bx+by*by+bz*bz));
        } else if (r > xc) {
          d0 = 1.4735;
          mx = 0.57538*d0;
          my = 0.047601*d0;
          mz = 0.24734*d0;
          by = 5.0074/std::sqrt(4.0*pi);
          bz = 1.5517/std::sqrt(4.0*pi);
          e0 = 1.9317/gm1 + 0.5*((mx*mx+my*my+mz*mz)/d0 + (bx*bx+by*by+bz*bz));
        } else if (r > xsm) {
          d0 = 1.6343;
          mx = 0.57538*d0;
          my = 0.047601*d0;
          mz = 0.24734*d0;
          by = 5.0074/std::sqrt(4.0*pi);
          bz = 1.5517/std::sqrt(4.0*pi);
          e0 = 1.9317/gm1 + 0.5*((mx*mx+my*my+mz*mz)/d0 + (bx*bx+by*by+bz*bz));
        } else if (r > xrm) {
          d0 = 1.4903;
          mx = 0.60588*d0;
          my = 0.22157*d0;
          mz = 0.30125*d0;
          by = 5.5713/std::sqrt(4.0*pi);
          bz = 1.7264/std::sqrt(4.0*pi);
          e0 = 1.6558/gm1 + 0.5*((mx*mx+my*my+mz*mz)/d0 + (bx*bx+by*by+bz*bz));
        } else if (r > xfm) {
          d0 = 1.4903;
          mx = 0.60588*d0;
          my = 0.11235*d0;
          mz = 0.55686*d0;
          by = 5.0987/std::sqrt(4.0*pi);
          bz = 2.8326/std::sqrt(4.0*pi);
          e0 = 1.6558/gm1 + 0.5*((mx*mx+my*my+mz*mz)/d0 + (bx*bx+by*by+bz*bz));
        } else {
          d0 = 1.08;
          mx = 1.2*d0;
          my = 0.01*d0;
          mz = 0.5*d0;
          by = 3.6/std::sqrt(4.0*pi);
          bz = 2.0/std::sqrt(4.0*pi);
          e0 = 0.95/gm1 + 0.5*((mx*mx+my*my+mz*mz)/d0 + (bx*bx+by*by+bz*bz));
        }

        err[IDN] += fabs(d0 - pmb->phydro->u(IDN,k,j,i));
        err[IM1] += fabs(mx - pmb->phydro->u(IM1,k,j,i));
        err[IM2] += fabs(my - pmb->phydro->u(IM2,k,j,i));
        err[IM3] += fabs(mz - pmb->phydro->u(IM3,k,j,i));
        err[IEN] += fabs(e0 - pmb->phydro->u(IEN,k,j,i));
        err[NHYDRO + IB1] += fabs(bx - pmb->pfield->bcc(IB1,k,j,i));
        err[NHYDRO + IB2] += fabs(by - pmb->pfield->bcc(IB2,k,j,i));
        err[NHYDRO + IB3] += fabs(bz - pmb->pfield->bcc(IB3,k,j,i));
      }
    }}

  // errors in sod solution
  } else {
    // positions of shock, contact, head and foot of rarefaction for sod test
    Real xs = 1.7522*tlim;
    Real xc = 0.92745*tlim;
    Real xf = -0.07027*tlim;
    Real xh = -1.1832*tlim;

    for (int k=pmb->ks; k<=pmb->ke; k++) {
    for (int j=pmb->js; j<=pmb->je; j++) {
      for (int i=pmb->is; i<=pmb->ie; i++) {
        Real r,d0,m0,e0;
        if (shk_dir == 1) r = pmb->pcoord->x1v(i);
        if (shk_dir == 2) r = pmb->pcoord->x2v(j);
        if (shk_dir == 3) r = pmb->pcoord->x3v(k);

        if (r > xs) {
          d0 = 0.125;
          m0 = 0.0;
          e0 = 0.25;
        } else if (r > xc) {
          d0 = 0.26557;
          m0 = 0.92745*d0;
          e0 = 0.87204;
        } else if (r > xf) {
          d0 = 0.42632;
          m0 = 0.92745*d0;
          e0 = 0.94118;
        } else if (r > xh) {
          Real v0 = 0.92745*(r-xh)/(xf-xh);
          d0 = 0.42632*pow((1.0+0.20046*(0.92745-v0)),5);
          m0 = v0*d0;
          e0 = (0.30313*pow((1.0+0.20046*(0.92745-v0)),7))/0.4 + 0.5*d0*v0*v0;
        } else {
          d0 = 1.0;
          m0 = 0.0;
          e0 = 2.5;
        }
        err[IDN] += fabs(d0  - pmb->phydro->u(IDN,k,j,i));
        err[IM1] += fabs(m0  - pmb->phydro->u(IM1,k,j,i));
        err[IM2] += fabs(0.0 - pmb->phydro->u(IM2,k,j,i));
        err[IM3] += fabs(0.0 - pmb->phydro->u(IM3,k,j,i));
        err[IEN] += fabs(e0  - pmb->phydro->u(IEN,k,j,i));
      }
    }}
  }

  // normalize errors by number of cells, compute rms
  for (int i=0; i<(NHYDRO+NFIELD); ++i) {
    err[i] = err[i]/static_cast<Real>(GetTotalCells());
  }
  Real rms_err = 0.0;
  for (int i=0; i<(NHYDRO+NFIELD); ++i) rms_err += SQR(err[i]);
  rms_err = std::sqrt(rms_err);

  // open output FILE and write out errors
  std::string fname;
  fname.assign("shock-errors.dat");
  std::stringstream msg;
  FILE *pfile;

  // the FILE exists -- reopen the FILE in append mode
  if ((pfile = fopen(fname.c_str(),"r")) != NULL) {
    if ((pfile = freopen(fname.c_str(),"a",pfile)) == NULL) {
      msg << "### fatal error in function [Mesh::UserWorkAfterLoop]"
          << std::endl << "error output FILE could not be opened" <<std::endl;
      throw std::runtime_error(msg.str().c_str());
    }

  // the FILE does not exist -- open the FILE in write mode and add headers
  } else {
    if ((pfile = fopen(fname.c_str(),"w")) == NULL) {
      msg << "### fatal error in function [Mesh::UserWorkAfterLoop]"
          << std::endl << "error output FILE could not be opened" <<std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    fprintf(pfile,"# nx1  nx2  nx3  ncycle  rms-error  d  m1  m2  m3  e");
    if (MAGNETIC_FIELDS_ENABLED) fprintf(pfile,"  b1c  b2c  b3c");
    fprintf(pfile,"\n");
  }

  // write errors
  fprintf(pfile,"%d  %d",pmb->block_size.nx1,pmb->block_size.nx2);
  fprintf(pfile,"  %d  %d  %e",pmb->block_size.nx3,ncycle,rms_err);
  fprintf(pfile,"  %e  %e  %e  %e  %e",err[IDN],err[IM1],err[IM2],err[IM3],err[IEN]);
  if (MAGNETIC_FIELDS_ENABLED) {
    fprintf(pfile,"  %e  %e  %e",err[NHYDRO+IB1],err[NHYDRO+IB2],err[NHYDRO+IB3]);
  }
  fprintf(pfile,"\n");
  fclose(pfile);

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief problem generator for the shock tube tests
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  std::stringstream msg;

  // parse shock direction: {1,2,3} -> {x1,x2,x3}
  int shk_dir = pin->GetInteger("problem","shock_dir");

  // parse shock location (must be inside grid)
  Real xshock = pin->GetReal("problem","xshock");
  if (shk_dir == 1 && (xshock < pmy_mesh->mesh_size.x1min ||
                       xshock > pmy_mesh->mesh_size.x1max)) {
    msg << "### fatal error in problem generator" << std::endl << "xshock="
        << xshock << " lies outside x1 domain for shkdir=" << shk_dir << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (shk_dir == 2 && (xshock < pmy_mesh->mesh_size.x2min ||
                       xshock > pmy_mesh->mesh_size.x2max)) {
    msg << "### fatal error in problem generator" << std::endl << "xshock="
        << xshock << " lies outside x2 domain for shkdir=" << shk_dir << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (shk_dir == 3 && (xshock < pmy_mesh->mesh_size.x3min ||
                       xshock > pmy_mesh->mesh_size.x3max)) {
    msg << "### fatal error in problem generator" << std::endl << "xshock="
        << xshock << " lies outside x3 domain for shkdir=" << shk_dir << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // parse left state read from input FILE: dl,ul,vl,wl,[pl]
  //Real wl[NHYDRO+NFIELD];
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

  // parse right state read from input FILE: dr,ur,vr,wr,[pr]
  //Real wr[NHYDRO+NFIELD];
  wr[IDN] = pin->GetReal("problem","dr");
  wr[IVX] = pin->GetReal("problem","ur");
  wr[IVY] = pin->GetReal("problem","vr");
  wr[IVZ] = pin->GetReal("problem","wr");
  if (NON_BAROTROPIC_EOS) wr[IPR] = pin->GetReal("problem","pr");
  if (MAGNETIC_FIELDS_ENABLED) {
    wr[NHYDRO  ] = pin->GetReal("problem","bxr");
    wr[NHYDRO+1] = pin->GetReal("problem","byr");
    wr[NHYDRO+2] = pin->GetReal("problem","bzr");
  }

// initialize the discontinuity in the hydro variables ---------------------------------

  switch(shk_dir) {

//--- shock in 1-direction
  case 1:
    d = wl[IDN];
    u = wl[IVX];
    p = wl[IPR];    
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        if (pcoord->x1v(i) < xshock) {
	  Real xcoord = pcoord->x1f(i);
          Real ycoord = pcoord->x2f(j);

          phydro->u(IDN,k,j,i) = return_d(xcoord, ycoord,1000,0.254031,0.5,d);
	  //phydro->u(IDN,k,j,i) = wl[IDN];
          phydro->u(IM1,k,j,i) = wl[IVX]*phydro->u(IDN,k,j,i);
          phydro->u(IM2,k,j,i) = wl[IVY]*wl[IDN];
          phydro->u(IM3,k,j,i) = wl[IVZ]*wl[IDN];
          if (NON_BAROTROPIC_EOS) phydro->u(IEN,k,j,i) =
            wl[IPR]/(peos->GetGamma() - 1.0)
            + 0.5*phydro->u(IDN,k,j,i)*(wl[IVX]*wl[IVX] + wl[IVY]*wl[IVY] + wl[IVZ]*wl[IVZ]);
				
					if (CLESS_ENABLED) {
						pcless->u(IDN ,k,j,i) = phydro->u(IDN,k,j,i); 
						pcless->u(IM1 ,k,j,i) = phydro->u(IM1,k,j,i);
						pcless->u(IM2 ,k,j,i) = phydro->u(IM2,k,j,i);
						pcless->u(IM3 ,k,j,i) = phydro->u(IM3,k,j,i);
						pcless->u(IE11,k,j,i) = wl[IPR] + wl[IDN]*wl[IVX]*wl[IVX];
						pcless->u(IE22,k,j,i) = wl[IPR] + wl[IDN]*wl[IVY]*wl[IVY];
						pcless->u(IE33,k,j,i) = wl[IPR] + wl[IDN]*wl[IVZ]*wl[IVZ];
						pcless->u(IE12,k,j,i) = wl[IPR] + wl[IDN]*wl[IVX]*wl[IVY];
						pcless->u(IE13,k,j,i) = wl[IPR] + wl[IDN]*wl[IVX]*wl[IVZ];
						pcless->u(IE23,k,j,i) = wl[IPR] + wl[IDN]*wl[IVY]*wl[IVZ];
					}
        } else {
          phydro->u(IDN,k,j,i) = wr[IDN];
          phydro->u(IM1,k,j,i) = wr[IVX]*wr[IDN];
          phydro->u(IM2,k,j,i) = wr[IVY]*wr[IDN];
          phydro->u(IM3,k,j,i) = wr[IVZ]*wr[IDN];
          if (NON_BAROTROPIC_EOS) phydro->u(IEN,k,j,i) =
            wr[IPR]/(peos->GetGamma() - 1.0)
            + 0.5*wr[IDN]*(wr[IVX]*wr[IVX] + wr[IVY]*wr[IVY] + wr[IVZ]*wr[IVZ]);
					
					if (CLESS_ENABLED) {
						pcless->u(IDN ,k,j,i) = phydro->u(IDN,k,j,i); 
						pcless->u(IM1 ,k,j,i) = phydro->u(IM1,k,j,i);
						pcless->u(IM2 ,k,j,i) = phydro->u(IM2,k,j,i);
						pcless->u(IM3 ,k,j,i) = phydro->u(IM3,k,j,i);
						pcless->u(IE11,k,j,i) = wr[IPR] + wr[IDN]*wr[IVX]*wr[IVX];
						pcless->u(IE22,k,j,i) = wr[IPR] + wr[IDN]*wr[IVY]*wr[IVY];
						pcless->u(IE33,k,j,i) = wr[IPR] + wr[IDN]*wr[IVZ]*wr[IVZ];
						pcless->u(IE12,k,j,i) = wr[IPR] + wr[IDN]*wr[IVX]*wr[IVY];
						pcless->u(IE13,k,j,i) = wr[IPR] + wr[IDN]*wr[IVX]*wr[IVZ];
						pcless->u(IE23,k,j,i) = wr[IPR] + wr[IDN]*wr[IVY]*wr[IVZ];
					}
        }
      }
    }}
    break;

//--- shock in 2-direction
  case 2:
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      if (pcoord->x2v(j) < xshock) {
        for (int i=is; i<=ie; ++i) {
          phydro->u(IDN,k,j,i) = wl[IDN];
          phydro->u(IM2,k,j,i) = wl[IVX]*wl[IDN];
          phydro->u(IM3,k,j,i) = wl[IVY]*wl[IDN];
          phydro->u(IM1,k,j,i) = wl[IVZ]*wl[IDN];
          if (NON_BAROTROPIC_EOS) phydro->u(IEN,k,j,i) =
            wl[IPR]/(peos->GetGamma() - 1.0)
            + 0.5*wl[IDN]*(wl[IVX]*wl[IVX] + wl[IVY]*wl[IVY] + wl[IVZ]*wl[IVZ]);
					if (CLESS_ENABLED) {
						pcless->u(IDN ,k,j,i) = phydro->u(IDN,k,j,i); 
						pcless->u(IM2 ,k,j,i) = phydro->u(IM2,k,j,i);
						pcless->u(IM3 ,k,j,i) = phydro->u(IM3,k,j,i);
						pcless->u(IM1 ,k,j,i) = phydro->u(IM1,k,j,i);
						pcless->u(IE22,k,j,i) = wl[IPR] + wl[IDN]*wl[IVX]*wl[IVX];
						pcless->u(IE33,k,j,i) = wl[IPR] + wl[IDN]*wl[IVY]*wl[IVY];
						pcless->u(IE11,k,j,i) = wl[IPR] + wl[IDN]*wl[IVZ]*wl[IVZ];
						pcless->u(IE23,k,j,i) = wl[IPR] + wl[IDN]*wl[IVX]*wl[IVY];
						pcless->u(IE12,k,j,i) = wl[IPR] + wl[IDN]*wl[IVX]*wl[IVZ];
						pcless->u(IE13,k,j,i) = wl[IPR] + wl[IDN]*wl[IVY]*wl[IVZ];
					}
				}

      } else {
        for (int i=is; i<=ie; ++i) {
          phydro->u(IDN,k,j,i) = wr[IDN];
          phydro->u(IM2,k,j,i) = wr[IVX]*wr[IDN];
          phydro->u(IM3,k,j,i) = wr[IVY]*wr[IDN];
          phydro->u(IM1,k,j,i) = wr[IVZ]*wr[IDN];
          if (NON_BAROTROPIC_EOS) phydro->u(IEN,k,j,i) =
            wr[IPR]/(peos->GetGamma() - 1.0)
            + 0.5*wr[IDN]*(wr[IVX]*wr[IVX] + wr[IVY]*wr[IVY] + wr[IVZ]*wr[IVZ]);
					if (CLESS_ENABLED) {
						pcless->u(IDN ,k,j,i) = phydro->u(IDN,k,j,i); 
						pcless->u(IM2 ,k,j,i) = phydro->u(IM2,k,j,i);
						pcless->u(IM3 ,k,j,i) = phydro->u(IM3,k,j,i);
						pcless->u(IM1 ,k,j,i) = phydro->u(IM1,k,j,i);
						pcless->u(IE22,k,j,i) = wr[IPR] + wr[IDN]*wr[IVX]*wr[IVX];
						pcless->u(IE33,k,j,i) = wr[IPR] + wr[IDN]*wr[IVY]*wr[IVY];
						pcless->u(IE11,k,j,i) = wr[IPR] + wr[IDN]*wr[IVZ]*wr[IVZ];
						pcless->u(IE23,k,j,i) = wr[IPR] + wr[IDN]*wr[IVX]*wr[IVY];
						pcless->u(IE12,k,j,i) = wr[IPR] + wr[IDN]*wr[IVX]*wr[IVZ];
						pcless->u(IE13,k,j,i) = wr[IPR] + wr[IDN]*wr[IVY]*wr[IVZ];
					}
        }
      }
    }}
    break;

//--- shock in 3-direction

  case 3:
    for (int k=ks; k<=ke; ++k) {
      if (pcoord->x3v(k) < xshock) {
        for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          phydro->u(IDN,k,j,i) = wl[IDN];
          phydro->u(IM3,k,j,i) = wl[IVX]*wl[IDN];
          phydro->u(IM1,k,j,i) = wl[IVY]*wl[IDN];
          phydro->u(IM2,k,j,i) = wl[IVZ]*wl[IDN];
          if (NON_BAROTROPIC_EOS) phydro->u(IEN,k,j,i) =
            wl[IPR]/(peos->GetGamma() - 1.0)
            + 0.5*wl[IDN]*(wl[IVX]*wl[IVX] + wl[IVY]*wl[IVY] + wl[IVZ]*wl[IVZ]);
        }}
      } else {
        for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          phydro->u(IDN,k,j,i) = wr[IDN];
          phydro->u(IM3,k,j,i) = wr[IVX]*wr[IDN];
          phydro->u(IM1,k,j,i) = wr[IVY]*wr[IDN];
          phydro->u(IM2,k,j,i) = wr[IVZ]*wr[IDN];
          if (NON_BAROTROPIC_EOS) phydro->u(IEN,k,j,i) =
            wr[IPR]/(peos->GetGamma() - 1.0)
            + 0.5*wr[IDN]*(wr[IVX]*wr[IVX] + wr[IVY]*wr[IVY] + wr[IVZ]*wr[IVZ]);
        }}
      }
    }
    break;

  default:
    msg << "### fatal error in problem generator" << std::endl
        << "shock_dir=" << shk_dir << " must be either 1,2, or 3" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

// now set face-centered (interface) magnetic fields -----------------------------------

  if (MAGNETIC_FIELDS_ENABLED) {
      
      for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {

 	  Real xcoord = pcoord->x1v(i);
          Real ycoord = pcoord->x1v(j);
     
          if(MAGNETIC_FIELDS_ENABLED) {
 	    ruser_meshblock_data[0](i, j) = b_potential(xcoord,ycoord,1000,0.254031,2.5,1.0);
        }

	}
    }

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
         if (shk_dir==1 && pcoord->x1v(i) < xshock) {
      	  
	  Real ycoord = pcoord->x1v(j);
	  Real xcoord = pcoord->x1v(i);
          int xdiff = (i < ie) ? 1 : -1;
          int ydiff = (j < je) ? 1 : -1;	
          Real ycoord2 = pcoord->x1v(j+ydiff);
          Real xcoord2 = pcoord->x1v(i+xdiff);
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
        
        } else if (shk_dir==2 && pcoord->x2v(j) < xshock) {
          pfield->b.x1f(k,j,i) = wl[NHYDRO+2];
          pfield->b.x2f(k,j,i) = wl[NHYDRO  ];
          pfield->b.x3f(k,j,i) = wl[NHYDRO+1];
        } else if (shk_dir==3 && pcoord->x3v(k) < xshock) {
          pfield->b.x1f(k,j,i) = wl[NHYDRO+1];
          pfield->b.x2f(k,j,i) = wl[NHYDRO+2];
          pfield->b.x3f(k,j,i) = wl[NHYDRO];
        }

        if (shk_dir==1 && pcoord->x1v(i) >= xshock) {
          pfield->b.x1f(k,j,i) = wr[NHYDRO  ];
          pfield->b.x2f(k,j,i) = wr[NHYDRO+1];
          pfield->b.x3f(k,j,i) = wr[NHYDRO+2];
        } else if (shk_dir==2 && pcoord->x2v(j) >= xshock) {
          pfield->b.x1f(k,j,i) = wr[NHYDRO+2];
          pfield->b.x2f(k,j,i) = wr[NHYDRO  ];
          pfield->b.x3f(k,j,i) = wr[NHYDRO+1];
        } else if (shk_dir==3 && pcoord->x3v(k) >= xshock)  {
          pfield->b.x1f(k,j,i) = wr[NHYDRO+1];
          pfield->b.x2f(k,j,i) = wr[NHYDRO+2];
          pfield->b.x3f(k,j,i) = wr[NHYDRO];
        }
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
  Real quot;
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
  //Real potential_values[ngh+2][je+2];
  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for (int i=1; i<=ngh; ++i) {
      Real xcoord = pmb->pcoord->x1f(is-i);
      Real ycoord = pmb->pcoord->x2f(j);
      prim(IDN,k,j,is-i) = return_d(xcoord,ycoord,1000,0.254031,0.5,wl[IDN]);
      //prim(IDN,k,j,is-i) = d;
      prim(IVX,k,j,is-i) = u;
      prim(IVY,k,j,is-i) = 0.0;
      prim(IVZ,k,j,is-i) = 0.0;
      prim(IPR,k,j,is-i) = p;
      /*
      if(MAGNETIC_FIELDS_ENABLED) {
        pmb->ruser_meshblock_data[1](i-1,j) = b_potential(xcoord,ycoord,1000,0.254031,2.5,1.0);
      }
      */
    }
  }}
  if(MAGNETIC_FIELDS_ENABLED) {
    for(int k=ks; k<=ke; ++k) {
    for(int j=js; j<=je; ++j) {
    for(int i=1; i<=ngh; ++i) {
/*
      Real xcoord = pmb->pcoord->x1v(is-i);
      Real ycoord = pmb->pcoord->x1v(j);
      int xdiff = (i < ngh) ? 1 : -1;
      int ydiff = (j < je) ? 1 : -1;	
      Real ycoord2 = pmb->pcoord->x1v(j+ydiff);
      Real xcoord2 = pmb->pcoord->x1v(is-i+xdiff);
      Real orig = pmb->ruser_meshblock_data[1](i-1, j);
      Real ychange = pmb->ruser_meshblock_data[1](i-1, j+ydiff); 
      Real xchange = pmb->ruser_meshblock_data[1](i+xdiff-1, j);
      //taking curl of magnetic potential function   	
      //discrete approxIMation of partial derivative 
      b.x1f(k,j,is-i) = (ychange - orig) / (ycoord2 - ycoord); 
      b.x2f(k,j,is-i) = (orig - xchange) / (xcoord2 - xcoord);
      b.x3f(k,j,is-i) = wl[NHYDRO+2];
*/
      b.x1f(k,j,i) = wl[NHYDRO];
      b.x2f(k,j,i) = wl[NHYDRO+1];
      b.x3f(k,j,i) = wl[NHYDRO+2];
        
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
#pragma omp sIMd
        for (int i=1; i<=ngh; ++i) {
          prim(IVX,k,j,ie+i) = -prim(IVX,k,j,(ie-i+1));  // reflect 1-velocity
        }
      }}
    } else {
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp sIMd
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
#pragma omp sIMd
      for (int i=1; i<=ngh; ++i) {
        b.x1f(k,j,(ie+i+1)) = -b.x1f(k,j,(ie-i+1  ));  // reflect 1-field
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma omp sIMd
      for (int i=1; i<=ngh; ++i) {
        b.x2f(k,j,(ie+i)) =  b.x2f(k,j,(ie-i+1));
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp sIMd
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


