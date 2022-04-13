// C++ headers
#include <algorithm>
#include <cmath>      // sqrt()
#include <iostream>
#include <limits>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../utils/utils.hpp"
#include "../utils/townsend_cooling.hpp"

// User defined boundary conditions 
void HydrostaticInnerX3(MeshBlock *pmb, Coordinates *pco,
                        AthenaArray<Real> &a,
                        FaceField &b, Real time, Real dt,
                        int il, int iu, int jl, int ju, 
                        int kl, int ku, int ngh);
void HydrostaticOuterX3(MeshBlock *pmb, Coordinates *pco,
                        AthenaArray<Real> &a,
                        FaceField &b, Real time, Real dt,
                        int il, int iu, int jl, int ju, 
                        int kl, int ku, int ngh);

// User defined source functions
void cooling(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim,
             const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc,
             AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar);

// User defined history functions
Real CalculateEmission(MeshBlock *pmb, int iout);
Real CalculateColdGasMass(MeshBlock *pmb, int iout);
Real CalculateColdGasVolume(MeshBlock *pmb, int iout);
Real GetFrameVel(MeshBlock *pmb, int iout);
Real GetBoxHeight(MeshBlock *pmb, int iout);

// Globals
static Real com_vel = 0.0;
static Real frame_vel = 0.0;
static Real box_height = 0.0;
static Real grav_acc = 0.0;

// Cooling
static Cooling cooler;

// Define code units in cgs
static Real const mh        = 1.6605e-24;   // atomic mass unit (g)
static Real const kb        = 1.380648e-16; // boltzmann constant (erg/K)
static Real const mu        = 0.6173;       // mean molecular weight (Solar)
                                            // X = 0.7; Z = 0.02

static Real const unit_len  = 3.086e20;     // 100 pc
static Real const unit_temp = 1.0e6;        // 10^6 K
static Real const unit_n    = 1.0e-4;       // cm^-3

static Real const unit_rho  = unit_n * mh * mu;
static Real const unit_pres = kb * unit_n * unit_temp;
static Real const unit_vel  = sqrt(unit_pres / unit_rho);
static Real const unit_time = unit_len / unit_vel;
static Real const unit_kap  = (unit_pres * unit_len * unit_len) /
                              (unit_time * unit_temp);
static Real const unit_gam  = (unit_pres/unit_time)/(unit_n*unit_n);

//====================================================================================
void Mesh::InitUserMeshData(ParameterInput *pin) {
  cooler = Cooling();

  // Enroll user-defined physical source terms
  EnrollUserExplicitSourceFunction(cooling);

  // Enroll user-defined boundary conditions
  EnrollUserBoundaryFunction(BoundaryFace::inner_x3, HydrostaticInnerX3);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x3, HydrostaticOuterX3);

  // Enroll user defined history outputs
  AllocateUserHistoryOutput(5);
  EnrollUserHistoryOutput(0,GetFrameVel,"frame_vel", UserHistoryOperation::max);
  EnrollUserHistoryOutput(1,GetBoxHeight,"box_height", UserHistoryOperation::max);
  EnrollUserHistoryOutput(2,CalculateColdGasMass,"cold_mass");
  EnrollUserHistoryOutput(3,CalculateColdGasVolume,"cold_vol");
  EnrollUserHistoryOutput(4,CalculateEmission,"emission");

  return;
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(1); // Emission in cell
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real x1, x2, x3;
  Real rd, rp, rvx, rvy, rvz;

  int64_t iseed = -1 - gid;

  Real gam = peos->GetGamma();
  Real T_bg = 1.0;
  Real P_bg = 1.0;
  Real cloud_overdensity = 100.0;
  Real cloud_radius = pin->GetReal("problem","r_cloud"); 

  grav_acc = phydro->hsrc.GetG3();

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        x1 = pcoord->x1v(i);
        x2 = pcoord->x2v(j);
        x3 = pcoord->x3v(k);

        rd = P_bg/T_bg;
        rvx = 0.0;
        rvy = 0.0;
        rvz = 0.0;
        rp = P_bg;

        phydro->u(IDN,k,j,i) = rd;
        if (sqrt(SQR(x1)+SQR(x2)+SQR(x3)) < cloud_radius) {
          // Create cloud with 1% perturbations
          phydro->u(IDN,k,j,i) *= cloud_overdensity*(1- 0.01*ran2(&iseed));
        }

        phydro->u(IM1,k,j,i) = rd*rvx;
        phydro->u(IM2,k,j,i) = rd*rvy;
        phydro->u(IM3,k,j,i) = rd*rvz;
        phydro->u(IEN,k,j,i) = rp/(gam-1.0) + 0.5*(SQR(phydro->u(IM1,k,j,i))
                                                   + SQR(phydro->u(IM2,k,j,i))
                                                   + SQR(phydro->u(IM3,k,j,i)))/rd;
      }
    }
  }

  return;
}

// Only called at the end of the full timestep. 
// Use if you want to do something fully operator split from the hydrodynamics 
// (as opposed to the semi-split explicit source function). However, this can 
// be a bit unsafe. In particular, this function is called after the last 
// conserved-to-primitive call, so if the user changes variables here, the 
// change must be done consistently to both the primitives and conserved variables. 
// Moreover, the function is called after the boundary functions, so the user must 
// update the ghost cells (both internal and external) here to reflect what the 
// boundary functions would have done had they been given the updated state.
void MeshBlock::UserWorkInLoop() {
  for (int k=ks-NGHOST; k<=ke+NGHOST; k++) {
    for (int j=js-NGHOST; j<=je+NGHOST; j++) {
      for (int i=is-NGHOST; i<=ie+NGHOST; i++) {
        Real u_d = phydro->u(IDN,k,j,i);
        Real vel_i = phydro->w(IVZ,k,j,i);

        // Boost to com frame
        Real vel_f = vel_i - com_vel;

        // Update both conserved and primitive variables
        phydro->u(IEN,k,j,i) -= 0.5*u_d*(vel_i*vel_i - vel_f*vel_f);
        phydro->u(IM3,k,j,i) -= u_d*com_vel;
        phydro->w(IVZ,k,j,i) = vel_f;

      }
    }
  }

  return;
}

// Calculate global properties for cloud tracking
void Mesh::UserWorkInLoop() { 
  // Update the distance box has travelled and it's current velocity
  box_height += frame_vel*dt;
  frame_vel += com_vel;

  Real sum_density = 0.0;
  Real sum_momentum = 0.0;

  // Loop over MeshBlocks
  for (int bn=0; bn<nblocal; ++bn) {
    MeshBlock *pmb = my_blocks(bn);
    Hydro *phyd = pmb->phydro;

    // Sum density and momentum over cells
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
    for (int i=pmb->is; i<=pmb->ie; ++i) {
      Real u_d  = phyd->u(IDN,k,j,i);
      Real u_mz = phyd->u(IM3,k,j,i);
      
      Real press = phyd->w(IPR,k,k,i);
      Real temp = press/u_d;

      if (temp < 0.02) { // only include contributions from cold gas
        sum_density += u_d;
        sum_momentum += u_mz;
      }
    }}}
  } // End loop over MeshBlocks

  #ifdef MPI_PARALLEL
  // Sum over all nodes
  MPI_Allreduce(MPI_IN_PLACE, &sum_density, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &sum_momentum, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  #endif

  // Compute center of mass velocity
  Real boost_vel = 0.0;
  if (sum_density > TINY_NUMBER) { // Avoid code crashing when no cold gas
    boost_vel = sum_momentum/sum_density; 
  }

  com_vel = boost_vel;
  
  return;
}

// Lower z boundary
// Exponentially extrapolate pressure and density in the ghost zones
// Zero-gradient extrapolation for velocities
void HydrostaticInnerX3(MeshBlock *pmb, Coordinates *pco,
                         AthenaArray<Real> &prim, FaceField &b,
                         Real time, Real dt,
                         int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int k=1; k<=ngh; k++) {
    for (int j=jl; j<=ju; j++) {
      for (int i=il; i<=iu; i++) {
        prim(IDN,kl-k,j,i) = 1.0;
        prim(IPR,kl-k,j,i) = 1.0;

        prim(IVX,kl-k,j,i) = 0.0; 
        prim(IVY,kl-k,j,i) = 0.0;
        prim(IVZ,kl-k,j,i) = -frame_vel; 
      }
    }
  }
  return;
}

// Upper z boundary
// Exponentially extrapolate pressure and density in the ghost zones
// Zero-gradient extrapolation for velocities
void HydrostaticOuterX3(MeshBlock *pmb, Coordinates *pco,
                         AthenaArray<Real> &prim,
                         FaceField &b, Real time, Real dt,
                         int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int k=1; k<=ngh; k++) {
    for (int j=jl; j<=ju; j++) {
      for (int i=il; i<=iu; i++) {
        prim(IDN,ku+k,j,i) = 1.0;
        prim(IPR,ku+k,j,i) = 1.0;

        prim(IVX,ku+k,j,i) = 0.0; 
        prim(IVY,ku+k,j,i) = 0.0;
        prim(IVZ,ku+k,j,i) = -frame_vel; 
      }
    }
  }
  return;
}

// User Defined Cooling Function
void cooling(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim,
             const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc,
             AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {
  Real g = pmb->peos->GetGamma();

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        // Need to take density and temperature at time step n from cons, not
        // prim because we do not need intermediate step to calculate the
        // cooling function
        Real rho = cons(IDN,k,j,i);
        Real eint = cons(IEN,k,j,i)
                    - 0.5 *(cons(IM1,k,j,i)*cons(IM1,k,j,i)
                          + cons(IM2,k,j,i)*cons(IM2,k,j,i)
                          + cons(IM3,k,j,i)*cons(IM3,k,j,i))/rho;        

        // T = P/rho
        Real temp = eint * (g-1.0)/rho;

        if (temp < 0) {
          std::cout << "Negative Temperature Detected..." << std::endl;
          continue;
        }

        Real temp_new = temp;

        if (temp < 0.5) {
          // Calculate new temperature using the Townsend Algorithm
          // The inputs should be given in cgs units
          Real temp_cgs = temp * unit_temp;
          Real rho_cgs  = rho  * unit_rho;
          Real dt_cgs   = dt   * unit_time;
          temp_new = cooler.townsend(temp_cgs,rho_cgs,dt_cgs)/unit_temp;
        }

        // Enforce temperature floor
        temp_new = std::max(temp_new, cooler.tfloor/unit_temp);

        // Update energy based on change in temperature
        cons(IEN,k,j,i) += (temp_new - temp) * (rho/(g-1.0));

        // Store the change in energy/time in a user defined output variable
        pmb->user_out_var(0,k,j,i) = (temp_new - temp) * (rho/(g-1.0)) / dt;

      }
    }
  }

  // Counter Gravity  
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real src = -dt*1.0*grav_acc;
        cons(IM3,k,j,i) += src;
        if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) += src*prim(IVZ,k,j,i);
      }
    }
  }

  return;
}

// User defined history functions
Real CalculateEmission(MeshBlock *pmb, int iout) {
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  AthenaArray<Real> vol(pmb->ncells1);

  Real sum = 0.0;

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
    pmb->pcoord->CellVolume(k,j,is,ie,vol);
    for (int i=is; i<=ie; i++) {
      sum += pmb->user_out_var(0,k,j,i)*vol(i);
    }
  }}

  return sum;
}

Real CalculateColdGasMass(MeshBlock *pmb, int iout) {
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  AthenaArray<Real> vol(pmb->ncells1);

  Real sum = 0.0;

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
    pmb->pcoord->CellVolume(k,j,is,ie,vol);
    for (int i=is; i<=ie; i++) {
      Real press = pmb->phydro->w(IPR,k,j,i);
      Real dens  = pmb->phydro->w(IDN,k,j,i);
      Real temp = press/dens;
      if (temp < 0.02) {sum += dens*vol(i);}
    }
  }}

  return sum;
}

Real CalculateColdGasVolume(MeshBlock *pmb, int iout) {
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  AthenaArray<Real> vol(pmb->ncells1);

  Real sum = 0.0;

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
    pmb->pcoord->CellVolume(k,j,is,ie,vol);
    for (int i=is; i<=ie; i++) {
      Real press = pmb->phydro->w(IPR,k,j,i);
      Real dens  = pmb->phydro->w(IDN,k,j,i);
      Real temp = press/dens;
      if (temp < 0.02) {sum += vol(i);}
    }
  }}

  return sum;
}

Real GetFrameVel(MeshBlock *pmb, int iout) {
  return frame_vel;
}

Real GetBoxHeight(MeshBlock *pmb, int iout) {
  return box_height;
}

/* Notes
Manipulating the data is safest with enrolling an explicit source function. 
This function is called at the end of every substep in the time integration. 
Its intention is to update the conserved variable array in the function input, 
which has already been updated according to the flux divergence for the substep. 
If you are imposing some constraint (e.g. a temperature ceiling) or doing something implicit, 
you should ignore the primitive input and deal only with the conserved values. 
If you are adding a proper physical source (e.g. explicit cooling) you should 
calculate the conserved change from the not-yet-updated primitive input and apply
it to the conserved values. After this function is called, the code performs 
conserved-to-primitive inversion to update the primitives to agree with the conserved values.

phydro->w is the primitives, consisting of density (IDN), gas pressure (IPR, not internal or any other energy),
and velocity (IVX, IVY, IVZ). phydro->u is the conserved variables, consisting of density 
(IDN, same density as before for Newtonian physics), total kinetic + internal + magnetic energy (IEN), and momentum (IM1, IM2, IM3).
*/

