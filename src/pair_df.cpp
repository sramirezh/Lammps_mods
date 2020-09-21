/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "pair_df.h"                                                                  // yukawa->Own
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Pairdf::Pairdf(LAMMPS *lmp) : Pair(lmp)                             //pairyukawa->pairOwn
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

Pairdf::~Pairdf()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

//    memory->destroy(rad);							// org
    memory->destroy(cut);
    memory->destroy(epsilon);
    memory->destroy(alpha);
    memory->destroy(offset);                                                    //use Uc
  }
}

/* ---------------------------------------------------------------------- */

void Pairdf::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,rinv,screening,forcedf,factor;                               //forceyukwa->force
  double cons1,var1,var2,r1inv,r2inv;
  double prefactor,p1,p2,f1,f2;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;                                       // 1-2,1-3,1-4 prefactor for lj ?
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor = special_lj[sbmask(j)];                                           // input units -> LJ units
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
//      delx = delx - int(delx);
//      dely = dely - int(dely);
//      delz = delz - int(delz);
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
      p1=(1+2.0*alpha[itype][jtype])/(2.0*alpha[itype][jtype]*(cutsq[itype][jtype]-1.0));
      p2= 2.0*alpha[itype][jtype]*cutsq[itype][jtype];
      prefactor=epsilon[itype][jtype]*p2*pow(p1,2*alpha[itype][jtype]+1);
        r = sqrt(rsq);
        r1inv = pow(1.0/r,2);
        r2inv= cutsq[itype][jtype]/rsq;
	var1= r1inv - 1.0 ;
       var2= r2inv - 1.0;
	//cons1 = 2.0*alpha[itype][jtype] - 1.0;
//        screening = exp(-kappa*r);
//        force = a[itype][jtype] * screening * (kappa + rinv);
	  forcedf=2.0*pow(var2,2*alpha[itype][jtype])
               +2.0*alpha[itype][jtype]*pow(var2,2*alpha[itype][jtype]-1)*2.0*cutsq[itype][jtype]*var1;
       
        fpair = prefactor*factor * forcedf * pow(r,-4);                                           

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
            evdwl = prefactor*var1*pow(var2,2.0*alpha[itype][jtype]);
//         evdwl = pow(var,alpha[itype][jtype]) / alpha[itype][jtype] - offset[itype][jtype];                                         // energy except coulobmic
          evdwl *= factor;                                                      //factor
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);                   //updata virial term
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void Pairdf::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

//  memory->create(rad,n+1,"pair:rad");
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(alpha,n+1,n+1,"pair:alpha");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void Pairdf::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");                //2->1

//  kappa = force->numeric(FLERR,arg[0]);
  cut_global = force->numeric(FLERR,arg[0]);                                    //1->0  set cut

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
---------------------------------------------------------------------- --- */

void Pairdf::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);                                   //if d(i) is not same, define two type and sigma,
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);                                   // look lj/cut.cpp

  double epsilon_one = force->numeric(FLERR,arg[2]);
  double alpha_one = force->numeric(FLERR,arg[3]);

  double cut_one = cut_global;
  if (narg == 5) cut_one = force->numeric(FLERR,arg[4]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      alpha[i][j] = alpha_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;                                                        // ?
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double Pairdf::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_distance(epsilon[i][i],epsilon[j][j]);
    alpha[i][j] = mix_distance(alpha[i][i],alpha[j][j]);                              // where?
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);                              // where?
  } 

  if (offset_flag) {                                                            // when use restart
    offset[i][j] = 0.0;
  } else offset[i][j] = 0.0;

  epsilon[j][i] = epsilon[i][j];
  alpha[j][i] = alpha[i][j];
  cut[j][i] = cut[i][j];
  offset[j][i] = offset[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */    // not use now

void Pairdf::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&alpha[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void Pairdf::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&epsilon[i][j],sizeof(double),1,fp);
          fread(&alpha[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&alpha[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void Pairdf::write_restart_settings(FILE *fp)
{
//  fwrite(&kappa,sizeof(double),1,fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void Pairdf::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
//    fread(&kappa,sizeof(double),1,fp);
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
//  MPI_Bcast(&kappa,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void Pairdf::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g\n",i,epsilon[i][i],alpha[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void Pairdf::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g\n",i,j,epsilon[i][j],alpha[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double Pairdf::single(int i, int j, int itype, int jtype, double rsq,
                          double factor_coul, double factor_lj,
                          double &fforce)
{
  double r1inv,r2inv,r,rinv,forcedf,phi;
  double var1,var2, p1,p2,f1,f2,prefactor,cons;

//  r2inv = 1.0/rsq;
 // r = sqrt(rsq);
  //rinv = 1.0/r;
  //var2 = 1.0 - r/sigma[itype][jtype];
 // cons = alpha[itype][jtype] - 1.0;
 // forceOwn = pow(var2,cons);
 // fforce = factor_lj*forceOwn * rinv;                                             //factor_lj ?

//  phi = pow(var2,alpha[itype][jtype]) / alpha[itype][jtype];

      p1=(1+2.0*alpha[itype][jtype])/(2.0*alpha[itype][jtype]*(cutsq[itype][jtype]-1.0));
      p2= 2.0*alpha[itype][jtype]*cutsq[itype][jtype];
      prefactor=epsilon[itype][jtype]*p2*pow(p1,2*alpha[itype][jtype]+1);
        r = sqrt(rsq);
        r1inv = pow(1.0/r,2);
        r2inv= cutsq[itype][jtype]/rsq;
	var1= r1inv - 1.0 ;
       var2= r2inv - 1.0;
	//cons1 = 2.0*alpha[itype][jtype] - 1.0;
//        screening = exp(-kappa*r);
//        force = a[itype][jtype] * screening * (kappa + rinv);
	  forcedf=2.0*pow(var2,2*alpha[itype][jtype])
          + 2.0*alpha[itype][jtype]*pow(var2,2*alpha[itype][jtype]-1)*2.0*cutsq[itype][jtype]*var1;
         
        fforce = prefactor*factor_lj * forcedf * pow(r,-4);                                           


        phi=  prefactor*var1*pow(var2,2.0*alpha[itype][jtype]);                                        
        return factor_lj*phi;                                                        // factor_lj ?
}
