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

/* ----------------------------------------------------------------------
   Contributing Author: Xipeng Wang, Simon Ramirez-Hinestrosa
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_wf.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairWF::PairWF(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairWF::~PairWF()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(nu);
    memory->destroy(mu);
    memory->destroy(nm);
    memory->destroy(e0nm);  //Alpha * epsilon
    memory->destroy(rcmu);
    memory->destroy(sigma_mu);
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairWF::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,factor_lj;
  double r,forcenm,rminv, rcn, rm, rn;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
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
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        r = sqrt(rsq);
        rminv = pow(r2inv,mu[itype][jtype]);
        rm = sigma_mu[itype][jtype]*rminv - 1.0;
        rn = rcmu[itype][jtype]*rminv - 1.0;

        forcenm = 2.0*mu[itype][jtype] *sigma_mu[itype][jtype]*pow(rn,2.0*nu[itype][jtype])
                + 4.0*nm[itype][jtype] *rcmu[itype][jtype]*rm*pow(rn,2.0*nu[itype][jtype]-1.0);
        fpair = factor_lj*e0nm[itype][jtype]*forcenm*pow(r2inv,mu[itype][jtype]+1.0);

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          evdwl = e0nm[itype][jtype] *
            (rm*pow(rn,2.0*nu[itype][jtype])) - offset[itype][jtype];
          evdwl *= factor_lj;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairWF::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(nu,n+1,n+1,"pair:nu");
  memory->create(mu,n+1,n+1,"pair:mu");
  memory->create(nm,n+1,n+1,"pair:nm");
  memory->create(e0nm,n+1,n+1,"pair:e0nm");
  memory->create(rcmu,n+1,n+1,"pair:rcmu");
  memory->create(sigma_mu,n+1,n+1,"pair:sigma_mu");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairWF::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairWF::coeff(int narg, char **arg)
{
  if (narg < 6 || narg > 7)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one = force->numeric(FLERR,arg[2]);
  double sigma_one = force->numeric(FLERR,arg[3]);
  double nu_one = force->numeric(FLERR,arg[4]);
  double mu_one = force->numeric(FLERR,arg[5]);

  double cut_one = cut_global;
  if (narg == 7) cut_one = force->numeric(FLERR,arg[6]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      nu[i][j] = nu_one;
      mu[i][j] = mu_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairWF::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  nm[i][j] = nu[i][j]*mu[i][j];
  e0nm[i][j] = epsilon[i][j]*2.0*nu[i][j]*pow(cut[i][j]/sigma[i][j],2.0*mu[i][j])
                       *pow((1+2.0*nu[i][j])/(2.0*nu[i][j])/(pow(cut[i][j]/sigma[i][j],2.0*mu[i][j])-1.0),
                              2.0*nu[i][j]+1.0);
  rcmu[i][j] = pow(cut[i][j],2.0*mu[i][j]);
  sigma_mu[i][j] = pow(sigma[i][j], 2.0*mu[i][j]);

  if (offset_flag && (cut[i][j] > 0.0)) {
    offset[i][j] = 0.0;
  } else offset[i][j] = 0.0;

  epsilon[j][i] = epsilon[i][j];
  nu[j][i] = nu[i][j];
  mu[j][i] = mu[i][j];
  nm[j][i] = nm[i][j];
  sigma[j][i] = sigma[i][j];
  e0nm[j][i] = e0nm[i][j];
  rcmu[j][i] = rcmu[i][j];
  sigma_mu[j][i] = sigma_mu[i][j];
  offset[j][i] = offset[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairWF::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&nu[i][j],sizeof(double),1,fp);
        fwrite(&mu[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairWF::read_restart(FILE *fp)
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
          fread(&sigma[i][j],sizeof(double),1,fp);
          fread(&nu[i][j],sizeof(double),1,fp);
          fread(&mu[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&nu[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&mu[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairWF::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
//  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairWF::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
//    fread(&tail_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
//  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairWF::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g\n",i,epsilon[i][i],sigma[i][i],nu[i][i],mu[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairWF::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g\n",i,j,
              epsilon[i][j],sigma[i][j],nu[i][j],mu[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairWF::single(int /*i*/, int /*j*/, int itype, int jtype,
                      double rsq, double /*factor_coul*/, double factor_lj,
                      double &fforce)
{
  double r2inv,r, rminv, rm, rn, forcenm,phinm;

  r2inv = 1.0/rsq;
  r = sqrt(rsq);
   
  r2inv = 1.0/rsq;
  r = sqrt(rsq);
  rminv = pow(r2inv,mu[itype][jtype]);
  rm = sigma_mu[itype][jtype]*rminv - 1.0;
  rn = rcmu[itype][jtype]*rminv - 1.0;
  forcenm = 2.0*mu[itype][jtype] *sigma_mu[itype][jtype]*pow(rn,2.0*nu[itype][jtype])
                + 4.0*nm[itype][jtype] *rcmu[itype][jtype]*rm*pow(rn,2.0*nu[itype][jtype]-1.0);
  fforce = factor_lj*e0nm[itype][jtype]*forcenm*pow(r2inv,mu[itype][jtype]+1.0);

  phinm = e0nm[itype][jtype] * rm*pow(rn,2.0*nu[itype][jtype]) -
    offset[itype][jtype];
  return factor_lj*phinm;
}

/* ---------------------------------------------------------------------- */

void *PairWF::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"epsilon") == 0) return (void *) epsilon;
  if (strcmp(str,"sigma") == 0) return (void *) sigma;
  if (strcmp(str,"nu") == 0) return (void *) nu;
  if (strcmp(str,"mu") == 0) return (void *) mu;
  return NULL;
}
