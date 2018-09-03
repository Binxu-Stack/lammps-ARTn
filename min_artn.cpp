/* -------------------------------------------------------------------------------------------------
 * ARTn: Activation Relaxation Technique nouveau
 * Bin Xu, xubinrun@gmail.com; Lingti Kong, konglt@gmail.com
------------------------------------------------------------------------------------------------- */

#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix_minimize.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "math.h"
#include "memory.h"
#include "min_artn.h"
#include "modify.h"
#include "output.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "timer.h"
#include "update.h"
#include "dump_atom.h"
#define MAXLINE 512
//#define DEBUG

#define ZERO  1.e-10

using namespace LAMMPS_NS;



/* -------------------------------------------------------------------------------------------------
 * lapack or MKL-lapack is used to evaluate the lowest eigenvalue of the matrix in Lanczos.
------------------------------------------------------------------------------------------------- */
#ifdef MKL
#include "mkl.h"
#define dstev_  dstev
#else
extern "C" {
extern void dstev_(char *, int*, double *, double *, double *, int *, double *, int *);
};
#endif

#define EPS_ENERGY 1.e-8

enum{MAXITER,MAXEVAL,ETOL,FTOL,DOWNHILL,ZEROALPHA,ZEROFORCE,ZEROQUAD};



/* -------------------------------------------------------------------------------------------------
 * Constructor of ARTn
------------------------------------------------------------------------------------------------- */
MinARTn::MinARTn(LAMMPS *lmp): MinLineSearch(lmp)
{
  random = NULL;
  pressure = NULL;
  dumpmin = dumpsad = dumpevent =  NULL;
  dumpmin_outside = dumpsad_outside = false;
  egvec = x0tmp = x00 = fperp = NULL;

  fp1 = fp2 = fp_sadlpress =  NULL;
  glist = NULL;
  groupname = flog = fevent = fconfg = c_fsadpress = fdump_direction = NULL;

  char *id_press = new char [13];
  strcpy(id_press,"thermo_press");
  int pressure_compute = modify->find_compute(id_press); delete [] id_press;

  if (pressure_compute < 0) error->all(FLERR,"Could not find compute pressure ID");

  pressure = modify->compute[pressure_compute];

  // set default control parameters
  set_defaults();

  MPI_Comm_rank(world, &me);
  MPI_Comm_size(world, &np);
return;
}

/* -------------------------------------------------------------------------------------------------
 * The main loops of ARTn
------------------------------------------------------------------------------------------------- */
int MinARTn::iterate(int maxevent)
{
  // read in control parameters
  read_control();
  artn_init();
  max_conv_steps = update->nsteps;

  // minimize before searching saddle points.
  if (me == 0) print_info(0);

  if(!min_fire)stop_condition = min_converge(max_conv_steps,0);
  else stop_condition = min_converge_fire(max_conv_steps); evalf += neval;
  eref = ecurrent;
  stopstr = stopstrings(stop_condition);
  ddum = sqrt(fnorm_sqr());

  if (me == 0) print_info(1);
  if (flag_press){
    ++update->ntimestep;
    double *vvec = atom->v[0];
    for (int i = 0; i < nvec; ++i) vvec[i] = 0.;
    pressure->addstep(update->ntimestep);
    energy_force(0); ++evalf; reset_x00();
    pressure->compute_vector();
    double * press = pressure->vector;
    --update->ntimestep;
    if (me == 0 && fp1){
      fprintf(fp1, "  - Pressure tensor           :");
      for (int ii = 0; ii < 6; ++ii) fprintf(fp1, " %g", press[ii]);
      fprintf(fp1, "\n");
    }
  }


  // dump the first stable configuration
  if (dumpmin){
    int idum = update->ntimestep;
    update->ntimestep = ref_id;
    dumpmin->write();
    update->ntimestep = idum;
  }
  if (dumpevent){
    int idum = update->ntimestep;
    update->ntimestep = sad_id*1000;
    dumpevent->write();
    update->ntimestep = idum;
  }

  // print header of event log file
  if (me == 0 && fp2){
    if (!flag_push_over){
      fprintf(fp2, "#  1       2        3       4       5     6      7      8\n");
      fprintf(fp2, "#Event   del-E   egv-sad   nsadl sad-dx sad-dy sad-dz sad-dr\n");
      fprintf(fp2, "#----------------------------------------------------------\n");

    } else {if (flag_press) print_info(2);
    else print_info(3);
    }
  }


  // main loop of ARTn
  int ievent = 0;

  // define the variable needed for events_per_atom
  int events_iterator = 0;
  iatom = 0;
  if (events_per_atom != 0) max_num_events = events_per_atom * ngroup;

  while (ievent < max_num_events){
    // activation
    stage = 1; ++nattempt;
    while ((flag_test?new_find_saddle():find_saddle()) == 0) {stage = 1; ++nattempt;}

    // confirm saddle
    ++sad_found;
    if (check_sad2min()) continue;

    if (flag_push_back) if (push_back_sad() == 0) continue;
    ++sad_id;

    if (flag_relax_sad) sad_converge(max_conv_steps);

    ++ievent; ++stage;
    if (me == 0 && fp2) fprintf(fp2, "%5d %9.6f %7.3f", ievent, delE, egval);

    if (dumpsad){
      int idum = update->ntimestep;
      update->ntimestep = sad_id;
      //domain->pbc();
      if(sad_id % dump_sad_every == 0) dumpsad->write();
      update->ntimestep = idum;
    }
    // analysis saddle point
    analysis_saddle();

    // Relaxation & Metropolis
    if(flag_push_over){
      push_down();
      metropolis();
    }else{
      if (me == 0 && fp2) fprintf(fp2, "\n");
      for (int i = 0; i < nvec; ++i) xvec[i] = x00[i];
      ecurrent = energy_force(1); ++evalf;
      artn_reset_vec();
    }
    fflush(fp2);
    if (++events_iterator >= events_per_atom){ ++iatom;events_iterator = 0;}
  }

  // finalize ARTn
  artn_final();

return MAXITER;
}

/* -------------------------------------------------------------------------------------------------
 * return 0 if distance between new saddle and original min is greater than initial kick & delE > 0.
------------------------------------------------------------------------------------------------- */
int MinARTn::check_sad2min()
{
  lanczos(flag_egvec, 1, num_lancz_vec_c);

  reset_x00();
  // check current center-of-mass
  group->xcm(groupall, masstot, com);
  double dxcm[3];
  dxcm[0] = com[0] - com0[0];
  dxcm[1] = com[1] - com0[1];
  dxcm[2] = com[2] - com0[2];

  // check the distance between new saddle and original min
  // fperp now stores the displacement vector, will be used by push_back and push_down
  double tmp_me[2], tmp_all[2];
  tmp_me[0] = tmp_me[1] = 0.;

  for (int i = 0; i < nvec; ++i){
    double dx = xvec[i] - x00[i] - dxcm[i%3];
    tmp_me[0] += dx * dx;
    tmp_me[1] += dx * egvec[i];
  }
  MPI_Allreduce(tmp_me, tmp_all, 2, MPI_DOUBLE, MPI_SUM, world);
  ddum = sqrt(tmp_all[0]);

  int status = 0;
  if (ddum >= disp_sad2min_thr && delE > 0.){
    if (me == 0) print_info(20);

    if (tmp_all[1] > 0.) for (int i = 0; i < nvec; ++i){ 
      fperp[i] = push_over_saddle > 0 ? egvec[i] * push_over_saddle : egvec[i] * (-push_over_saddle * ddum);
    }
    else for (int i = 0; i < nvec; ++i) fperp[i] = -egvec[i] * push_over_saddle;

  } else {
    if (me == 0) print_info(21);
    
    for (int i = 0; i < nvec; ++i) xvec[i] = x00[i];
    status = 1;
  }

return status;
}

/* -------------------------------------------------------------------------------------------------
 * return 1, saddle connect with starting minimum, else return 0
------------------------------------------------------------------------------------------------- */
int MinARTn::push_back_sad()
{
  ++stage;
  for (int i = 0; i < nvec; ++i) x0tmp[i] = xvec[i];  // x0tmp stores the saddle position

  // push back the saddle; fperp carries the direction vector, set by check_sad2min
  for (int i = 0; i < nvec; ++i) xvec[i] -= fperp[i];

  if (me == 0) print_info(30);

  // minimization using CG
  neval = 0;
  if(!min_fire){
    stop_condition = min_converge(max_conv_steps,2); 
  }else{
    stop_condition = min_converge_fire(max_conv_steps);
  }
  evalf += neval;
  stopstr = stopstrings(stop_condition); artn_reset_vec();
  reset_coords();
  ddum = sqrt(fnorm_sqr());

  // output minimization information
  if (me == 0) print_info(31);
  ddum = fabs(ecurrent - eref);

  if (ddum > max_ener_tol) {
    if (me == 0) print_info(32);

    for (int i = 0; i < nvec; ++i) xvec[i] = x00[i];
    return 0;
  }

  // check current center-of-mass
  double dxcm[3];
  group->xcm(groupall, masstot, com);
  dxcm[0] = com[0] - com0[0];
  dxcm[1] = com[1] - com0[1];
  dxcm[2] = com[2] - com0[2];

  double dr = 0., drall;
  double drmaxone = 0., drmaxall, dx, dy, dz;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; ++i){
    dx = xvec[i*3] - x00[i*3] - dxcm[0];
    dy = xvec[i*3+1] - x00[i*3+1] - dxcm[1];
    dz = xvec[i*3+2] - x00[i*3+2] - dxcm[2];
    // get one atom displacement
    double tmp = dx*dx + dy*dy + dz*dz;
    if (tmp > drmaxone) drmaxone = tmp;
    dr += tmp;
  }
  MPI_Allreduce(&dr,&drall,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&drmaxone,&drmaxall,1,MPI_DOUBLE,MPI_MAX,world);

  ddum = sqrt(drall);
  tmpddum = drmaxall = sqrt(drmaxall);

  if (max_disp_tol > 0){
    if (ddum < max_disp_tol ) {
      if (me == 0) print_info(33);
      for (int i = 0; i < nvec; ++i) xvec[i] = x0tmp[i];
      return 1;
    } else {
      if (me == 0) print_info(34);
      for (int i = 0; i < nvec; ++i) xvec[i] = x00[i];
    }
  } else {
    if (drmaxall < fabs(max_disp_tol)) {
      if (me == 0) print_info(33);
      for (int i = 0; i < nvec; ++i) xvec[i] = x0tmp[i];
      return 1;
    }else {
      if (me == 0) print_info(34);
      for (int i = 0; i < nvec; ++i) xvec[i] = x00[i];
    }
  }


return 0;
}

/* -------------------------------------------------------------------------------------------------
 * analysis saddle point
------------------------------------------------------------------------------------------------- */
void MinARTn::analysis_saddle(){
  // calculate the moved atoms at saddle point
  group->xcm(groupall, masstot, com);
  double **x = atom->x;
  double dx, dy, dz;
  double dxcm[3];
  dxcm[0] = com[0] - com0[0];
  dxcm[1] = com[1] - com0[1];
  dxcm[2] = com[2] - com0[2];

  int n_moved = 0, n_movedall, n = 0;
  double tmp, disp_thr2 = atom_disp_thr*atom_disp_thr;
  double temp[3]={0}, tempall[3]={0};
  for (int i = 0; i < atom->nlocal; ++i) {
    dx = x[i][0] - x00[n]   - dxcm[0];
    dy = x[i][1] - x00[n+1] - dxcm[1];
    dz = x[i][2] - x00[n+2] - dxcm[2];

    temp[0] += dx*dx;
    temp[1] += dy*dy;
    temp[2] += dz*dz;
    tmp = dx*dx + dy*dy + dz*dz;
    n += 3;

    if (tmp > disp_thr2) ++n_moved;
  }
  MPI_Reduce(&n_moved, &n_movedall, 1, MPI_INT, MPI_SUM, 0, world);
  MPI_Reduce(temp,tempall,3,MPI_DOUBLE,MPI_SUM,0,world);
  if (me == 0 && fp2) fprintf(fp2, " %5d %7.2f %7.2f %7.2f %7.2f", n_movedall, sqrt(tempall[0]),sqrt(tempall[1])
      ,sqrt(tempall[2]),sqrt(tempall[0]+tempall[1]+tempall[2]));
  if (flag_sadl_press){
    double *vvec = atom->v[0];
    for (int i = 0; i < nvec; ++i) vvec[i] = 0.;

    ++update->ntimestep;
    pressure->addstep(update->ntimestep);
    energy_force(0); ++evalf; reset_x00();
    pressure->compute_vector();
    double * press = pressure->vector;
    --update->ntimestep;


    if (me == 0 && fp_sadlpress){
      fprintf(fp_sadlpress, "%i",sad_id);
      for (int i = 0; i < 6; ++i) fprintf(fp_sadlpress, " %10g", press[i]);
      fprintf(fp_sadlpress, "\n");
    }
    fflush(fp_sadlpress);
  }

}

/* -------------------------------------------------------------------------------------------------
 * push the configuration push_down
------------------------------------------------------------------------------------------------- */
void MinARTn::push_down()
{
  // push down the saddle; fperp carries the direction vector, set by check_sad2min

  for (int i = 0; i < nvec; ++i) xvec[i] += fperp[i];

  if (dumpevent){
    int idum = update->ntimestep;
    update->ntimestep = 1000*sad_id+500;
    dumpevent->write();
    update->ntimestep = idum;
  }

  ecurrent = energy_force(1); ++evalf;
  if (me == 0) print_info(50);

  // minimization using CG or FIRE
  if(!min_fire){SD_min_converge(SD_steps,1); evalf += neval;
  stop_condition = min_converge(max_conv_steps,1); evalf += neval;
  }else{stop_condition = min_converge_fire(max_conv_steps); evalf += neval;}
  stopstr = stopstrings(stop_condition);
  artn_reset_vec();
  reset_x00();
  ddum = sqrt(fnorm_sqr());

  // output minimization information
  if (me == 0) print_info(51);

  // store min configuration
  ++min_id;
  if (dumpmin){
    int idum = update->ntimestep;
    update->ntimestep = min_id;
    //domain->pbc();
    if(min_id % dump_min_every == 0) dumpmin->write();
    update->ntimestep = idum;
  }

  return;
}

/* -------------------------------------------------------------------------------------------------
 * decide whether reject or accept the new configuration
------------------------------------------------------------------------------------------------- */
void MinARTn::metropolis()
{
  // output reference energy ,current energy and pressure.
  if (me == 0 && fp2) fprintf(fp2, " %5d %4d %5d %7d %10.3f %10.3f", ref_id, sad_id, min_id, that, eref, ecurrent);

  double dr = 0., drall;
  double **x = atom->x;
  double dx,dy,dz;
  double disp_me[2], disp_all[2];
  disp_all[0] = disp_all[1] = 0.;

  // check current center-of-mass
  group->xcm(groupall, masstot, com);
  double dxcm[3];
  dxcm[0] = com[0] - com0[0];
  dxcm[1] = com[1] - com0[1];
  dxcm[2] = com[2] - com0[2];

  // calculate displacment w.r.t. ref
  double tmp, disp_thr2 = atom_disp_thr*atom_disp_thr;
  int n_moved = 0, n_movedall, n = 0;
  for (int i = 0; i < atom->nlocal; ++i) {
    dx = x[i][0] - x00[n]   - dxcm[0];
    dy = x[i][1] - x00[n+1] - dxcm[1];
    dz = x[i][2] - x00[n+2] - dxcm[2];

    tmp = dx*dx + dy*dy + dz*dz;
    dr += tmp; n += 3;

    if (tmp > disp_thr2) ++n_moved;
  }
  disp_me[0] = dr; disp_me[1] = double(n_moved);
  MPI_Reduce(disp_me, disp_all, 2, MPI_DOUBLE, MPI_SUM, 0, world);
  drall = sqrt(disp_all[0]); n_movedall = int(disp_all[1]);

  if (me == 0 && fp2) fprintf(fp2, " %5d", n_movedall);

  // set v = 0 to calculate pressure
  if (flag_press){
    double *vvec = atom->v[0];
    for (int i = 0; i < nvec; ++i) vvec[i] = 0.;

    ++update->ntimestep;
    pressure->addstep(update->ntimestep);
    energy_force(0); ++evalf; reset_x00();
    pressure->compute_vector();
    --update->ntimestep;
    double * press = pressure->vector;

    if (me == 0 && fp2) for (int i = 0; i < 6; ++i) fprintf(fp2, " %10g", press[i]);
  }

  // Metropolis
  int acc = 0;
  if (me == 0){
    ddum = drall;
    print_info(60);

    // temporary modification
    if (temperature > 0. && (ecurrent < eref || random->uniform() < exp((eref - ecurrent)/temperature)) && ecurrent - eref > -5) acc = 1;
  }
  MPI_Bcast(&acc, 1, MPI_INT, 0, world);

  if (acc){
    if (me == 0) print_info(61);
    ref_id = min_id; eref = ecurrent;

  } else {
    if (me == 0) print_info(62);

    for (int i = 0; i < nvec; ++i) xvec[i] = x00[i];
    ecurrent = energy_force(1); ++evalf;
    artn_reset_vec();
  }

  if (me == 0 && fp2){
    for (int i = 0; i < 3; ++i) dxcm[i] *= double(atom->natoms);
    fprintf(fp2, " %12.5f %2d %8.5f\n", ecurrent, acc, drall);
    fflush(fp2);
  }

return;
}

/* -------------------------------------------------------------------------------------------------
 * setup default parameters, some values come from Norman's code
------------------------------------------------------------------------------------------------- */
void MinARTn::set_defaults()
{
  // global
  seed = 12345;
  nattempt = ref_id = min_id = sad_id = 0;
  sad_found = 0;
  max_num_events   = 1000;
  flag_press       = 0;
  flag_sadl_press  = 0;
  min_fire         = 0;
  events_per_atom  = 0;

  // activation, harmonic well escape
  cluster_radius   = 5.0;
  fire_lanczos_every = 5;
  init_step_size   = 0.1;
  basin_factor     = 2.5;
  max_perp_move_h  = 20;
  max_iter_basin   = 30;
  min_num_ksteps   = 0;		
  increment_size   = 0.09;
  force_th_perp_h  = 0.5;
  eigen_th_well    = -0.01;
  flag_dump_direction = 0;
  flag_deformation_gradient = 0;
  flag_random_deformation_gradient = 0;
  for (int i = 0; i < 9; ++i) deformation_gradient[i] = 0.0;
  random_kick_factor = 0.0;

  // activation, converge to saddle
  max_activat_iter  = 100;
  flag_test         = 0;
  para_factor       = 1.;
  fire_output_every = 10;
  use_fire          = 0;
  force_th_saddle   = 0.005;
  eigen_th_fail     = 0.000;
  conv_perp_inc     = 40;
  max_perp_moves_c  = 15;
  force_th_perp_sad = 0.001;

  // confirmatom of new saddle
  disp_sad2min_thr = -1.;
  flag_push_back   = 0;
  max_disp_tol     = 0.2;
  max_ener_tol     = 0.001;
  flag_relax_sad   = 0;

  // convergence to new minimum
  flag_push_over   = 1;
  push_over_saddle = 0.2;
  atom_disp_thr    = 0.2;
  temperature      = 0.1;
  SD_steps         = 5;

  // for lanczos
  num_lancz_vec_h  = 30;
  num_lancz_vec_c  = 20;
  del_disp_lancz   = 0.001;
  eigen_th_lancz   = 0.01;

  // output
  log_level        = 1;
  print_freq       = 1;
  dump_min_every   = 1;
  dump_sad_every   = 1;
  dump_event_every = 0;

return;
}

/* -------------------------------------------------------------------------------------------------
 * read ARTn control parameters from file "artn.control"
------------------------------------------------------------------------------------------------- */
void MinARTn::read_control()
{
  char oneline[MAXLINE], str[MAXLINE], *token1, *token2;
  FILE *fp = fopen("artn.control", "r");
  char *fmin, *fsad, *fproc; fmin = fsad = fproc = NULL;
  if (fp == NULL){
    error->warning(FLERR, "Cannot open ARTn control parameter file. Default parameters will be used.\n");

  } else {
    while ( 1 ) {
      fgets(oneline, MAXLINE, fp);
      if (feof(fp)) break;

      if (token1 = strchr(oneline,'#')) *token1 = '\0';

      token1 = strtok(oneline," \t\n\r\f");
      if (token1 == NULL) continue;
      token2 = strtok(NULL," \t\n\r\f");
      if (token2 == NULL){
	sprintf(str, "Insufficient parameter for %s of ARTn!", token1);
	error->all(FLERR, str);
      }

      if (strcmp(token1, "random_seed") == 0){
	seed = force->inumeric(FLERR, token2);
	if (seed < 1) error->all(FLERR, "ARTn: seed must be greater than 0");

      } else if (strcmp(token1, "temperature") == 0){
	temperature = force->numeric(FLERR, token2);

      } else if (strcmp(token1, "max_num_events") == 0){
	max_num_events = force->inumeric(FLERR, token2);
	if (max_num_events < 1) error->all(FLERR, "ARTn: max_num_events must be greater than 0");

      } else if (strcmp(token1, "max_activat_iter") == 0){
	max_activat_iter = force->inumeric(FLERR, token2);
	if (max_activat_iter < 1) error->all(FLERR, "ARTn: max_activat_iter must be greater than 0");

      } else if (strcmp(token1, "increment_size") == 0){
	increment_size = force->numeric(FLERR, token2);
	if (increment_size <= 0.) error->all(FLERR, "ARTn: increment_size must be greater than 0.");

      } else if (!strcmp(token1, "cluster_radius")){
	cluster_radius = force->numeric(FLERR, token2);

      } else if (strcmp(token1, "group_4_activat") == 0){
	if (groupname) delete [] groupname;
	groupname = new char [strlen(token2)+1];
	strcpy(groupname, token2);

      } else if (strcmp(token1, "init_step_size") == 0){
	init_step_size = force->numeric(FLERR, token2);
	if (init_step_size <= 0.) error->all(FLERR, "ARTn: init_step_size must be greater than 0.");

      } else if (strcmp(token1, "basin_factor") == 0){
	basin_factor = force->numeric(FLERR, token2);
	if (basin_factor <= 0.) error->all(FLERR, "ARTn: basin_factor must be greater than 0.");

      } else if (strcmp(token1, "max_perp_move_h") == 0){
	max_perp_move_h = force->inumeric(FLERR, token2);
	if (max_perp_move_h < 1) error->all(FLERR, "ARTn: max_perp_move_h must be greater than 0.");

      } else if (strcmp(token1, "min_num_ksteps") == 0){
	min_num_ksteps = force->inumeric(FLERR, token2);
	if (min_num_ksteps < 1) error->all(FLERR, "ARTn: min_num_ksteps must be greater than 0");

      } else if (strcmp(token1, "eigen_th_well") == 0){
	eigen_th_well = force->numeric(FLERR, token2);
	if (eigen_th_well > 0.) error->all(FLERR, "ARTn: eigen_th_well must be less than 0.");

      } else if (strcmp(token1, "max_iter_basin") == 0){
	max_iter_basin = force->inumeric(FLERR, token2);
	if (max_iter_basin < 1) error->all(FLERR, "ARTn: max_iter_basin must be greater than 0");

      } else if (strcmp(token1, "force_th_perp_h") == 0){
	force_th_perp_h = force->numeric(FLERR, token2);
	if (force_th_perp_h <= 0.) error->all(FLERR, "ARTn: force_th_perp_h must be greater than 0.");

      } else if (strcmp(token1, "num_lancz_vec_h") == 0){
	num_lancz_vec_h = force->inumeric(FLERR, token2);
	if (num_lancz_vec_h < 1) error->all(FLERR, "ARTn: num_lancz_vec_h must be greater than 0");

      } else if (strcmp(token1, "num_lancz_vec_c") == 0){
	num_lancz_vec_c = force->inumeric(FLERR, token2);
	if (num_lancz_vec_c < 1) error->all(FLERR, "ARTn: num_lancz_vec_c must be greater than 0");

      } else if (strcmp(token1, "del_disp_lancz") == 0){
	del_disp_lancz = force->numeric(FLERR, token2);
	if (del_disp_lancz  <=  0.) error->all(FLERR, "ARTn: del_disp_lancz must be greater than 0.");

      } else if (strcmp(token1, "eigen_th_lancz") == 0){
	eigen_th_lancz = force->numeric(FLERR, token2);
	if (eigen_th_lancz <=  0.) error->all(FLERR, "ARTn: eigen_th_lancz must be greater than 0.");

      } else if (strcmp(token1, "force_th_saddle") == 0){
	force_th_saddle = force->numeric(FLERR, token2);
	if (force_th_saddle <=  0.) error->all(FLERR, "ARTn: force_th_saddle must be greater than 0.");

      } else if (strcmp(token1, "disp_sad2min_thr") == 0){
	disp_sad2min_thr = force->numeric(FLERR, token2);
	if (disp_sad2min_thr <=  0.) error->all(FLERR, "ARTn: disp_sad2min_thr must be greater than 0.");

      } else if (strcmp(token1, "push_over_saddle") == 0){
	push_over_saddle = force->numeric(FLERR, token2);
	//if (push_over_saddle <=  0.) error->all(FLERR, "ARTn: push_over_saddle must be greater than 0.");

      } else if (strcmp(token1, "eigen_th_fail") == 0){
	eigen_th_fail = force->numeric(FLERR, token2);

      } else if (!strcmp(token1, "atom_disp_thr")){
	atom_disp_thr = force->numeric(FLERR, token2);
	if (atom_disp_thr <= 0.) error->all(FLERR, "ARTn: atom_disp_thr must be greater than 0.");

      } else if (strcmp(token1, "max_perp_moves_c") == 0){
	max_perp_moves_c = force->inumeric(FLERR, token2);
	if (max_perp_moves_c < 1) error->all(FLERR, "ARTn: max_perp_moves_c must be greater than 0.");

      } else if (strcmp(token1, "force_th_perp_sad") == 0){
	force_th_perp_sad = force->numeric(FLERR, token2);
	if (force_th_perp_sad <= 0.) error->all(FLERR, "ARTn: force_th_perp_sad must be greater than 0.");

      } else if (strcmp(token1, "use_fire") == 0){
	use_fire = force->inumeric(FLERR, token2);

      } else if (!strcmp(token1, "flag_push_back")){
	flag_push_back = force->inumeric(FLERR, token2);

      } else if (!strcmp(token1, "flag_relax_sad")){
	flag_relax_sad = force->inumeric(FLERR, token2);

      } else if (!strcmp(token1, "max_disp_tol")){
	max_disp_tol = force->numeric(FLERR, token2);

      } else if (!strcmp(token1, "max_ener_tol")){
	max_ener_tol = force->numeric(FLERR, token2);
	if (max_ener_tol <= 0.) error->all(FLERR, "ARTn: max_ener_tol must be greater than 0.");

      } else if (!strcmp(token1, "flag_press")){
	flag_press = force->inumeric(FLERR, token2);

      } else if (!strcmp(token1, "log_file")){
	if (flog) delete []flog;
	flog = new char [strlen(token2)+1];
	strcpy(flog, token2);

      } else if (!strcmp(token1, "log_level")){
	log_level = force->inumeric(FLERR, token2);

      } else if (!strcmp(token1, "print_freq")){
	print_freq = force->inumeric(FLERR, token2);

      } else if (strcmp(token1, "event_list_file") == 0){
	if (fevent) delete [] fevent;
	fevent = new char [strlen(token2)+1];
	strcpy(fevent, token2);

      } else if (strcmp(token1, "init_config_id") == 0){
	ref_id = force->inumeric(FLERR, token2);

      } else if (strcmp(token1, "dump_min_config") == 0){
	if (fmin) delete []fmin;
	fmin = new char [strlen(token2)+1];
	strcpy(fmin, token2);

      } else if (strcmp(token1, "dump_sad_config") == 0){
	if (fsad) delete []fsad;
	fsad = new char [strlen(token2)+1];
	strcpy(fsad, token2);

      } else if (strcmp(token1, "dump_event_config") == 0){
	if (fproc) delete []fproc;
	fproc = new char [strlen(token2)+1];
	strcpy(fproc, token2);

      } else if (strcmp(token1, "conv_perp_inc") == 0){
	conv_perp_inc = force->inumeric(FLERR, token2);

      } else if (strcmp(token1, "SD_steps") == 0){
	SD_steps = force->inumeric(FLERR, token2);

      } else if (strcmp(token1, "dump_min_every") == 0){
	dump_min_every = force->inumeric(FLERR, token2);

      } else if (strcmp(token1, "dump_sad_every") == 0){
	dump_sad_every = force->inumeric(FLERR, token2);

      } else if (strcmp(token1, "dump_event_every") == 0){
	dump_event_every = force->inumeric(FLERR, token2);

      } else if (strcmp(token1, "min_fire") == 0){
	min_fire = force->inumeric(FLERR, token2);

      } else if (strcmp(token1, "flag_push_over") == 0){
	flag_push_over = force->inumeric(FLERR, token2);

      } else if (strcmp(token1, "para_factor") == 0){
	para_factor = force->numeric(FLERR, token2);

      } else if (strcmp(token1, "fire_lanczos_every") == 0){
	fire_lanczos_every = force->inumeric(FLERR, token2);

      } else if (strcmp(token1, "fire_output_every") == 0){
	fire_output_every = force->inumeric(FLERR, token2);

      } else if (strcmp(token1, "flag_test") == 0){
	flag_test = force->inumeric(FLERR, token2);

      } else if (strcmp(token1, "events_per_atom") == 0){
	events_per_atom = force->inumeric(FLERR, token2);

      } else if (strcmp(token1, "flag_sadl_press") == 0){
	flag_sadl_press = force->inumeric(FLERR, token2);

      } else if (strcmp(token1, "sadl_press_file") == 0){
	if (c_fsadpress) delete []c_fsadpress;
	c_fsadpress = new char [strlen(token2)+1];
	strcpy(c_fsadpress, token2);

      } else if (strcmp(token1, "fdump_direction") == 0){
	if (fdump_direction) delete []fdump_direction;
	fdump_direction = new char [strlen(token2)+1];
	strcpy(fdump_direction, token2);

      } else if (strcmp(token1, "flag_dump_direction") == 0){
	flag_dump_direction = force->inumeric(FLERR, token2);

	// here to parase this command just for historic reason.
      } else if (strcmp(token1, "dump_direction_random_factor") == 0){
	random_kick_factor = force->numeric(FLERR, token2);

      } else if (strcmp(token1, "random_kick_factor") == 0){
	random_kick_factor = force->numeric(FLERR, token2);

      } else if (strcmp(token1, "deformation_gradient") == 0){
	flag_deformation_gradient = 1;
	if (strcmp(token2, "random") == 0){
	  flag_random_deformation_gradient = 1;
	}else {
	  for (int i =0; i < 9; ++i){
	    deformation_gradient[i] = force->numeric(FLERR, token2);
	    token2 = strtok(NULL," \t\n\r\f");
	  }
	}

      } else {
	sprintf(str, "Unknown control parameter for ARTn: %s", token1);
	error->all(FLERR, str);
      } 

    }
    fclose(fp);
  }

  if (flag_deformation_gradient && cluster_radius <= ZERO ) {
    sprintf(str, "Cluster_radius should be a positive value when using deformation gradient.");
    error->all(FLERR, str);
  }

  if (flag_deformation_gradient && flag_dump_direction){
    sprintf(str, "Deformation gradient and dump_direction should not be used at the same time.");
    error->all(FLERR, str);
  }

  // set default output file names
  if (flog == NULL){
    flog = new char [9];
    strcpy(flog, "log.artn");
  }
  if (fevent == NULL){
    fevent = new char [10];
    strcpy(fevent, "log.event");
  }
  if (fmin == NULL){
    fmin = new char [14];
    strcpy(fmin, "min.lammpstrj");
  }
  if (fsad == NULL){
    fsad = new char [14];
    strcpy(fsad, "sad.lammpstrj");
  }
  if (fproc == NULL){
    fproc = new char [14];
    strcpy(fproc, "eve.lammpstrj");
  }
  if (c_fsadpress == NULL){
    c_fsadpress = new char [15];
    strcpy(c_fsadpress, "sadl_press.dat");
  }
  min_id = ref_0 = ref_id;

  // if disp_sad2min_thr not set, set as twice init_step_size
  //if (disp_sad2min_thr <= 0.) disp_sad2min_thr = init_step_size + init_step_size;

  // default group name is all
  if (groupname == NULL) {groupname = new char [4]; strcpy(groupname, "all");}

  int igroup = group->find(groupname);

  if (igroup == -1){
    sprintf(str, "Cannot find ARTn group: %s", groupname);
    error->all(FLERR, str);
  }
  groupbit = group->bitmask[igroup];
  ngroup = group->count(igroup);
  if (ngroup < 1) error->all(FLERR, "No atom is found in your desired group for activation!");

  igroup = group->find("FREEZE");
  if (igroup == -1) {
    flag_freeze = 0;
  }else {
    flag_freeze = 1;
    freeze_group_bit = group->bitmask[igroup];
  }

  // group info for all
  groupall = group->find("all");
  masstot = group->mass(groupall);

  // open log file and output control parameter info
  if (me == 0 && flag_sadl_press && strcmp(c_fsadpress, "NULL") != 0){
    fp_sadlpress = fopen(c_fsadpress, "w");
    if(fp_sadlpress == NULL){
      sprintf(str, "Cannot open ARTn sadlpress file: %s for writing", c_fsadpress);
      error->one(FLERR,str);
    }
  }
  if (me == 0 && strcmp(flog, "NULL") != 0){
    fp1 = fopen(flog, "w");
    if (fp1 == NULL){
      sprintf(str, "Cannot open ARTn log file: %s for writing", flog);
      error->one(FLERR, str);
    }

    fprintf(fp1, "\n#===================================== ARTn based on LAMMPS ========================================\n");
    fprintf(fp1, "# global control parameters\n");
    fprintf(fp1, "max_num_events      %-18d  # %s\n", max_num_events,"Max number of events");
    fprintf(fp1, "min_fire            %-18d  # %s\n", min_fire, "use FIRE to do minimization both in push back & push forward");
    fprintf(fp1, "flag_press          %-18d  # %s\n", flag_press, "Flag whether the pressure info will be monitored");
    fprintf(fp1, "flag_sadl_press     %-18d  # %s\n", flag_press, "Flag whether the sadlle point pressure info will be monitored");
    fprintf(fp1, "random_seed         %-18d  # %s\n", seed, "Seed for random generator");
    fprintf(fp1, "init_config_id      %-18d  # %s\n", min_id, "ID of the initial stable configuration");
    fprintf(fp1, "flag_push_over      %-18d  # %s\n", flag_push_over, "Flag whether to push over saddle to find another minimum");
    fprintf(fp1, "events_per_atom     %-18d  # %s\n", events_per_atom, "Find designed events per atom, set to 0 to shutoff this method");
    fprintf(fp1, "\n# activation, harmonic well escape\n");
    fprintf(fp1, "group_4_activat     %-18s  # %s\n", groupname, "The lammps group ID of the atoms that can be activated");
    fprintf(fp1, "flag_dump_direction %-18d  # %s\n", flag_dump_direction, "Use dump direction file as the initial kick direction");
    fprintf(fp1, "random_kick_factor %-18g  # %s\n", random_kick_factor, "Norm of the initial displacement (activation)");
    fprintf(fp1, "cluster_radius      %-18g  # %s\n", cluster_radius, "The radius of the cluster that will be activated");
    fprintf(fp1, "init_step_size      %-18g  # %s\n", init_step_size, "Norm of the initial displacement (activation)");
    fprintf(fp1, "basin_factor        %-18g  # %s\n", basin_factor, "Factor multiplying Increment_Size for leaving the basin");
    fprintf(fp1, "min_num_ksteps      %-18d  # %s\n", min_num_ksteps, "Min # of k-steps before calling Lanczos");
    fprintf(fp1, "max_perp_move_h     %-18d  # %s\n", max_perp_move_h, "Max # of perpendicular steps leaving basin");
    fprintf(fp1, "max_iter_basin      %-18d  # %s\n", max_iter_basin, "Maximum # of iteration for leaving the basin");
    fprintf(fp1, "increment_size      %-18g  # %s\n", increment_size, "Overall scale for the increment moves");
    fprintf(fp1, "force_th_perp_h     %-18g  # %s\n", force_th_perp_h, "Perpendicular force threshold in harmonic well");
    fprintf(fp1, "eigen_th_well       %-18g  # %s\n", eigen_th_well, "Eigenvalue threshold for leaving basin");
    fprintf(fp1, "\n# activation, converging to saddle\n");
    fprintf(fp1, "max_activat_iter    %-18d  # %s\n", max_activat_iter, "Maximum # of iteraction to approach the saddle");
    fprintf(fp1, "fire_lanczos_every  %-18d  # %s\n", fire_lanczos_every, "call lanczos every # fire steps approaching saddle");
    fprintf(fp1, "use_fire            %-18d  # %s\n", use_fire, "Use FIRE for perpendicular steps approaching the saddle");
    fprintf(fp1, "para_factor         %-18g  # %s\n", para_factor, "Factor use to make parallel force nagative when approaching saddle");
    fprintf(fp1, "eigen_th_fail       %-18g  # %s\n", eigen_th_fail, "Eigen threshold for failure in searching the saddle");
    fprintf(fp1, "force_th_saddle     %-18g  # %s\n", force_th_saddle, "Force threshold for convergence at saddle point");
    fprintf(fp1, "conv_perp_inc       %-18d  # %s\n", conv_perp_inc, "The basic steps of max # of perpendicular steps approaching the saddle");
    fprintf(fp1, "max_perp_moves_c    %-18d  # %s\n", max_perp_moves_c, "Maximum # of perpendicular steps approaching the saddle");
    fprintf(fp1, "force_th_perp_sad   %-18g  # %s\n", force_th_perp_sad, "Perpendicular force threshold approaching saddle point");
    fprintf(fp1, "\n# confirmation of new found saddle\n");
    fprintf(fp1, "disp_sad2min_thr    %-18g  # %s\n", disp_sad2min_thr, "Minimum distance between saddle and original minimum");
    fprintf(fp1, "flag_push_back      %-18d  # %s\n", flag_push_back, "Push back the saddle to check its linkage to the start");
    fprintf(fp1, "push_over_saddle    %-18g  # %s\n", push_over_saddle, "Scale of displacement when pushing over the saddle");
    fprintf(fp1, "max_disp_tol        %-18g  # %s\n", max_disp_tol, "Tolerance displacement to claim the saddle is linked");
    fprintf(fp1, "max_ener_tol        %-18g  # %s\n", max_ener_tol, "Tolerance displacement to claim the saddle is linked");
    fprintf(fp1, "flag_relax_sad      %-18d  # %s\n", flag_relax_sad, "Further relax the newly found saddle via SD algorithm");
    fprintf(fp1, "SD_steps            %-18d  # %s\n", SD_steps, "Steepest Descent steps before CG minimizationm");
    fprintf(fp1, "\n# Lanczos related parameters\n");
    fprintf(fp1, "num_lancz_vec_h     %-18d  # %s\n", num_lancz_vec_h, "Num of vectors included in Lanczos for escaping well");
    fprintf(fp1, "num_lancz_vec_c     %-18d  # %s\n", num_lancz_vec_c, "Num of vectors included in Lanczos for convergence");
    fprintf(fp1, "del_disp_lancz      %-18g  # %s\n", del_disp_lancz, "Step of the numerical derivative of forces in Lanczos");
    fprintf(fp1, "eigen_th_lancz      %-18g  # %s\n", eigen_th_lancz, "Eigenvalue threshold for Lanczos convergence");
    fprintf(fp1, "\n# Metropolis\n");
    fprintf(fp1, "temperature         %-18g  # %s\n", temperature, "Temperature for Metropolis algorithm, in eV");
    fprintf(fp1, "atom_disp_thr       %-18g  # %s\n", atom_disp_thr, "Displacement threshold to identify an atom as displaced");
    fprintf(fp1, "\n# Output related parameters\n");
    fprintf(fp1, "log_file            %-18s  # %s\n", flog, "File to write ARTn log info; NULL to skip");
    fprintf(fp1, "log_level           %-18d  # %s\n", log_level, "Level of ARTn log ouput: 1, high; 0, low.");
    fprintf(fp1, "sadl_press_file     %-18s  # %s\n", c_fsadpress, "File to write ARTn sadl pressure info; NULL to skip");
    fprintf(fp1, "print_freq          %-18d  # %s\n", print_freq, "Print ARTn log ouput frequency, if log_level is 1.");
    fprintf(fp1, "event_list_file     %-18s  # %s\n", fevent, "File to record the event info; NULL to skip");
    fprintf(fp1, "dump_min_config     %-18s  # %s\n", fmin, "File for atomic dump of stable configurations; NULL to skip");
    fprintf(fp1, "dump_sad_config     %-18s  # %s\n", fsad, "file for atomic dump of saddle configurations; NULL to skip");
    fprintf(fp1, "dump_event_config   %-18s  # %s\n", fproc, "file for atomic dump of event configurations; NULL to skip");
    fprintf(fp1, "dump_min_every      %-18d  # %s\n", dump_min_every, "Dump min configuration every # step (if 0, no dump)");
    fprintf(fp1, "dump_sad_every      %-18d  # %s\n", dump_sad_every, "Dump sad configuration every # step (if 0, no dump)");
    fprintf(fp1, "dump_event_every    %-18d  # %s\n", dump_event_every, "Dump event process configuration every # step (if 0, no dump)");
    fprintf(fp1, "#====================================================================================================\n");
  }

  if (me == 0 && strcmp(fevent, "NULL") != 0){
    fp2 = fopen(fevent, "w");
    if (fp2 == NULL){
      sprintf(str, "Cannot open file: %s for writing", fevent);
      error->one(FLERR, str);
    }
  }

  // open dump files
  char **tmp;
  memory->create(tmp, 5, MAX(MAX(10,strlen(fmin)+1),MAX(strlen(fsad)+1,strlen(fproc)+1)), "ARTn string");

  char **format;
  memory->create(format, 3, 30, "ARTn format");
  strcpy(format[0],"format");
  strcpy(format[1],"line");
  strcpy(format[2],"%i %i %20.16g %20.16g %20.16g");
  int idump;
  for (idump = 0; idump < output-> ndump; idump++){
    if (strcmp("ARTnmin", output->dump[idump]->id) == 0) break;
  }
  if (idump == output->ndump) {
    if (strcmp(fmin, "NULL") != 0){
      strcpy(tmp[0],"ARTnmin");
      strcpy(tmp[1],"all");
      strcpy(tmp[2],"atom");
      strcpy(tmp[3],"1");
      strcpy(tmp[4],fmin);
      if (dump_min_every) {
	dumpmin = new DumpAtom(lmp, 5, tmp);
	dumpmin->modify_params(3, format);
      }
    }
  }else{
    // dump_min already defined in the lammps and use that one
    if(dump_min_every) {
      dumpmin = output->dump[idump];
      dumpmin_outside = true;
    }
  }

  for (idump = 0; idump < output-> ndump; idump++){
    if (strcmp("ARTnsad", output->dump[idump]->id) == 0) break;
  }
  if (idump == output->ndump) {
    if (strcmp(fsad, "NULL") != 0){
      strcpy(tmp[0],"ARTnsad");
      strcpy(tmp[1],"all");
      strcpy(tmp[2],"atom");
      strcpy(tmp[3],"1");
      strcpy(tmp[4],fsad);
      if(dump_sad_every) {
	dumpsad = new DumpAtom(lmp, 5, tmp);
        dumpsad->modify_params(3, format);
      }
    }
  }else{
    // dump_sad already defined in the lammps and use that one
    if(dump_sad_every) {
      dumpsad = output->dump[idump];
      dumpsad_outside = true;
    }
  }

  if (strcmp(fproc, "NULL") != 0){
    strcpy(tmp[0],"ARTneve");
    strcpy(tmp[1],"all");
    strcpy(tmp[2],"atom");
    strcpy(tmp[3],"1");
    strcpy(tmp[4],fproc);
    if(dump_event_every) dumpevent = new DumpAtom(lmp, 5, tmp);
  }

  memory->destroy(tmp);
  memory->destroy(format);

  delete []fmin;
  delete []fsad;
  delete []fproc;

return;
}

/* -------------------------------------------------------------------------------------------------
 * initializing ARTn
------------------------------------------------------------------------------------------------- */
void MinARTn::artn_init()
{
  random = new RanPark(lmp, seed+me);

  evalf = 0;
  flag_egvec = 0;

  // peratom vector I use
  fix_minimize->store_box();
  fix_minimize->add_vector(3);			//3
  fix_minimize->add_vector(3);			//4
  fix_minimize->add_vector(3);			//5
  fix_minimize->add_vector(3);			//6
  x0tmp = fix_minimize->request_vector(3);  //3
  egvec = fix_minimize->request_vector(4);  //4
  x00   = fix_minimize->request_vector(5);  //5
  fperp = fix_minimize->request_vector(6);  //6

  // group list
  int *tag  = atom->tag;
  int *mask = atom->mask;
  int *llist; memory->create(llist, MAX(1, atom->nlocal), "llist");

  int n = 0;
  for (int i = 0; i < atom->nlocal; ++i) if (mask[i] & groupbit) llist[n++] = tag[i];

  int nsingle = n, nall;
  MPI_Allreduce(&nsingle,&nall,1,MPI_INT,MPI_SUM,world);

  if (nall != ngroup) error->all(FLERR, "# of atoms in group mismatch!");
  if (me == 0) memory->create(glist, ngroup, "glist");

  int *disp = new int [np];
  int *recv = new int [np];
  for (int i = 0; i < np; ++i) disp[i] = recv[i] = 0;
  MPI_Gather(&nsingle,1,MPI_INT,recv,1,MPI_INT,0,world);
  for (int i = 1; i < np; ++i) disp[i] = disp[i-1] + recv[i-1];

  MPI_Gatherv(llist,nsingle,MPI_INT,glist,recv,disp,MPI_INT,0,world);
  delete [] disp;
  delete [] recv;
  memory->destroy(llist);


  // sort glist in the acending order
  if(me == 0) {
    for (int i = 0; i < ngroup; ++i) {
      for (int j = i + 1; j < ngroup; ++j) {
	if (glist[i] > glist[j]){
	  int tmp = glist[i];
	  glist[i] = glist[j];
	  glist[j] = tmp;
	}
      }
    }
  }

  if (dumpmin && !dumpmin_outside) dumpmin->init();
  if (dumpsad && !dumpsad_outside) dumpsad->init();
  if (dumpevent) dumpevent->init();

return;
}

/* -------------------------------------------------------------------------------------------------
 * reset vectors
------------------------------------------------------------------------------------------------- */
void MinARTn::artn_reset_vec()
{
  x0tmp = fix_minimize->request_vector(3);
  egvec = fix_minimize->request_vector(4);
  x00   = fix_minimize->request_vector(5);
  fperp = fix_minimize->request_vector(6);

return;
}

/* -------------------------------------------------------------------------------------------------
 * Try to find saddle point. If failed, return 0; else return 1
------------------------------------------------------------------------------------------------- */
int MinARTn::find_saddle( )
{
  
  int flag = 0;
  int nlanc = 0;
  double ftot = 0., ftotall, fpar2 = 0.,fpar2all = 0., fperp2 = 0., fperp2all = 0.,  delr;
  double fdoth = 0., fdothall = 0.;
  double step, preenergy;
  double tmp;
  int m_perp = 0, trial = 0, nfail = 0;
  double tmp_me[3], tmp_all[3];

  double force_thh_p2 = force_th_perp_h * force_th_perp_h;
  double force_thc_p2 = force_th_perp_sad * force_th_perp_sad;

  // record center-of-mass
  group->xcm(groupall, masstot, com0);

  egval = 0.;
  flag_egvec = 0;

  // record the original atomic positions
  for (int i = 0; i < nvec; ++i) x0[i] = x00[i] = xvec[i];

  // randomly displace the desired atoms: activation
  random_kick();
  if (dumpevent){
    int idum = update->ntimestep;
    update->ntimestep = 1000*sad_id+900;
    dumpevent->write();
    update->ntimestep = idum;
  }
  if (me == 0) print_info(10);

  int nmax_perp = max_perp_move_h;

  // try to leave harmonic well
  for (int it = 0; it < max_iter_basin; ++it){
    // minimizing perpendicularly by using SD method
    ecurrent = energy_force(1); ++evalf;
    artn_reset_vec(); reset_x00();

    m_perp = nfail = trial = 0;
    step = increment_size * 0.4;
    preenergy = ecurrent;
    while ( 1 ){
      //preenergy = ecurrent;
      fdoth = 0.;
      for (int i = 0; i < nvec; ++i) fdoth += fvec[i] * h[i];
      MPI_Allreduce(&fdoth, &fdothall,1,MPI_DOUBLE,MPI_SUM,world);
      fperp2 = 0.;
      
      for (int i = 0; i < nvec; ++i){ 
        fperp[i] = fvec[i] - fdothall * h[i];
        fperp2 += fperp[i] * fperp[i];
      }
      MPI_Allreduce(&fperp2, &fperp2all,1,MPI_DOUBLE,MPI_SUM,world);
      if (fperp2all < force_thh_p2 || m_perp > nmax_perp || nfail > 5) break; // condition to break
      
      for (int i = 0; i < nvec; ++i){
        x0tmp[i] = xvec[i];
        xvec[i] += MIN(dmax, MAX(-dmax, step * fperp[i]));
      }
      ecurrent = energy_force(1); ++evalf;
      artn_reset_vec(); reset_coords();
      
      if (ecurrent < preenergy){
        step *= 1.2;
        ++m_perp; nfail = 0;
        preenergy = ecurrent;

      } else {
        for (int i = 0; i < nvec; ++i) xvec[i] = x0tmp[i];
        step *= 0.5; ++nfail;
        ecurrent = energy_force(1); ++evalf;
        artn_reset_vec(); reset_x00();
      }
      ++trial;
    }
    
    ftot = delr = 0.;
    for (int i = 0; i < nvec; ++i) {
      ftot += fvec[i] * fvec[i];
      tmp = xvec[i] - x0[i];
      delr += tmp * tmp;
    }
    tmp_me[0] = ftot; tmp_me[1] = delr;
    MPI_Reduce(tmp_me, tmp_all,2,MPI_DOUBLE,MPI_SUM,0,world);
    ftotall = tmp_all[0]; delr = tmp_all[1];
    
    if (it > min_num_ksteps) nlanc = lanczos(flag_egvec, 1, num_lancz_vec_h);
    
    fpar2 = 0.;
    for (int i = 0; i < nvec; ++i) fpar2 += fvec[i] * h[i];
    MPI_Allreduce(&fpar2, &fpar2all,1,MPI_DOUBLE,MPI_SUM,world);
    
    delE = ecurrent-eref;
    if (dumpevent && (it%dump_event_every == 0)){
      int idum = update->ntimestep;
      update->ntimestep = 1000*sad_id+it;
      dumpevent->write();
      update->ntimestep = idum;
    }

    if (me == 0){
      fperp2 = sqrt(fperp2all);
      ftot   = sqrt(ftotall);
      delr   = sqrt(delr);
      if (fp1 && log_level && it%print_freq == 0) fprintf(fp1, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %10.5f %10.5f " BIGINT_FORMAT "\n", it,
      delE, m_perp, trial, nlanc, ftot, fpar2all, fperp2, egval, delr, evalf);
      if (screen && it%print_freq == 0) fprintf(screen, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %10.5f %10.5f " BIGINT_FORMAT "\n", it,
      delE, m_perp, trial, nlanc, ftot, fpar2all, fperp2, egval, delr, evalf);
    }
    
    if (it > min_num_ksteps && egval < eigen_th_well){
      idum = it;
      if (me == 0){
        if (fp1 && log_level && it%print_freq) fprintf(fp1, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %10.5f %10.5f " BIGINT_FORMAT "\n", it,
        delE, m_perp, trial, nlanc, ftot, fpar2all, fperp2, egval, delr, evalf);
        if (screen && it%print_freq) fprintf(screen, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %10.5f %10.5f " BIGINT_FORMAT "\n", it,
        delE, m_perp, trial, nlanc, ftot, fpar2all, fperp2, egval, delr, evalf);
      
        print_info(11);
      }

      flag = 1;
      break;
    }

    // push along the search direction
    step = basin_factor * increment_size;
    for(int i = 0; i < nvec; ++i) xvec[i] += MIN(dmax, MAX(-dmax, step * h[i]));
  }

  delE = ecurrent-eref;
  if (flag == 0){
    if (me == 0){
      if (fp1 && log_level && (max_iter_basin-1)%print_freq) fprintf(fp1, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %10.5f %10.5f " BIGINT_FORMAT "\n", (max_iter_basin-1),
      delE, m_perp, trial, nlanc, ftot, fpar2all, fperp2, egval, delr, evalf);
      if (screen && (max_iter_basin-1)%print_freq) fprintf(screen, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %10.5f %10.5f " BIGINT_FORMAT "\n", (max_iter_basin-1),
      delE, m_perp, trial, nlanc, ftot, fpar2all, fperp2, egval, delr, evalf);

      print_info(12);
    }
    reset_x00();
    for (int i = 0; i < nvec; ++i) xvec[i] = x00[i];

    return 0;
  }

  flag = 0; ++stage;
  if (me == 0) print_info(13);
  double hdot, hdotall;

  // now try to move close to the saddle point according to the egvec.
  int inc = conv_perp_inc;
  for (int it_s = 0; it_s < max_activat_iter; ++it_s){

    // g store old h
    for (int i = 0; i < nvec; ++i) g[i] = h[i];

    // caculate egvec use lanczos
    nlanc = lanczos(flag_egvec, 1, num_lancz_vec_c);
    if (flag_freeze){
      for (int i = 0; i < atom->nlocal; ++i){
	if (freeze_group_bit & atom->mask[i]){
	  int n = 3*i;
	  egvec[n] = 0.0;
	  egvec[n+1] = 0.0;
	  egvec[n+2] = 0.0;
	}
      }
    }
    for (int i = 0; i < nvec; ++i) h[i] = egvec[i];

    hdot = hdotall = 0.;
    for(int i = 0; i < nvec; ++i) hdot += h[i] * g[i];
    MPI_Reduce(&hdot,&hdotall,1,MPI_DOUBLE,MPI_SUM,0,world);

    // do minimizing perpendicular use SD or FIRE
    if (use_fire) {
      m_perp = trial = min_perp_fire(MIN(max_perp_moves_c, it_s + inc));
    } else {
      m_perp = trial = nfail = 0;
      step = increment_size * 0.25;
      int max_perp = MIN(max_perp_moves_c,it_s + inc);
      preenergy = ecurrent;
      while ( 1 ){
        fdoth = 0.;
        for (int i = 0; i < nvec; ++i) fdoth += fvec[i] * h[i];
        MPI_Allreduce(&fdoth, &fdothall,1,MPI_DOUBLE,MPI_SUM,world);
        fperp2 = 0.;
        for (int i = 0; i < nvec; ++i){
          fperp[i] = fvec[i] - fdothall * h[i];
          fperp2 += fperp[i] * fperp[i];
        }
        MPI_Allreduce(&fperp2, &fperp2all,1,MPI_DOUBLE,MPI_SUM,world);

        // condition to break
        if (fperp2all < force_thc_p2 || m_perp > max_perp || nfail > 5) break; 
        
        for (int i = 0; i < nvec; ++i){
          x0tmp[i] = xvec[i];
          xvec[i] += MIN(dmax, MAX(-dmax, step * fperp[i]));
        }
        ecurrent = energy_force(1); ++evalf;
        artn_reset_vec(); reset_coords(); 
        
        if (ecurrent < preenergy){
          step *= 1.2;
          ++m_perp; nfail = 0;
          preenergy = ecurrent;
        
        } else {
          for (int i = 0; i < nvec; ++i) xvec[i] = x0tmp[i];
          step *= 0.6; ++nfail;
          ecurrent = energy_force(1); ++evalf;
          artn_reset_vec(); reset_x00();
        }
        ++trial;
      }
    }
        
    tmp_me[0] = tmp_me[1] = tmp_me[2] = 0.;
    for (int i = 0; i < nvec; ++i) {
      tmp_me[0] += fvec[i] * fvec[i];
      delr = xvec[i] - x0[i];
      tmp_me[1] += delr * delr;
      tmp_me[2] += fvec[i] * h[i];
    }
    MPI_Allreduce(tmp_me, tmp_all, 3, MPI_DOUBLE, MPI_SUM, world);
    ftotall = sqrt(tmp_all[0]); delr = sqrt(tmp_all[1]); fpar2all = tmp_all[2];
   
    // output information
    fperp2 = 0.;
    for (int i = 0; i < nvec; ++i){
      tmp = fvec[i] - fpar2all * h[i];
      fperp2 +=  tmp * tmp;
    }
    MPI_Reduce(&fperp2, &fperp2all,1,MPI_DOUBLE,MPI_SUM,0,world);
  
    delE = ecurrent - eref;
    if (dumpevent && (it_s%dump_event_every == 0)){
      int idum = update->ntimestep;
      update->ntimestep = 1000*sad_id+100+it_s;
      dumpevent->write();
      update->ntimestep = idum;
    }

    if (me == 0){
      fperp2 = sqrt(fperp2all);
      if (fp1 && log_level && it_s%print_freq==0) fprintf(fp1, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %8.4f %8.4f %6.3f " BIGINT_FORMAT "\n",
      it_s, delE, m_perp, trial, nlanc, ftotall, fpar2all, fperp2, egval, delr, hdotall, evalf);
      if (screen && it_s%print_freq==0) fprintf(screen, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %8.4f %8.4f %6.3f " BIGINT_FORMAT "\n",
      it_s, delE, m_perp, trial, nlanc, ftotall, fpar2all, fperp2, egval, delr, hdotall, evalf);
    }
   
    if (egval > eigen_th_fail){
      if (me == 0){
        if (fp1 && log_level && it_s%print_freq) fprintf(fp1, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %8.4f %8.4f %6.3f " BIGINT_FORMAT "\n",
        it_s, delE, m_perp, trial, nlanc, ftotall, fpar2all, fperp2, egval, delr, hdotall, evalf);
        if (screen && it_s%print_freq) fprintf(screen, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %8.4f %8.4f %6.3f " BIGINT_FORMAT "\n",
        it_s, delE, m_perp, trial, nlanc, ftotall, fpar2all, fperp2, egval, delr, hdotall, evalf);
 
        print_info(14);
      }

      reset_x00();
      for(int i = 0; i < nvec; ++i) xvec[i] = x00[i];
  
      return 0;
    }

    if (delr < disp_sad2min_thr){
      if (me == 0){
        if (fp1 && log_level && it_s%print_freq) fprintf(fp1, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %8.4f %8.4f %6.3f " BIGINT_FORMAT "\n",
        it_s, delE, m_perp, trial, nlanc, ftotall, fpar2all, fperp2, egval, delr, hdotall, evalf);
        if (screen && it_s%print_freq) fprintf(screen, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %8.4f %8.4f %6.3f " BIGINT_FORMAT "\n",
        it_s, delE, m_perp, trial, nlanc, ftotall, fpar2all, fperp2, egval, delr, hdotall, evalf);
 
        ddum = delr;
        print_info(15);
      }

      reset_x00();
      for(int i = 0; i < nvec; ++i) xvec[i] = x00[i];
  
      return 0;
    }
  
    if (ftotall < force_th_saddle){
      if (me == 0){
        if (fp1 && log_level && it_s%print_freq) fprintf(fp1, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %8.4f %8.4f %6.3f " BIGINT_FORMAT "\n",
        it_s, delE, m_perp, trial, nlanc, ftotall, fpar2all, fperp2, egval, delr, hdotall, evalf);
        if (screen && it_s%print_freq) fprintf(screen, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %8.4f %8.4f %6.3f " BIGINT_FORMAT "\n",
        it_s, delE, m_perp, trial, nlanc, ftotall, fpar2all, fperp2, egval, delr, hdotall, evalf);

        idum = it_s;
        print_info(16);
      }

      return 1;
    }
  
    // caculate egvec use lanczos
    double sign = 1.;
    nlanc = lanczos(flag_egvec, 1, num_lancz_vec_c);
    if (flag_freeze){
      for (int i = 0; i < atom->nlocal; ++i){
	if (freeze_group_bit & atom->mask[i]){
	  int n = 3*i;
	  egvec[n] = 0.0;
	  egvec[n+1] = 0.0;
	  egvec[n+2] = 0.0;
	}
      }
    }
 
    double tmpsum = 0., tmpsumall;
    for (int i = 0; i < nvec; ++i) tmpsum += egvec[i] * fvec[i];
    MPI_Allreduce(&tmpsum,&tmpsumall,1,MPI_DOUBLE,MPI_SUM,world);

    if (tmpsumall > 0.) sign = -1.;

#define MinEGV 0.5 // was 0.5
    // push along the search direction; E. Cances, et al. JCP, 130, 114711 (2009)
    double factor = sign * MIN(2.*increment_size, fabs(fpar2all)/MAX(fabs(egval), MinEGV));
    for (int i = 0; i < nvec; ++i) xvec[i] += MIN(dmax, MAX(-dmax, factor * egvec[i]));
    ecurrent = energy_force(1); ++evalf;
    artn_reset_vec(); reset_x00();
  }

  delE = ecurrent - eref;
  if (me == 0){
    if (fp1 && log_level && (max_activat_iter-1)%print_freq) fprintf(fp1, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %8.4f %8.4f %6.3f " BIGINT_FORMAT "\n",
    (max_activat_iter-1), delE, m_perp, trial, nlanc, ftotall, fpar2all, fperp2, egval, delr, hdotall, evalf);
    if (screen && (max_activat_iter-1)%print_freq) fprintf(screen, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %8.4f %8.4f %6.3f " BIGINT_FORMAT "\n",
    (max_activat_iter-1), delE, m_perp, trial, nlanc, ftotall, fpar2all, fperp2, egval, delr, hdotall, evalf);

    print_info(17);
  }

  reset_x00();
  for (int i = 0; i < nvec; ++i) xvec[i] = x00[i];

return 0;
}
/* -------------------------------------------------------------------------------------------------
 * Try to find saddle point. If failed, return 0; else return 1
------------------------------------------------------------------------------------------------- */
int MinARTn::new_find_saddle( )
{
  int flag = 0;
  int nlanc = 0;
  double ftot = 0., ftotall, fpar2 = 0.,fpar2all = 0., fperp2 = 0., fperp2all = 0.,  delr;
  double fdoth = 0., fdothall = 0.;
  double step, preenergy;
  double tmp;
  int m_perp = 0, trial = 0, nfail = 0;
  double tmp_me[3], tmp_all[3];

  double force_thh_p2 = force_th_perp_h * force_th_perp_h;
  double force_thc_p2 = force_th_perp_sad * force_th_perp_sad;

  // record center-of-mass
  group->xcm(groupall, masstot, com0);

  egval = 0.;
  flag_egvec = 0;

  // record the original atomic positions
  for (int i = 0; i < nvec; ++i) x0[i] = x00[i] = xvec[i];

  // randomly displace the desired atoms: activation
  random_kick();

  if (me == 0) print_info(10);

  int nmax_perp = max_perp_move_h;

  // try to leave harmonic well
  for (int it = 0; it < max_iter_basin; ++it){
    // minimizing perpendicularly by using SD method
    ecurrent = energy_force(1); ++evalf;
    artn_reset_vec(); reset_x00();

    m_perp = nfail = trial = 0;
    step = increment_size * 0.4;
    preenergy = ecurrent;
    while ( 1 ){
      //preenergy = ecurrent;
      fdoth = 0.;
      for (int i = 0; i < nvec; ++i) fdoth += fvec[i] * h[i];
      MPI_Allreduce(&fdoth, &fdothall,1,MPI_DOUBLE,MPI_SUM,world);
      fperp2 = 0.;

      for (int i = 0; i < nvec; ++i){ 
	fperp[i] = fvec[i] - fdothall * h[i];
	fperp2 += fperp[i] * fperp[i];
      }
      MPI_Allreduce(&fperp2, &fperp2all,1,MPI_DOUBLE,MPI_SUM,world);
      if (fperp2all < force_thh_p2 || m_perp > nmax_perp || nfail > 5) break; // condition to break

      for (int i = 0; i < nvec; ++i){
	x0tmp[i] = xvec[i];
	xvec[i] += MIN(dmax, MAX(-dmax, step * fperp[i]));
      }
      ecurrent = energy_force(1); ++evalf;
      artn_reset_vec(); reset_coords();

      if (ecurrent < preenergy){
	step *= 1.2;
	++m_perp; nfail = 0;
	preenergy = ecurrent;

      } else {
	for (int i = 0; i < nvec; ++i) xvec[i] = x0tmp[i];
	step *= 0.5; ++nfail;
	ecurrent = energy_force(1); ++evalf;
	artn_reset_vec(); reset_x00();
      }
      ++trial;
    }

    ftot = delr = 0.;
    for (int i = 0; i < nvec; ++i) {
      ftot += fvec[i] * fvec[i];
      tmp = xvec[i] - x0[i];
      delr += tmp * tmp;
    }
    tmp_me[0] = ftot; tmp_me[1] = delr;
    MPI_Reduce(tmp_me, tmp_all,2,MPI_DOUBLE,MPI_SUM,0,world);
    ftotall = tmp_all[0]; delr = tmp_all[1];

    if (it > min_num_ksteps) nlanc = lanczos(flag_egvec, 1, num_lancz_vec_h);

    fpar2 = 0.;
    for (int i = 0; i < nvec; ++i) fpar2 += fvec[i] * h[i];
    MPI_Allreduce(&fpar2, &fpar2all,1,MPI_DOUBLE,MPI_SUM,world);

    delE = ecurrent-eref;
    if (me == 0){
      fperp2 = sqrt(fperp2all);
      ftot   = sqrt(ftotall);
      delr   = sqrt(delr);
      if (fp1 && log_level && it%print_freq == 0) fprintf(fp1, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %10.5f %10.5f " BIGINT_FORMAT "\n", it,
	  delE, m_perp, trial, nlanc, ftot, fpar2all, fperp2, egval, delr, evalf);
      if (screen && it%print_freq == 0) fprintf(screen, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %10.5f %10.5f " BIGINT_FORMAT "\n", it,
	  delE, m_perp, trial, nlanc, ftot, fpar2all, fperp2, egval, delr, evalf);
    }

    if (it > min_num_ksteps && egval < eigen_th_well){
      idum = it;
      if (me == 0){
	if (fp1 && log_level && it%print_freq) fprintf(fp1, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %10.5f %10.5f " BIGINT_FORMAT "\n", it,
	    delE, m_perp, trial, nlanc, ftot, fpar2all, fperp2, egval, delr, evalf);
	if (screen && it%print_freq) fprintf(screen, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %10.5f %10.5f " BIGINT_FORMAT "\n", it,
	    delE, m_perp, trial, nlanc, ftot, fpar2all, fperp2, egval, delr, evalf);

	print_info(11);
      }

      flag = 1;
      break;
    }

    // push along the search direction
    step = basin_factor * increment_size;
    for(int i = 0; i < nvec; ++i) xvec[i] += MIN(dmax, MAX(-dmax, step * h[i]));
  }

  delE = ecurrent-eref;
  if (flag == 0){
    if (me == 0){
      if (fp1 && log_level && (max_iter_basin-1)%print_freq) fprintf(fp1, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %10.5f %10.5f " BIGINT_FORMAT "\n", (max_iter_basin-1),
	  delE, m_perp, trial, nlanc, ftot, fpar2all, fperp2, egval, delr, evalf);
      if (screen && (max_iter_basin-1)%print_freq) fprintf(screen, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %10.5f %10.5f " BIGINT_FORMAT "\n", (max_iter_basin-1),
	  delE, m_perp, trial, nlanc, ftot, fpar2all, fperp2, egval, delr, evalf);

      print_info(12);
    }
    reset_x00();
    for (int i = 0; i < nvec; ++i) xvec[i] = x00[i];

    return 0;
  }

  flag = 0; ++stage;
  if (me == 0) print_info(13);
  double hdot, hdotall;

  // now try to move close to the saddle point according to the egvec.
  int inc = conv_perp_inc;
  // g store old h
  for (int i = 0; i < nvec; ++i) g[i] = h[i];

  // caculate egvec use lanczos
  nlanc = lanczos(flag_egvec, 1, num_lancz_vec_c);
  for (int i = 0; i < nvec; ++i) h[i] = egvec[i];

  hdot = hdotall = 0.;
  for(int i = 0; i < nvec; ++i) hdot += h[i] * g[i];
  MPI_Reduce(&hdot,&hdotall,1,MPI_DOUBLE,MPI_SUM,0,world);

  // do minimizing perpendicular use  FIRE
  stop_condition = new_min_perp_fire(max_conv_steps);
  stopstr = stopstrings(stop_condition);
  // output minimization information
  if(stop_condition == FTOL) return 1;else return 0;
}
/* -----------------------------------------------------------------------------
 * reset coordinates x0tmp
 * ---------------------------------------------------------------------------*/
void MinARTn::reset_coords()
{
  domain->set_global_box();

  double **x = atom->x;
  int nlocal = atom->nlocal;
  double dx,dy,dz,dx0,dy0,dz0;

  int n = 0;
  for (int i = 0; i < nlocal; ++i) {
    dx = dx0 = x[i][0] - x0tmp[n];
    dy = dy0 = x[i][1] - x0tmp[n+1];
    dz = dz0 = x[i][2] - x0tmp[n+2];
    domain->minimum_image(dx,dy,dz);

    if (dx != dx0) x0tmp[n]   = x[i][0] - dx;
    if (dy != dy0) x0tmp[n+1] = x[i][1] - dy;
    if (dz != dz0) x0tmp[n+2] = x[i][2] - dz;

    n += 3;
  }

  domain->set_global_box();

  reset_x00();

return;
}

/* -----------------------------------------------------------------------------
 * reset coordinates x00
 * ---------------------------------------------------------------------------*/
void MinARTn::reset_x00()
{
  domain->set_global_box();

  double **x = atom->x;
  int nlocal = atom->nlocal;
  double dx,dy,dz,dx0,dy0,dz0;

  int n = 0;
  for (int i = 0; i < nlocal; ++i) {
    dx = dx0 = x[i][0] - x00[n];
    dy = dy0 = x[i][1] - x00[n+1];
    dz = dz0 = x[i][2] - x00[n+2];
    domain->minimum_image(dx,dy,dz);

    if (dx != dx0) x00[n]   = x[i][0] - dx;
    if (dy != dy0) x00[n+1] = x[i][1] - dy;
    if (dz != dz0) x00[n+2] = x[i][2] - dz;

    n += 3;
  }

  domain->set_global_box();

return;
}
/* ------------------------------------------------------------------------------
 * random kick atoms in the defined group and radius
 * ----------------------------------------------------------------------------*/
void MinARTn::random_kick()
{
  // define the central atom that will be activated
  if (me == 0){
    if(events_per_atom == 0){
      int index = int(random->uniform()*double(ngroup))%ngroup;
      that = glist[index];
    }else{
      that = glist[iatom];
    }
    print_info(18);
  }
  MPI_Bcast(&that, 1, MPI_INT, 0, world);

  double cord[3];

  double *delpos = fix_minimize->request_vector(4);
  double *tmpdelpos = fix_minimize->request_vector(6);
  for (int i = 0; i < nvec; ++i) delpos[i] = 0.;
  int nlocal = atom->nlocal;
  int natoms = atom->natoms;
  int *tag   = atom->tag;
  int nhit = 0;


  if (flag_dump_direction){
    read_dump_direction(fdump_direction,tmpdelpos);
  }
  if (fabs(cluster_radius) < ZERO){ // only the cord atom will be kicked
    for (int i = 0; i < nlocal; ++i){
      if (tag[i] == that){
        int n = 3*i;
        delpos[n  ] = 0.5 - random->uniform();
        delpos[n+1] = 0.5 - random->uniform();
        if(domain->dimension == 3)delpos[n+2] = 0.5 - random->uniform();
        else delpos[n+2] = 0.;

        ++nhit; break;
      }
    }

  } else if (cluster_radius < 0.){ // all atoms in group will be kicked
    for (int i = 0; i < nlocal; ++i){
      if (groupbit & atom->mask[i]){
        int n = 3*i;
        delpos[n  ] = 0.5 - random->uniform();
        delpos[n+1] = 0.5 - random->uniform();
        if(domain->dimension == 3)delpos[n+2] = 0.5 - random->uniform();
        else delpos[n+2] = 0.;
        ++nhit;
      }
    }

  } else { // only atoms within a radius to the central atom will be kicked
    double one[3]; one[0] = one[1] = one[2] = 0.;
    for (int i = 0; i < nlocal; ++i){
      if (tag[i] == that){
        one[0] = atom->x[i][0];
        one[1] = atom->x[i][1];
        one[2] = atom->x[i][2];
      }
    }
    MPI_Allreduce(&one[0], &cord[0], 3, MPI_DOUBLE, MPI_SUM, world);
    if (flag_deformation_gradient && flag_random_deformation_gradient){
      // creat a random orthogonal tensor
      if (me == 0){
	double * tmp, sum, dotsum, dotsum1, dotsum2;
	int i;
	tmp = new double[9];
	// a random tensor
	for (i = 0; i < 9; ++i){
	  tmp[i] = random->uniform();
	}
	// normolize the first vector
	sum = 0.;
	for (i = 0; i < 3; ++i){
	  sum = sum + tmp[i] * tmp[i];
	}
	sum = 1. / sqrt(sum);
	for (i = 0; i < 3; ++i){
	  tmp[i] = tmp[i] * sum;
	}
	// calculate the second vector
	dotsum = 0.;
	for (i = 3; i < 6; ++i){
	  dotsum = dotsum + tmp[i] * tmp[i-3];
	}
	for (i = 3; i < 6; ++i){
	  tmp[i] = tmp[i] - dotsum * tmp[i-3];
	}
	sum = 0.;
	for (i = 3; i < 6; ++i){
	  sum = sum + tmp[i] * tmp[i];
	}
	sum = 1./ sqrt(sum);
	for (i = 3; i < 6; ++i){
	  tmp[i] = tmp[i] * sum;
	}
	// calculate the third vector
	dotsum1 = 0.; dotsum2 = 0.;
	for (i = 6; i < 9; ++i){
	  dotsum1 = dotsum1 + tmp[i] * tmp[i-3];
	  dotsum2 = dotsum2 + tmp[i] * tmp[i-6];
	}
	for (i = 6; i < 9; ++i){
	  tmp[i] = tmp[i] - dotsum1 * tmp[i-3] - dotsum2 * tmp[i-6];
	}
	sum = 0.;
	for (i = 6; i < 9; ++i){
	  sum = sum + tmp[i] * tmp[i];
	}
	sum = 1./ sqrt(sum);
	for (i = 6; i < 9; ++i){
	  tmp[i] = tmp[i] * sum;
	}
	for (i = 0; i < 9; ++i){
	  deformation_gradient[i] = tmp[i];
	}
	delete []tmp;
      }
      // broadcast the tensor
      MPI_Bcast(deformation_gradient, 9, MPI_DOUBLE, 0, world);
    }
    double rcut2 = cluster_radius * cluster_radius;
    for (int i = 0; i < nlocal; ++i){
      //if (groupbit & atom->mask[i]){
      double dx = atom->x[i][0] - cord[0];
      double dy = atom->x[i][1] - cord[1];
      double dz = atom->x[i][2] - cord[2];
      domain->minimum_image(dx, dy, dz);
      double r2 = dx*dx + dy*dy + dz*dz;
      if (r2 <= rcut2){
        if (flag_deformation_gradient){
          double tmpx = deformation_gradient[0] * dx + deformation_gradient[1] * dy + deformation_gradient[2] * dz;
          double tmpy = deformation_gradient[3] * dx + deformation_gradient[4] * dy + deformation_gradient[5] * dz;
          double tmpz = deformation_gradient[6] * dx + deformation_gradient[7] * dy + deformation_gradient[8] * dz;
	  tmpdelpos[3*i] = tmpx - dx;
	  tmpdelpos[3*i+1] = tmpy - dy;
	  tmpdelpos[3*i+2] = tmpz - dz;
        }
        int n = 3*i;
        delpos[n  ] = 0.5 - random->uniform();
        delpos[n+1] = 0.5 - random->uniform();
        if(domain->dimension == 3)delpos[n+2] = 0.5 - random->uniform();
        else delpos[n+2] = 0.;
        ++nhit;
      }
      //}
    }
  }
  double norm = 0., normall;
  double norm2 = 0., normall2;
  for (int i = 0; i < nvec; ++i){
    norm += delpos[i] * delpos[i];
    norm2 += tmpdelpos[i] * tmpdelpos[i];
  }
  MPI_Allreduce(&norm,&normall,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&norm2,&normall2,1,MPI_DOUBLE,MPI_SUM,world);
  double norm_i = 1./sqrt(normall);
  double norm2_i = 1./sqrt(normall2);
  for (int i = 0; i < nvec; ++i){
    delpos[i] = delpos[i] * norm_i;
    tmpdelpos[i] = tmpdelpos[i] * norm2_i;
  }

  if (flag_dump_direction || flag_deformation_gradient){
    for (int i = 0; i < 3*nlocal; ++i){
      delpos[i] = delpos[i] * random_kick_factor + ( 1 - random_kick_factor) * tmpdelpos[i];
    }
  }

  // minus x,y,z drift
  //double dx, dy, dz;
  //double tmp[3],tmpall[3];
  //dx = dy = dz = 0.0;
  //for (int i = 0; i < nlocal; ++i){
  //  dx += delpos[i*3];
  //  dy += delpos[i*3+1];
  //  dz += delpos[i*3+2];
  //}
  //tmp[0] = dx; tmp[1] = dy; tmp[2] = dz;
  //MPI_Allreduce(tmp, tmpall, 3, MPI_DOUBLE, MPI_SUM, world);
  ////if(me == 0)fprintf(screen, "\n dx = %f, dy = %f, dz = %f \n", tmpall[0], tmpall[1] ,tmpall[2]);
  //dx = tmpall[0] / natoms;
  //dy = tmpall[1] / natoms;
  //dz = tmpall[2] / natoms;
  //for (int i = 0; i < nlocal; ++i){
  //  delpos[i*3] -= dx;
  //  delpos[i*3+1] -= dy;
  //  delpos[i*3+2] -= dz;
  //}

  MPI_Reduce(&nhit,&idum,1,MPI_INT,MPI_SUM,0,world);
  if (me == 0 && idum < 1) error->one(FLERR, "No atom to kick!");
  if (me == 0) print_info(19);

  // now normalize and apply the kick to the selected atom(s)
  norm = 0.;
  for (int i = 0; i < nvec; ++i) norm += delpos[i] * delpos[i];
  MPI_Allreduce(&norm,&normall,1,MPI_DOUBLE,MPI_SUM,world);

  norm_i = 1./sqrt(normall);
  for (int i = 0; i < nvec; ++i){
    h[i] = delpos[i] * norm_i;
    xvec[i] += MIN(dmax, MAX(-dmax, init_step_size * h[i]));
  }

return;
}

/* -----------------------------------------------------------------------------
 * converge to minimum, here I use conjugate gradient method.
 * ---------------------------------------------------------------------------*/
int MinARTn::min_converge(int maxiter, const int flag)
{
  neval = 0;
  int i,fail;
  int ntimestep;
  double beta,gg,dot[2],dotall[2];

  // nlimit = max # of CG iterations before restarting
  // set to ndoftotal unless too big
  int nlimit = static_cast<int> (MIN(MAXSMALLINT,ndoftotal));

  // initialize working vectors

  eprevious = ecurrent = energy_force(1); ++neval;
  for (i = 0; i < nvec; ++i) h[i] = g[i] = fvec[i];

  gg = fnorm_sqr();
  double ftol_sq = update->ftol*update->ftol;

  niter = 0;
  for (int iter = 0; iter < maxiter; ++iter) {
    ntimestep = ++update->ntimestep;
    ++niter;

    // line minimization along direction h from current atom->x
    fail = (this->*linemin)(ecurrent,alpha_final);
    if (fail) return fail;

    // function evaluation criterion
    if (neval >= update->max_eval) return MAXEVAL;

    // energy tolerance criterion
    if (fabs(ecurrent-eprevious) < update->etol * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS_ENERGY)) return ETOL;

    // force tolerance criterion

    dot[0] = dot[1] = 0.;
    for (i = 0; i < nvec; ++i) {
      dot[0] += fvec[i]*fvec[i];
      dot[1] += fvec[i]*g[i];
    }
    MPI_Allreduce(dot,dotall,2,MPI_DOUBLE,MPI_SUM,world);

    if (dotall[0] < ftol_sq) return FTOL;

    // update new search direction h from new f
    // = -Grad(x) and old g
    // this is Polak-Ribieri formulation
    // beta = dotall[0]/gg would be
    // Fletcher-Reeves
    // reinitialize CG every ndof iterations by
    // setting beta = 0.0

    beta = MAX(0.,(dotall[0] - dotall[1])/gg);
    if ((niter+1) % nlimit == 0) beta = 0.;
    gg = dotall[0];

    for (i = 0; i < nvec; ++i) {
      g[i] = fvec[i];
      h[i] = g[i] + beta*h[i];
    }

    // reinitialize CG
    // if new search
    // direction h is
    // not push_down

    dot[0] = 0.;
    for (i = 0; i < nvec; ++i) dot[0] += g[i]*h[i];
    MPI_Allreduce(dot,dotall,1,MPI_DOUBLE,MPI_SUM,world);

    if (dotall[0] <= 0.) for (i = 0; i < nvec; ++i) h[i] = g[i];
    if (flag == 2) reset_coords();
    else if (flag == 1) reset_x00();
    // output for thermo, dump, restart files
    //if (output->next == ntimestep) {
    //  timer->stamp();
    //  output->write(ntimestep);
    //  timer->stamp(Timer::OUTPUT);
    //}
  }

return MAXITER;
}

/* -----------------------------------------------------------------------------
 * The lanczos method to get the lowest eigenvalue and corresponding eigenvector
 * Ref: R.A. Olsen , G. J. Kroes, Comparison of methods for 
 *   fanding saddle points without knowledge of final states,
 *   121, 20(2004)
 * ---------------------------------------------------------------------------*/
int MinARTn::lanczos(bool egvec_exist, int flag, int maxvec){
  FixMinimize * fix_lanczos;
  char **fixarg = new char*[3];
  fixarg[0] = (char *) "artn_lanczos";
  fixarg[1] = (char *) "all";
  fixarg[2] = (char *) "MINIMIZE";
  modify->add_fix(3,fixarg);
  delete [] fixarg;
  fix_lanczos = (FixMinimize *) modify->fix[modify->nfix-1];

  fix_lanczos->add_vector(3);		// 0, for r_k_1
  fix_lanczos->add_vector(3);		// 1, for q_k_1
  fix_lanczos->add_vector(3);		// 2, for q_k
  fix_lanczos->add_vector(3);		// 3, for u_k
  fix_lanczos->add_vector(3);		// 4, for r_k
  fix_lanczos->add_vector(3);		// 5, for f0;
  double *r_k_1 = fix_lanczos->request_vector(0);
  double *q_k_1 = fix_lanczos->request_vector(1);
  for (int i = 0; i < nvec; ++i) q_k_1[i] = 0.;
  double *q_k = fix_lanczos->request_vector(2);
  double *u_k = fix_lanczos->request_vector(3);
  double *r_k = fix_lanczos->request_vector(4);
  double *f0  = fix_lanczos->request_vector(5);
  double *d = new double [maxvec];
  double *e = new double [maxvec];
  double *d_bak = new double [maxvec];
  double *e_bak = new double [maxvec];
  double tmp;
  double beta_k_1 = 0., alpha_k, beta_k;
  int con_flag = 0;
  double ** lanc = new double *[maxvec];
  for (int i = 0; i < maxvec; ++i){
    fix_lanczos->add_vector(3);
    lanc[i] = fix_lanczos->request_vector(i+6);
  }
  // DEL_LANCZOS is the step we use in the finite differece approximation.
  const double DEL_LANCZOS = del_disp_lancz;
  const double IDEL_LANCZOS = 1.0 / DEL_LANCZOS;

  // set r(k-1) according to egvec or random vector
  int nlocal = atom->nlocal;
  if (egvec_exist) for (int i = 0; i < nvec; ++i) r_k_1[i] = egvec[i];
  else for (int i = 0; i < nlocal; ++i) {
    int n = i*3;
    r_k_1[n]   = 0.5 - random->uniform();
    r_k_1[n+1] = 0.5 - random->uniform();
    if(domain->dimension == 3)r_k_1[n+2] = 0.5 - random->uniform();
    else r_k_1[n+2] = 0.;
  }
  for (int i =0; i < nvec; ++i) beta_k_1 += r_k_1[i] * r_k_1[i];

  MPI_Allreduce(&beta_k_1,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  beta_k_1 = sqrt(tmp);

  double eigen1 = 0., eigen2 = 0.;
  char jobs = 'V';
  double *work, *z;
  int ldz = maxvec, info;
  z = new double [ldz*maxvec];
  work = new double [2*maxvec];

  // store origin configuration and force
  for (int i = 0; i < nvec; ++i){
    x0tmp[i] = xvec[i];
    f0[i] = fvec[i];
  }
  int n;
  for (n = 1; n <= maxvec; ++n){
    for (int i = 0; i < nvec; ++i){
      q_k[i] = r_k_1[i] / beta_k_1;
      lanc[n-1][i] = q_k[i];
    }

    //reset_coords();
    // random move to caculate u(k) with the finite difference approximation
    for (int i = 0; i < nvec; ++i) xvec[i] = x0tmp[i] + q_k[i] * DEL_LANCZOS;

    energy_force(1); ++evalf;
    reset_coords();

    r_k_1 = fix_lanczos->request_vector(0);
    q_k_1 = fix_lanczos->request_vector(1);
    q_k = fix_lanczos->request_vector(2);
    u_k = fix_lanczos->request_vector(3);
    r_k = fix_lanczos->request_vector(4);
    f0  = fix_lanczos->request_vector(5);
    artn_reset_vec();
    for (int i = 0; i < maxvec; ++i){
      lanc[i] = fix_lanczos->request_vector(i+6);
    }
    alpha_k = 0.;

    for (int i = 0; i < nvec; ++i){
      u_k[i] = (f0[i] - fvec[i]) * IDEL_LANCZOS;
      r_k[i] = u_k[i] - beta_k_1 * q_k_1[i];
      alpha_k += q_k[i] * r_k[i];
    }
    MPI_Allreduce(&alpha_k,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    alpha_k = tmp;
    beta_k = 0.;

    for (int i = 0; i < nvec; ++i){
      r_k[i] -= alpha_k * q_k[i];
      beta_k += r_k[i] * r_k[i];
    }
    MPI_Allreduce(&beta_k,&tmp,1,MPI_DOUBLE,MPI_SUM,world);

    beta_k = sqrt(tmp);
    d[n-1] = alpha_k;
    e[n-1] = beta_k;
    for (int i = 0; i != maxvec; ++i){
      d_bak[i] = d[i];
      e_bak[i] = e[i];
    }
    if (n >= 2){
      dstev_(&jobs, &n, d_bak, e_bak, z, &ldz, work, &info);

      if (info != 0){
	char str[MAXLINE];
        sprintf(str, "ARTn: dstev_ error in Lanczos subroute. Error Info = %i. \n(Info < 0: the i-th argument had an illegal value; Info > 0: i off-diagonal elements did not converged to zero.", info);
       	error->all(FLERR,str);
      }

      eigen1 = eigen2; eigen2 = d_bak[0];
    }
    if (n >= 3 && fabs((eigen2-eigen1)/eigen1) < eigen_th_lancz) {
      con_flag = 1;
      for (int i = 0; i < nvec; ++i){
	     xvec[i] = x0tmp[i];
	     fvec[i] = f0[i];
      }
      egval = eigen2;
      if (flag > 0){
        flag_egvec = 1;
        for (int i = 0; i < nvec; ++i) egvec[i] = 0.;
        for (int i = 0; i < nvec; ++i)
        for (int j = 0; j < n; ++j) egvec[i] += z[j] * lanc[j][i];

        // normalize egvec.
        double sum = 0., sumall;
        for (int i = 0; i < nvec; ++i) sum += egvec[i] * egvec[i];
        
        MPI_Allreduce(&sum, &sumall,1,MPI_DOUBLE,MPI_SUM,world);
        sumall = 1. / sqrt(sumall);
        for (int i = 0; i < nvec; ++i) egvec[i] *= sumall;
      }
      break;
    }
    for (int i = 0; i < nvec; ++i){
      r_k_1[i] = r_k[i];
      q_k_1[i] = q_k[i];
    }
    beta_k_1 = beta_k;
  }

  if (con_flag == 0){
    //outlog("WARNING, LNACZOS method not converged!\n");
    egval = eigen2;
    if (flag > 0){
      flag_egvec = 1;
      for (int i = 0; i < nvec; ++i) egvec[i] = 0.;
      for (int i = 0; i < nvec; ++i)
      for (int j = 0; j < n-1; ++j) egvec[i] += z[j] * lanc[j][i];

      // normalize egvec.
      double sum = 0., sumall;
      for (int i = 0; i < nvec; ++i) sum += egvec[i] * egvec[i];

      MPI_Allreduce(&sum, &sumall,1,MPI_DOUBLE,MPI_SUM,world);
      sumall = 1. / sqrt(sumall);
      for (int i = 0; i < nvec; ++i) egvec[i] *= sumall;
    }

    for (int i = 0; i < nvec; ++i){
      xvec[i] = x0tmp[i];
      fvec[i] = f0[i];
    }
    energy_force(1); ++evalf;
    reset_coords();
    artn_reset_vec();
  }

  //if (fabs(egval - 0.0) < 1e-4) {
  //  flag_egvec = 0;
  //}

  delete []d;
  delete []e;
  delete []d_bak;
  delete []e_bak;
  delete []z;
  delete []work;
  delete []lanc;

  modify->delete_fix("artn_lanczos");

return MIN(int(n),maxvec);
}

/* ---------------------------------------------------------------------------
 *  FIRE: fast interial relaxation engine, return iteration number
 * -------------------------------------------------------------------------*/
int MinARTn::min_perp_fire(int maxiter)
{
  double dt = update->dt;
  const int n_min = 5;
  const double f_inc = 1.1;
  const double f_dec = 0.5;
  const double alpha_start = 0.1;
  const double f_alpha = 0.99;
  const double  TMAX = 10.;
  const double dtmax = TMAX * dt;
  double vdotf, vdotfall;
  double vdotvall;
  double fdotfall;
  double fdoth, fdothall;
  double scale1, scale2;
  double alpha;
  int last_negative = 0;

  double force_thr2 = force_th_perp_sad*force_th_perp_sad;
  double *vvec = atom->v[0];

  for (int i = 0; i < nvec; ++i) vvec[i] = 0.;

  alpha = alpha_start;
  for (int iter = 0; iter < maxiter; ++iter){

    fdoth = 0.;
    for (int i = 0; i < nvec; ++i) fdoth += fvec[i] * h[i];
    MPI_Allreduce(&fdoth,&fdothall,1,MPI_DOUBLE,MPI_SUM,world);

    for (int i = 0; i < nvec; ++i) fperp[i] = fvec[i] - fdothall * h[i];

    vdotf = 0.;
    vvec = atom->v[0];
    for (int i = 0; i < nvec; ++i) vdotf += vvec[i] * fperp[i];
    MPI_Allreduce(&vdotf,&vdotfall,1,MPI_DOUBLE,MPI_SUM,world);
    
    // if (v dot f) > 0:
    // v = (1-alpha) v + alpha |v| Fhat
    // |v| = length of v, Fhat = unit f
    // if more than DELAYSTEP since v dot f was negative:
    // increase timestep and decrease alpha
    if (vdotfall > 0.) {
      double tmp_me[2], tmp_all[2];
      scale1 = 1. - alpha;
      tmp_me[0] = tmp_me[1] = 0.;
      for (int i = 0; i < nvec; ++i) {
        tmp_me[0] += vvec[i] * vvec[i];
        tmp_me[1] += fperp[i] * fperp[i];
      }
      MPI_Allreduce(tmp_me, tmp_all,2,MPI_DOUBLE,MPI_SUM,world);
      vdotvall = tmp_all[0]; fdotfall = tmp_all[1];
      
      if (fdotfall < force_thr2) return iter;

      if (fdotfall == 0.) scale2 = 0.;
      else scale2 = alpha * sqrt(vdotvall/fdotfall);
      for (int i = 0; i < nvec; ++i) vvec[i] = scale1 * vvec[i] + scale2 * fperp[i];
      
      if ((iter - last_negative) > n_min) {
        dt = MIN(dt * f_inc, dtmax);
        alpha = alpha * f_alpha;
      }
      
    } else {
      last_negative = iter;
      dt = dt * f_dec;
      for (int i = 0; i < nvec; ++i) vvec[i] = 0.;
      alpha = alpha_start;
    }
    double **x = atom->x;
    double **v = atom->v;
    double *rmass = atom->rmass;
    double *mass = atom->mass;
    int *type = atom->type;
    double dtfm;
    // limit timestep so no particle moves further than dmax
    double dtvone = dt;
    double vmax = 0.;
    double dtv;
    for (int i = 0; i < atom->nlocal; i++) {
      vmax = MAX(fabs(v[i][0]),fabs(v[i][1]));
      vmax = MAX(vmax,fabs(v[i][2]));
      if (dtvone*vmax > dmax) dtvone = dmax/vmax;
    }
    MPI_Allreduce(&dtvone,&dtv,1,MPI_DOUBLE,MPI_MIN,world);

    double dtf = dtv * force->ftm2v;
    int n = 0;
    if (rmass) {
      for (int i = 0; i < atom->nlocal; ++i) {
        dtfm = dtf / rmass[i];
        x[i][0] += dt * v[i][0];
        x[i][1] += dt * v[i][1];
        x[i][2] += dt * v[i][2];
        v[i][0] += dtfm * fperp[n++];
        v[i][1] += dtfm * fperp[n++];
        v[i][2] += dtfm * fperp[n++];
      }
    } else {
      for (int i = 0; i < atom->nlocal; ++i) {
        dtfm = dtf / mass[type[i]];
        x[i][0] += dt * v[i][0];
        x[i][1] += dt * v[i][1];
        x[i][2] += dt * v[i][2];
        v[i][0] += dtfm * fperp[n++];
        v[i][1] += dtfm * fperp[n++];
        v[i][2] += dtfm * fperp[n++];
      }
    }
    eprevious = ecurrent;
    ecurrent = energy_force(1); ++evalf;
    artn_reset_vec();reset_x00();
  }

return maxiter;
} 

/* ---------------------------------------------------------------------------
 *  FIRE: fast interial relaxation engine, return iteration number
 * -------------------------------------------------------------------------*/
int MinARTn::new_min_perp_fire(int maxiter)
{
  double dt = update->dt;
  const int n_min = 5;
  const double f_inc = 1.1;
  const double f_dec = 0.5;
  const double alpha_start = 0.1;
  const double f_alpha = 0.99;
  const double  TMAX = 10.;
  const double dtmax = TMAX * dt;
  double vdotf, vdotfall;
  double vdotvall;
  double fdotfall;
  double fdoth, fdothall;
  double scale1, scale2;
  double alpha;
  int last_negative = 0;
  int nlanc = 0;
  double hdot = 0., hdotall = 0.;
  double tmp_me[3]={0.}, tmp_all[3]={0.};
  double tmp = 0.;
  double delr, ftotall, fpar2all;
  double fperp2, fperp2all;

  const double force_thr2 = force_th_saddle*force_th_saddle;

  double *vvec = atom->v[0];
  for (int i = 0; i < nvec; ++i) vvec[i] = 0.;

  alpha = alpha_start;
  for (int iter = 0; iter < maxiter; ++iter){
    
    //ntimestep = ++update->ntimestep;

    if(iter % fire_lanczos_every == 0){
      // g store old h
      for (int i = 0; i < nvec; ++i) g[i] = h[i];

      // caculate egvec use lanczos
      nlanc = lanczos(flag_egvec, 1, num_lancz_vec_c);
      for (int i = 0; i < nvec; ++i) h[i] = egvec[i];

      hdot = hdotall = 0.;
      for(int i = 0; i < nvec; ++i) hdot += h[i] * g[i];
      MPI_Reduce(&hdot,&hdotall,1,MPI_DOUBLE,MPI_SUM,0,world);
    }
    if(iter % fire_output_every == 0){
      tmp_me[0] = tmp_me[1] = tmp_me[2] = 0.;
      for (int i = 0; i < nvec; ++i) {
	tmp_me[0] += fvec[i] * fvec[i];
	delr = xvec[i] - x0[i];
	tmp_me[1] += delr * delr;
	tmp_me[2] += fvec[i] * h[i];
      }
      MPI_Allreduce(tmp_me, tmp_all, 3, MPI_DOUBLE, MPI_SUM, world);
      ftotall = sqrt(tmp_all[0]); delr = sqrt(tmp_all[1]); fpar2all = tmp_all[2];
      fperp2 = 0.;
      for (int i = 0; i < nvec; ++i){
	tmp = fvec[i] - fpar2all * h[i];
	fperp2 +=  tmp * tmp;
      }
      MPI_Reduce(&fperp2, &fperp2all,1,MPI_DOUBLE,MPI_SUM,0,world);
      delE = ecurrent - eref;
      int m_perp= 0, trial = 0;
      if (me == 0){
	fperp2 = sqrt(fperp2all);
	if (fp1 && log_level ) fprintf(fp1, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %8.4f %8.4f %6.3f " BIGINT_FORMAT "\n",
	    iter, delE, m_perp, trial, nlanc, ftotall, fpar2all, fperp2, egval, delr, hdotall, evalf);
	if (screen) fprintf(screen, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %8.4f %8.4f %6.3f " BIGINT_FORMAT "\n",
	    iter, delE, m_perp, trial, nlanc, ftotall, fpar2all, fperp2, egval, delr, hdotall, evalf);
      }
      if (egval > eigen_th_fail){
	if (me == 0){
	  if (fp1 && log_level) fprintf(fp1, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %8.4f %8.4f %6.3f " BIGINT_FORMAT "\n",
	      iter, delE, m_perp, trial, nlanc, ftotall, fpar2all, fperp2, egval, delr, hdotall, evalf);
	  if (screen) fprintf(screen, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %8.4f %8.4f %6.3f " BIGINT_FORMAT "\n",
	      iter, delE, m_perp, trial, nlanc, ftotall, fpar2all, fperp2, egval, delr, hdotall, evalf);

	  print_info(14);
	}

	reset_x00();
	for(int i = 0; i < nvec; ++i) xvec[i] = x00[i];

	return ETOL;
      }
    }
    fdoth = 0.;
    for (int i = 0; i < nvec; ++i) fdoth += fvec[i] * h[i];
    MPI_Allreduce(&fdoth,&fdothall,1,MPI_DOUBLE,MPI_SUM,world);

    for (int i = 0; i < nvec; ++i) fperp[i] = fvec[i] - (1 + para_factor) * fdothall * h[i];

    vdotf = 0.;
    vvec = atom->v[0];
    for (int i = 0; i < nvec; ++i) vdotf += vvec[i] * fperp[i];
    MPI_Allreduce(&vdotf,&vdotfall,1,MPI_DOUBLE,MPI_SUM,world);
    
    // if (v dot f) > 0:
    // v = (1-alpha) v + alpha |v| Fhat
    // |v| = length of v, Fhat = unit f
    // if more than DELAYSTEP since v dot f was negative:
    // increase timestep and decrease alpha
    if (vdotfall > 0.) {
      double tmp_me[3], tmp_all[3];
      scale1 = 1. - alpha;
      tmp_me[0] = tmp_me[1] = 0.;
      for (int i = 0; i < nvec; ++i) {
        tmp_me[0] += vvec[i] * vvec[i];
        tmp_me[1] += fperp[i] * fperp[i];
      }
      MPI_Allreduce(tmp_me, tmp_all,2,MPI_DOUBLE,MPI_SUM,world);
      vdotvall = tmp_all[0]; fdotfall = tmp_all[1];
      
      if (fdotfall < force_thr2){
	tmp_me[0] = tmp_me[1] = tmp_me[2] = 0.;
	for (int i = 0; i < nvec; ++i) {
	  tmp_me[0] += fvec[i] * fvec[i];
	  delr = xvec[i] - x0[i];
	  tmp_me[1] += delr * delr;
	  tmp_me[2] += fvec[i] * h[i];
	}
	MPI_Allreduce(tmp_me, tmp_all, 3, MPI_DOUBLE, MPI_SUM, world);
	ftotall = sqrt(tmp_all[0]); delr = sqrt(tmp_all[1]); fpar2all = tmp_all[2];
	fperp2 = 0.;
	for (int i = 0; i < nvec; ++i){
	  tmp = fvec[i] - fpar2all * h[i];
	  fperp2 +=  tmp * tmp;
	}
	MPI_Reduce(&fperp2, &fperp2all,1,MPI_DOUBLE,MPI_SUM,0,world);
	delE = ecurrent - eref;
	int m_perp= 0, trial = 0;
	if (me == 0){
	  if (fp1 && log_level) fprintf(fp1, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %8.4f %8.4f %6.3f " BIGINT_FORMAT "\n",
	      iter, delE, m_perp, trial, nlanc, ftotall, fpar2all, fperp2, egval, delr, hdotall, evalf);
	  if (screen ) fprintf(screen, "%8d %10.5f %3d %3d %5d %10.5f %10.5f %10.5f %8.4f %8.4f %6.3f " BIGINT_FORMAT "\n",
	      iter, delE, m_perp, trial, nlanc, ftotall, fpar2all, fperp2, egval, delr, hdotall, evalf);

	  idum = iter;
	  print_info(16);
	}
	return FTOL;
      }

      if (fdotfall == 0.) scale2 = 0.;
      else scale2 = alpha * sqrt(vdotvall/fdotfall);
      for (int i = 0; i < nvec; ++i) vvec[i] = scale1 * vvec[i] + scale2 * fperp[i];
      
      if ((iter - last_negative) > n_min) {
        dt = MIN(dt * f_inc, dtmax);
        alpha = alpha * f_alpha;
      }
      
    } else {
      last_negative = iter;
      dt = dt * f_dec;
      for (int i = 0; i < nvec; ++i) vvec[i] = 0.;
      alpha = alpha_start;
    }
    double **x = atom->x;
    double **v = atom->v;
    double *rmass = atom->rmass;
    double *mass = atom->mass;
    int *type = atom->type;
    double dtfm;
    // limit timestep so no particle moves further than dmax
    double dtvone = dt;
    double vmax = 0.;
    double dtv;
    for (int i = 0; i < atom->nlocal; i++) {
      vmax = MAX(fabs(v[i][0]),fabs(v[i][1]));
      vmax = MAX(vmax,fabs(v[i][2]));
      if (dtvone*vmax > dmax) dtvone = dmax/vmax;
    }
    MPI_Allreduce(&dtvone,&dtv,1,MPI_DOUBLE,MPI_MIN,world);

    double dtf = dtv * force->ftm2v;
    int n = 0;
    if (rmass) {
      for (int i = 0; i < atom->nlocal; ++i) {
        dtfm = dtf / rmass[i];
        x[i][0] += dt * v[i][0];
        x[i][1] += dt * v[i][1];
        x[i][2] += dt * v[i][2];
        v[i][0] += dtfm * fperp[n++];
        v[i][1] += dtfm * fperp[n++];
        v[i][2] += dtfm * fperp[n++];
      }
    } else {
      for (int i = 0; i < atom->nlocal; ++i) {
        dtfm = dtf / mass[type[i]];
        x[i][0] += dt * v[i][0];
        x[i][1] += dt * v[i][1];
        x[i][2] += dt * v[i][2];
        v[i][0] += dtfm * fperp[n++];
        v[i][1] += dtfm * fperp[n++];
        v[i][2] += dtfm * fperp[n++];
      }
    }
    eprevious = ecurrent;
    ecurrent = energy_force(1); ++evalf;
    artn_reset_vec();reset_x00();

    
  }

return MAXITER;
}

/* ---------------------------------------------------------------------------
 *  A little bit statistics on exit
 * -------------------------------------------------------------------------*/
void MinARTn::artn_final()
{
  if (me == 0){
    if (fp1){
      fprintf(fp1, "\n");
      fprintf(fp1, "#=========================================================================================\n");
      fprintf(fp1, "# Total number of ARTn attempts : %d\n", nattempt);
      fprintf(fp1, "# Number of new found saddle    : %d (%4.1f%% success)\n", sad_found, double(sad_found)/double(MAX(1,nattempt))*100.);
      fprintf(fp1, "# Number of accepted new saddle : %d (%4.1f%% acceptance)\n", sad_id, double(sad_id)/double(MAX(1,sad_found))*100.);
      fprintf(fp1, "# Overall mission success rate  : %4.1f%%\n", double(sad_id)/double(MAX(1,nattempt))*100.);
      fprintf(fp1, "# Number of new minimumi found  : %d\n", min_id-ref_0);
      fprintf(fp1, "# Number of accepted minima     : %d (%g%% acceptance)\n", ref_id-ref_0, double(ref_id-ref_0)/double(MAX(1,min_id-ref_0))*100.);
      fprintf(fp1, "# Number of force evaluation    : " BIGINT_FORMAT "\n", evalf);
      fprintf(fp1, "#=========================================================================================\n");
      fclose(fp1); fp1 = NULL;
    }
   
    if (fp2) fclose(fp2);

    if (screen){
      fprintf(screen, "\n");
      fprintf(screen, "#=========================================================================================\n");
      fprintf(screen, "# Total number of ARTn attempts : %d\n", nattempt);
      fprintf(screen, "# Number of new found saddle    : %d (%.1f%% success)\n", sad_found, double(sad_found)/double(MAX(1,nattempt))*100.);
      fprintf(screen, "# Number of accepted new saddle : %d (%.1f%% acceptance)\n", sad_id, double(sad_id)/double(MAX(1,sad_found))*100.);
      fprintf(screen, "# Overall mission success rate  : %.1f%%\n", double(sad_id)/double(MAX(1,nattempt))*100.);
      fprintf(screen, "# Number of new minimumi found  : %d\n", min_id-ref_0);
      fprintf(screen, "# Number of accepted minima     : %d (%g%% acceptance)\n", ref_id-ref_0, double(ref_id-ref_0)/double(MAX(1,min_id-ref_0))*100.);
      fprintf(screen, "# Number of force evaluation    : " BIGINT_FORMAT "\n", evalf);
      fprintf(screen, "#=========================================================================================\n");
    }
  }

  if (fp_sadlpress) fclose(fp_sadlpress);

  if (flog)   delete [] flog;
  if (fevent) delete [] fevent;
  if (fconfg) delete [] fconfg;
  if (groupname) delete [] groupname;

  if (glist)  delete [] glist;

  if (random)  delete random;
  if (dumpmin && !dumpmin_outside) delete dumpmin;
  if (dumpsad && !dumpsad_outside) delete dumpsad;
  if (dumpevent) delete dumpevent;

return;
}

/* -------------------------------------------------------------------------------------------------
 *  Write out related info
------------------------------------------------------------------------------------------------- */
void MinARTn::print_info(const int flag)
{
  if (flag == 0){
    if (fp1) fprintf(fp1, "\nMinimizing the initial configuration, id = %d ....\n", ref_id);
    if (screen) fprintf(screen, "\nMinimizing the initial configuration, id = %d ....\n", ref_id);

  } else if (flag == 1){
    if (fp1){
      if (log_level) fprintf(fp1, "  - Minimizer stop condition  : %s\n",  stopstr);
      fprintf(fp1, "  - Current (ref) energy (eV) : %.6f\n", ecurrent);
      fprintf(fp1, "  - Temperature   (eV)        : %.6f\n", temperature);
      fprintf(fp1, "  - Force two-norm final      : %.6f\n", ddum);
    }
    if (screen){
      if (log_level) fprintf(screen, "  - Minimizer stop condition  : %s\n",  stopstr);
      fprintf(screen, "  - Current (ref) energy (eV) : %.6f\n", ecurrent);
      fprintf(screen, "  - Temperature   (eV)        : %.6f\n", temperature);
      fprintf(screen, "  - Force two-norm final      : %.6f\n", ddum);
    }

  } else if (flag == 2){
      fprintf(fp2, "#  1       2        3       4      5      6      7       8        9   10    11      12      13        14        15         16         17       18     19        20        21        22       23      24\n");
      fprintf(fp2, "#Event   del-E   egv-sad   nsadl sad-dx sad-dy sad-dz  sad-dr    ref  sad   min   center   Eref      Emin     nMove        pxx        pyy      pzz    pxy      pxz       pyz     Efinal   status    dr\n");
      fprintf(fp2, "#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");

  } else if (flag == 3){
      fprintf(fp2, "#  1       2        3       4      5      6       7       8       9    10   11      12     13        14        15      16       17    18\n");
      fprintf(fp2, "#Event   del-E   egv-sad   nsadl sad-dx sad-dy  sad-dz  sad-dr   ref  sad   min   center   Eref      Emin     nMove   Efinal    status dr\n");
      fprintf(fp2, "#----------------------------------------------------------------------------------------------------------------------------------------\n");

  } else if (flag == 10){
    if (fp1){
      fprintf(fp1, "  Stage %d, search for the saddle from configuration %d\n", stage, ref_id);
      if (log_level){
        fprintf(fp1, "  ----------------------------------------------------------------------------------------------------\n");
        fprintf(fp1, "    Steps  E-Eref m_perp trial nlanc ftot       fpar        fperp     eigen       delr   evalf\n");
        fprintf(fp1, "  ----------------------------------------------------------------------------------------------------\n");
      }
      fflush(fp1);
    }
    if (screen){
      fprintf(screen, "  Stage %d, search for the saddle from configuration %d\n", stage, ref_id);
      fprintf(screen, "  ----------------------------------------------------------------------------------------------------\n");
      fprintf(screen, "    Steps  E-Eref m_perp trial nlanc ftot       fpar        fperp     eigen       delr   evalf\n");
      fprintf(screen, "  ----------------------------------------------------------------------------------------------------\n");
    }

  } else if (flag == 11){
    if (fp1){
      if (log_level) fprintf(fp1, "  ----------------------------------------------------------------------------------------------------\n");
      fprintf(fp1, "  Stage %d succeeded after %d iterations, continue searching based on eigen-vector.\n", stage, idum);
    }
    if (screen){
      fprintf(screen, "  ----------------------------------------------------------------------------------------------------\n");
      fprintf(screen, "  Stage %d succeeded after %d iterations, continue searching based on eigen-vector.\n", stage, idum);
    }

  } else if (flag == 12){
    if (fp1){
      if (log_level) fprintf(fp1, "  ----------------------------------------------------------------------------------------------------\n");
      fprintf(fp1, "  Stage %d failed, cannot get out of the harmonic well after %d steps.\n\n", stage, max_iter_basin);
    }
    if (screen){
      fprintf(screen, "  ----------------------------------------------------------------------------------------------------\n");
      fprintf(screen, "  Stage %d failed, cannot get out of the harmonic well after %d steps.\n\n", stage, max_iter_basin);
    }

  } else if (flag == 13){
    if (fp1){
      fprintf(fp1, "  Stage %d, converge to the saddle by using Lanczos\n", stage);
      if (log_level){
        fprintf(fp1, "  ----------------------------------------------------------------------------------------------------\n");
        fprintf(fp1, "    Iter   E-Eref m_perp trial nlanc ftot      fpar        fperp     eigen     delr    h.h' evalf\n");
        fprintf(fp1, "  ----------------------------------------------------------------------------------------------------\n");
      }
    }
    if (screen){
      fprintf(screen, "  Stage %d, converge to the saddle by using Lanczos\n", stage);
      fprintf(screen, "  ----------------------------------------------------------------------------------------------------\n");
      fprintf(screen, "    Iter   E-Eref m_perp trial nlanc ftot      fpar        fperp     eigen     delr    h.h' evalf\n");
      fprintf(screen, "  ----------------------------------------------------------------------------------------------------\n");
    }

  } else if (flag == 14){
    if (fp1){
      if (log_level) fprintf(fp1, "  ----------------------------------------------------------------------------------------------------\n");
      fprintf(fp1, "  Stage %d failed, the smallest eigen value is %g > %g\n", stage, egval, eigen_th_fail);
    }
    if (screen){
      fprintf(screen, "  ----------------------------------------------------------------------------------------------------\n");
      fprintf(screen, "  Stage %d failed, the smallest eigen value is %g > %g\n", stage, egval, eigen_th_fail);
    }

  } else if (flag == 15){
    if (fp1){
      if (log_level) fprintf(fp1, "  ----------------------------------------------------------------------------------------------------\n");
      fprintf(fp1, "  Stage %d failed, the distance to min-%d is %g < %g\n", stage, ref_id, ddum, disp_sad2min_thr);
    }
    if (screen){
      fprintf(screen, "  ----------------------------------------------------------------------------------------------------\n");
      fprintf(screen, "  Stage %d failed, the distance to min-%d is %g < %g\n", stage, ref_id, ddum, disp_sad2min_thr);
    }

  } else if (flag == 16){
    if (fp1){
      if (log_level) fprintf(fp1, "  ----------------------------------------------------------------------------------------------------\n");
      fprintf(fp1, "  Stage %d converged at a new saddle after %d iterations, dE = %g\n", stage, idum, ecurrent-eref);
    }
    if (screen){
      fprintf(screen, "  ----------------------------------------------------------------------------------------------------\n");
      fprintf(screen, "  Stage %d converged at a new saddle after %d iterations, dE = %g\n", stage, idum, ecurrent-eref);
    }

  } else if (flag == 17){
    if (fp1){
      if (log_level) fprintf(fp1, "  ----------------------------------------------------------------------------------------------------\n");
      fprintf(fp1, "  Stage %d failed, the max Lanczos steps %d reached.\n", stage, max_activat_iter);
    }
    if (screen){
      fprintf(screen, "  ----------------------------------------------------------------------------------------------------\n");
      fprintf(screen, "  Stage %d failed, the max Lanczos steps %d reached.\n", stage, max_activat_iter);
    }

  } else if (flag == 18){
    if (fp1) fprintf(fp1, "\nAttempt %d, new activation centered on atom %d", nattempt, that);
    if (screen) fprintf(screen, "\nAttempt %d, new activation centered on atom %d", nattempt, that);

  } else if (flag == 19){
    if (fp1) fprintf(fp1, " with total %d atoms. %d success till now.\n", idum, sad_id);
    if (screen) fprintf(screen, " with total %d atoms. %d success till now.\n", idum, sad_id);

  } else if (flag == 20){
    if (fp1) fprintf(fp1, "    The distance between new found saddle and min-%d is %g, > %g, acceptable.\n", ref_id, ddum, disp_sad2min_thr);
    if (screen) fprintf(screen, "    The distance between new found saddle and min-%d is %g, > %g, acceptable.\n", ref_id, ddum, disp_sad2min_thr);

  } else if (flag == 21){
    if (fp1) fprintf(fp1, "    The distance between new saddle and min-%d is %g, < %g, rejected.\n", ref_id, ddum, disp_sad2min_thr);
    if (screen) fprintf(screen, "    The distance between new saddle and min-%d is %g, < %g, rejected.\n", ref_id, ddum, disp_sad2min_thr);

  } else if (flag == 30){
    if (fp1) fprintf(fp1, "  Stage %d, push back the saddle to confirm if it is linked with min-%d\n", stage, ref_id);
    if (screen) fprintf(screen, "  Stage %d, push back the saddle to confirm if it is linked with min-%d\n", stage, ref_id);

  } else if (flag == 31){
    if (fp1) {
      fprintf(fp1, "    - Current   energy (eV)     : %.6f\n", ecurrent);
      fprintf(fp1, "    - Reference energy (eV)     : %.6f\n", eref);
      fprintf(fp1, "    - Force two-norm final      : %.6f\n", ddum);
      if (log_level) fprintf(fp1, "    - Minimizer stop condition  : %s\n",  stopstr);
      if (log_level) fprintf(fp1, "    - # of force evaluations    : %d\n", neval);
    }
    if (screen) {
      fprintf(screen, "    - Current   energy (eV)     : %.6f\n", ecurrent);
      fprintf(screen, "    - Reference energy (eV)     : %.6f\n", eref);
      fprintf(screen, "    - Force two-norm final      : %.6f\n", ddum);
      if (log_level) fprintf(screen, "    - Minimizer stop condition  : %s\n",  stopstr);
      if (log_level) fprintf(screen, "    - # of force evaluations    : %d\n", neval);
    }

  } else if (flag == 32){
    if (fp1) fprintf(fp1, "  Stage %d failed, |Ecurrent - Eref| = %g > %g, reject the new saddle.\n", stage, ddum, max_ener_tol);
    if (screen) fprintf(screen, "  Stage %d failed, |Ecurrent - Eref| = %g > %g, reject the new saddle.\n", stage, ddum, max_ener_tol);

  } else if (flag == 33){
    if (fp1) fprintf(fp1, "  Stage %d succeeded, dr = %g, drmax_atom = %g, tolerance = %g, accept the new saddle.\n", stage, ddum, tmpddum, max_disp_tol);
    if (screen) fprintf(screen, "  Stage %d succeeded, dr = %g, drmax_atom = %g, tolerance = %g, accept the new saddle.\n", stage, ddum, tmpddum, max_disp_tol);

  } else if (flag == 34){
    if (fp1) fprintf(fp1, "  Stage %d failed, dr = %g, drmax_atom = %g, tolerance = %g, reject the new saddle.\n", stage, ddum, tmpddum, max_disp_tol);
    if (screen) fprintf(screen, "  Stage %d failed, dr = %g, drmax_atom = %g, tolerance = %g, reject the new saddle.\n", stage, ddum, tmpddum, max_disp_tol);

  } else if (flag == 40){
    if (fp1) fprintf(fp1, "  Stage %d, further relax the newly found sad-%d ...\n", stage, sad_id);
    if (screen) fprintf(screen, "  Stage %d, further relax the newly found sad-%d ...\n", stage, sad_id);

  } else if (flag == 41){
    if (fp1){
      fprintf(fp1, "    The new sad-%d is now converged as:\n", sad_id);
      fprintf(fp1, "      - Current energy  (eV)      : %.6f\n", ecurrent);
      fprintf(fp1, "      - Energy  barrier (eV)      : %.6f\n", delE);
      if (log_level) fprintf(fp1, "      - Norm2  of total force     : %lg\n", sqrt(ddum));
      if (log_level) fprintf(fp1, "      - Minimizer stop condition  : %s\n",  stopstr);
      if (log_level) fprintf(fp1, "      - # of force evaluations    : %d\n", neval);
    }
    if (screen){
      fprintf(screen, "    The new sad-%d is now converged as:\n", sad_id);
      fprintf(screen, "      - Current energy  (eV)      : %.6f\n", ecurrent);
      fprintf(screen, "      - Energy  barrier (eV)      : %.6f\n", delE);
      if (log_level) fprintf(screen, "      - Norm2  of total force     : %lg\n", sqrt(ddum));
      if (log_level) fprintf(screen, "      - Minimizer stop condition  : %s\n",  stopstr);
      if (log_level) fprintf(screen, "      - # of force evaluations    : %d\n", neval);
    }

  } else if (flag == 50){
    if (fp1) fprintf(fp1, "  Stage %d, push over the new saddle, Enew: %g; Eref= %g. Relaxing...\n", stage, ecurrent, eref);
    if (screen) fprintf(screen, "  Stage %d, push over the new saddle, Enew: %g; Eref= %g. Relaxing...\n", stage, ecurrent, eref);

  } else if (flag == 51){
    if (fp1){
      fprintf(fp1, "    Relaxed to a nearby minimum to sad-%d\n", sad_id);
      fprintf(fp1, "    - Current  min  energy (eV) : %.6f\n", ecurrent);
      fprintf(fp1, "    - Force two-norm final      : %.6f\n", ddum);
      if (log_level)fprintf(fp1, "      - Reference     energy (eV) : %.6f\n", eref);
      if (log_level) fprintf(fp1, "      - Minimizer stop condition  : %s\n",  stopstr);
    }
    if (screen){
      fprintf(screen, "    Relaxed to a nearby minimum to sad-%d\n", sad_id);
      fprintf(screen, "    - Current  min  energy (eV) : %.6f\n", ecurrent);
      fprintf(screen, "    - Force two-norm final      : %.6f\n", ddum);
      if (log_level) fprintf(screen, "      - Reference     energy (eV) : %.6f\n", eref);
      if (log_level) fprintf(screen, "      - Minimizer stop condition  : %s\n",  stopstr);
    }

  } else if (flag == 60){
    if (fp1 && log_level) fprintf(fp1, "      - Distance to min-%8.8d  : %g\n", ref_id, ddum);
    if (screen && log_level) fprintf(screen, "      - Distance to min-%8.8d  : %g\n", ref_id, ddum);

  } else if (flag == 61){
    if (fp1) fprintf(fp1, "  Stage %d done, the new min (E= %g) of ID  %d is accepted by Metropolis.\n", stage, ecurrent, min_id);
    if (screen) fprintf(screen, "  Stage %d done, the new min (E= %g) of ID %d is accepted by Metropolis.\n", stage, ecurrent, min_id);

  } else if (flag == 62){
    if (fp1) fprintf(fp1, "  Stage %d done, the new min (E= %g) of ID %d is rejected by Metropolis.\n", stage, ecurrent, min_id);
    if (screen) fprintf(screen, "  Stage %d done, the new min (E= %g) of ID %d is rejected by Metropolis.\n", stage, ecurrent, min_id);

  }
return;
}

/* -------------------------------------------------------------------------------------------------
 * Converge the saddle point by using SD method; the force parallel to the eigenvector corresponding
 * to the smallest egval is reversed.
------------------------------------------------------------------------------------------------- */
void MinARTn::sad_converge(int maxiter)
{
  ++stage;
  if (me == 0) print_info(40);

  neval = 0;
  int i,fail;
  double edf, edf_all;

  lanczos(flag_egvec, 1, num_lancz_vec_c);
  // initialize working vectors
  edf = 0.;
  for (i =0; i < nvec; ++i) edf += egvec[i] * fvec[i];
  MPI_Allreduce(&edf, &edf_all,1,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < nvec; ++i) h[i] = fvec[i] - 2.*edf_all*egvec[i];

  stop_condition = MAXITER;

  for (int iter = 0; iter < maxiter; ++iter) {
    // line minimization along h from current position x
    // h = downhill gradient direction
    eprevious = ecurrent;
    fail = (this->*linemin)(ecurrent,alpha_final);
    if (fail) {stop_condition = fail; break;}

    // function evaluation criterion
    if (neval >= update->max_eval) {stop_condition = MAXEVAL; break;}

    // energy tolerance criterion
    if (fabs(ecurrent-eprevious) < update->etol * 0.5*(fabs(ecurrent)
    + fabs(eprevious) + EPS_ENERGY)) {stop_condition = ETOL; break;}

    // force tolerance criterion
    double fdotf = fnorm_sqr();
    if (fdotf < update->ftol*update->ftol) {stop_condition = FTOL; break;}

    // set new search direction h to f = -Grad(x)
    lanczos(flag_egvec, 1, num_lancz_vec_c);
    edf = 0.;
    for (i = 0; i < nvec; ++i) edf += egvec[i] * fvec[i];
    MPI_Allreduce(&edf, &edf_all,1,MPI_DOUBLE,MPI_SUM,world);

    for (i = 0; i < nvec; ++i) h[i] = fvec[i] - 2.*edf_all*egvec[i];
  }

  evalf += neval;
  stopstr = stopstrings(stop_condition);

  // output minimization information
  delE = ecurrent-eref;
  ddum = fnorm_sqr();
  if (me == 0) print_info(41);

return;
}
/* ---------------------------------------------------------------------------------------------- */
int MinARTn::SD_min_converge(int maxiter, const int flag)
{
  neval = 0;
  int i,m,n,fail,ntimestep;
  double fdotf;
  double *fatom,*hatom;

  // initialize working vectors

  for (i = 0; i < nvec; i++) h[i] = fvec[i];
  if (nextra_atom)
    for (m = 0; m < nextra_atom; m++) {
      fatom = fextra_atom[m];
      hatom = hextra_atom[m];
      n = extra_nlen[m];
      for (i = 0; i < n; i++) hatom[i] = fatom[i];
    }
  if (nextra_global)
    for (i = 0; i < nextra_global; i++) hextra[i] = fextra[i];

  for (int iter = 0; iter < maxiter; iter++) {
    //ntimestep = ++update->ntimestep;
    niter++;

    // line minimization along h from current position x
    // h = downhill gradient direction

    eprevious = ecurrent = energy_force(1); ++neval;
    fail = (this->*linemin)(ecurrent,alpha_final);
    if (fail) return fail;

    // function evaluation criterion

    if (neval >= update->max_eval) return MAXEVAL;

    // energy tolerance criterion

    if (fabs(ecurrent-eprevious) <
        update->etol * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS_ENERGY))
      return ETOL;

    // force tolerance criterion

    fdotf = fnorm_sqr();
    if (fdotf < update->ftol*update->ftol) return FTOL;

    // set new search direction h to f = -Grad(x)

    for (i = 0; i < nvec; i++) h[i] = fvec[i];
    if (nextra_atom)
      for (m = 0; m < nextra_atom; m++) {
        fatom = fextra_atom[m];
        hatom = hextra_atom[m];
        n = extra_nlen[m];
        for (i = 0; i < n; i++) hatom[i] = fatom[i];
      }
    if (nextra_global)
      for (i = 0; i < nextra_global; i++) hextra[i] = fextra[i];
    if (flag == 2) reset_coords();
    else if (flag == 1) reset_x00();

    // output for thermo, dump, restart files
    //if (output->next == ntimestep) {
    //  timer->stamp();
    //  output->write(ntimestep);
    //  timer->stamp(Timer::OUTPUT);
    //}
  }
  return MAXITER;
   
}
int MinARTn::min_converge_fire(int maxiter){
  double dt = update->dt;
  const int n_min = 5;
  const double f_inc = 1.1;
  const double f_dec = 0.5;
  const double alpha_start = 0.1;
  const double f_alpha = 0.99;
  const double  TMAX = 10.;
  const double dtmax = TMAX * dt;
  double vdotf, vdotfall;
  double vdotvall;
  double fdotfall;
  double fdoth, fdothall;
  double scale1, scale2;
  double alpha;
  int last_negative = 0;
  int ntimestep;
  neval = 0;

  double force_thr2 = force_th_perp_sad*force_th_perp_sad;

  double *vvec = atom->v[0];
  for (int i = 0; i < nvec; ++i) vvec[i] = 0.;

  alpha = alpha_start;
  for (int iter = 0; iter < maxiter; ++iter){
    ntimestep = ++update->ntimestep;
    vdotf = 0.;
    vvec = atom->v[0];
    for (int i = 0; i < nvec; ++i) vdotf += vvec[i] * fvec[i];
    MPI_Allreduce(&vdotf,&vdotfall,1,MPI_DOUBLE,MPI_SUM,world);

    // if (v dot f) > 0:
    // v = (1-alpha) v + alpha |v| Fhat
    // |v| = length of v, Fhat = unit f
    // if more than DELAYSTEP since v dot f was negative:
    // increase timestep and decrease alpha
    if (vdotfall > 0.) {
      double tmp_me[2], tmp_all[2];
      scale1 = 1. - alpha;
      tmp_me[0] = tmp_me[1] = 0.;
      for (int i = 0; i < nvec; ++i) {
	tmp_me[0] += vvec[i] * vvec[i];
	tmp_me[1] += fvec[i] * fvec[i];
      }
      MPI_Allreduce(tmp_me, tmp_all,2,MPI_DOUBLE,MPI_SUM,world);
      vdotvall = tmp_all[0]; fdotfall = tmp_all[1];

      if (fdotfall == 0.) scale2 = 0.;
      else scale2 = alpha * sqrt(vdotvall/fdotfall);
      for (int i = 0; i < nvec; ++i) vvec[i] = scale1 * vvec[i] + scale2 * fvec[i];

      if ((iter - last_negative) > n_min) {
	dt = MIN(dt * f_inc, dtmax);
	alpha = alpha * f_alpha;
      }

    } else {
      last_negative = iter;
      dt = dt * f_dec;
      for (int i = 0; i < nvec; ++i) vvec[i] = 0.;
      alpha = alpha_start;
    }
    double **x = atom->x;
    double **v = atom->v;
    double *rmass = atom->rmass;
    double *mass = atom->mass;
    int *type = atom->type;
    double dtfm;
    // limit timestep so no particle moves further than dmax
    double dtvone = dt;
    double vmax = 0.;
    double dtv;
    for (int i = 0; i < atom->nlocal; i++) {
      vmax = MAX(fabs(v[i][0]),fabs(v[i][1]));
      vmax = MAX(vmax,fabs(v[i][2]));
      if (dtvone*vmax > dmax) dtvone = dmax/vmax;
    }
    MPI_Allreduce(&dtvone,&dtv,1,MPI_DOUBLE,MPI_MIN,world);

    double dtf = dtv * force->ftm2v;
    int n = 0;
    if (rmass) {
      for (int i = 0; i < atom->nlocal; ++i) {
	dtfm = dtf / rmass[i];
	x[i][0] += dt * v[i][0];
	x[i][1] += dt * v[i][1];
	x[i][2] += dt * v[i][2];
	v[i][0] += dtfm * fvec[n++];
	v[i][1] += dtfm * fvec[n++];
	v[i][2] += dtfm * fvec[n++];
      }
    } else {
      for (int i = 0; i < atom->nlocal; ++i) {
	dtfm = dtf / mass[type[i]];
	x[i][0] += dt * v[i][0];
	x[i][1] += dt * v[i][1];
	x[i][2] += dt * v[i][2];
	v[i][0] += dtfm * fvec[n++];
	v[i][1] += dtfm * fvec[n++];
	v[i][2] += dtfm * fvec[n++];
      }
    }
    eprevious = ecurrent;
    ecurrent = energy_force(1); ++neval;
    artn_reset_vec();reset_x00();
    // energy tolerance criterion
    // only check after DELAYSTEP elapsed since velocties reset to 0

    if (update->etol > 0.0 && iter-last_negative > n_min) {
      if (fabs(ecurrent-eprevious) <
	  update->etol * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS_ENERGY))
	return ETOL;
    }
    // force tolerance criterion

    if (update->ftol > 0.0) {
      double fdotf = fnorm_sqr();
      if (fdotf < update->ftol*update->ftol) return FTOL;
    }

    //if (output->next == ntimestep) {
    //  timer->stamp();
    //  output->write(ntimestep);
    //  timer->stamp(Timer::OUTPUT);
    //}


  }

  return MAXITER;

}

/*-------------------------------------------------------------------
 *  read displacement for initial kick
 *  @file: the dump direction file to open, native and orthogonal 
 *  coordinate format suported.
 *  @delpos: displacement vector (3*nlocal element)
 *  Important: atom_modify map hash/array shoud be used in lammps input
 *  script.
 * ----------------------------------------------------------------*/
void MinARTn::read_dump_direction(char * file, double * delpos){
  if (me == 0 && screen){
    fprintf(screen,"\nNOTE: atom_modify map hash or atom_modify map array should be used in the LAMMPS input file when using ARTn fdump feature\n");
  }
  FILE * fp;
  char str[MAXLINE], oneline[MAXLINE], *token;
  int flag_scale = 0;
  int flag_tilt = 0;
  int id,type,ilocal;
  bigint natoms; 
  double * dumppos;
  double lox,hix,loy,hiy,loz,hiz,lx,ly,lz;
  double xy,xz,yz; 
  lox = hix = loy = hiy = loz = lx = ly = lz = 0.0;
  double loxb,hixb,loyb,hiyb,lozb,hizb;
  loxb = hixb = loyb = hiyb = lozb = hizb = 0.0;
  if (file == NULL){
    error->one(FLERR,"Dump direction file not set.");
    return;
  }
  if (me == 0){
    fp = fopen(file,"r"); 	// open file
    if (fp == NULL){
       sprintf(str, "Cannot open ARTn dump direction file: %s for reading", file);
       error->one(FLERR,str);
    }
    for(int i = 0;i < 3;++i)fgets(oneline,MAXLINE,fp);
    fgets(oneline,MAXLINE,fp);
    sscanf(oneline,BIGINT_FORMAT,&natoms);
  }
  MPI_Bcast(&natoms, 1, MPI_INT,0,world);
  // 3N + 1 to set id to position x[id*3], e.g. id = 1 are stored in dumppos[1*3]
  dumppos = new double [3*natoms+3];
  if (me == 0){
    fgets(oneline,MAXLINE,fp);
    if(oneline[strlen("ITEM: BOX BOUNDS p")] == 'p'){
      fgets(oneline,MAXLINE,fp);
      sscanf(oneline,"%lg %lg",&lox,&hix);
      lx = hix - lox;
      fgets(oneline,MAXLINE,fp);
      sscanf(oneline,"%lg %lg",&loy,&hiy);
      ly = hiy - loy;
      fgets(oneline,MAXLINE,fp);
      sscanf(oneline,"%lg %lg",&loz,&hiz);
      lz = hiz - loz;
    }else{
      flag_tilt = 1;
      fgets(oneline,MAXLINE,fp);
      sscanf(oneline,"%lg %lg %lg",&loxb,&hixb,&xy);
      fgets(oneline,MAXLINE,fp);
      sscanf(oneline,"%lg %lg %lg",&loyb,&hiyb,&xz);
      fgets(oneline,MAXLINE,fp);
      sscanf(oneline,"%lg %lg %lg",&lozb,&hizb,&yz);
      lox = loxb - MIN(0.0, MIN(xy, MIN(xz, xy + xz))); 
      hix = hixb - MAX(0.0, MAX(xy, MAX(xz, xy + xz)));
      loy = loyb - MIN(0.0, yz);
      hiy = hiyb - MAX(0.0, yz);
      loz = lozb;
      hiz = hizb;
      lx = hix - lox;
      ly = hiy - loy;
      lz = hiz - loz;
    }
    fgets(oneline,MAXLINE,fp);
    if (oneline[strlen("ITEM: ATOMS id type x")] == 's') flag_scale = 1;
    for (bigint i = 0; i < natoms; ++i){
      fgets(oneline,MAXLINE,fp);
      sscanf(oneline,"%i %i", &id, &type);
      sscanf(oneline,"%i %i %lg %lg %lg", &id, &type, dumppos+(3*id), dumppos+(3*id)+1,dumppos+(3*id)+2);
      if (flag_scale){
	if(!flag_tilt){
	  dumppos[3*id]   = dumppos[3*id] * lx + lox;
	  dumppos[3*id+1] = dumppos[3*id+1] * ly + loy;
	  dumppos[3*id+2] = dumppos[3*id+2] * lz + loz;
	}else{
	  dumppos[3*id]   = dumppos[3*id] * lx + lox + xy * dumppos[3*id+1] + xz * dumppos[3*id+2];
	  dumppos[3*id+1] = dumppos[3*id+1] * ly + loy + yz * dumppos[3*id+2];
	  dumppos[3*id+2] = dumppos[3*id+2] * lz + loz;
	}
      }
    }
  }
  MPI_Bcast(dumppos,3*natoms+3, MPI_DOUBLE,0,world);
  double dx,dy,dz;
  int nlocal = atom->nlocal;
  xvec = atom->x[0];
  tagint tag;
  for (bigint i = 1; i <= natoms; ++i){
    tag = static_cast<tagint>(i);
    ilocal = atom->map(tag);
    if (ilocal >= 0 && ilocal < nlocal){
      int n = ilocal*3;
      dx = dumppos[i*3] - xvec[n];
      dy = dumppos[i*3+1] - xvec[n+1];
      dz = dumppos[i*3+2] - xvec[n+2];
      domain->minimum_image(dx, dy, dz);
      delpos[n] = dx;
      delpos[n+1] = dy;
      delpos[n+2] = dz;
    }
  }
  if (dumppos) delete []dumppos;
  if (me == 0) fclose(fp);
}
void MinARTn::check_freeze(){
  int nlocal = atom->nlocal;
  xvec = atom->x[0];
  tagint tag;
  double tmp;
  double sum = 0.;
  for (int i = 0; i < nlocal; ++i){
    if(!(groupbit & atom->mask[i])){
      int n = 3*i;
      int id = atom->tag[i];
      tmp = xvec[n]-x00[n];
      if (tmp > ZERO){
	fprintf(screen, "\nAtom %i displacement %f > ZERO\n", id, sqrt(tmp));
      }
      sum += tmp*tmp;
      tmp = xvec[n+1]-x00[n+1];
      sum += tmp*tmp;
      tmp = xvec[n+2]-x00[n+2];
      sum += tmp*tmp;
    }
  } 
  double sumall;
  MPI_Allreduce(&sum,&sumall,1,MPI_DOUBLE, MPI_SUM, world);
  sumall = sqrt(sumall);
  if(me == 0) fprintf(screen,"\nDisplacemnt: %f\n", sumall);
  //return sumall;
}
