#ifdef MINIMIZE_CLASS

MinimizeStyle(artn,MinARTn)

#else

#ifndef ARTN_H
#define ARTN_H

#include "dump_atom.h"
#include "random_park.h"
#include "min_linesearch.h"

using namespace std;

namespace LAMMPS_NS {

class MinARTn: public MinLineSearch {
public:
  MinARTn(class LAMMPS *);
  int iterate(int);

private:
  void set_defaults();			
  void artn_init();
  void read_control();

  void random_kick();
  void read_dump_direction(char *dumpfile, double *delpos); 	// allocate direction outside the function

  int find_saddle();
  int new_find_saddle();
  int check_sad2min();
  void analysis_saddle();
  int min_converge(int, const int);
  int min_converge_fire(int);
  int SD_min_converge(int, const int);
  void sad_converge(int);

  int push_back_sad();
  void push_down();
  void metropolis();

  void artn_reset_vec();
  void reset_coords();
  void reset_x00();
  int lanczos(bool, int, int);
  int min_perp_fire(int);
  int new_min_perp_fire(int);

  int sad_found;
  int ref_id, min_id, sad_id, ref_0;
  DumpAtom * dumpmin;
  DumpAtom * dumpsad;
  DumpAtom * dumpevent;
  Compute *pressure;

  int me, np;
  //int ntimestep;
  int iatom;               // index of current atom to find saddlepoint when using events_per_atom
  bigint evalf;
  bool flag_egvec;

  double eref, delE;

  double egval;
  double *egvec;
  double *x0tmp;
  double *x00;
  double *fperp;

  int seed;
  RanPark *random;

  // global control
  int max_conv_steps;
  double temperature;      // Fictive temperature, if negative always reject the event

  // for art
  int nattempt, stage;
  int fire_lanczos_every;
  int max_num_events;      // Maximum number of events
  int max_activat_iter;    // Maximum number of iteractions for reaching the saddle point
  double increment_size;   // Overall scale for the increment moves
  int use_fire;            // use FIRE to do minimization in the perpendicular direction
  int min_fire;            // use FIRE to do minimization both in push back & push forward
  int events_per_atom;     // Find designed number of ARTn events on each atom, set to 0 to shutoff this method
  int flag_push_back;      // if 1, will push back saddle point to check if it connect with the minimum
  int flag_push_over;      // if 1, wiil push over saddle point to reach another minimum
  int flag_relax_sad;      // further relax to the newly found saddle
  double max_disp_tol;     // tolerance displacement between ref and push-back that can claim saddle is indeed linked to ref (negative value corresponds to the max atom displacement)
  double max_ener_tol;     // energy tolerance 
  int flag_press;          // Pressure will be calculated.
  int flag_sadl_press;     // Saddle point's pressure will be calculated.
  double atom_disp_thr;    // threshold to identify whether an atom is displaced or not

  int groupbit, ngroup;    // group bit & # of atoms for initial kick
  int that, *glist;        // ID and list for kick
  char *groupname;         // group name for initial kick
  double cluster_radius;   // radius for kick; <0, all; ==0, single; >0, cluster

  // for harmonic well
  int flag_dump_direction; // use dump direction file as the initial kick direction
  int flag_deformation_gradient; // use deformation gradient as the initial kick direction
  double deformation_gradient[9]; // deformation gradient, in the component order 11 12 13 21 22 23 31 32 33
  double random_kick_factor; // add a random disturbation to the dump direction.
  double init_step_size;   // Size of initial displacement
  double basin_factor;     // Factor multiplying Increment_Size for leaving the basin
  int max_perp_move_h;     // Maximum number of perpendicular steps leaving basin
  int min_num_ksteps;      // Min. number of ksteps before calling lanczos
  double eigen_th_well;    // Eigenvalue threshold for leaving basin
  int max_iter_basin;      // Maximum number of iteraction for leaving the basin (kter)
  double force_th_perp_h;  // Perpendicular force threhold in harmonic well

  // for lanczos
  int num_lancz_vec_h;     // Number of vectors included in lanczos procedure in the Harmonic well
  int num_lancz_vec_c;     // Number of vectors included in lanczos procedure in convergence
  double del_disp_lancz;   // Step of the numerical derivative of forces in lanczos
  double eigen_th_lancz;   // Eigen_threhold for lanczos convergence

  // for convergence
  int    SD_steps;
  int    flag_test;
  double para_factor;
  int fire_output_every;
  double force_th_saddle;  // Threshold for convergence at saddle point
  double disp_sad2min_thr; // minimum distance between saddle and original min
  double push_over_saddle; // Fraction of displacement over the saddle
  double eigen_th_fail;    // the eigen cutoff for failing in searching saddle point
  int    max_perp_moves_c; // Maximum number of perpendicular steps approaching saddle point
  double force_th_perp_sad;// Perpendicular force threhold approaching saddle point
  int    conv_perp_inc;

  // for output
  FILE *fp1, *fp2, *fp_sadlpress;
  char *flog, *fevent, *fconfg, *c_fsadpress, *fdump_direction;
  int log_level;           // 1, all; 0, main
  int print_freq;          // default 1
  int dump_min_every;      // dump min configuration every # step
  int dump_sad_every;      // dump sadl configuration every # step
  int dump_event_every;
  int idum;
  double ddum;
  double tmpddum;

  void artn_final();
  void print_info(const int);

  // center of mass info
  int groupall;
  double com0[3], com[3];
  double masstot;
};
}
#endif
#endif
