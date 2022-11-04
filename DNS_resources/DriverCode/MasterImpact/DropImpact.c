// Main developer: Radu Cimpeanu
// Contributor: Thomas Sykes
// Date: 26/09/2022

// Mollify density and viscosity jumps
#define FILTERED
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))

#include "axi.h"                     // axisymmetric geometry
#include "navier-stokes/centered.h"  // solve NS equations (centered formulation)
#include "two-phase.h"               // two fluid phases
#include "tension.h"                 // include surf tension between phases
// Curvature is computed via a geometric implementation ~ PLIC (see curvature.h)
#include "vof.h"                     // solve using VoF method
#include "fractions.h"               // initially define fraction of fluid phases
#include "view.h"                    // need to make the animations
#include "tag.h"                     // helps track droplet properties
// Tag indexes connected fluid regions (i.e. droplets), used to track satellite
// droplet properties, and remove satellites near the free slip boundary
#include "draw.h"                    // visualisation helper
#include "tracer.h"                  // passive scalar tracer
#include "maxruntime.h"              // kill simulation before ARC/HPC does


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 INITIALISATION
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Initialise the passive scalar:
scalar dTracer[];
scalar * tracers = {dTracer};

// Dimensional quantities (fluids at 20 degrees C):
#define rhoLiquid   785.0              // Liquid phase density (kg/m^3)
#define rhoGas      1.21               // Gas phase density (kg/m^3)

#define muLiquid    1.10e-3            // Liquid dynamic viscosity (kg/ms)
#define muGas       1.81e-5            // Gas dynamic viscosity(kg/ms)

#define sig         0.0222             // Surface tension (N/m)

#define g_accel     9.81               // Gravitational acceleration (m/s^2)

#define dRadius     1.02e-3            // Drop radius (m)

// Impact velocity (controlled by the Weber number):
#define v_init      (pow((sig*ND_Weber)/(rhoLiquid*2*dRadius),0.5))

// Dimensionless quantities (key groupings defined in main below):
#define rho_ratio   (rhoGas/rhoLiquid) // Density ratio
#define mu_ratio    (muGas/muLiquid)   // Viscosity ratio

// Computational box size (in diameters - remember axisymmetric):
#define domainSize  1.5

// The acceleration term a in the centered form of the NS equations, which has
// been divided through by density (i.e. the body force term):
face vector av[];

// Initialise file pointers (files themselves will be created below):
FILE * fp_stats;
FILE * fp_droplets;

// Initialise the dimensionless numbers. Those that are not populated by/
// dependent on arguments are set here. Otherwise just initialise empty:
double ND_Weber;        // Weber number (non-dimensionalises surface tension)
double ND_Reynolds;     // Reynolds number (non-dimensionalises viscosity)
double ND_Froude;       // Froude number (non-dimensionalises gravity)

double ND_thickness;    // Pool depth (needed for initial condition)
double filmHeight;      // Position of free surface in centred coordinate system

int minLevel = 6;       // Minimum refinement level (2^n across the box)
int maxLevel;           // Maximum refinement level (2^n across the box)

double tEnd;            // Simulation end time (total time to simulate)


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 BOUNDARY CONDITIONS
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
The simulation environment is orientated such that the axis of symmetry is on the
bottom (y=0), whilst the pool base is on the left. The coordinate system is
shifted below to ensure that x=0
*/

// Bottom of Pool = LEFT of domain
// No slip, no permeability
u.n[left]=dirichlet(0.);
u.t[left]=dirichlet(0.); 
p[left]=neumann(0.);
pf[left]=neumann(0.);

// Side of Pool = TOP of domain
// Free slip, no permeability
u.n[top]=dirichlet(0.); // Impermeable
u.t[top]=neumann(0.);   // Slip

// Above the Pool = RIGHT of domain
// Outflow
u.n[right]=neumann(0.);
p[right]=dirichlet(0.);
pf[right]=dirichlet(0.);


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 MAIN FUNCTION
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* Everything within main is either defined above, or is defined in the included
files. For example, rho* and mu* are defned in */

int main(int argc, char * argv[]) {
  // argc: number of arguments passed in from the execution environment
  // argv: pointer to the first element of an array of pointers to those arguments

  // Function to call for maximum runtime, which is an optional arguement (-m)
  maxruntime (&argc,argv);

  // Pull in the arguments (atof for float; atoi for integer)
  ND_Weber     = atof(argv[1]); // We number
  ND_thickness = atof(argv[2]); // Film thickness (relative to drop diameter)
  maxLevel     = atoi(argv[3]); // Maximum refinement level
  tEnd         = atof(argv[4]); // Maximum simulation run time (dimensionless!)
  
  // Print all arguments to the terminal:
  fprintf(stdout, "Arguments:\n"); fflush(stdout);
  fprintf(stdout, "We = %0.1f \n", ND_Weber); fflush(stdout);
  fprintf(stdout, "h = %0.2f \n", ND_thickness); fflush(stdout);
  fprintf(stdout, "Max. Level = %d \n", maxLevel); fflush(stdout);
  fprintf(stdout, "tend = %0.3f \n", tEnd); fflush(stdout);
 
  // Determine other dimensionless numbers from the given arguments
  ND_Reynolds = (rhoLiquid*v_init*2.0*dRadius)/muLiquid;
  ND_Froude   = v_init/pow(2.0*dRadius*g_accel,0.5);
  
  // Initialise the grid at level 2^7 (note: *1 << n* is C for $2^n$)
  // Rule of thumb: Use 4-5 levels below the maximum refinement
  init_grid(1 << 7);
  
  // Set the size of the simulation box, as initialised above
  size(domainSize);

  // Shift the coordinate system so that the origin is in the middle of the
  // bottom of the domain (i.e. middle of the axis of symmetry)
  origin(-0.5*domainSize, 0.0);

  // Make folders to store the outputs
  mkdir("Slices", 0700);
  mkdir("Animations", 0700);
  mkdir("Interfaces", 0700);

  // Print dimensionless numbers to the terminal
  fprintf(stdout, "Derived values:\n"); fflush(stdout);
  fprintf(stdout, "Reynolds number = %0.6f \n", ND_Reynolds); fflush(stdout);
  fprintf(stdout, "Weber number = %0.6f \n", ND_Weber); fflush(stdout);
  fprintf(stdout, "Froude number = %0.6f \n", ND_Froude); fflush(stdout);

  // Print impact velocity to the terminal
  fprintf(stdout, "impact velocity = %0.3f \n", v_init); fflush(stdout);

  // Set dimensionless densities and viscosities
  
  rho1 = 1.;                     // = muLiquid/muLiquid
  rho2 = rho_ratio;              // = rhoGas/rhoLiquid
  
  mu1  = 1./ND_Reynolds;         // = muLiquid/(rhoLiquid*v_init*2.0*dRadius)
  mu2  = mu_ratio*mu1;           // = (muGas/muLiquid)*mu1 = muGas/(...)

  // Print fluid properties to the terminal
  fprintf(stdout, "rho1 = %0.9f \n", rho1); fflush(stdout);
  fprintf(stdout, "rho2 = %0.9f \n", rho2); fflush(stdout);
  fprintf(stdout, "mu1 = %0.9f \n", mu1); fflush(stdout);
  fprintf(stdout, "mu2 = %0.9f \n", mu2); fflush(stdout);
  
  // Acceleration term, on the cell face (initialised above):
  a = av;

  // Dimensionless surface tension:
  f.sigma = 1./ND_Weber;

  // Pointer of the file to save stats (appends, generates file if required):
  {
    char name[200];
    sprintf(name, "logstats.dat");
    fp_stats = fopen(name, "a");
  }

  // Pointer of the file to droplet info (appends, generates file if required):
  {
    char name[200];
    sprintf(name, "logdroplets.dat");
    fp_droplets = fopen(name, "a");
  }

  // Linear (multigrid) solver settings (system/fvsolutions equiv. of OpenFOAM):
  DT        = 1e-3;  // Max timestep
  NITERMIN  = 1;     // Min number of linear solver iterations (default 1)
  NITERMAX  = 200;   // Max number of linear solver iterations (default 100)
  TOLERANCE = 1e-4;  // Linear solver tolerance (default 1e-3)
  
  // Main execution command
  run();

  // Close the files opened for writing just above
  fclose(fp_stats);
  fclose(fp_droplets);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 EVENTS
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Initialisation for the events
scalar omega[], viewingfield[], mylevel[], myrho[], mymu[];

// INCLUDE GRAVITY The Froude number provides the non-dimensionalisation:
event acceleration (i++) {
  foreach_face(x)  
    av.x[] -= 1./pow(ND_Froude,2.0);
  foreach_face(y)  
    av.y[] += 0.0;
}

// INITIAL CONDITION
event init (t = 0.0) {

  // If statement to enable restarting of the simulation using the dump file "restart"
  if (!restore(file = "restart")) {

    // Determine filmHeight, such that y=filmHeight is the flat free surface of the
    // pool in the centered coordinate system (see the main function)
    filmHeight = -domainSize/2. + ND_thickness;

    // Strong refinement around the interfacial regions
    refine (((sq(x - (filmHeight + 0.5 + 0.1)) + sq(y) < sq(0.5*1.005) && sq(x - (filmHeight + 0.5 + 0.1)) + sq(y) > sq(0.5*0.995)) || fabs(x - filmHeight) <= 0.005) && level < maxLevel);
  
    // Create active liquid phase (f=1) as a union between the droplet and film
    fraction (f, union(- x + filmHeight, sq(0.5) - sq(x - (filmHeight + 0.5 + 0.1)) - sq(y)));
  
    // Add passive scalar inside droplet (for visualisation)
    fraction (dTracer, sq(0.5) - sq(x - (filmHeight + 0.5 + 0.1)) - sq(y));
  
    /* Initialise uniform velocity inside droplet (cell centres). dTracer=1
    inside the droplet and zero elsewhere. We do not use the volume fraction,
    because it is also non zero in the pool. Remember: the left boundary of the
    domain is the base of the pool, so the droplet should "fall" in the negative
    x direction. We have non-dimensionalised wrt the initial droplet velocity
    above */
    foreach()
    {
      u.x[]   = -1.0*dTracer[];
      u.y[]   = 0.0;
      p[]     = 0.0;
      omega[] = 0.0;
    }
  }
  
}

// DYNAMIC MESH
event adapt (i++) {

  // Refine only with respect to interfacial shape location and velocity component
  // magnitude
  adapt_wavelet ((scalar *){f, dTracer, u}, (double[]){1e-4, 1e-2, 5e-3, 5e-3}, maxLevel, minLevel);


  // remove droplets in boundary region for both active VOF and tracer
  foreach() {
    if (f[] > 0 && ((y > domainSize - 0.1) && ((x > filmHeight + 0.1)))) {
      f[] = 0;
    }

    if (dTracer[] > 0 && ((y > domainSize - 0.1) && ((x > filmHeight + 0.1)))) {
      dTracer[] = 0;
    }
  }

}

/* REMOVE SMALL DROPLETS AND BUBBLES
Remove_droplets is defined in tag.h (included above). Droplets/bubbles
constituting less than 4 cells are removed.
Fourth argument: true = bubble, false = droplet. */
event filterDroplets (i++)
{
    remove_droplets(f, 4, 1e-3);              // VOF droplet
    remove_droplets(f, 4, 1e-3, true);        // VOF bubble
    remove_droplets(dTracer, 4, 1e-3);        // Passive scalar droplet
    remove_droplets(dTracer, 4, 1e-3, true);  // Passive scalar bubble
}

// SAVE SLICES with which we can image process
event gfsview (t = 0.0; t += 0.01; t <= tEnd) {
    char name_gfs[200];
    sprintf(name_gfs,"Slices/DropImpact-%0.2f.gfs",t);

    FILE* fp_gfs = fopen (name_gfs, "w");
    output_gfs(fp_gfs);
    fclose(fp_gfs);
}

// MAKE VIDEOS
event movies (t += 0.005){

  char timestring[100];
 
  foreach(){
      omega[] = (u.y[1,0] - u.y[-1,0])/(2.*Delta) - (u.x[0,1] - u.x[0,-1])/(2.*Delta);
      viewingfield[] = 1.0 - f[] - dTracer[];
      mylevel[] = level;
      myrho[]   = rho(f[]);       // Add mixture density to the slice information
      mymu[]    = mu(f[]);        // Add mixture viscosity to the slice information
  }
  
  view(width=1900, height=1050, fov=10.0, ty = 0.35, quat = {0,0,-0.707,0.707});
	
  clear();
  draw_vof("f", lw=2);
  squares("viewingfield", map = cool_warm, min = -0.5, max = 2.5);
  mirror({0,1}) {
  	draw_vof("f", lw=2);	
	cells(lw=0.5);
	squares("mylevel", map = cool_warm, min = minLevel, max = maxLevel);
  } 

  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  
  save ("Animations/ImpactSummaryClose.mp4");

  view(width=1900, height=1050, fov=10.0, ty = 0.35, quat = {0,0,-0.707,0.707});
  clear();
  
  draw_vof("f", lw=2);
  squares("u.x", map = cool_warm, min = -1., max = 1,);
  mirror({0,1}) {
  	draw_vof("f", lw=2);	
	squares("u.y", map = cool_warm, min = -0.5, max = 2.);
  } 

  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  
  save ("Animations/ImpactVelocitiesClose.mp4");

  view(width=1900, height=1050, fov=10.0, ty = 0.35, quat = {0,0,-0.707,0.707});
  clear();
  
  draw_vof("f", lw=2);
  squares("omega", map = cool_warm, min = -100., max = 100.);
  mirror({0,1}) {
  	draw_vof("f", lw=2);	
	squares("p", map = cool_warm, min = -0.2, max = 0.5);
  } 

  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  
  save ("Animations/ImpactPVortClose.mp4");

  view(width=1900, height=1050, fov=22.5, ty = 0.01, quat = {0,0,-0.707,0.707});
  clear();
  
  draw_vof("f", lw=2);
  squares("viewingfield", map = cool_warm, min = -0.5, max = 2.5);
  mirror({0,1}) {
  	draw_vof("f", lw=2);	
	cells(lw=0.5);
	squares("mylevel", map = cool_warm, min = minLevel, max = maxLevel);
  } 

  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  
  save ("Animations/ImpactSummaryFar.mp4");

  view(width=1900, height=1050, fov=22.5, ty = 0.01, quat = {0,0,-0.707,0.707});
  clear();
  
  draw_vof("f", lw=2);
  squares("u.x", map = cool_warm, min = -1., max = 1,);
  mirror({0,1}) {
  	draw_vof("f", lw=2);	
	squares("u.y", map = cool_warm, min = -0.5, max = 2.);
  } 

  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  
  save ("Animations/ImpactVelocitiesFar.mp4");

  view(width=1900, height=1050, fov=22.5, ty = 0.01, quat = {0,0,-0.707,0.707});
  clear();
  
  draw_vof("f", lw=2);
  squares("omega", map = cool_warm, min = -100., max = 100.);
  mirror({0,1}) {
  	draw_vof("f", lw=2);	
	squares("p", map = cool_warm, min = -0.2, max = 0.5);
  } 

  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  
  save ("Animations/ImpactPVortFar.mp4");

}

// SAVE DROPLET AND POOL INTERFACE POSITIONS
event saveInterfaces (t += 0.001) {

    char nameInterfaces1[200];

    sprintf(nameInterfaces1,"Interfaces/interfaces-%0.3f.dat",t);

    FILE * fp1 = fopen(nameInterfaces1, "w");
    output_facets (f, fp1);	
    fclose(fp1);

    char nameInterfaces2[200];

    sprintf(nameInterfaces2,"Interfaces/interfacesDrop-%0.3f.dat",t);

    FILE * fp2 = fopen(nameInterfaces2, "w");
    output_facets (dTracer, fp2);	
    fclose(fp2);
}

// OBTAIN DATA ON DROPLETS (SEPARATE FLUID REGIONS)
event droplets (t += 0.01)
{
  scalar m[];
  foreach()
    m[] = f[] > 1e-2;
  int n = tag (m);

  double v[n];
  coord b[n];
  for (int j = 0; j < n; j++)
    v[j] = b[j].x = b[j].y = b[j].z = 0.;
  foreach_leaf()
    if (m[] > 0) {
      int j = m[] - 1;
      v[j] += dv()*f[];
      coord p = {x,y,z};
      foreach_dimension()
	b[j].x += dv()*f[]*p.x;
    }

  #if _MPI
    MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #endif
  for (int j = 0; j < n; j++)
    fprintf (fp_droplets, "%d %g %d %g %g %g\n", i, t,
	     j, v[j], b[j].x/v[j], b[j].y/v[j]);
  fflush (fp_droplets);
}

// WRITE DATA TO LOG FILE
event logstats (t += 0.001) {

    timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
 
    // i, timestep, no of cells, real time elapsed, cpu time
    fprintf(fp_stats, "i: %i t: %g dt: %g #Cells: %ld Wall clock time (s): %g CPU time (s): %g \n", i, t, dt, grid->n, perf.t, s.cpu);
    fflush(fp_stats);
}

