/*! Header for the standard file format.
 */
extern struct io_header
{
  int npart[6];     /*!< number of particles of each type in this file */
  double mass[6];   /*!< mass of particles of each type. If 0, then the masses are explicitly
                                stored in the mass-block of the snapshot file, otherwise they are omitted */
  double time;      /*!< time of snapshot file */
  double redshift;    /*!< redshift of snapshot file */
  int flag_sfr;     /*!< flags whether the simulation was including star formation */
  int flag_feedback;    /*!< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[6]; /*!< total number of particles of each type in this snapshot. This can be
           different from npart if one is dealing with a multi-file snapshot. */
  int flag_cooling;   /*!< flags whether cooling was included  */
  int num_files;    /*!< number of files in multi-file snapshot */
  double BoxSize;   /*!< box-size of simulation in case periodic boundaries were used */
  double Omega0;    /*!< matter density in units of critical density */
  double OmegaLambda;   /*!< cosmological constant parameter */
  double HubbleParam;   /*!< Hubble parameter in units of 100 km/sec/Mpc */
  int flag_stellarage;    /*!< flags whether the file contains formation times of star particles */
  int flag_metals;    /*!< flags whether the file contains metallicity values for gas and star particles */
  unsigned int npartTotalHighWord[6]; /*!< High word of the total number of particles of each type (needed to combine with npartTotal to allow >2^31 particles of a given type) */
  int flag_doubleprecision; /*!< flags that snapshot contains double-precision instead of single precision */

  int flag_ic_info;             /*!< flag to inform whether IC files are generated with ordinary Zeldovich approximation,
                                     or whether they ocontains 2nd order lagrangian perturbation theory initial conditions.
                                     For snapshots files, the value informs whether the simulation was evolved from
                                     Zeldoch or 2lpt ICs. Encoding is as follows:
                                        FLAG_ZELDOVICH_ICS     (1)   - IC file based on Zeldovich
                                        FLAG_SECOND_ORDER_ICS  (2)   - Special IC-file containing 2lpt masses
                                        FLAG_EVOLVED_ZELDOVICH (3)   - snapshot evolved from Zeldovich ICs
                                        FLAG_EVOLVED_2LPT      (4)   - snapshot evolved from 2lpt ICs
                                        FLAG_NORMALICS_2LPT    (5)   - standard gadget file format with 2lpt ICs
                                     All other values, including 0 are interpreted as "don't know" for backwards compatability.
                                 */
  float lpt_scalingfactor;      /*!< scaling factor for 2lpt initial conditions */

  char fill[18];    /*!< fills to 256 Bytes */
  char names[15][2];
}
header;       /*!< holds header for snapshot files */

// Struct for files to be read
struct fileArray {
  char **fileArray;  // Array of file names
  int len; // number of files
}files;

// Struct for data to be returned by readsnap
struct dataArray {
  double ***data; 
  int len[2]; // dimensions of data array
};

struct dataStruct {
  double **coordinates;
  double **velocities;
  double **metallicity;
  double *masses;
  double *density;
  int *pids;
  int *cids;
  double *ElectronAbundance;
  double *NeutralHydrogenAbundance;
  double *smoothingLengths;
  double *internalEnergy;
  int len[2];
};

int get_values_per_blockelement(char *name);

void read_header_attributes_in_hdf5(char *fname, struct io_header *header);

void getFileNames(char *fname_base, int minSnapNum, int maxSnapNum, int snapStep);

//struct dataArray readsnap(char *fileName, int ptype, char **params, int num_params) ;
struct dataStruct readsnap(char *fileName, int ptype, char **params, int num_params) ;


// These functions do postprocessing on snapshot data

int shrinking_sphere_parallel(double *gas_densities, double **gas_positions, int Ngas, double *masses, double **positions, int Nstars, int numFilesPerSnap, double rShrinkSphere, double shrinkFactor, double rMinSphere);

double* find_disk_orientation(double *hydrogen_densities, double **gas_positions, double **gas_masses, double **gas_velocities, double *gas_temperatures, int Ngas, double pos_center[]);

double* calcHydogenNumberDensity(double** gas_metallicities, double** gas_densities, int Ngas);

double* calcH1Abundance(double** gas_masses, double** neutral_hydrogen_densities, double** kernalLengths, double** gas_densities, double** gas_metallicities, int Ngas);

struct dataArray calcTemp();
