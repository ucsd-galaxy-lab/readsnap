#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <hdf5.h>
#include <mpi.h>
#include <stdbool.h>
#include <unistd.h>
#include "allvars.h"


int main( int argc, char *argv[])
{
  int numtasks, rank, rc;
  //char *file_base = "/oasis/tscc/scratch/ctrapp/m12m_res7100/output/snapdir_600/snapshot";
  char *file_base = "/oasis/scratch/comet/ctrapp/temp_project/m12m_res7100/snapdir_600/snapshot";


  int ptype;
  //struct dataArray dataArray_stars,dataArray_gas;
  struct dataStruct dataGas,dataStars;
  //double ***data_gas;
  //double ***data_stars; // data to be loaded from snapshots
  int Ngas,Nstars;
  char *gas_params[] = {
      "Coordinates",
      "Density",
      "Velocities",
      "Metallicity"
      };
  char *star_params[] = {
      "Coordinates",
      "Masses"
      };
  int num_gas_params = sizeof(gas_params)/sizeof(gas_params[0]);
  int num_star_params = sizeof(star_params)/sizeof(star_params[0]);

  int minSnapNum = 600;
  int maxSnapNum = 600;
  int snapStep = 1;

  rc = MPI_Init(&argc, &argv);
  if (rc != MPI_SUCCESS) {
      printf("Error starting MPI program. Terminating.\n");
      MPI_Abort(MPI_COMM_WORLD, rc);
  }

  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  printf("Number of tasks= %d My rank %d\n",numtasks,rank);


  // Create the array of file names 
  getFileNames(file_base, minSnapNum, maxSnapNum, snapStep);
  printf("Number of files: %d \n", files.len);
  fflush(stdout);

  int fileCount = rank;
  int numFiles = files.len;

  // Now go through files taking strides of numtasks until all files have been read
  while (fileCount<numFiles)
  {

    //Load Gas
      ptype=0;
      printf("Rank %d is loading GAS: %s\n",rank,files.fileArray[fileCount]);
      dataGas = readsnap(files.fileArray[fileCount], ptype, gas_params, num_gas_params);
      //data_gas = dataArray_gas.data;
      Ngas = dataGas.len[1];
      //printf("Rank %d: Metallicity entries are: %f, %f\n",rank,data_gas[2][0][0],data_gas[2][0][1]);


    //Load Stars
      ptype=1;
      printf("Rank %d is loading STARS: %s\n",rank,files.fileArray[fileCount]);
      dataStars = readsnap(files.fileArray[fileCount], ptype, star_params, num_star_params);
      //data_stars = dataArray_stars.data;
      Nstars = dataStars.len[1];

      fileCount+=numtasks;

      //printf("Dimensions of data: params x num_particles x data_elems\n");
      //printf("%d x %d x (depends on data type)\n", dataGas.len[0],dataGas.len[1]);
  }
  printf("Rank %d: Trying to print a thing!!\n",rank);
  printf("Rank %d: gasMetal[0][0] is %f\n",rank,dataGas.metallicity[0][0]);
  printf("Rank %d: gasMetal[0][1] is %f\n",rank,dataGas.metallicity[0][1]);

  printf("Rank %d. Out of files to load\n",rank);
  printf("fileCount: %d  numFiles: %d\n",fileCount,numFiles);
  fflush(stdout);

  nH = calcHydrogenNumberDensity(dataGas.metallicity,dataGas.density,Ngas);
  NH1 = calcH1Abundance(dataGas.masses,dataGas.NeutralHydrogenAbundance,dataGas.smoothingLengths,dataGas.density,dataGas.metallicity,Ngas);
  gasTemp = calcTemperatures(dataGas.internalEnergy);
  double rShrinkSphere = 5000.0;
  double shrinkFactor = 0.7;
  double rMinSphere = 10.0;
  int numFilesPerSnap = 4;
  shrinking_sphere_parallel(dataGas.density, dataGas.coordinates, Ngas, dataStars.masses, dataStars.coordinates , Nstars, numFilesPerSnap, rShrinkSphere, shrinkFactor, rMinSphere);

  //free(dataArray_gas.data);
  //free(dataArray_stars.data);
  MPI_Finalize();

  return 0;

}
