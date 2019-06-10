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
  char *dirc = "/oasis/scratch/comet/ctrapp/temp_project/m12m_res7100/";
  char *file_base = "snapshot";

  int ptype;
  struct dataStruct dataGas,dataStars;
  int Ngas,Nstars;

  char *gas_params[] = {
      "Coordinates",
      "Density",
      "Masses",
      "Velocities",
      "Metallicity",
      "NeutralHydrogenAbundance",
      "ElectronAbundance",
      "SmoothingLength",
      "InternalEnergy"
      };
  char *star_params[] = {
      "Coordinates",
      "Masses",
      "Velocities"
      };
  int num_gas_params = sizeof(gas_params)/sizeof(gas_params[0]);
  int num_star_params = sizeof(star_params)/sizeof(star_params[0]);

  int minSnapNum = 581;
  int maxSnapNum = 600;
  int snapStep = 1;


  double start, end;
  rc = MPI_Init(&argc, &argv);

  // Get the start time
  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();


  if (rc != MPI_SUCCESS) {
      printf("Error starting MPI program. Terminating.\n");
      MPI_Abort(MPI_COMM_WORLD, rc);
  }

  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  printf("Number of tasks= %d My rank %d\n",numtasks,rank);


  // Create the array of file names 
  getFileNames(dirc, file_base, minSnapNum, maxSnapNum, snapStep);
  printf("Number of files to read: %d \n", files.len);
  fflush(stdout);

  int fileCount = rank;
  int numFiles = files.len;
  int numFilesPerSnap = files.filesPerSnap;


  // Check that we have the right number of nodes for the given number of files per snapshot
  MPI_Barrier(MPI_COMM_WORLD);
  if (numtasks%numFilesPerSnap != 0){
    if (rank == 0) {
      printf("Need to make sure number of MPI tasks is divisible by number of files per snapshot!\n");
      printf("Number of tasks = %d Number of files per snapshot = %d \n", numtasks, numFilesPerSnap);
      fflush(stdout);
    }
    MPI_Finalize();
    return 0;
  }

  // Needed for shrinking sphere and disk orientation calculation
  double* nH;
  double* NH1;
  double* gasTemp;
  double rShrinkSphere = 5000.0;
  double shrinkFactor = 0.7;
  double rMinSphere = 10.0;
  double* pos_center;
  double* vel_center;
  double* Lhat;

  // Now go through files taking strides of numtasks until all files have been read and ata processed
  while (fileCount<numFiles)
  {

    //Load Gas
    ptype=0;
    printf("Rank %d is loading GAS: %s\n",rank,files.fileArray[fileCount]);
    dataGas = readsnap(files.fileArray[fileCount], ptype, gas_params, num_gas_params);
    Ngas = dataGas.len[1];


    //Load Stars
    ptype=1;
    printf("Rank %d is loading STARS: %s\n",rank,files.fileArray[fileCount]);
    dataStars = readsnap(files.fileArray[fileCount], ptype, star_params, num_star_params);
    Nstars = dataStars.len[1];

    fileCount+=numtasks;
        printf("Rank %d: calculating parameters\n",rank);
    fflush(stdout);
    

    nH = calcHydrogenNumberDensity(dataGas.metallicity,dataGas.density,Ngas);
    NH1 = calcH1Abundance(dataGas.masses,dataGas.NeutralHydrogenAbundance,dataGas.SmoothingLength,dataGas.density,dataGas.metallicity,Ngas);
    gasTemp = calcTemperatures(dataGas.InternalEnergy,dataGas.ElectronAbundance,dataGas.metallicity,Ngas);

    printf("Rank %d: finding galactic disk\n",rank);
    fflush(stdout);


    pos_center = shrinking_sphere_parallel(dataGas.density, dataGas.coordinates, Ngas, dataStars.masses, dataStars.coordinates , Nstars, numFilesPerSnap, rShrinkSphere, shrinkFactor, rMinSphere);

    vel_center = center_of_mass_velocity(dataStars.masses, dataStars.velocities, dataStars.coordinates, pos_center, Nstars, numFilesPerSnap);

    Lhat = find_disk_orientation(nH, dataGas.coordinates, dataGas.masses, dataGas.velocities, gasTemp, Ngas, vel_center, pos_center, numFilesPerSnap);
 


    // free up all the memory we had returned
    free(dataGas.coordinates[0]);
    free(dataGas.coordinates);
    free(dataGas.velocities[0]);
    free(dataGas.velocities);
    free(dataGas.metallicity[0]);
    free(dataGas.metallicity);
    free(dataGas.masses);
    free(dataGas.density);
    //free(dataGas.pids);
    //free(dataGas.cids);
    free(dataGas.ElectronAbundance);
    free(dataGas.NeutralHydrogenAbundance);
    free(dataGas.SmoothingLength);
    free(dataGas.InternalEnergy);
    //free(dataGas.len);
    //free(dataGas);
    free(dataStars.coordinates[0]);
    free(dataStars.coordinates);
    free(dataStars.velocities[0]);
    free(dataStars.velocities);
    //free(dataStars.metallicity);
    free(dataStars.masses);
    //free(dataStars.len);
    //free(dataStars);
    free(nH);
    free(NH1);
    free(gasTemp);
    //free(pos_center);
    //free(vel_center);
    //free(Lhat);
  }

  printf("Rank %d. Out of files to load\n",rank);
  printf("fileCount: %d  numFiles: %d\n",fileCount,numFiles);
  fflush(stdout);


  

  //free(dataArray_gas.data);
  //free(dataArray_stars.data);

  // get end time
  MPI_Barrier(MPI_COMM_WORLD); 
  end = MPI_Wtime();

  MPI_Finalize();

  if (rank == 0) { /* use time on master node */
    printf("Runtime = %f\n", end-start);
  }

  return 0;

}
