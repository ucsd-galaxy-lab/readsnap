#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <hdf5.h>
#include <mpi.h>
#include <stdbool.h>
#include <unistd.h>

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


int get_values_per_blockelement(char *name, int flag_metals)
{
    int values = 0;
    
    if (strcmp(name, "Coordinates") == 0 || strcmp(name, "Velocities") == 0)
      values = 3;
    else if (strcmp(name, "Masses") == 0 || strcmp(name, "ParticleIDs") == 0 || strcmp(name, "Density") == 0 || strcmp(name, "ElectronAbundance") == 0 || strcmp(name, "NeutralHydrogenAbundance") == 0)
      values = 1;
    else if (strcmp(name, "Metallicity") == 0){
      values = flag_metals; //Only take total metallicity and helium for now
    }
    else
      values = 0;

    return values;
}

/* Read header attributes for file */
void read_header_attributes_in_hdf5(char *fname, struct io_header *header)
{
    hid_t hdf5_file, hdf5_headergrp, hdf5_attribute;
    
    hdf5_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    hdf5_headergrp = H5Gopen1(hdf5_file, "/Header");

    int temp;
    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_ThisFile");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, header->npart);
    H5Aclose(hdf5_attribute);
    
    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_Total");
    H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header->npartTotal);
    H5Aclose(hdf5_attribute);
    
    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_Total_HighWord");
    H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header->npartTotalHighWord);
    H5Aclose(hdf5_attribute);
    
    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "MassTable");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, header->mass);
    H5Aclose(hdf5_attribute);
    
    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Time");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header->time);
    H5Aclose(hdf5_attribute);
    
    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumFilesPerSnapshot");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header->num_files);
    H5Aclose(hdf5_attribute);
    
    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_DoublePrecision");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header->flag_doubleprecision);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_Metals");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header->flag_metals);
    H5Aclose(hdf5_attribute);
    
    H5Gclose(hdf5_headergrp);
    H5Fclose(hdf5_file);
}

struct fileArray {
  char **fileArray;  // Array of file names
  int len; // number of files
  int filesPerSnap; // files per snapshot
}files;

struct dataArray {
  double ***data;
  int len[2];
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
  double *SmoothingLength;
  double *InternalEnergy;
  int len[2];
};

/* Gets list of snapshot files given the base name of the snapshots, the min and max snapshot number, and the set size 
 * between each snapshot.
 */
void getFileNames(char *dirc, char *fname_base, int minSnapNum, int maxSnapNum, int snapStep) {
  int i;
  bool isMultiPartFile,needToLoadMore,startingNewFile;
  char fnamePart[200],toAppend[200];
  int filePartIdx=0;
  int numFiles=0;
  int fileIdx = maxSnapNum;
  int fileNum = (maxSnapNum-minSnapNum)/snapStep+1;
  char **fileArray; //initialize scope



  printf("Getting array of files to open...\n");
  fflush(stdout);

  /*Generate Array of Files that need to be loaded*/

  strcpy(fnamePart,dirc);
  strcat(fnamePart,fname_base);
  sprintf(toAppend,"_%d.hdf5",maxSnapNum);
  strcat(fnamePart,toAppend);

  if(  access(  fnamePart, F_OK) != -1 ) { //Check if the given file exists
    isMultiPartFile=false;
    /* Allocate array of pointers to rows.*/
    fileArray = (char **) malloc (fileNum * sizeof (char *));
    /*Allocate space for floating point data.*/
    fileArray[0] = (char *) malloc (fileNum * 200 * sizeof (char));
    filePartIdx = 1;
    for (i=0;i<fileNum;i++) {
      strcpy(fnamePart,dirc);
      strcat(fnamePart,fname_base);
      sprintf(toAppend,"_%d.hdf5",maxSnapNum-i);
      strcat(fnamePart,toAppend);
      strcpy(fileArray[i],fnamePart);
    }
    files.filesPerSnap = 1;
  } else {
    isMultiPartFile=true;
    needToLoadMore=true;
    startingNewFile=false;

    while (needToLoadMore) {
      sprintf(toAppend,"_%d/",fileIdx);
      strcpy(fnamePart,dirc);
      strcat(fnamePart,"snapdir");
      strcat(fnamePart,toAppend);
      strcat(fnamePart,fname_base);
      sprintf(toAppend,"_%d.%d.hdf5",fileIdx,filePartIdx);
      strcat(fnamePart,toAppend);
      if( access( fnamePart, F_OK) != -1) {
        filePartIdx = filePartIdx+1;
        numFiles = numFiles+1;
      } else {
        if ( startingNewFile )
        {
          fileIdx = fileIdx-snapStep; //Go to next file
          filePartIdx=0;
        }
        if (fileIdx<minSnapNum) needToLoadMore=false;
        startingNewFile=true;
      }
    }

    files.filesPerSnap = numFiles/fileNum;

    /* Allocate array of pointers to rows.*/
    fileArray = (char **) malloc (numFiles * sizeof (char *));
    /*Allocate space for floating point data.*/
    fileArray[0] = (char *) malloc (numFiles * 200 * sizeof (char));
    /* Set the rest of the pointers to rows to the correct addresses.*/

    for (i=1; i<numFiles; i++) {
      fileArray[i] = fileArray[0] + i * 200;
    }       

    filePartIdx=0;
    fileIdx=maxSnapNum;
    for (i=0;i<numFiles;i++) {
      sprintf(toAppend,"_%d/",fileIdx);
      strcpy(fnamePart,dirc);
      strcat(fnamePart,"snapdir");
      strcat(fnamePart,toAppend);
      strcat(fnamePart,fname_base);
      sprintf(toAppend,"_%d.%d.hdf5",fileIdx,filePartIdx);
      strcat(fnamePart,toAppend);
      if( access( fnamePart, F_OK) != -1) {
        strcpy(fileArray[i],fnamePart);
        filePartIdx=filePartIdx+1;
      } else {
        fileIdx=fileIdx-1;
        filePartIdx=0;
        sprintf(toAppend,"_%d/",fileIdx);
        strcpy(fnamePart,dirc);
        strcat(fnamePart,"snapdir");
        strcat(fnamePart,toAppend);
        strcat(fnamePart,fname_base);
        sprintf(toAppend,"_%d.%d.hdf5",fileIdx,filePartIdx);
        strcat(fnamePart,toAppend);
        strcpy(fileArray[i],fnamePart);
        filePartIdx=filePartIdx+1;
      }
    }
  }

  // Assign the file array and size to the files struct
  files.len = numFiles;
  files.fileArray = fileArray;
  
  printf("Done getting files array.\n");
  fflush(stdout);

}


/* Read position for ptype particles from file named fname */
struct dataStruct readsnap(char *fileName, int ptype, char **params, int num_params) {

  hid_t          fid,group_id,dset_id,dtype;
  hid_t          dataspace,memspace;
  herr_t         status;

  hsize_t        memDim[2],dsetDim[2];
  hsize_t        memCount[2],dsetCount[2];
  hsize_t        memOffset[2],dsetOffset[2];

  struct dataStruct dataStruct;
  char *name;

  int rank,numtasks;

  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  struct io_header header;

  int i;
  char ptype_str[200],dset_str[200];
   

  /*  Get info from header*/
  read_header_attributes_in_hdf5(fileName, &header);
  int N = header.npart[ptype];
  /*Open the HDF5 File*/
  fid = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);


  int p_index; // parameter index
  int num_elems; // number of elements for given parameter


  if(header.npart[ptype] > 0)
  {
    /*Define the group id string*/
    sprintf(ptype_str, "/PartType%d",ptype);
    /*Open the particle type group*/
    group_id=H5Gopen2(fid,ptype_str,H5P_DEFAULT);

    for (p_index=0;p_index<num_params;p_index++)
    {
      /*Set the dataset id*/
      sprintf(dset_str,params[p_index]);

      /*Set the dataset type*/
      dtype = H5Tcopy(H5T_NATIVE_DOUBLE); //Can just define a library for various keywords/steal from gizmo?
      /*Open the dataset*/
      dset_id=H5Dopen2(group_id,dset_str,H5P_DEFAULT);


      num_elems = get_values_per_blockelement(params[p_index],header.flag_metals);

      /*Define Hyperslab in dataset. Not necesarry for non parallel version*/
      dsetOffset[0]=0;
      dsetOffset[1]=0;
      dsetCount[0]= N; //Just take everything for non parallel version
      dsetCount[1]= num_elems;
      /* Create the memory allocation in the variable space*/
      dataspace = H5Screate_simple(2,dsetCount,NULL);
      status = H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,dsetOffset,NULL,dsetCount,NULL);

      /*Define memory dataspace*/
      memDim[0]=N;
      memDim[1]=num_elems;
      memspace = H5Screate_simple(2,memDim,NULL); //Rank 2 for vector

      /*Define memory hyperslab*/
      memOffset[0]=0;
      memOffset[1]=0;
      memCount[0]=N;
      memCount[1]=num_elems;
      status = H5Sselect_hyperslab(memspace,H5S_SELECT_SET,memOffset,NULL,memCount,NULL);

      /*Initialize variables to allocate memory for the data*/
      hsize_t dims[2] = {num_elems,N};
      hid_t space;
      int ndims;
      space = H5Dget_space(dset_id);
      ndims = H5Sget_simple_extent_dims(space,dims,NULL);
      /*****************************************************/

      name = params[p_index];

      if (strcmp(name, "Coordinates")==0){
        dataStruct.coordinates = (double **) malloc (dims[0] * sizeof (double*));
        dataStruct.coordinates[0] = (double *) malloc (dims[0] * dims[1] *sizeof (double));
        for (i=1;i<dims[0];i++){dataStruct.coordinates[i] = dataStruct.coordinates[0]+i*dims[1];}
        status = H5Dread(dset_id,dtype,memspace,dataspace,H5P_DEFAULT,dataStruct.coordinates[0]);
      }

      if (strcmp(name, "Velocities")==0){
        dataStruct.velocities = (double **) malloc (dims[0] * sizeof (double*));
        dataStruct.velocities[0] = (double *) malloc (dims[0] * dims[1] *sizeof (double));
        for (i=1;i<dims[0];i++){dataStruct.velocities[i] = dataStruct.velocities[0]+i*dims[1];}
        status = H5Dread(dset_id,dtype,memspace,dataspace,H5P_DEFAULT,dataStruct.velocities[0]);
      }

      if (strcmp(name, "Density")==0){
        dataStruct.density = (double *) malloc (dims[1] * sizeof (double));
        status = H5Dread(dset_id,dtype,H5S_ALL,H5S_ALL,H5P_DEFAULT,dataStruct.density);
        
      }

      if (strcmp(name, "Masses")==0){
        dataStruct.masses = (double *) malloc (dims[1] * sizeof (double));
        status = H5Dread(dset_id,dtype,H5S_ALL,H5S_ALL,H5P_DEFAULT,dataStruct.masses);
      }

      if (strcmp(name, "Metallicity")==0){
        dataStruct.metallicity = (double **) malloc (dims[0] * sizeof (double*));
        dataStruct.metallicity[0] = (double *) malloc (dims[0] * dims[1] *sizeof (double));
        for (i=1;i<dims[0];i++){dataStruct.metallicity[i] = dataStruct.metallicity[0]+i*dims[1];}
        status = H5Dread(dset_id,dtype,memspace,dataspace,H5P_DEFAULT,dataStruct.metallicity[0]);
      }

      if (strcmp(name, "NeutralHydrogenAbundance")==0){
        dataStruct.NeutralHydrogenAbundance = (double *) malloc (dims[1] * sizeof (double));
        status = H5Dread(dset_id,dtype,H5S_ALL,H5S_ALL,H5P_DEFAULT,dataStruct.NeutralHydrogenAbundance);
      }

      if (strcmp(name, "ElectronAbundance")==0){
        dataStruct.ElectronAbundance = (double *) malloc (dims[1] * sizeof (double));
        status = H5Dread(dset_id,dtype,H5S_ALL,H5S_ALL,H5P_DEFAULT,dataStruct.ElectronAbundance);
      }

      if (strcmp(name, "SmoothingLength")==0){
        dataStruct.SmoothingLength = (double *) malloc (dims[1] * sizeof (double));
        status = H5Dread(dset_id,dtype,H5S_ALL,H5S_ALL,H5P_DEFAULT,dataStruct.SmoothingLength);
      }

      if (strcmp(name, "InternalEnergy")==0){
        dataStruct.InternalEnergy = (double *) malloc (dims[1] * sizeof (double));
        status = H5Dread(dset_id,dtype,H5S_ALL,H5S_ALL,H5P_DEFAULT,dataStruct.InternalEnergy);
      }

      H5Tclose(dtype);
      H5Sclose(memspace);
      H5Sclose(dataspace);
      H5Dclose(dset_id);

    }
  }
  dataStruct.len[0]=num_params; dataStruct.len[1]=N;

  return dataStruct;
}
