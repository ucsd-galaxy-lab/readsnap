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


int get_values_per_blockelement(char *name)
{
    int values = 0;
    
    if (strcmp(name, "Coordinates") == 0 || strcmp(name, "Velocities") == 0)
      values = 3;
    else if (strcmp(name, "Masses") == 0 || strcmp(name, "ParticleIDs") == 0 || strcmp(name, "Density") == 0 || strcmp(name, "ElectronAbundance") == 0 || strcmp(name, "NeutralHydrogenAbundance") == 0)
      values = 1;
    else if (strcmp(name, "metallicity") == 0)
      values = 15;
    else
      values = 0;

    return values;
}

/* Read header attributes for file */
void read_header_attributes_in_hdf5(char *fname, struct io_header *header)
{

    printf("File name: %s\n", fname);
    fflush(stdout);
    hid_t hdf5_file, hdf5_headergrp, hdf5_attribute;
    
    hdf5_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    hdf5_headergrp = H5Gopen1(hdf5_file, "/Header");

    //struct io_header *retVal = malloc(sizeof(struct io_header))

    int temp;
    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_ThisFile");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, header->npart);
    printf("npart[0]: %d\n", (*header).npart[0]);
    fflush(stdout);
    H5Aclose(hdf5_attribute);
    //retVal->npart = temp;
    //*header = retVal;
    
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
    
    H5Gclose(hdf5_headergrp);
    H5Fclose(hdf5_file);
}

struct fileArray {
  char **fileArray;  // Array of file names
  int len; // number of files
}files;

struct dataArray {
  double ***data;
  int len[2];
};

/* Gets list of snapshot files given the base name of the snapshots, the min and max snapshot number, and the set size 
 * between each snapshot.
 */
void getFileNames(char *fname_base, int minSnapNum, int maxSnapNum, int snapStep) {
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

  strcpy(fnamePart,fname_base);
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
       strcpy(fnamePart,fname_base);
       sprintf(toAppend,"_%d.hdf5",maxSnapNum-i);
       strcat(fnamePart,toAppend);
       strcpy(fileArray[i],fnamePart);
    }      
  } else {
    isMultiPartFile=true;
    needToLoadMore=true;
    startingNewFile=false;

    while (needToLoadMore) {
      sprintf(toAppend,"_%d.%d.hdf5",fileIdx,filePartIdx);
      strcpy(fnamePart,fname_base);
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
    printf("maxSnapNum is now: %d\n",maxSnapNum);
    for (i=0;i<numFiles;i++) {
      sprintf(toAppend,"_%d.%d.hdf5",fileIdx,filePartIdx);
      strcpy(fnamePart,fname_base);
      strcat(fnamePart,toAppend);
      if( access( fnamePart, F_OK) != -1) {
        strcpy(fileArray[i],fnamePart);
        filePartIdx=filePartIdx+1;
      } else {
        fileIdx=fileIdx-1;
        filePartIdx=0;
        sprintf(toAppend,"_%d.%d.hdf5",fileIdx,filePartIdx);
        strcpy(fnamePart,fname_base);
        strcat(fnamePart,toAppend);
        strcpy(fileArray[i],fnamePart);
        filePartIdx=filePartIdx+1;
      }
    }
  }

  // Assign the file array and size to the files struct
  files.len = numFiles;
  files.fileArray = fileArray;
}


/* Read position for ptype particles from file named fname */
struct dataArray readsnap(char *fileName, int ptype, char **params, int num_params) {
  printf("Running Readsnap...\n");

  hid_t          fid,group_id,dset_id,dtype;
  hid_t          dataspace,memspace;
  herr_t         status;

  hsize_t        memDim[2],dsetDim[2];
  hsize_t        memCount[2],dsetCount[2];
  hsize_t        memOffset[2],dsetOffset[2];

  struct dataArray dataArray;
  double ***data;
  data = (double ***) malloc (num_params * sizeof (double **));

  int rank,numtasks;

  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  printf("Rank %d: Running Readsnap...\n",rank);

  struct io_header header;

  int i;
  char ptype_str[200],dset_str[200];
   

  /*  Get info from header*/
  read_header_attributes_in_hdf5(fileName, &header);
  int N = header.npart[ptype];
  printf("Rank %d: Loading %d particles.\n",rank, N);
  fflush(stdout);

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
      printf("P_Index: %d\n",p_index);
      fflush(stdout);
      /*Set the dataset id*/
      sprintf(dset_str,params[p_index]);

      /*Set the dataset type*/
      dtype = H5Tcopy(H5T_NATIVE_DOUBLE); //Can just define a library for various keywords/steal from gizmo?
      /*Open the dataset*/
      dset_id=H5Dopen2(group_id,dset_str,H5P_DEFAULT);


      num_elems = get_values_per_blockelement(params[p_index]);

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

      /* Allocate array of pointers to rows.*/
      data[p_index] = (double **) malloc (dims[0] * sizeof (double *));

     /*Allocate space for floating point data.*/
      data[p_index][0] = (double *) malloc (dims[0] * dims[1] * sizeof (double));
      /* Set the rest of the pointers to rows to the correct addresses.*/
      for (i=1; i<dims[0]; i++)
         data[p_index][i] = data[p_index][0] + i * dims[1];



      /*Push the dataset into the position vector*/
      if (num_elems==1){ // In case the data is a 1D array
        double *buffer = (double *) malloc (N * sizeof (double));
        fflush(stdout);
        status = H5Dread(dset_id,dtype,H5S_ALL,H5S_ALL,H5P_DEFAULT,buffer); //memspace/dataspace not needed for nonparallel version.
        printf("Rank %d: buffer[0] element: %f\n",rank, buffer[0]);
        fflush(stdout);
        data[p_index][0] = buffer;
        printf("Rank %d: data[0][0] element: %f\n",rank, data[p_index][0][0]);
        fflush(stdout);


      }
      else status = H5Dread(dset_id,dtype,H5S_ALL,H5S_ALL,H5P_DEFAULT,data[p_index][0]); //memspace/dataspace not needed for nonparallel version.

      H5Tclose(dtype);
      H5Sclose(memspace);
      H5Sclose(dataspace);
      H5Dclose(dset_id);

    }
  }

  dataArray.data = data;
  dataArray.len[0]=num_params; dataArray.len[1]=N;

  return dataArray;
}


int main( int argc, char *argv[])
{
  int numtasks, rank, rc;
  char *file_base = "/home/cchoban/test/test";
  struct dataArray dataArray;
  double ***data; // data to be loaded from snapshots

  // List of parameters you want to be loaded from each snapshot
  char *params[] = {
  "Coordinates",
  "Masses"
  };
  int num_params = sizeof(params)/sizeof(params[0]);

  int minSnapNum = 598;
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
    printf("Rank %d is loading: %s\n",rank,files.fileArray[fileCount]);
    int ptype = 0;
    dataArray = readsnap(files.fileArray[fileCount], ptype, params, num_params);
    data = dataArray.data;
    fileCount+=numtasks;

    printf("Dimensions of data: params x num_particles x data_elems\n");
    printf("%d x %d x (depends on data type)\n", dataArray.len[0],dataArray.len[1]);

    int p_index;
    for (p_index=0;p_index<dataArray.len[0];p_index++)
    {
      printf("Rank %d: [0][1] element: %f\n",rank, data[p_index][0][1]);
      fflush(stdout);
    } 

  }

  printf("Rank %d. Out of files to load\n",rank);
  printf("fileCount: %d  numFiles: %d\n",fileCount,numFiles);
  fflush(stdout);

  MPI_Finalize();

}
