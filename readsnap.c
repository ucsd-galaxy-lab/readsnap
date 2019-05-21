#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <hdf5.h>
#include <mpi.h>

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


/* Read position for ptype particles from file named fname */
int readsnap(char *fname, int ptype) {
   hid_t          fid,group_id,dset_id,dtype;
   hid_t          dataspace,memspace;
   herr_t         status;

   hsize_t        memDim[2],dsetDim[2];
   hsize_t        memCount[2],dsetCount[2];
   hsize_t        memOffset[2],dsetOffset[2];
 
   double **pos;
   int i;

   char ptype_str[200],dset_str[200];

   struct io_header header;



   /*  Get info from header*/
   read_header_attributes_in_hdf5(fname, &header);
   int N = header.npart[ptype];
   printf("Loading %d particles.\n", N);
   fflush(stdout);

   /*Open the HDF5 File*/
   fid = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

   printf("Ptype: %d\n", ptype);
   fflush(stdout);

   if(header.npart[ptype] > 0)
   {
     /*Define the group id string*/
     sprintf(ptype_str, "/PartType%d",ptype);
     /*Open the particle type group*/
     group_id=H5Gopen2(fid,ptype_str,H5P_DEFAULT);

     /*Set the dataset id*/
     sprintf(dset_str,"Coordinates");
     /*Set the dataset type*/
     dtype = H5Tcopy(H5T_NATIVE_DOUBLE); //Can just define a library for various keywords/steal from gizmo?
     /*Open the dataset*/
     dset_id=H5Dopen2(group_id,dset_str,H5P_DEFAULT);




     /*Define Hyperslab in dataset. Not necesarry for non parallel version*/
     dsetOffset[0]=0;
     dsetOffset[1]=0;
     dsetCount[0]= N; //Just take everything for non parallel version
     dsetCount[1]= 3;
     /* Create the memory allocation in the variable space*/
     dataspace = H5Screate_simple(2,dsetCount,NULL);
     status = H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,dsetOffset,NULL,dsetCount,NULL);
   
     /*Define memory dataspace*/
     memDim[0]=N;
     memDim[1]=3;
     memspace = H5Screate_simple(2,memDim,NULL); //Rank 2 for vector
   
     /*Define memory hyperslab*/
     memOffset[0]=0;
     memOffset[1]=0;
     memCount[0]=N;
     memCount[1]=3;
     status = H5Sselect_hyperslab(memspace,H5S_SELECT_SET,memOffset,NULL,memCount,NULL);

  
     /*Initialize variables to allocate memory for the data*/
     hsize_t dims[2] = {3,N};
     hid_t space;
     int ndims;
     space = H5Dget_space(dset_id);
     ndims = H5Sget_simple_extent_dims(space,dims,NULL);
     /*****************************************************/


 
     /* Allocate array of pointers to rows.*/
    pos = (double **) malloc (dims[0] * sizeof (double *));

    /*Allocate space for floating point data.*/
    pos[0] = (double *) malloc (dims[0] * dims[1] * sizeof (double));

 
     /* Set the rest of the pointers to rows to the correct addresses.*/
    for (i=1; i<dims[0]; i++)
        pos[i] = pos[0] + i * dims[1];



     /*Push the dataset into the position vector*/
     status = H5Dread(dset_id,dtype,memspace,dataspace,H5P_DEFAULT,pos[0]); //memspace/dataspace not needed for nonparallel version.

     H5Tclose(dtype);
     H5Sclose(memspace);
     H5Sclose(dataspace);
     H5Dclose(dset_id);

    

   }
   
   /*Print some elements to check if reasonable*/
   printf("[0][0] element: %f\n", pos[0][0]);
   printf("[1][0] element: %f\n", pos[1][0]);
   printf("[2][0] element: %f\n\n", pos[2][0]);
   printf("[0][1] element: %f\n", pos[0][1]);
   printf("[1][1] element: %f\n", pos[1][1]);
   printf("[2][1] element: %f\n\n", pos[2][1]);

   fflush(stdout);
   

   return 0;
}


int main( int argc, char *argv[])
{
    int numtasks, rank, rc;
    char *file_name = "snapshot_600.hdf5";
    rc = MPI_Init(&argc, &argv);
    if (rc != MPI_SUCCESS) {
        printf("Error starting MPI program. Terminating.\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    printf("Number of tasks= %d My rank %d\n",numtasks,rank);
    
    
    int ptype = 0;
    printf("Ptype: %d\n", ptype);
    printf("File name: %s\n", file_name);
    fflush(stdout);
    readsnap(file_name, ptype);
    

    MPI_Finalize();

}
