#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <hdf5.h>
#include <mpi.h>
#include <stdbool.h>
#include <unistd.h>
#include <math.h>

double[3] find_disk_orientation(double *hydrogen_densities, double **gas_positions, double **gas_masses, double **gas_velocities, double *gas_temperatures, int Ngas, double[3] pos_center)
{

    int rank,group_number,group_rank,numtasks;
    int i;
    int group_ranks[numFilesPerSnap];
    int maxrank,maxindex;



    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

    printf("Rank %d:Starting Disk Orientation Finder...\n",rank);

    MPI_Group orig_group, new_group;
    MPI_Comm new_comm;

    MPI_Comm_group(MPI_COMM_WORLD, &orig_group);
    group_number = rank / numFilesPerSnap;
    group_rank = rank % numFilesPerSnap;

    //Generate the array of ranks assigned to a single snapshot
    for (i=0;i < numFilesPerSnap;i++){
        group_ranks[i] = group_number*numFilesPerSnap + i;
    }

    printf("Rank %d: Defining new subgroup\n",rank);
    fflush(stdout);
    //Define subgroup of only tasks working on a single snapshot
    MPI_Group_incl(orig_group, numFilesPerSnap, group_ranks, &new_group);
    printf("Rank %d: Defining new comm\n",rank);
    fflush(stdout);
    MPI_Comm_create(MPI_COMM_WORLD, new_group, &new_comm);
    

    double Lsum[3],Lranksum[3],Lhat[3];
    double r2,dx,dy,dz;

    for (i=0;i<Ngas;i++) {
        dx = gas_positions[i][0]-pos_center[0];
        dy = gas_positions[i][1]-pos_center[1];
        dz = gas_positions[i][2]-pos_center[2];
        r2 = dx*dx+dy*dy+dz*dz;
      //  h_massfrac = 1 - (gas_metallicities[i][0]+gas_metallicities[i][1]);
      //  nH = (gas_densities[0][i]*h_massfrac) / proton_mass;
        if (r2 < (rCore*rCore) & hydrogen_densities[i] > 1 & gas_temperatures[i] < 8000) {
            Lsum[0] += gas_masses[i]*(
                       gas_positions[i][1]*gas_velocities[i][2]
                       - gas_positions[i][2]*gas_velocities[i][1]);
            Lsum[1] += gas_masses[i]*(
                       gas_positions[i][2]*gas_velocities[i][0]
                       - gas_positions[i][0]*gas_velocities[i][2]);
            Lsum[2] += gas_masses[i]*(
                       gas_positions[i][0]*gas_velocities[i][1]
                       - gas_positions[i][1]*gas_velocities[i][0]);
    }
    MPI_Reduce( &Lsum, &Lranksum, 3, MPI_DOUBLE, MPI_SUM, 0, new_comm);


    if (rank==0) {
        Lmag = sqrt( Lranksum[0]*Lranksum[0]+Lranksum[1]*Lranksum[1]+Lranksum[2]*Lranksum[2] );
        Lhat[0] = Lranksum[0]/Lmag;
        Lhat[1] = Lranksum[1]/Lmag;
        Lhat[2] = Lranksum[2]/Lmag;
        printf("Lhat is (%f,%f,%f)\n",Lhat[0],Lhat[1],Lhat[2]);
    }
    MPI_Bcast(&Lhat,3,MPI_DOUBLE,0,new_comm);
    MPI_Barrier(new_comm);
    return Lhat;
}



