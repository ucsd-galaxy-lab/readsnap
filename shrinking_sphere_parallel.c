#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <hdf5.h>
#include <mpi.h>
#include <stdbool.h>
#include <unistd.h>

int shrinking_sphere_parallel(double **gas_densities, double **gas_positions, int Ngas, double **masses, double **positions, int nStars, int numFilesPerSnap, float rShrinkSphere, float shrinkFactor, float rMinSphere)
{
//Find the maximum density as first guess for center
    printf("Starting Shrinking Sphere Function...\n");
    int rank,group_number,group_rank,numtasks;
    int i;
    int group_ranks[numFilesPerSnap];
    int maxrank,maxindex;
    int pos_center[3];

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

    MPI_Group orig_group, new_group;
    MPI_Comm new_comm;

    printf("ping\n");
    fflush(stdout);
    MPI_Comm_group(MPI_COMM_WORLD, &orig_group);
    group_number = rank / numFilesPerSnap;
    group_rank = rank % numFilesPerSnap;

    //Generate the array of ranks assigned to a single snapshot
    for (i=0;i < numFilesPerSnap;i++){
        group_ranks[i] = group_number*numFilesPerSnap + i;
    }

    printf("Defining new subgroup\n");
    fflush(stdout);
    //Define subgroup of only tasks working on a single snapshot
    MPI_Group_incl(orig_group, numFilesPerSnap, group_ranks, &new_group);
    printf("Defining new comm\n");
    fflush(stdout);
    MPI_Comm_create(MPI_COMM_WORLD, new_group, &new_comm);
    

    printf("Defining struct for MPI_MAXLOC\n");
    fflush(stdout);
    struct {
        double value;
        int idx_rank;
    } in, out;
    in.value = gas_densities[0][0];
    in.idx_rank = group_rank;
/*Find the maximum value on this rank*/
    for (i=1; i < Ngas; i++){
        if (in.value < gas_densities[i][0]) {
            in.value = gas_densities[i][0];
            maxindex = i;
        }
    }
    //in.index=rank*Ngas+in.index; //Define a global index
/*Pass Maximum Value to Reduction Algorithm*/
    printf("Running MPI_MAXLOC\n");
    fflush(stdout);
    MPI_Reduce( &in, &out , 1, MPI_DOUBLE, MPI_MAXLOC, 0, new_comm);

/*Maximum resides in root*/

    if (group_rank==0) {
        printf("Finding maximum density location\n");
        fflush(stdout);
        maxrank = out.idx_rank;
        //maxindex = out.index % Ngas;
        /*Sends the rank and index to all ranks in group*/
        printf("Broadcasting maximum density location\n");
        fflush(stdout);
        MPI_Bcast(&maxrank,1,MPI_INT,0,new_comm);
        //MPI_Bcast(&maxindex,1,MPI_INT,0,new_comm);
    }
    MPI_Barrier(new_comm); /*Make sure everything is synced up Needed here? */

    /*Find the initial guess for central position*/
    if (group_rank==maxrank) {
        printf("Identifying first guess pos_center: ");
        fflush(stdout);
        pos_center[0] = gas_positions[maxindex][0];
        pos_center[1] = gas_positions[maxindex][1];
        pos_center[2] = gas_positions[maxindex][2];
        MPI_Bcast(&pos_center,1,MPI_DOUBLE,maxrank,new_comm);
        printf("(%f, %f, %f)\n",pos_center[0],pos_center[1],pos_center[2]);
        fflush(stdout);
    }
    
    double deltaPos, mw_pos_sum[3], mass_sum;
    double mass_ranksum, mw_pos_ranksum;
    double old_pos[3];
    double r2;

    double deltaMin = 0.000001;
    int wiggles = 0;
    double dx,dy,dz;
    printf("Rank %d: Beginning Shrinking Sphere Algorithm...\n",rank);
    fflush(stdout);
    while (rShrinkSphere > rMinSphere) {
        printf("Rank %d: Looping over particles...\n",rank);
        fflush(stdout);
        for (i=0;i<Nstars;i++){
            dx = positions[i][0]-pos_center[0];
            dy = positions[i][1]-pos_center[1];
            dz = positions[i][2]-pos_center[2];
            r2 = dx*dx+dy*dy+dz*dz;

            if (r2 < (rShrinkSphere*rShrinkSphere)){
                mass_sum += masses[i][0];
                mw_pos_sum[0] += masses[i][0]*positions[i][0]; 
                mw_pos_sum[1] += masses[i][0]*positions[i][1]; 
                mw_pos_sum[2] += masses[i][0]*positions[i][2]; 
            }
        }
        printf("Rank %d: Summing Masses...\n",rank);
        fflush(stdout);
        MPI_Reduce( &mass_sum, &mass_ranksum, 1, MPI_DOUBLE, MPI_SUM, 0, new_comm);
        printf("Rank %d: Summing Positions...",rank);
        fflush(stdout);
        MPI_Reduce( &mw_pos_sum, &mw_pos_ranksum, 3, MPI_DOUBLE, MPI_SUM, 0, new_comm);
 
        if (group_rank==0){
            printf("Rank %d: Calculating Average\n",rank);
            fflush(stdout);
            pos_center[0] = mw_pos_sum[0] / mass_sum; //Calculate center of mass
            pos_center[1] = mw_pos_sum[1] / mass_sum; //Calculate center of mass
            pos_center[2] = mw_pos_sum[2] / mass_sum; //Calculate center of mass
            MPI_Bcast(&pos_center,1,MPI_DOUBLE,0,new_comm);
        }
        MPI_Barrier(new_comm); //Make sure all ranks have new position
        rShrinkSphere = shrinkFactor*rShrinkSphere;
        printf("Rank %d: Shrinking Sphere...\n",rank);
        fflush(stdout);  
    }

    printf("Rank %d: Beginning sphere wiggling.\n",rank);
    fflush(stdout);
    rShrinkSphere = rMinSphere;
    deltaPos = 1000.0;
    /* Wiggle the sphere without changing the radius to ensure convergence */
    while ((deltaPos > deltaMin) & (wiggles < 10)){
        printf("Rank %d: Wiggling Sphere...\n",rank);
        fflush(stdout);
        wiggles+=1;
        for (i=0;i<Nstars;i++){
            dx = positions[i][0]-pos_center[0];
            dy = positions[i][1]-pos_center[1];
            dz = positions[i][2]-pos_center[2];
            r2 = dx*dx+dy*dy+dz*dz;

            if (r2 < (rShrinkSphere*rShrinkSphere)){
                mass_sum += masses[i][0];
                mw_pos_sum[0] += masses[i][0]*positions[i][0]; 
                mw_pos_sum[1] += masses[i][0]*positions[i][1]; 
                mw_pos_sum[2] += masses[i][0]*positions[i][2]; 
            }
        }

        MPI_Reduce( &mass_sum, &mass_ranksum, 1, MPI_DOUBLE, MPI_SUM, 0, new_comm);
        MPI_Reduce( &mw_pos_sum, &mw_pos_ranksum, 3, MPI_DOUBLE, MPI_SUM, 0, new_comm);
 
        old_pos[0] = pos_center[0];
        old_pos[1] = pos_center[1];
        old_pos[2] = pos_center[2];
        if (group_rank==0){
            pos_center[0] = mw_pos_sum[0] / mass_sum; //Calculate center of mass
            pos_center[1] = mw_pos_sum[1] / mass_sum; //Calculate center of mass
            pos_center[2] = mw_pos_sum[2] / mass_sum; //Calculate center of mass
            MPI_Bcast(&pos_center,1,MPI_DOUBLE,0,new_comm);
        }
        MPI_Barrier(new_comm);        
        dx = pos_center[0]-old_pos[0];
        dy = pos_center[1]-old_pos[1];
        dz = pos_center[2]-old_pos[2];
        deltaPos = dx*dx+dy*dy+dz*dz;
    }


    printf("Final Central Position is: (%f, %f, %f)\n",pos_center[0],pos_center[1],pos_center[2]);
    fflush(stdout);
    return 0;
}


