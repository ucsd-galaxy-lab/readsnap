#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <mpi.h>
#include <stdbool.h>
#include <unistd.h>
#include <math.h>

double* shrinking_sphere_parallel(double *gas_densities, double **gas_positions, int Ngas, double *masses, double **positions, int Nstars, int numFilesPerSnap, double rShrinkSphere, double shrinkFactor, double rMinSphere)
{
//Find the maximum density as first guess for center
    int rank,group_number,group_rank,numtasks;
    int i;
    int group_ranks[numFilesPerSnap];
    int maxrank,maxindex;
    double pos_center[3];
    double* to_return;
    to_return = (double *) malloc(3*sizeof(double));

    MPI_Status status;
   
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

    printf("Rank %d:Starting Shrinking Sphere Function...\n",rank);

    MPI_Group orig_group, new_group;
    MPI_Comm new_comm;

    MPI_Comm_group(MPI_COMM_WORLD, &orig_group);
    group_number = rank / numFilesPerSnap;
    group_rank = rank % numFilesPerSnap;

    //Generate the array of ranks assigned to a single snapshot
    for (i=0;i < numFilesPerSnap;i++){
        group_ranks[i] = group_number*numFilesPerSnap + i;
    }

    //Define subgroup of only tasks working on a single snapshot
    MPI_Group_incl(orig_group, numFilesPerSnap, group_ranks, &new_group);
    MPI_Comm_create(MPI_COMM_WORLD, new_group, &new_comm);
    
    struct {
        double value;
        int idx_rank;
    } in, out;

    in.value = gas_densities[0];
    in.idx_rank = group_rank;



    /*Find the maximum value on this rank*/
    for (i=1; i < Ngas; i++){
        if (in.value < gas_densities[i]) {
            in.value = gas_densities[i];
            maxindex = i;
        }
    }
    printf("Rank %d: looped over in/out\n",rank);
    fflush(stdout);

    /*Pass Maximum Value to Reduction Algorithm*/
    MPI_Reduce( &in, &out , 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, new_comm);
    MPI_Barrier(new_comm);
    /*Maximum resides in root*/

    if (group_rank==0) {
        printf("Finding maximum density location\n");
        fflush(stdout);
        maxrank = out.idx_rank;
        //maxindex = out.index % Ngas;
        /*Sends the rank and index to all ranks in group*/
        printf("Broadcasting maximum density location\n");
        fflush(stdout);
    }
    MPI_Bcast(&maxrank,1,MPI_INT,0,new_comm);
    MPI_Barrier(new_comm);

    MPI_Barrier(new_comm); /*Make sure everything is synced up Needed here? */
    /*Find the initial guess for central position*/
    if (group_rank==maxrank) {
        pos_center[0] = gas_positions[maxindex][0];
        pos_center[1] = gas_positions[maxindex][1];
        pos_center[2] = gas_positions[maxindex][2];
    }
    MPI_Bcast(&pos_center,3,MPI_DOUBLE,maxrank,new_comm);
    MPI_Barrier(new_comm);

    printf("Rank %d: Identifying first guess pos_center: (%f, %f, %f)\n",rank,pos_center[0],pos_center[1],pos_center[2]);
    fflush(stdout);
    double deltaPos, mw_pos_sum[3], mass_sum;
    double mass_ranksum, mw_pos_ranksum[3];
    double old_pos[3];
    double r2;

    double deltaMin = 0.0001;
    int wiggles = 0;
    double dx,dy,dz;
    printf("Rank %d: Beginning Shrinking Sphere Algorithm...\n",rank);
    fflush(stdout);

    while (rShrinkSphere > rMinSphere) {
        fflush(stdout);
        mass_sum=0;
        mw_pos_sum[0]=0;
        mw_pos_sum[1]=0;
        mw_pos_sum[2]=0;
        for (i=0;i<Nstars;i++){
            dx = positions[i][0]-pos_center[0];
            dy = positions[i][1]-pos_center[1];
            dz = positions[i][2]-pos_center[2];
            r2 = dx*dx+dy*dy+dz*dz;


            if (r2 < (rShrinkSphere*rShrinkSphere)){
                mass_sum += masses[i];
                mw_pos_sum[0] += masses[i]*positions[i][0]; 
                mw_pos_sum[1] += masses[i]*positions[i][1]; 
                mw_pos_sum[2] += masses[i]*positions[i][2];
            }
        }

        MPI_Barrier(new_comm); /*Make sure everything is synced up Needed here? */
        MPI_Reduce( &mass_sum, &mass_ranksum, 1, MPI_DOUBLE, MPI_SUM, 0, new_comm);
        MPI_Reduce( &mw_pos_sum, &mw_pos_ranksum, 3, MPI_DOUBLE, MPI_SUM, 0, new_comm);
        MPI_Barrier(new_comm); /*Make sure everything is synced up Needed here? */

        if (group_rank==0){
            pos_center[0] = mw_pos_ranksum[0] / mass_ranksum; //Calculate center of mass
            pos_center[1] = mw_pos_ranksum[1] / mass_ranksum; //Calculate center of mass
            pos_center[2] = mw_pos_ranksum[2] / mass_ranksum; //Calculate center of mass
            printf("Mass_sum is %f, new pos_center is (%f,%f,%f)\n",mass_ranksum,pos_center[0],pos_center[1],pos_center[2]);
        }

        MPI_Bcast(&pos_center,3,MPI_DOUBLE,0,new_comm);
        MPI_Barrier(new_comm); //Make sure all ranks have new position
        rShrinkSphere = shrinkFactor*rShrinkSphere;
        printf("Rank %d: Shrinking Sphere to %f\n",rank,rShrinkSphere);
        fflush(stdout);  
    }

    printf("Rank %d: Beginning sphere wiggling.\n",rank);
    fflush(stdout);
    rShrinkSphere = rMinSphere;
    deltaPos = 1000.0;
    /* Wiggle the sphere without changing the radius to ensure convergence */
    while ((deltaPos > deltaMin*deltaMin) & (wiggles < 10)){
        printf("Rank %d: Wiggling Sphere...\n",rank);
        wiggles+=1;
        mass_sum=0;
        mw_pos_sum[0]=0;
        mw_pos_sum[1]=0;
        mw_pos_sum[2]=0;
        for (i=0;i<Nstars;i++){
            dx = positions[i][0]-pos_center[0];
            dy = positions[i][1]-pos_center[1];
            dz = positions[i][2]-pos_center[2];
            r2 = dx*dx+dy*dy+dz*dz;

            if (r2 < (rShrinkSphere*rShrinkSphere)){
                mass_sum += masses[i];
                mw_pos_sum[0] += masses[i]*positions[i][0]; 
                mw_pos_sum[1] += masses[i]*positions[i][1]; 
                mw_pos_sum[2] += masses[i]*positions[i][2]; 
            }
        }

        MPI_Reduce( &mass_sum, &mass_ranksum, 1, MPI_DOUBLE, MPI_SUM, 0, new_comm);
        MPI_Reduce( &mw_pos_sum, &mw_pos_ranksum, 3, MPI_DOUBLE, MPI_SUM, 0, new_comm);

        old_pos[0] = pos_center[0];
        old_pos[1] = pos_center[1];
        old_pos[2] = pos_center[2];
        printf("Rank %d: old_pos is (%f,%f,%f)\n",rank,old_pos[0],old_pos[1],old_pos[2]);
        if (group_rank==0){
            pos_center[0] = mw_pos_ranksum[0] / mass_ranksum; //Calculate center of mass
            pos_center[1] = mw_pos_ranksum[1] / mass_ranksum; //Calculate center of mass
            pos_center[2] = mw_pos_ranksum[2] / mass_ranksum; //Calculate center of mass
        }
        MPI_Bcast(&pos_center,3,MPI_DOUBLE,0,new_comm);
        MPI_Barrier(new_comm);
        dx = pos_center[0]-old_pos[0];
        dy = pos_center[1]-old_pos[1];
        dz = pos_center[2]-old_pos[2];
        deltaPos = dx*dx+dy*dy+dz*dz;
    }

    printf("Rank %d: Wiggled sphere %d times!\n",rank,wiggles);
    printf("Final Central Position is: (%f, %f, %f)\n",pos_center[0],pos_center[1],pos_center[2]);
    fflush(stdout);
    /*Do Not Change.*/
    to_return[0] = pos_center[0];
    to_return[1] = pos_center[1];
    to_return[2] = pos_center[2];
    /*               */

    //Star positions correct here...
    return to_return;
}

double* center_of_mass_velocity(double *masses, double **velocities, double **positions, double *pos_center, int Npart, int numFilesPerSnap){
    int i;
    double mass_sum = 0;
    double mass_vel[3];
    mass_vel[0] = 0;
    mass_vel[1] = 0;
    mass_vel[2] = 0;
    double rCore = 15;
    double dx,dy,dz,r2;
    double vel_center[3],mass_ranksum,mass_rankvel[3];
    double* toReturn;
    toReturn = (double *) malloc (3 * sizeof (double));

    int rank,group_number,group_rank,numtasks;
    int group_ranks[numFilesPerSnap];

    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

    printf("Rank %d:Starting Velocity Finder...\n",rank);

    MPI_Group orig_group, new_group;
    MPI_Comm new_comm;

    MPI_Comm_group(MPI_COMM_WORLD, &orig_group);
    group_number = rank / numFilesPerSnap;
    group_rank = rank % numFilesPerSnap;

    //Generate the array of ranks assigned to a single snapshot
    for (i=0;i < numFilesPerSnap;i++){
        group_ranks[i] = group_number*numFilesPerSnap + i;
    }

    //Define subgroup of only tasks working on a single snapshot
    MPI_Group_incl(orig_group, numFilesPerSnap, group_ranks, &new_group);

    MPI_Comm_create(MPI_COMM_WORLD, new_group, &new_comm);

    for (i=0;i<Npart;i++){
      dx = positions[i][0]-pos_center[0];
      dy = positions[i][1]-pos_center[1];
      dz = positions[i][2]-pos_center[2];
      r2 = dx*dx+dy*dy+dz*dz;

      if (r2 < (rCore*rCore)){
        mass_sum += masses[i];
        mass_vel[0] += masses[i]*velocities[i][0];
        mass_vel[1] += masses[i]*velocities[i][1];
        mass_vel[2] += masses[i]*velocities[i][2];
      }
    }
    
    MPI_Reduce( &mass_sum, &mass_ranksum, 1, MPI_DOUBLE, MPI_SUM, 0, new_comm);
    MPI_Reduce( &mass_vel, &mass_rankvel, 3, MPI_DOUBLE, MPI_SUM, 0, new_comm);

    if (group_rank==0){
        vel_center[0] = mass_rankvel[0] / mass_ranksum;
        vel_center[1] = mass_rankvel[1] / mass_ranksum;
        vel_center[2] = mass_rankvel[2] / mass_ranksum;
        printf("Vel_center is (%f, %f, %f)\n",vel_center[0],vel_center[1],vel_center[2]);
    }

    MPI_Bcast(&vel_center,3,MPI_DOUBLE,0,new_comm);
    MPI_Barrier(new_comm);

    /*Do Not Change.*/
    toReturn[0] = vel_center[0];
    toReturn[1] = vel_center[1];
    toReturn[2] = vel_center[2];
    /*              */
    return toReturn;
}

double* find_disk_orientation(double *hydrogen_densities, double **gas_positions, double *gas_masses, double **gas_velocities, double *gas_temperatures, int Ngas, double* pos_center, double* vel_center, int numFilesPerSnap)
{

    int rank,group_number,group_rank,numtasks;
    int i;
    int group_ranks[numFilesPerSnap];
    int maxrank,maxindex;
    double unit_L = 3.086*pow(10,21);



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

    //Define subgroup of only tasks working on a single snapshot
    MPI_Group_incl(orig_group, numFilesPerSnap, group_ranks, &new_group);

    MPI_Comm_create(MPI_COMM_WORLD, new_group, &new_comm);
    

    double Lsum[3],Lranksum[3],Lmag;
    double r2,dx,dy,dz,vx,vy,vz;
    double rCore=10;
    double Lhat[3],Lavg[3];
    double* toReturn;
    toReturn = (double *) malloc (3 * sizeof (double));

    for (i=0;i<Ngas;i++) {
        dx = gas_positions[i][0]-pos_center[0];
        dy = gas_positions[i][1]-pos_center[1];
        dz = gas_positions[i][2]-pos_center[2];
        r2 = dx*dx+dy*dy+dz*dz;

        vx = gas_velocities[i][0]-vel_center[0];
        vy = gas_velocities[i][1]-vel_center[1];
        vz = gas_velocities[i][2]-vel_center[2];
      //  h_massfrac = 1 - (gas_metallicities[i][0]+gas_metallicities[i][1]);
      //  nH = (gas_densities[0][i]*h_massfrac) / proton_mass;
        if (r2 < (rCore*rCore) & hydrogen_densities[i] > pow(unit_L,3) & gas_temperatures[i] < 8000) {
            Lsum[0] += gas_masses[i]*(dy*vz - dz*vy);
            Lsum[1] += gas_masses[i]*(dz*vx - dx*vz);
            Lsum[2] += gas_masses[i]*(dx*vy - dy*vx);
        }
        
    }
    MPI_Reduce( &Lsum, &Lranksum, 3, MPI_DOUBLE, MPI_SUM, 0, new_comm);


    if (group_rank==0) {
        Lavg[0] = Lranksum[0];
        Lavg[1] = Lranksum[1];
        Lavg[2] = Lranksum[2];
        Lmag = sqrt( Lavg[0]*Lavg[0]+Lavg[1]*Lavg[1]+Lavg[2]*Lavg[2] );
        Lhat[0] = Lavg[0]/Lmag;
        Lhat[1] = Lavg[1]/Lmag;
        Lhat[2] = Lavg[2]/Lmag;
        printf("Lmag is %f, Lhat is (%f,%f,%f)\n",Lmag,Lhat[0],Lhat[1],Lhat[2]);
    }
    MPI_Bcast(&Lhat,3,MPI_DOUBLE,0,new_comm);
    MPI_Barrier(new_comm);
    toReturn[0]=Lhat[0];
    toReturn[1]=Lhat[1];
    toReturn[2]=Lhat[2];
    return toReturn;
}

