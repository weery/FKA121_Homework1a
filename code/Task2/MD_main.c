/*
 MD_main.c

 Created by Anders Lindman on 2013-10-31.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "initfcc.h"
#include "alpotential.h"
#define nbr_of_particles 256
#define nbr_of_timesteps 1000
#define nbr_of_dimensions 3


/* Main program */
int main()
{
    srand(time(NULL));

    /* Simulation parameters */
    double m_AL;
    double cell_length; // Called something else

    double lattice_spacing;
    double initial_displacement;

    int n_unit_cell;
    double lattice_param;

    double timestep;
    int equilibrium_time;
    double temperature_eq;
    double current_time;

    FILE *file1;
    FILE *file2;
    FILE *file3;


    /* Current displacement, velocities, and acceleratons */
    double q[nbr_of_particles][nbr_of_dimensions] = { 0 }; // Displacements
    double v[nbr_of_particles][nbr_of_dimensions] = { 0 }; // Velocities
    double f[nbr_of_particles][nbr_of_dimensions] = { 0 }; // Forces

    /* Allocate memory for large vectors */
    /* Create 3 dimensional data by placing iniitalizeing a 1-dimensional array*/
    #define qq(i,j,k) (disp_arr[nbr_of_particles*nbr_of_dimensions*i+nbr_of_dimensions*j+k])
    double* disp_arr =(double*)malloc(nbr_of_timesteps*nbr_of_particles*nbr_of_dimensions*sizeof(double));

    double* energy =(double*)malloc(nbr_of_timesteps*sizeof(double));
    double* energy_kin = (double*)malloc(nbr_of_timesteps*sizeof(double));
    double* virial =(double*)malloc(nbr_of_timesteps*sizeof(double));


    equilibrium_time = 1000;


    double* energy_equilibrium =(double*)malloc(equilibrium_time*sizeof(double));
    double* energy_kin_equilibrium = (double*)malloc(equilibrium_time*sizeof(double));
    double* virial_equilibrium =(double*)malloc(equilibrium_time*sizeof(double));
    double* temperature_equilibrium=(double*)malloc(equilibrium_time*sizeof(double));
    double* pressure_equilibrium=(double*)malloc(equilibrium_time*sizeof(double));
    #define qq_equilibrium(i,j,k) (disp_arr_equilibrium[nbr_of_particles*nbr_of_dimensions*i+nbr_of_dimensions*j+k])
    double* disp_arr_equilibrium =(double*)malloc(equilibrium_time*nbr_of_particles*nbr_of_dimensions*sizeof(double));
    #define vv_equilibrium(i,j,k) (vel_arr_equilibrium[nbr_of_particles*nbr_of_dimensions*i+nbr_of_dimensions*j+k])
    double* vel_arr_equilibrium =(double*)malloc(equilibrium_time*nbr_of_particles*nbr_of_dimensions*sizeof(double));
    double* alpha_equilibrium =(double*)malloc(equilibrium_time*sizeof(double));


    //TODO go over parameters again
    /* Initialize parameters*/
    initial_displacement = 0.05;
    lattice_param = 4.046; // For aluminium (Å)
    lattice_spacing = lattice_param/sqrt(2.0);
    timestep = 0.01; // 0.1 Bad, 0.01 Seems decent
    m_AL = 0.0027964; // In ASU
    cell_length = 4*lattice_param; // Check this, or 8?
    temperature_eq =500.0;


    // Initialize all displacements, for all times, as 0
    for (int i  = 0; i < nbr_of_timesteps; i++){
        for (int j = 0; j < nbr_of_particles; j++){
            for (int k = 0; k < nbr_of_dimensions; k++){
                qq(i,j,k) = 0;
                vv_equilibrium(i,j,k)=0;
                qq_equilibrium(i,j,k)=0;
            }
        }
    }

    /* Initial conditions */

    /* Put atoms on lattice */
    init_fcc(q, 4, lattice_param);

    for (int i = 0; i < nbr_of_particles; i++){
        for (int j = 0; j < nbr_of_dimensions; j++){

            // Initial perturbation from equilibrium
            q[i][j] +=lattice_param* initial_displacement
                * ((double)rand()/(double)RAND_MAX);
            v[i][j]=0;
        }
        energy_equilibrium[i]=0;
        virial_equilibrium[i]=0;
        energy_kin_equilibrium[i]=0;
        temperature_equilibrium[i]=0;
        pressure_equilibrium[i]=0;
        alpha_equilibrium[i]=1;
    }

    energy_equilibrium[0]=get_energy_AL(q,cell_length,nbr_of_particles);
    virial_equilibrium[0]=get_virial_AL(q,cell_length,nbr_of_particles);
    energy_kin_equilibrium[0]=get_kinetic_AL(v,nbr_of_dimensions,nbr_of_particles,m_AL);
    temperature_equilibrium[0]=instantaneous_temperature(energy_kin_equilibrium[0], nbr_of_particles);
    pressure_equilibrium[0]= instantaneous_pressure(virial[0], temperature_equilibrium[0],  nbr_of_particles, cell_length*cell_length*cell_length);
    alpha_equilibrium[0]=1;

    get_forces_AL(f,q,cell_length,nbr_of_particles);
    for (int i = 0; i < nbr_of_particles; i++){
        for (int j = 0; j < nbr_of_dimensions; j++){
            qq_equilibrium(0,i,j)=q[i][j];
            vv_equilibrium(0,i,j)=v[i][j];
        }
    }

    double alpha_current=1;
    double alpha_last=1;
    /* Simulation */
    /* Simulate to equilibrium*/
    for (int i = 1; i < equilibrium_time; i++)
    {
        /** Verlet algorithm **/
        /* Half step for velocity */
        for (int j = 0; j < nbr_of_particles; j++){
            for (int k = 0; k < nbr_of_dimensions; k++){
                v[j][k] += sqrt(alpha_current)*timestep * 0.5 * f[j][k]/m_AL;
            }
        }

        /* Update displacement*/
        for (int j = 0; j < nbr_of_particles; j++){
            for (int k = 0; k < nbr_of_dimensions; k++){
                q[j][k] += timestep * v[j][k];
            }
        }

        /* Forces */
        get_forces_AL(f,q,cell_length,nbr_of_particles);

        /* Final velocity*/
        for (int j = 0; j < nbr_of_particles; j++){
            for (int k = 0; k < nbr_of_dimensions; k++){
                v[j][k] += sqrt(alpha_current)*timestep * 0.5 * f[j][k]/m_AL;
            }
        }

        energy_equilibrium[i]=get_energy_AL(q,cell_length,nbr_of_particles);
        virial_equilibrium[i]=get_virial_AL(q,cell_length,nbr_of_particles);
        energy_kin_equilibrium[i]=get_kinetic_AL(v,nbr_of_dimensions,nbr_of_particles,m_AL);
        temperature_equilibrium[i]=instantaneous_temperature(energy_kin_equilibrium[i], nbr_of_particles);
        pressure_equilibrium[i]= instantaneous_pressure(virial[i], temperature_equilibrium[i],  nbr_of_particles, cell_length*cell_length*cell_length);
        double tmp_temp = temperature_equilibrium[i];
        for (int j = 0; i < nbr_of_particles; i++){
            for (int k = 0; j < nbr_of_dimensions; j++){
                qq_equilibrium(i,j,k)=q[j][k];
                vv_equilibrium(i,j,k)=v[j][k];
            }
        }
        alpha_last = alpha_current;
        alpha_equilibrium[i]=alpha_last;
        alpha_current = 1+ (temperature_eq-tmp_temp)/(equilibrium_time*tmp_temp);



    }

    file1 = fopen("displacement1.dat","w");

    for (int i = 0; i < equilibrium_time; i ++)
    {
        current_time = i*timestep;
        fprintf(file1, "%.4f \t", current_time );
        for (int j = 0; j < nbr_of_particles; j++)
        {
            for (int k = 0; k < nbr_of_dimensions; k++)
            {
                fprintf(file1, "%.4f \t", qq_equilibrium(i,j,k));
            }
        }
        fprintf(file1, "\n");
    }
    fclose(file1);

    file1 = fopen("velocities1.dat","w");

    for (int i = 0; i < equilibrium_time; i ++)
    {
        current_time = i*timestep;
        fprintf(file1, "%.4f \t", current_time );
        for (int j = 0; j < nbr_of_particles; j++)
        {
            for (int k = 0; k < nbr_of_dimensions; k++)
            {
                fprintf(file1, "%.4f \t", vv_equilibrium(i,j,k));
            }
        }
        fprintf(file1, "\n");
    }
    fclose(file1);

    file2 = fopen("alpha1.dat","w");

    for (int i = 0; i < equilibrium_time; i ++)
    {
        current_time = i*timestep;
        fprintf(file2, "%.4f \t", current_time);
        fprintf(file2, "%.4f \n", alpha_equilibrium[i]);
    }
    fclose(file2);

    /* Save energies to file */
    file2 = fopen("energy1.dat","w");

    for (int i = 0; i < equilibrium_time; i ++)
    {
        current_time = i*timestep;
        fprintf(file2, "%.4f \t", current_time);
        fprintf(file2, "%.4f \t", energy_equilibrium[i]);
        fprintf(file2, "%.4f \n", energy_kin_equilibrium[i]);
    }
    fclose(file2);

    file3 = fopen("temppres1.dat","w");
    for (int i = 0; i < equilibrium_time; i ++)
    {
        current_time = i*timestep;
        fprintf(file3, "%.4f \t", current_time);
        fprintf(file3, "%.4f \t", temperature_equilibrium[i]);
        fprintf(file3, "%.4f \n", pressure_equilibrium[i]);
    }
    fclose(file3);


    free(energy_equilibrium); energy_equilibrium=NULL;
    free(virial_equilibrium); virial_equilibrium=NULL;
    free(energy_kin_equilibrium); energy_kin_equilibrium=NULL;
    free(temperature_equilibrium); temperature_equilibrium=NULL;
    free(pressure_equilibrium); pressure_equilibrium=NULL;
    free(vel_arr_equilibrium); vel_arr_equilibrium=NULL;
    free(alpha_equilibrium);alpha_equilibrium=NULL;


    /* Start saving displacements, energies and virials after equilibrium is done*/
    for (int i = 0; i < nbr_of_particles; i++){
        for (int j = 0; j < nbr_of_dimensions; j++){
            qq(0,i,j)=q[i][j];
        }
    }



    energy[0]=get_energy_AL(q,cell_length,nbr_of_particles);
    virial[0]=get_virial_AL(q,cell_length,nbr_of_particles);
    energy_kin[0]=get_kinetic_AL(v,nbr_of_dimensions,nbr_of_particles,m_AL);

    /* After equilibrium time*/
    for (int i = 1; i < nbr_of_timesteps; i++)
    {
        /** Verlet algorithm **/
        /* Half step for velocity */
        for (int j = 0; j < nbr_of_particles; j++){
            for (int k = 0; k < nbr_of_dimensions; k++){
                v[j][k] += timestep * 0.5 * f[j][k]/m_AL;
            }
        }

        /* Update displacement*/
        for (int j = 0; j < nbr_of_particles; j++){
            for (int k = 0; k < nbr_of_dimensions; k++){
                q[j][k] += timestep * v[j][k];
            }
        }

        /* Forces */
        get_forces_AL(f,q,cell_length,nbr_of_particles);

        /* Final velocity*/
        for (int j = 0; j < nbr_of_particles; j++){
            for (int k = 0; k < nbr_of_dimensions; k++){
                v[j][k] += timestep * 0.5 * f[j][k]/m_AL;
            }
        }

        /* Calculate energy */
        // Potential energy
        energy[i] = get_energy_AL(q,cell_length,nbr_of_particles);
        // Kinetic energy
        energy_kin[i] = get_kinetic_AL(v,nbr_of_dimensions,nbr_of_particles,m_AL);

        virial[i]=get_virial_AL(q,cell_length,nbr_of_particles);

        /* Save current displacements to array*/
        for (int j = 0; j < nbr_of_particles; j++){
            for (int k = 0; k < nbr_of_dimensions; k++){
                qq(i,j,k)=q[j][k];
            }
        }
    }

    /* Save data to file*/
    file1 = fopen("displacement.dat","w");


    for (int i = 0; i < nbr_of_timesteps; i ++)
    {
        current_time = i*timestep;
        fprintf(file1, "%.4f \t", current_time );
        for (int j = 0; j < nbr_of_particles; j++)
        {
            for (int k = 0; k < nbr_of_dimensions; k++)
            {
                fprintf(file1, "%.4f \t", qq(i,j,k));
            }
        }
        fprintf(file1, "\n");
    }
    fclose(file1);

    /* Save energies to file */
    file2 = fopen("energy.dat","w");

    for (int i = 0; i < nbr_of_timesteps; i ++)
    {
        current_time = i*timestep;
        fprintf(file2, "%.4f \t", current_time);
        fprintf(file2, "%.4f \t", energy[i]);
        fprintf(file2, "%.4f \n", energy_kin[i]);
    }
    fclose(file2);

    /* Save energies to file */
    file3 = fopen("virial.dat","w");

    for (int i = 0; i < nbr_of_timesteps; i ++)
    {
        current_time = i*timestep;
        fprintf(file3, "%.4f \t", current_time);
        fprintf(file3, "%.4f \n", virial[i]);
    }
    fclose(file3);

    free(energy_kin); energy_kin=NULL;
    free(energy); energy=NULL;
    free(disp_arr); disp_arr=NULL;
	free(virial); virial=NULL;

    return 0;
}
