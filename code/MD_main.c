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
    double* virial =(double*)malloc(nbr_of_timesteps*sizeof(double));

    //TODO go over parameters again
    /* Initialize parameters*/
    initial_displacement = 0.05;
    lattice_param = 4.046; // For aluminium
    lattice_spacing = lattice_param/sqrt(2.0);
    timestep = 0.01;
    m_AL = 0.0027964; // In ASU
    cell_length = 4*lattice_param; // Check this



    // Initialize all displacements, for all times, as 0
    for (int i  = 0; i < nbr_of_timesteps; i++){
        for (int j = 0; j < nbr_of_particles; j++){
            for (int k = 0; k < nbr_of_dimensions; k++){
                qq(i,j,k) = 0;
            }
        }
    }

    /* Initial conditions */
    for (int i = 0; i < nbr_of_particles; i++){
        for (int j = 0; j < nbr_of_dimensions; j++){

            // Initial perturbation from equilibrium
            q[i][j] = lattice_spacing * initial_displacement
                * ((double)rand()/(double)RAND_MAX);
            qq(0,i,j)=q[i][j];
        }
    }

    init_fcc(q, 4, lattice_param);

    energy[0]=get_energy_AL(q,cell_length,nbr_of_particles);
    virial[0]=get_virial_AL(q,cell_length,nbr_of_particles);


    get_forces_AL(f,q,cell_length,nbr_of_particles);


    /* Simulation */
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
        for (int j = 0; j < nbr_of_particles; j++) {
            for (int k = 0; k < nbr_of_dimensions; k++) {
                energy[i] += m_AL * pow(v[j][k], 2) / 2;
            }
        }

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

    double current_time;
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
        fprintf(file2, "%.4f \n", energy[i]);
    }
    fclose(file2);

    /* Save energies to file */
    file3 = fopen("virial.dat","w");

    for (int i = 0; i < nbr_of_timesteps; i ++)
    {
        current_time = i*timestep;
        fprintf(file3, "%.4f \t", current_time);
        fprintf(file3, "%.4f \n", energy[i]);
    }
    fclose(file3);

    free(disp_arr);

    return 0;
}
