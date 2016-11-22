/*
 MD_main.c

 Created by Anders Lindman on 2013-10-31.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "initfcc.h"
#include "alpotential.h"
#define nbr_of_particles 256
#define nbr_of_dimensions 3

/* Main program */
int main()
{
    srand(time(NULL));

    /* Simulation parameters */
    double m_AL; // Mass of atom
    double cell_length; // Side length of supercell

    double lattice_spacing; // Smallest length between atoms
    double initial_displacement;    // Initial displacement of the atoms from their
                                    // lattice positions

    double lattice_param;   // Lattice parameter, length of each side in the
                            // unit cell
    double timesteps[8];

    FILE *file;

    int time_length = 10;

    /* Initialize parameters*/
    initial_displacement = 0.05;
    lattice_param = 4.046; // For aluminium (Å)
    lattice_spacing = lattice_param/sqrt(2.0);
    //timestep = 0.01; // 0.1 Bad, 0.01 Seems decent

    m_AL = 0.0027964; // In ASU
    cell_length = 4*lattice_param;  // Side of the supercell: The 256 atoms are
                                    // structured in a block of 4x4x4 unit cells

    // Test different timestep with 0.01 Å difference
    for (int i = 0; i < 8; i++)
        timesteps[i]=0.005*(i+1);

    for (int t = 0; t < 8; t++)
    {
        // Current timestep and number of timesteps
        double timestep = timesteps[t];
        int nbr_of_timesteps = (int)(time_length/timestep);

        /* Current displacement, velocities, and acceleratons */
        double q[nbr_of_particles][nbr_of_dimensions] = { 0 }; // Displacements
        double v[nbr_of_particles][nbr_of_dimensions] = { 0 }; // Velocities
        double f[nbr_of_particles][nbr_of_dimensions] = { 0 }; // Forces



        /* Allocate memory for large vectors */
        /* Simulate 3 dimensional data by placing iniitalizeing a 1-dimensional array*/
        double* energy_pot =(double*)malloc(nbr_of_timesteps*sizeof(double));
        double* energy_kin = (double*)malloc(nbr_of_timesteps*sizeof(double));

        /* Put atoms on lattice */
        init_fcc(q, 4, lattice_param);

        /* Initial conditions */

        for (int i = 0; i < nbr_of_particles; i++){
            for (int j = 0; j < nbr_of_dimensions; j++){

                // Initial perturbation from equilibrium
                q[i][j] +=lattice_spacing* initial_displacement
                    * ((double)rand()/(double)RAND_MAX);
            }
        }
        energy_pot[0]=get_energy_AL(q,cell_length,nbr_of_particles);
        energy_kin[0]=get_kinetic_AL(v,nbr_of_dimensions,nbr_of_particles,m_AL);

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
            energy_pot[i] = get_energy_AL(q,cell_length,nbr_of_particles);
            // Kinetic energy
            energy_kin[i] = get_kinetic_AL(v,nbr_of_dimensions,nbr_of_particles,m_AL);
        }
        char str[80];
        char S[3];
        sprintf(S, "%.3f", timestep);

          strcpy (str,"data/energy");
          strcat (str,S);
          strcat (str,".dat");


        /* Save energies to file */
        file = fopen(str,"w");

        double current_time;
        for (int i = 0; i < nbr_of_timesteps; i ++)
        {
            current_time = i*timestep;
            fprintf(file, "%.4f \t", current_time);
            fprintf(file, "%.4f \t", energy_pot[i]);
            fprintf(file, "%.4f \n", energy_kin[i]);
        }
        fclose(file);



        free(energy_kin); energy_kin=NULL;
        free(energy_pot); energy_pot=NULL;
    }

    return 0;
}
