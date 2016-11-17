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
#define nbr_of_timesteps 1e4
#define nbr_of_timesteps_eq 4000
#define nbr_of_dimensions 3

double boundary_condition(double,double);



/* Main program */
int main()
{
    srand(time(NULL));

    /* Simulation parameters */
    double m_AL; // Mass of atom
    double cell_length; // Side length of supercell
    double volume;
    double lattice_spacing; // Smallest length between atoms
    double initial_displacement;    // Initial displacement of the atoms from their
                                    // lattice positions
    double lattice_param;   // Lattice parameter, length of each side in the
                            // unit cell
    double timestep;
    double temperature_eq[] = { 1000.0+273.15, 700.0+273.15 };
    double delta_temperature[] = { -5.0, 5.0 };
    double pressure_eq = 101325e-11/1.602; // 1 atm in ASU

    FILE *file;


    /* Current displacement, velocities, and acceleratons */
    double q[nbr_of_particles][nbr_of_dimensions] = { 0 }; // Displacements
    double v[nbr_of_particles][nbr_of_dimensions] = { 0 }; // Velocities
    double f[nbr_of_particles][nbr_of_dimensions] = { 0 }; // Forces

    double heat_capacity;
    double energy_avg[2] = { 0 };
    double temperature_avg[2] = { 0 };


    /* Allocate memory for large vectors */

    double* energy_pot 		= (double*) malloc(nbr_of_timesteps * sizeof(double));
    double* energy_kin 		= (double*) malloc(nbr_of_timesteps * sizeof(double));


    /* Initialize parameters*/
    initial_displacement 	= 0.05;
    lattice_param 			= 4.046; // For aluminium (Å)
    lattice_spacing 		= lattice_param/sqrt(2.0);
    timestep 				= 0.001; // 0.1 Bad, 0.01 Seems decent
    m_AL 					= 0.0027964; // In ASU
    cell_length 			= 4*lattice_param;  // Side of the supercell: The 256 atoms are
                                    			// structured in a block of 4x4x4 unit cells
    volume 					= pow(cell_length, 3);

    

    /* Put atoms on lattice */
    init_fcc(q, 4, lattice_param);


    /* Initial conditions */
    for (int i = 0; i < nbr_of_particles; i++){
        for (int j = 0; j < nbr_of_dimensions; j++){

            // Initial perturbation from equilibrium
            q[i][j] += lattice_spacing * initial_displacement
                * ((double)rand()/(double)RAND_MAX);

        }
    }


    get_forces_AL(f, q, cell_length, nbr_of_particles);

    /* Simulation */
    /* Equilibrium stage */

    double inst_temperature_eq;
    double inst_pressure_eq;
    double alpha_T = 1.0;
    double alpha_P = 1.0;
    double energy_kin_eq = get_kinetic_AL(v,nbr_of_dimensions,nbr_of_particles,m_AL);
    double virial_eq = get_virial_AL(q,cell_length,nbr_of_particles);

    for (int d = 0; d < 2; d++) {

        for (int equil = 0; equil < 2; equil++) {

            double target_temp = temperature_eq[equil] + delta_temperature[d];

            for (int i = 1; i < nbr_of_timesteps_eq; i++)
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
                        v[j][k] += timestep * 0.5* f[j][k]/m_AL;
                    }
                }

                /* Calculate energy */
                // Kinetic energy
                energy_kin_eq = get_kinetic_AL(v, nbr_of_dimensions, nbr_of_particles, m_AL);

                virial_eq = get_virial_AL(q, cell_length, nbr_of_particles);


                inst_temperature_eq = instantaneous_temperature(energy_kin_eq, nbr_of_particles);
                inst_pressure_eq = instantaneous_pressure(virial_eq, inst_temperature_eq,
                    nbr_of_particles, volume);


                // Update alhpas
                alpha_T = 1.0 + 0.01*(target_temp-inst_temperature_eq)/inst_temperature_eq;
                alpha_P = 1.0 - 0.01*(pressure_eq - inst_pressure_eq);

                // Scale velocities
                for (int j = 0; j < nbr_of_particles; j++){
                    for (int k = 0; k < nbr_of_dimensions; k++){
                        v[j][k] *= sqrt(alpha_T);
                    }
                }

                // Scale positions and volume
                cell_length *= pow(alpha_P, 1.0/3.0);
                volume = pow(cell_length, 3);
                for (int j = 0; j < nbr_of_particles; j++) {
                    for (int k = 0; k < nbr_of_dimensions; k++) {
                        q[j][k] *= pow(alpha_P, 1.0/3.0);
                    }
                }

            }
        }


        // Compute energies, temperature etc. at equilibrium
        energy_pot[0] = get_energy_AL(q, cell_length, nbr_of_particles);
        energy_kin[0] = get_kinetic_AL(v, nbr_of_dimensions, nbr_of_particles, m_AL);


        /* Simulation after equilibrium*/
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

            /* Update Forces */
            get_forces_AL(f, q, cell_length, nbr_of_particles);

            /* Final velocity*/
            for (int j = 0; j < nbr_of_particles; j++){
                for (int k = 0; k < nbr_of_dimensions; k++){
                    v[j][k] += timestep * 0.5 * f[j][k]/m_AL;
                }
            }

            /* Calculate energy */
            // Potential energy
            energy_pot[i] = get_energy_AL(q, cell_length, nbr_of_particles);
            // Kinetic energy
            energy_kin[i] = get_kinetic_AL(v, nbr_of_dimensions, nbr_of_particles, m_AL);

            
        } // equilibration/simulation

        // Compute heat capacity
        temperature_avg[d] = averaged_temperature(energy_kin, nbr_of_particles, nbr_of_timesteps-1);
        // Compute average total energy
        for (int i = 0; i < nbr_of_timesteps; i++)
            energy_avg[d] += energy_pot[i] + energy_kin[i];
        energy_avg[d] /= nbr_of_timesteps;

        printf("Temp: %f\nAverage total energy: %.10f\n", temperature_avg[d], energy_avg[d]);

    }

    // Compute heat capacity
    heat_capacity = (energy_avg[1]-energy_avg[0])/(temperature_avg[1]-temperature_avg[0]);

    printf("heat capacity: %f\n", heat_capacity);

    // Save results to file
    file = fopen("heat_capacity.dat", "w");
    fprintf(file, "%.2f\t%e\n", temperature_eq[1], heat_capacity);
    fclose(file);



    free(energy_kin);		energy_kin = NULL;
    free(energy_pot);		energy_pot = NULL;

    return 0;
}

