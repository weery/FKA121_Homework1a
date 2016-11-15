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
#define nbr_of_timesteps 10000
#define nbr_of_timesteps_eq 5000
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
    double temperature_eq = 500.0+273.15;

    FILE *file1;
    FILE *file2;
    FILE *file3;


    /* Current displacement, velocities, and acceleratons */
    double q[nbr_of_particles][nbr_of_dimensions] = { 0 }; // Displacements
    double v[nbr_of_particles][nbr_of_dimensions] = { 0 }; // Velocities
    double f[nbr_of_particles][nbr_of_dimensions] = { 0 }; // Forces

    /* Allocate memory for large vectors */
    /* Simulate 3 dimensional data by placing iniitalizeing a 1-dimensional array*/
    #define qq(i,j,k) (disp_arr[nbr_of_particles*nbr_of_dimensions*i+nbr_of_dimensions*j+k])
    double* disp_arr = (double*)malloc(nbr_of_timesteps*nbr_of_particles*nbr_of_dimensions*sizeof(double));

    double* energy 			= (double*) malloc(nbr_of_timesteps * sizeof(double));
    double* energy_kin 		= (double*) malloc(nbr_of_timesteps * sizeof(double));
    double* virial 			= (double*) malloc(nbr_of_timesteps * sizeof(double));
    double* temperature_avg = (double*) malloc(nbr_of_timesteps * sizeof(double));
    double* pressure_avg 	= (double*) malloc(nbr_of_timesteps * sizeof(double));

    //TODO go over parameters again
    /* Initialize parameters*/
    initial_displacement 	= 0.05;
    lattice_param 			= 4.046; // For aluminium (Å)
    lattice_spacing 		= lattice_param/sqrt(2.0);
    timestep 				= 0.0001; // 0.1 Bad, 0.01 Seems decent
    m_AL 					= 0.0027964; // In ASU
    cell_length 			= 4*lattice_param;  // Side of the supercell: The 256 atoms are
                                    			// structured in a block of 4x4x4 unit cells
    volume 					= pow(cell_length, 3);

    // Initialize all displacements, for all times, as 0
    for (int i  = 0; i < nbr_of_timesteps; i++){
        for (int j = 0; j < nbr_of_particles; j++){
            for (int k = 0; k < nbr_of_dimensions; k++){
                qq(i,j,k) = 0;
            }
        }
    }

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


    for (int i = 0; i < nbr_of_particles; i++){
        for (int j = 0; j < nbr_of_dimensions; j++){
            qq(0,i,j)=q[i][j];
        }
    }


    get_forces_AL(f,q,cell_length,nbr_of_particles);

    /* Simulation */
    /* Equilibrium stage */
    double energy_eq = get_energy_AL(q,cell_length,nbr_of_particles);
    // Kinetic energy
    double energy_kin_eq = get_kinetic_AL(v,nbr_of_dimensions,nbr_of_particles,m_AL);
    double virial_eq = get_virial_AL(q,cell_length,nbr_of_particles);
    double inst_temperature_eq;
    double alpha = 1;
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
                q[j][k] = boundary_condition(q[j][k],cell_length);
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
        // Potential energy
        energy_eq = get_energy_AL(q,cell_length,nbr_of_particles);
        // Kinetic energy
        energy_kin_eq = get_kinetic_AL(v,nbr_of_dimensions,nbr_of_particles,m_AL);

        virial_eq = get_virial_AL(q,cell_length,nbr_of_particles);


        inst_temperature_eq = instantaneous_temperature(energy_kin_eq, nbr_of_particles);


        alpha = 1+1.0/100.0*(temperature_eq-inst_temperature_eq)/inst_temperature_eq;
        //alpha = 1+1/(double)(nbr_of_timesteps)*(temperature_eq-inst_temperature_eq)/inst_temperature_eq;


        //printf("Alpha: %.4f\n",alpha );
        for (int j = 0; j < nbr_of_particles; j++){
            for (int k = 0; k < nbr_of_dimensions; k++){
                v[j][k] *= sqrt(alpha);
            }
        }
    }

    // Compute energies, temperature etc. at equilibrium
    energy[0] = get_energy_AL(q,cell_length,nbr_of_particles);
    virial[0] = get_virial_AL(q,cell_length,nbr_of_particles);
    energy_kin[0] = get_kinetic_AL(v,nbr_of_dimensions,nbr_of_particles,m_AL);
    temperature_avg[0] = instantaneous_temperature(energy_kin[0], nbr_of_particles);
    pressure_avg[0] = instantaneous_pressure(virial[0], temperature_avg[0],
    	nbr_of_particles, volume);
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
                q[j][k] = boundary_condition(q[j][k],cell_length);
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

        virial[i] = get_virial_AL(q,cell_length,nbr_of_particles);

		// Temperature
        temperature_avg[i] = averaged_temperature(energy_kin, nbr_of_particles,
        	timestep, i);


        // Pressure
        pressure_avg[i] = averaged_pressure(virial, energy_kin, volume,
        	timestep, i);


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

    // Save temperature to file
    file3 = fopen("temperature.dat", "w");
    for (int i = 0; i < nbr_of_timesteps; i++)
    {

    	current_time = i*timestep;
    	fprintf(file3, "%.3f \t %e\n", current_time, temperature_avg[i]);
    }
    fclose(file3);

    // Save pressure to file
    file3 = fopen("pressure.dat", "w");
    for (int i = 0; i < nbr_of_timesteps; i++)
    {
    	current_time = i*timestep;
    	fprintf(file3, "%.3f \t%e \n", current_time, pressure_avg[i]);
    }
    fclose(file3);

    free(energy_kin); energy_kin=NULL;
    free(energy); energy=NULL;
    free(disp_arr); disp_arr=NULL;
	free(virial); virial=NULL;
	free(temperature_avg); temperature_avg = NULL;
	free(pressure_avg); pressure_avg = NULL;

    return 0;
}


double boundary_condition(double u, double L)
{
    return fmod(u,L);
}
