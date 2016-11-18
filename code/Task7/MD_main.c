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
#define nbr_of_timesteps_eq 4000
#define nbr_of_dimensions 3

#define PI 3.141592653589
int get_bin(double , double , double , double );

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
    double temperature_eq[] = { 1500.0+273.15, 700.0+273.15 };
    double pressure_eq = 101325e-11/1.602; // 1 atm in ASU
    double isothermal_compressibility = 1.0; //0.8645443196; // 1.385e-11 m^2/N = 1.385/1.602 Å^3/eV

    FILE *file;


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
    double* temperature     = (double*) malloc((2 * nbr_of_timesteps_eq + nbr_of_timesteps) * sizeof(double));
    double* pressure        = (double*) malloc((2 * nbr_of_timesteps_eq + nbr_of_timesteps) * sizeof(double));


    int k_bins  = 250;

    //TODO go over parameters again
    /* Initialize parameters*/
    initial_displacement 	= 0.05;
    lattice_param 			= 4.046; // For aluminium (Å)
    lattice_spacing 		= lattice_param/sqrt(2.0);
    timestep 				= 0.01; // 0.1 Bad, 0.01 Seems decent
    m_AL 					= 0.0027964; // In ASU
    cell_length 			= 4*lattice_param;  // Side of the supercell: The 256 atoms are
                                    			// structured in a block of 4x4x4 unit cells
    volume 					= pow(cell_length, 3);

    // Initialize all displacements, for all times, as 0
    for (int i  = 0; i < nbr_of_timesteps; i++) {
        for (int j = 0; j < nbr_of_particles; j++) {
            for (int k = 0; k < nbr_of_dimensions; k++) {
                qq(i,j,k) = 0;
            }
        }
    }

    /* Put atoms on lattice */
    init_fcc(q, 4, lattice_param);


    /* Initial conditions */
    for (int i = 0; i < nbr_of_particles; i++) {
        for (int j = 0; j < nbr_of_dimensions; j++) {

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

    temperature[0]  = instantaneous_temperature(energy_kin_eq, nbr_of_particles);
    pressure[0]     = instantaneous_pressure(virial_eq, temperature[0], nbr_of_particles, volume);

    for (int equil = 0; equil < 2; equil++) {
        for (int i = 1; i < nbr_of_timesteps_eq; i++) {

            /** Verlet algorithm **/
            /* Half step for velocity */
            for (int j = 0; j < nbr_of_particles; j++) {
            	for (int k = 0; k < nbr_of_dimensions; k++) {
            		v[j][k] += timestep * 0.5 * f[j][k]/m_AL;
                }
            }

            /* Update displacement*/
            for (int j = 0; j < nbr_of_particles; j++) {
                for (int k = 0; k < nbr_of_dimensions; k++) {
                    q[j][k] += timestep * v[j][k];
                }
            }

            /* Forces */
            get_forces_AL(f,q,cell_length,nbr_of_particles);

            /* Final velocity*/
            for (int j = 0; j < nbr_of_particles; j++) {
                for (int k = 0; k < nbr_of_dimensions; k++) {
                    v[j][k] += timestep * 0.5* f[j][k]/m_AL;
                }
            }

            /* Calculate energy */
            // Kinetic energy
            energy_kin_eq = get_kinetic_AL(v, nbr_of_dimensions, nbr_of_particles, m_AL);

            virial_eq = get_virial_AL(q, cell_length, nbr_of_particles);


            inst_temperature_eq = instantaneous_temperature(energy_kin_eq, nbr_of_particles);
            temperature[equil*(nbr_of_timesteps_eq-1) + i] = inst_temperature_eq;
            inst_pressure_eq = instantaneous_pressure(virial_eq, inst_temperature_eq,
                nbr_of_particles, volume);
            pressure[equil*(nbr_of_timesteps_eq-1) + i] = inst_pressure_eq;


            // Update alhpas
            alpha_T = 1.0 + 0.01*(temperature_eq[equil]-inst_temperature_eq)/inst_temperature_eq;
            alpha_P = 1.0 - 0.01*isothermal_compressibility*(pressure_eq - inst_pressure_eq);


            // Scale velocities
            for (int j = 0; j < nbr_of_particles; j++) {
                for (int k = 0; k < nbr_of_dimensions; k++) {
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

    for (int i = 0; i < nbr_of_particles; i++) {
        for (int j = 0; j < nbr_of_dimensions; j++) {
            qq(0,i,j)=q[i][j];
        }
    }


    // Compute energies, temperature etc. at equilibrium
    double min = 0.0;
    double max = sqrt(3*cell_length*cell_length);
    double d_r = (max-min)/(1.0*k_bins);
    int bins[k_bins];
    int* bins2 = (int*) malloc(k_bins * sizeof(int));

    for (int i = 0; i < k_bins; i++) {
        bins[i]=0;
        bins2[i]=0;
    }

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

        virial[i]=get_virial_AL(q,cell_length,nbr_of_particles);

        /* Save current displacements to array*/
        for (int j = 0; j < nbr_of_particles; j++){
            for (int k = 0; k < nbr_of_dimensions; k++){
                qq(i,j,k)=q[j][k];
            }
        }
    }


    // COPY OVER THIS TODO
    // Create Histogram

    printf("Debug");

    for (int i = 1; i < nbr_of_timesteps; i++)
    {
        for (int j = 1 ; j < nbr_of_particles; j++) {
            for (int k = j+1 ; k < nbr_of_particles; k++) {
                double sum =  0;
                for (int d = 0; d < nbr_of_dimensions; d++) {
                    double q1 = qq(i,j,d);
                    double q2 = qq(i,k,d);
                    q1=boundary_condition(q1,cell_length);
                    q2=boundary_condition(q2,cell_length);

                    sum += pow(q1-q2,2);
                }
                sum = sqrt(sum);
                int bin = get_bin(sum,min,max,d_r);
                if (bin >= k_bins)
                    printf("Fuck");
                bins2[bin]++;
            }
        }
    }
    double Nideal[k_bins];
    double factor =((double)(nbr_of_particles-1.0))/volume * 4.0*PI/3.0;
    for (int i = 0; i < k_bins; i++) {
        Nideal[i] = factor*(3.0*i*i-3.0*i+1.0)*d_r*d_r*d_r;
    }




    /* Save data to file*/
    file = fopen("histogram.dat","w");
    for (int i = 0; i < k_bins; i ++) {
        fprintf(file, "%e \t %i \t %i \t %e \n",d_r*(i-0.5), bins[i],bins2[i], Nideal[i]);
    }
    fclose(file);
    // TO THIS ISH TODO

    printf("Saved files\n");

    free(energy_kin);		energy_kin = NULL;
    printf("Här är jag1\n");
    free(energy); 			energy = NULL;
    printf("Här är jag2\n");
    free(disp_arr); 		disp_arr = NULL;
    printf("Här är jag3\n");
	free(virial); 			virial = NULL;
    printf("Här är jag4\n");
	free(temperature_avg); 	temperature_avg = NULL;
    printf("Här är jag5\n");
	free(pressure_avg);		pressure_avg = NULL;

    return 0;
}

int get_bin(double val , double min , double max , double  d_r)
{
    int bin = 0;
    double current = min;
    while (current <= val)
    {
        current += d_r;
        bin++;
        //printf("current: %.8f \t Val: %.8f, d_r: %.8f \n", current, val, d_r);
    }
    if (current > max)
        return --bin;
    return bin;
}

double boundary_condition(double u, double L)
{
    double sum = u;
    while (sum > 0)
    {
        sum -= L;
    }
    return sum +L;
}
