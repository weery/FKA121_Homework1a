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


    int k_bins  = 100;

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
            temperature[equil*(nbr_of_timesteps_eq-1) + i] = inst_temperature_eq;
            inst_pressure_eq = instantaneous_pressure(virial_eq, inst_temperature_eq,
                nbr_of_particles, volume);
            pressure[equil*(nbr_of_timesteps_eq-1) + i] = inst_pressure_eq;


            // Update alhpas
            alpha_T = 1.0 + 0.01*(temperature_eq[equil]-inst_temperature_eq)/inst_temperature_eq;
            alpha_P = 1.0 - 0.01*isothermal_compressibility*(pressure_eq - inst_pressure_eq);

            // DEBUG:alpha
            //printf("%.8f \t %.8f \n", alpha_T, alpha_P);

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

    for (int i = 0; i < nbr_of_particles; i++){
        for (int j = 0; j < nbr_of_dimensions; j++){
            qq(0,i,j)=q[i][j];
        }
    }

    // Compute energies, temperature etc. at equilibrium
    energy[0] = get_energy_AL(q, cell_length, nbr_of_particles);
    virial[0] = get_virial_AL(q, cell_length, nbr_of_particles);
    energy_kin[0] = get_kinetic_AL(v, nbr_of_dimensions, nbr_of_particles, m_AL);
    temperature_avg[0] = instantaneous_temperature(energy_kin[0], nbr_of_particles);
    pressure_avg[0] = instantaneous_pressure(virial[0], temperature_avg[0],
    	nbr_of_particles, volume);
    double min = 0;
    double max = sqrt(3*cell_length*cell_length);
    double d_r = (max-min)/(1.0*k_bins);

    #define dists(i,j,k) (dists_arr[nbr_of_particles*nbr_of_particles*i+nbr_of_particles*j+k])
    double* dists_arr = (double*)malloc(nbr_of_timesteps*nbr_of_particles*nbr_of_particles*sizeof(double));

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
        energy[i] = get_energy_AL(q, cell_length, nbr_of_particles);
        // Kinetic energy
        energy_kin[i] = get_kinetic_AL(v, nbr_of_dimensions, nbr_of_particles, m_AL);

        virial[i] = get_virial_AL(q, cell_length, nbr_of_particles);

		// Temperature
        temperature_avg[i] = averaged_temperature(energy_kin, nbr_of_particles, i);
        temperature[2*(nbr_of_timesteps_eq-1) + i] = instantaneous_temperature(energy_kin[i],
            nbr_of_particles);


        // Pressure
        pressure_avg[i] = averaged_pressure(virial, energy_kin, volume, i);
        pressure[2*(nbr_of_timesteps_eq-1) + i] = instantaneous_pressure(virial[i],
            temperature[2*(nbr_of_timesteps_eq-1) + i],
            nbr_of_particles, volume);


        /* Save current displacements to array*/
        for (int j = 0; j < nbr_of_particles; j++){
            for (int k = 0; k < nbr_of_dimensions; k++){
                qq(i,j,k)=q[j][k];
            }
        }

        double distances[nbr_of_particles][nbr_of_particles] = {0};
        for (int l =0 ; l < nbr_of_particles; l++)
        {
            for (int j =0 ; j < nbr_of_particles; j++)
            {
                for (int d = 0; d < nbr_of_dimensions; d++)
                {
                    double q1 = q[l][d];
                    double q2 = q[j][d];
                    q1=boundary_condition(q1,cell_length);
                    q2=boundary_condition(q2,cell_length);

                    distances[l][j] += pow(q1-q2,2);
                }
                distances[l][j] = sqrt(distances[l][j]);
                dists(i,l,j)= distances[l][j];
            }
        }


    } // equilibration/simulation

    printf("Här?");

    // COPY OVER THIS TODO
    // Create Histogram

    double distances[nbr_of_particles][nbr_of_particles] = {0};
    for (int i =0 ; i < nbr_of_particles; i++)
    {
        for (int j =0 ; j < nbr_of_particles; j++)
        {
            for (int d = 0; d < nbr_of_dimensions; d++)
            {
                distances[i][j] += pow(q[i][d]-q[j][d],2);
            }
            distances[i][j] = sqrt(distances[i][j]);
        }
    }

    /*
    // Get min distance between two atoms
    max=0
    min=1e10
    double dist = 0;
    for (int i = 1; i < nbr_of_particles; i++)
        for (int j = i+1; j< nbr_of_particles; j++)
        {
            dist = distances[i][j];
            if (dist < min)
            {
                min = dist;
                printf("%i, %i, %.8f \n", i,j,dist );
            }
        }
    // Get max distance between two atoms
    dist = 0;
    for (int i = 1; i < nbr_of_particles; i++)
        for (int j = i+1; j< nbr_of_particles; j++)
        {
            dist = distances[i][j];
            if (dist > max)
            {
                max = dist;
                printf("%i, %i, %.8f \n", i,j,dist );
            }
        }
        min = 0;
        max = sqrt(3*cell_length*cell_length);
    */

    // Take tke max and min distances to be zero and the diagonal of the cube



// AVG OF distances
double dists_avg[nbr_of_particles][nbr_of_particles];
for (int i = 0; i < nbr_of_particles; i++)
    for (int j = 0; j < nbr_of_particles; j++){
        for (int t = 0; t < nbr_of_timesteps; t++)
        {
            dists_avg[i][j] += dists(t,i,j);
        }
        dists_avg[i][j]/=nbr_of_timesteps;
    }


    int bins[k_bins];
    int bins2[k_bins];
    double Nideal[k_bins];
    for (int i = 0; i < k_bins; i++)
    {
        bins[i]=0;
        bins2[i]=0;
    }

    double factor =((double)(nbr_of_particles-1.0))/volume * 4.0*PI/3.0;
    for (int i = 0; i < k_bins; i++)
    {
        Nideal[i] = factor*(3.0*i*i-3.0*i+1.0)*d_r*d_r*d_r;
    }

    // Use only upper triangle of matrix
    for (int i = 1; i < nbr_of_particles; i++)
        for (int j = 1+i; j < nbr_of_particles; j++)
        {
            int bin = get_bin(dists_avg[i][j],min,max,d_r);
            bins[bin]++;
        }

    for (int i =1 ; i < nbr_of_particles; i++)
    {
        for (int j =i+1 ; j < nbr_of_particles; j++)
        {
            int bin = get_bin(distances[i][j],min,max,d_r);
            bins2[bin]++;
        }
    }

    /* Save data to file*/
    file = fopen("histogram.dat","w");
    for (int i = 0; i < k_bins; i ++)
    {
        fprintf(file, "%e \t %i \t %i \t %e \n",d_r*(i-0.5) ,bins[i], bins2[i], Nideal[i]);
    }
    fclose(file);
    // TO THIS ISH TODO


    printf("Saved files");

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
    printf("Här är jag6\n");
    free(dists_arr);        dists_arr = NULL;

    return 0;
}

int get_bin(double val , double min , double max , double  d_r)
{
    int bin =0;
    double current=min;
    while (current <= val)
    {
        current += d_r;
        bin++;
    }
    return bin;
}

double boundary_condition(double u, double L)
{

    double f = fmod(u,L);
    if (f < 0)
        return -f;
    else
        return f;
}
