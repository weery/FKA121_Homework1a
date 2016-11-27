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
#include <complex.h>
#define nbr_of_particles 256
#define nbr_of_timesteps 1e4
#define nbr_of_timesteps_eq 6000
#define nbr_of_dimensions 3

#define PI 3.141592653589
int get_bin(double , double , double , double );
double get_length(double*, int);

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


    FILE *file;


    /* Current displacement, velocities, and acceleratons */
    double q[nbr_of_particles][nbr_of_dimensions] = { 0 }; // Displacements
    double v[nbr_of_particles][nbr_of_dimensions] = { 0 }; // Velocities
    double f[nbr_of_particles][nbr_of_dimensions] = { 0 }; // Forces

    /* Allocate memory for large vectors */
    /* Simulate 3 dimensional data by placing iniitalizeing a 1-dimensional array */
    #define qq(i,j,k) (disp_arr[nbr_of_particles*nbr_of_dimensions*i+nbr_of_dimensions*j+k])
    double* disp_arr = (double*)malloc(nbr_of_timesteps*nbr_of_particles*nbr_of_dimensions*sizeof(double));

    double* energy 			= (double*) malloc(nbr_of_timesteps * sizeof(double));
    double* energy_kin 		= (double*) malloc(nbr_of_timesteps * sizeof(double));
    double* virial 			= (double*) malloc(nbr_of_timesteps * sizeof(double));
    double* temperature_avg = (double*) malloc(nbr_of_timesteps * sizeof(double));
    double* pressure_avg 	= (double*) malloc(nbr_of_timesteps * sizeof(double));
    double* temperature     = (double*) malloc((2 * nbr_of_timesteps_eq + nbr_of_timesteps) * sizeof(double));
    double* pressure        = (double*) malloc((2 * nbr_of_timesteps_eq + nbr_of_timesteps) * sizeof(double));

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

            /* Update displacement */
            for (int j = 0; j < nbr_of_particles; j++){
                for (int k = 0; k < nbr_of_dimensions; k++){
                    q[j][k] += timestep * v[j][k];
                }
            }

            /* Forces */
            get_forces_AL(f,q,cell_length,nbr_of_particles);

            /* Final velocity */
            for (int j = 0; j < nbr_of_particles; j++){
                for (int k = 0; k < nbr_of_dimensions; k++){
                    v[j][k] += timestep * 0.5* f[j][k]/m_AL;
                }
            }

            /* Calculate energy */
            /* Kinetic energy */
            energy_kin_eq = get_kinetic_AL(v, nbr_of_dimensions, nbr_of_particles, m_AL);

            virial_eq = get_virial_AL(q, cell_length, nbr_of_particles);


            inst_temperature_eq = instantaneous_temperature(energy_kin_eq, nbr_of_particles);
            temperature[equil*(nbr_of_timesteps_eq-1) + i] = inst_temperature_eq;
            inst_pressure_eq = instantaneous_pressure(virial_eq, inst_temperature_eq,
                nbr_of_particles, volume);
            pressure[equil*(nbr_of_timesteps_eq-1) + i] = inst_pressure_eq;


            // Update alphas
            alpha_T = 1.0 + 0.01*(temperature_eq[equil]-inst_temperature_eq)/inst_temperature_eq;
            //alpha_P = 1.0 - 0.01*isothermal_compressibility*(pressure_eq - inst_pressure_eq);
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
    } // equilibration/simulation

    int n_x = 40;
    int n_y = 40;
    int n_z = 40;

    double factor = PI*2.0/cell_length;

    #define qs(i,d) (qs_arr[(nbr_of_dimensions)*i+d])
    double* qs_arr = (double*)malloc((2*n_x+1)*(2*n_y+1)*(2*n_z+1)*(nbr_of_dimensions)*sizeof(double));

    int nk = 0;

    for (int i = -n_x; i <= n_x; i++)
        for (int j = -n_y; j<= n_y; j++)
            for (int k = -n_z; k <= n_z ; k++)
            {
                qs(nk,0)=i*factor;
                qs(nk,1)=j*factor;
                qs(nk,2)=k*factor;
                nk++;
            }

    #define s(i,q) (s_arr[3*i+q])
    double* s_arr = (double*)malloc(nk*(3)*sizeof(double));



    double max =0;
    for (int i = 0; i < nk; i++)
    {
        printf("Current nk: %i out of: %i\n", i, nk );
        double len_sq = 0;
        for (int d = 0; d < nbr_of_dimensions; d++)
        {
            len_sq += qs(i,d)*qs(i,d);
        }
        len_sq = sqrt(len_sq);
        if (len_sq> max)
        {
            max = len_sq;
            printf("Current max: %e\n",max );
        }
        s(i,0)=0;
        s(i,1)=0;
        for (int t = 9000; t <nbr_of_timesteps;t++ )
        {
            double sum1 = 0;
            double sum2 = 0;
            for (int r=0; r < nbr_of_particles; r++)
            {
                double expo=0;
                for (int d = 0; d < nbr_of_dimensions; d++)
                {
                    double ri = qq(t,r,d);
                    //ri=boundary_condition(ri,cell_length);
                    expo+= qs(i,d)*ri;
                }
                sum1+= cos(expo);
                sum2+= sin(expo);
            }
            sum1=abs(sum1*sum1);
            sum2=abs(sum2*sum2);

            s(i,0)+=sum1;
            s(i,1)+=sum2;
        }
        s(i,0)/=nbr_of_particles;
        s(i,1)/=nbr_of_particles;
        s(i,0)/=1000;
        s(i,1)/=1000;
        s(i,2) = len_sq;
    }


    file = fopen("sq.dat","w");
    for (int i = 0; i < nk; i ++) {
        fprintf(file, "%e \t %e \t %e \n", s(i,0), s(i,1),s(i,2));
    }
    fclose(file);

    printf("Data saved \n");

    double min =0;
    int k_bins = 100;
    //int* bins = (double*)malloc(k_bins*sizeof(double));
    double* bins = (double*)malloc(k_bins*sizeof(double));
    double d_r = min + (max-min)/k_bins;
    for (int i = 1; i < nbr_of_timesteps; i++)
    {
        double dist = s(i,2);
        int bin = get_bin(dist,min,max,d_r);
        if (bin < k_bins)
        {
            //bins[bin]++;
            bins[bin]+= s(i,1)+s(i,0);
        }
    }

    file = fopen("sq_bin.dat","w");
    for (int i = 0; i < k_bins; i++ )
    {
        fprintf(file, "%e \t %e \n", d_r*(i-0.5),bins[i]);
    }


    free(energy_kin);		energy_kin = NULL;
    free(energy); 			energy = NULL;
    free(disp_arr); 		disp_arr = NULL;
	free(virial); 			virial = NULL;
	free(temperature_avg); 	temperature_avg = NULL;
	free(pressure_avg);		pressure_avg = NULL;
    free(s_arr);            s_arr = NULL;
    free(qs_arr);           qs_arr= NULL;

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
    double f = u/L;
    f-= floor(f);
    return f*L;
}

double get_length(double* arr, int N)
{
    double len = 0;
    for (int i = 0; i < N;  i++)
        len += arr[i]*arr[i];
    return sqrt(len);
}
