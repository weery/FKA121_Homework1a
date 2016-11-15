/*
alpotential.h

Created by Anders Lindman on 2013-03-15.
*/

#ifndef _alpotential_h
#define _alpotential_h

extern void get_forces_AL(double[][3] , double[][3], double, int);
extern double get_energy_AL(double[][3], double, int);
extern double get_kinetic_AL(double[][3], int, int, double);
extern double get_virial_AL(double[][3], double, int);
extern double instantaneous_temperature(double,int);
extern double averaged_temperature(double*,int,double,int);
extern double instantaneous_pressure(double,double,int,double);
extern double averaged_pressure(double*,double*,double,double,int);

#endif
