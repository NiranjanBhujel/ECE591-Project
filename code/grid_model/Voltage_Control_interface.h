/*
Filename:       Voltage_Control_interface.h
Written by:     Niranjan Bhujel
Date:           2022-11-30
*/


#ifndef Voltage_Control_interface_h
#define Voltage_Control_interface_h

#include "sim_code/Voltage_Control.h"

void initialize(double *action, double *rand_param, double *t_step, double *rand_input, double *obs, double *y_meas, double *y_true);

void one_step(double *action, double *rand_param, double *t_step, double *rand_input, double *obs, double *y_meas, double *y_true);

#endif