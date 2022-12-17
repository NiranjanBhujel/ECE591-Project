/*
Filename:       Voltage_Control_interface.c
Written by:     Niranjan Bhujel
Date:           2022-11-30
*/


#include <stdio.h>
#include "sim_code/Voltage_Control.h"
#include "Voltage_Control_interface.h"

#if defined(_WIN32) || defined(_WIN64)
#define EXPORT __declspec(dllexport)
#define IMPORT __declspec(dllimport)
#elif defined(__linux__)
#define EXPORT __attribute__((visibility("default")))
#define IMPORT
#else
#define EXPORT
#define IMPORT
#pragma warning Unknown dynamic link import/export semantics.
#endif

EXPORT void initialize(double *action, double *rand_param, double *t_step, double *rand_input, double *obs, double *y_meas, double *y_true)
{
    // freopen ("nul", "w", stdout);
    
    for (int i = 0; i < 2; i++)
    {
        Voltage_Control_U.action[i] = action[i];
    }
    
    for (int i = 0; i < 3; i++)
    {
        Voltage_Control_U.rand_param[i] = rand_param[i];
    }
    
    Voltage_Control_U.t_step = t_step[0];
    
    for (int i = 0; i < 4; i++)
    {
        Voltage_Control_U.rand_input[i] = rand_input[i];
    }
    
    Voltage_Control_initialize();
    Voltage_Control_output();
    
    for (int i = 0; i < 175; i++)
    {
        obs[i] = Voltage_Control_Y.obs[i];
    }
    for (int i = 0; i < 4; i++)
    {
        y_meas[i] = Voltage_Control_Y.y_meas[i];
    }
    for (int i = 0; i < 4; i++)
    {
        y_true[i] = Voltage_Control_Y.y_true[i];
    }
}

EXPORT void one_step(double *action, double *rand_param, double *t_step, double *rand_input, double *obs, double *y_meas, double *y_true)
{
    // freopen ("nul", "w", stdout);
    
    for (int i = 0; i < 2; i++)
    {
        Voltage_Control_U.action[i] = action[i];
    }
    
    for (int i = 0; i < 3; i++)
    {
        Voltage_Control_U.rand_param[i] = rand_param[i];
    }
    
    Voltage_Control_U.t_step = t_step[0];
    
    for (int i = 0; i < 4; i++)
    {
        Voltage_Control_U.rand_input[i] = rand_input[i];
    }
    
    Voltage_Control_update();
    Voltage_Control_output();
    
    for (int i = 0; i < 175; i++)
    {
        obs[i] = Voltage_Control_Y.obs[i];
    }
    
    for (int i = 0; i < 4; i++)
    {
        y_meas[i] = Voltage_Control_Y.y_meas[i];
    }
    
    for (int i = 0; i < 4; i++)
    {
        y_true[i] = Voltage_Control_Y.y_true[i];
    }
    
}