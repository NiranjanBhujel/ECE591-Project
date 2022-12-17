/*
 * Voltage_Control_private.h
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "Voltage_Control".
 *
 * Model version              : 1.284
 * Simulink Coder version : 9.3 (R2020a) 18-Nov-2019
 * C source code generated on : Wed Nov 30 17:24:56 2022
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objective: Execution efficiency
 * Validation result: Not run
 */

#ifndef RTW_HEADER_Voltage_Control_private_h_
#define RTW_HEADER_Voltage_Control_private_h_
#include "rtwtypes.h"
#include "multiword_types.h"
#include <math.h>
#include <stdlib.h>
#include "Voltage_Control.h"

/* Private macros used by the generated code to access rtModel */
#ifndef rtmSetFirstInitCond
# define rtmSetFirstInitCond(rtm, val) ((rtm)->Timing.firstInitCondFlag = (val))
#endif

#ifndef rtmIsFirstInitCond
# define rtmIsFirstInitCond(rtm)       ((rtm)->Timing.firstInitCondFlag)
#endif

#ifndef rtmIsMajorTimeStep
# define rtmIsMajorTimeStep(rtm)       (((rtm)->Timing.simTimeStep) == MAJOR_TIME_STEP)
#endif

#ifndef rtmIsMinorTimeStep
# define rtmIsMinorTimeStep(rtm)       (((rtm)->Timing.simTimeStep) == MINOR_TIME_STEP)
#endif

#ifndef rtmSetTPtr
# define rtmSetTPtr(rtm, val)          ((rtm)->Timing.t = (val))
#endif

#ifndef CodeFormat
#define CodeFormat                     S-Function
#else
#undef CodeFormat
#define CodeFormat                     S-Function
#endif

#ifndef S_FUNCTION_NAME
#define S_FUNCTION_NAME                simulink_only_sfcn
#else
#undef S_FUNCTION_NAME
#define S_FUNCTION_NAME                simulink_only_sfcn
#endif

#ifndef S_FUNCTION_LEVEL
#define S_FUNCTION_LEVEL               2
#else
#undef S_FUNCTION_LEVEL
#define S_FUNCTION_LEVEL               2
#endif

#ifndef RTW_GENERATED_S_FUNCTION
#define RTW_GENERATED_S_FUNCTION
#endif

#ifndef rtmGetDataMapInfo
# define rtmGetDataMapInfo(rtm)        NULL
#endif

#ifndef rtmSetDataMapInfo
# define rtmSetDataMapInfo(rtm, val)
#endif

#if !defined(RTW_SFUNCTION_DEFINES)
#define RTW_SFUNCTION_DEFINES
#ifndef _RTW_COMMON_DEFINES_
#define _RTW_COMMON_DEFINES_
#endif
#endif

extern real_T rt_hypotd_snf(real_T u0, real_T u1);
real_T rt_VTDelayfindtDInterpolate(
  real_T x,real_T* tBuf,real_T* uBuf,real_T* xBuf,int_T bufSz,int_T head,int_T
  tail,int_T* pLast,real_T t,real_T tStart,boolean_T discrete,boolean_T
  minorStepAndTAtLastMajorOutput,real_T initOutput,real_T* appliedDelay);
extern void Voltage__Subsystempi2delay_Init(B_Subsystempi2delay_Voltage_C_T
  *localB, P_Subsystempi2delay_Voltage_C_T *localP);
extern void Volta_Subsystempi2delay_Disable(DW_Subsystempi2delay_Voltage__T
  *localDW);
extern void Voltage_Contr_Subsystempi2delay(RT_MODEL_Voltage_Control_T * const
  Voltage_Control_M, uint8_T rtu_Enable, const real_T rtu_alpha_beta[2], real_T
  rtu_wt, B_Subsystempi2delay_Voltage_C_T *localB,
  DW_Subsystempi2delay_Voltage__T *localDW);
extern void Voltage_Control_Subsystem1_Init(B_Subsystem1_Voltage_Control_T
  *localB, P_Subsystem1_Voltage_Control_T *localP);
extern void Voltage_Cont_Subsystem1_Disable(DW_Subsystem1_Voltage_Control_T
  *localDW);
extern void Voltage_Control_Subsystem1(RT_MODEL_Voltage_Control_T * const
  Voltage_Control_M, uint8_T rtu_Enable, const real_T rtu_alpha_beta[2], real_T
  rtu_wt, B_Subsystem1_Voltage_Control_T *localB,
  DW_Subsystem1_Voltage_Control_T *localDW);

/* private model entry point functions */
extern void Voltage_Control_derivatives(void);

#endif                               /* RTW_HEADER_Voltage_Control_private_h_ */
