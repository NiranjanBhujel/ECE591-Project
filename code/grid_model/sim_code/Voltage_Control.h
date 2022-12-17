/*
 * Voltage_Control.h
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

#ifndef RTW_HEADER_Voltage_Control_h_
#define RTW_HEADER_Voltage_Control_h_
#include <string.h>
#include <math.h>
#ifndef Voltage_Control_COMMON_INCLUDES_
# define Voltage_Control_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#endif                                 /* Voltage_Control_COMMON_INCLUDES_ */

#include "Voltage_Control_types.h"

/* Shared type includes */
#include "multiword_types.h"
#include "rtGetInf.h"
#include "rt_nonfinite.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetContStateDisabled
# define rtmGetContStateDisabled(rtm)  ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
# define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
# define rtmGetContStates(rtm)         ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
# define rtmSetContStates(rtm, val)    ((rtm)->contStates = (val))
#endif

#ifndef rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag
# define rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm) ((rtm)->CTOutputIncnstWithState)
#endif

#ifndef rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag
# define rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm, val) ((rtm)->CTOutputIncnstWithState = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
# define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
# define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetIntgData
# define rtmGetIntgData(rtm)           ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
# define rtmSetIntgData(rtm, val)      ((rtm)->intgData = (val))
#endif

#ifndef rtmGetOdeF
# define rtmGetOdeF(rtm)               ((rtm)->odeF)
#endif

#ifndef rtmSetOdeF
# define rtmSetOdeF(rtm, val)          ((rtm)->odeF = (val))
#endif

#ifndef rtmGetOdeY
# define rtmGetOdeY(rtm)               ((rtm)->odeY)
#endif

#ifndef rtmSetOdeY
# define rtmSetOdeY(rtm, val)          ((rtm)->odeY = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
# define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
# define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
# define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
# define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetZCCacheNeedsReset
# define rtmGetZCCacheNeedsReset(rtm)  ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
# define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetdX
# define rtmGetdX(rtm)                 ((rtm)->derivs)
#endif

#ifndef rtmSetdX
# define rtmSetdX(rtm, val)            ((rtm)->derivs = (val))
#endif

#ifndef rtmGetErrorStatus
# define rtmGetErrorStatus(rtm)        ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
# define rtmSetErrorStatus(rtm, val)   ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
# define rtmGetStopRequested(rtm)      ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
# define rtmSetStopRequested(rtm, val) ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
# define rtmGetStopRequestedPtr(rtm)   (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
# define rtmGetT(rtm)                  (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTPtr
# define rtmGetTPtr(rtm)               ((rtm)->Timing.t)
#endif

/* Block signals for system '<S49>/Subsystem - pi//2 delay' */
typedef struct {
  real_T Fcn;                          /* '<S53>/Fcn' */
  real_T Fcn1;                         /* '<S53>/Fcn1' */
} B_Subsystempi2delay_Voltage_C_T;

/* Block states (default storage) for system '<S49>/Subsystem - pi//2 delay' */
typedef struct {
  boolean_T Subsystempi2delay_MODE;    /* '<S49>/Subsystem - pi//2 delay' */
} DW_Subsystempi2delay_Voltage__T;

/* Block signals for system '<S49>/Subsystem1' */
typedef struct {
  real_T Fcn;                          /* '<S54>/Fcn' */
  real_T Fcn1;                         /* '<S54>/Fcn1' */
} B_Subsystem1_Voltage_Control_T;

/* Block states (default storage) for system '<S49>/Subsystem1' */
typedef struct {
  boolean_T Subsystem1_MODE;           /* '<S49>/Subsystem1' */
} DW_Subsystem1_Voltage_Control_T;

/* Block signals (default storage) */
typedef struct {
  real_T Product[3];                   /* '<S4>/Product' */
  real_T Initial;                      /* '<S37>/Initial' */
  real_T Integrator;                   /* '<S37>/Integrator' */
  real_T Gain3[3];                     /* '<S133>/Gain3' */
  real_T UnitDelay;                    /* '<S67>/Unit Delay' */
  real_T UnitDelay_b;                  /* '<S68>/Unit Delay' */
  real_T UnitDelay_i;                  /* '<S69>/Unit Delay' */
  real_T StateSpace[24];               /* '<S139>/State-Space' */
  real_T Gain1[3];                     /* '<S122>/Gain1' */
  real_T Kv1[3];                       /* '<S10>/Kv1' */
  real_T Gain1_l[3];                   /* '<S128>/Gain1' */
  real_T Sum[4];                       /* '<Root>/Sum' */
  real_T donotdeletethisgain;          /* '<S71>/do not delete this gain' */
  real_T Divide;                       /* '<S67>/Divide' */
  real_T donotdeletethisgain_a;        /* '<S74>/do not delete this gain' */
  real_T Divide_b;                     /* '<S68>/Divide' */
  real_T donotdeletethisgain_p;        /* '<S77>/do not delete this gain' */
  real_T Divide_n;                     /* '<S69>/Divide' */
  real_T Memory;                       /* '<S37>/Memory' */
  real_T integrator;                   /* '<S59>/integrator' */
  real_T Memory_h;                     /* '<S59>/Memory' */
  real_T Switch;                       /* '<S59>/Switch' */
  real_T Kp5;                          /* '<S39>/Kp5' */
  real_T Kp6;                          /* '<S39>/Kp6' */
  real_T Saturation2;                  /* '<S39>/Saturation2' */
  real_T RateLimiter;                  /* '<S37>/Rate Limiter' */
  real_T x1;                           /* '<S55>/A*x1+ B* u1' */
  real_T x2;                           /* '<S55>/A*x2+ B*u2' */
  real_T y;                            /* '<S55>/C*x + D*u' */
  real_T period;                       /* '<S59>/period' */
  real_T Switch_o[2];                  /* '<S60>/Switch' */
  real_T Gain1_a;                      /* '<S4>/Gain1' */
  real_T Gain3_a;                      /* '<S4>/Gain3' */
  real_T Fcn;                          /* '<S138>/Fcn' */
  real_T Fcn1;                         /* '<S138>/Fcn1' */
  real_T Fcn_b;                        /* '<S137>/Fcn' */
  real_T Fcn1_b;                       /* '<S137>/Fcn1' */
  real_T Switch_ow[2];                 /* '<S49>/Switch' */
  real_T integrator_p;                 /* '<S47>/integrator' */
  real_T Memory_p;                     /* '<S47>/Memory' */
  real_T Switch_j;                     /* '<S47>/Switch' */
  real_T integrator_m;                 /* '<S48>/integrator' */
  real_T Memory_j;                     /* '<S48>/Memory' */
  real_T Switch_i;                     /* '<S48>/Switch' */
  real_T MathFunction;                 /* '<S38>/Math Function' */
  real_T period_p;                     /* '<S47>/period' */
  real_T period_a;                     /* '<S48>/period' */
  uint8_T Compare;                     /* '<S135>/Compare' */
  uint8_T Compare_p;                   /* '<S136>/Compare' */
  uint8_T Compare_m;                   /* '<S123>/Compare' */
  uint8_T Compare_a;                   /* '<S124>/Compare' */
  uint8_T Compare_f;                   /* '<S129>/Compare' */
  uint8_T Compare_a5;                  /* '<S130>/Compare' */
  uint8_T Compare_o;                   /* '<S62>/Compare' */
  uint8_T Compare_n;                   /* '<S63>/Compare' */
  uint8_T Compare_l;                   /* '<S52>/Compare' */
  uint8_T Compare_pa;                  /* '<S51>/Compare' */
  boolean_T RelationalOperator;        /* '<S37>/Relational Operator' */
  B_Subsystem1_Voltage_Control_T Subsystem1_ex;/* '<S127>/Subsystem1' */
  B_Subsystempi2delay_Voltage_C_T Subsystempi2delay_n;/* '<S127>/Subsystem - pi//2 delay' */
  B_Subsystem1_Voltage_Control_T Subsystem1_e;/* '<S121>/Subsystem1' */
  B_Subsystempi2delay_Voltage_C_T Subsystempi2delay_b;/* '<S121>/Subsystem - pi//2 delay' */
  B_Subsystem1_Voltage_Control_T Subsystem1;/* '<S60>/Subsystem1' */
  B_Subsystempi2delay_Voltage_C_T Subsystempi2delay;/* '<S60>/Subsystem - pi//2 delay' */
  B_Subsystem1_Voltage_Control_T Subsystem1_h;/* '<S49>/Subsystem1' */
  B_Subsystempi2delay_Voltage_C_T Subsystempi2delay_d;/* '<S49>/Subsystem - pi//2 delay' */
} B_Voltage_Control_T;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  real_T UnitDelay_DSTATE;             /* '<S67>/Unit Delay' */
  real_T UnitDelay_DSTATE_h;           /* '<S68>/Unit Delay' */
  real_T UnitDelay_DSTATE_i;           /* '<S69>/Unit Delay' */
  real_T StateSpace_DSTATE[6];         /* '<S139>/State-Space' */
  real_T UnitDelay_DSTATE_n[2];        /* '<Root>/Unit Delay' */
  real_T Initial_FirstOutputTime;      /* '<S37>/Initial' */
  real_T Memory_PreviousInput;         /* '<S37>/Memory' */
  real_T Memory_PreviousInput_f;       /* '<S59>/Memory' */
  real_T PrevY;                        /* '<S37>/Rate Limiter' */
  real_T LastMajorTime;                /* '<S37>/Rate Limiter' */
  real_T Memory_PreviousInput_d;       /* '<S47>/Memory' */
  real_T Memory_PreviousInput_o;       /* '<S48>/Memory' */
  real_T yData[125];                   /* '<Root>/MATLAB Function1' */
  real_T uData[50];                    /* '<Root>/MATLAB Function1' */
  struct {
    real_T modelTStart;
    real_T TUbufferArea[24576];
  } VariableTransportDelay_RWORK;      /* '<S59>/Variable Transport Delay' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[24576];
  } VariableTransportDelay_RWORK_k;    /* '<S47>/Variable Transport Delay' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[24576];
  } VariableTransportDelay_RWORK_c;    /* '<S48>/Variable Transport Delay' */

  struct {
    void *AS;
    void *BS;
    void *CS;
    void *DS;
    void *DX_COL;
    void *BD_COL;
    void *TMP1;
    void *TMP2;
    void *XTMP;
    void *SWITCH_STATUS;
    void *SWITCH_STATUS_INIT;
    void *SW_CHG;
    void *G_STATE;
    void *USWLAST;
    void *XKM12;
    void *XKP12;
    void *XLAST;
    void *ULAST;
    void *IDX_SW_CHG;
    void *Y_SWITCH;
    void *SWITCH_TYPES;
    void *IDX_OUT_SW;
    void *SWITCH_TOPO_SAVED_IDX;
    void *SWITCH_MAP;
  } StateSpace_PWORK;                  /* '<S139>/State-Space' */

  struct {
    void *TUbufferPtrs[3];
  } VariableTransportDelay_PWORK;      /* '<S59>/Variable Transport Delay' */

  struct {
    void *TUbufferPtrs[3];
  } VariableTransportDelay_PWORK_p;    /* '<S47>/Variable Transport Delay' */

  struct {
    void *TUbufferPtrs[3];
  } VariableTransportDelay_PWORK_f;    /* '<S48>/Variable Transport Delay' */

  int_T Integrator_IWORK;              /* '<S37>/Integrator' */
  int_T StateSpace_IWORK[11];          /* '<S139>/State-Space' */
  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } VariableTransportDelay_IWORK;      /* '<S59>/Variable Transport Delay' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } VariableTransportDelay_IWORK_l;    /* '<S47>/Variable Transport Delay' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } VariableTransportDelay_IWORK_g;    /* '<S48>/Variable Transport Delay' */

  boolean_T Subsystem1_MODE;           /* '<S134>/Subsystem1' */
  boolean_T Subsystempi2delay_MODE;    /* '<S134>/Subsystem - pi//2 delay' */
  boolean_T AutomaticGainControl_MODE; /* '<S37>/Automatic Gain Control' */
  DW_Subsystem1_Voltage_Control_T Subsystem1_ex;/* '<S127>/Subsystem1' */
  DW_Subsystempi2delay_Voltage__T Subsystempi2delay_n;/* '<S127>/Subsystem - pi//2 delay' */
  DW_Subsystem1_Voltage_Control_T Subsystem1_e;/* '<S121>/Subsystem1' */
  DW_Subsystempi2delay_Voltage__T Subsystempi2delay_b;/* '<S121>/Subsystem - pi//2 delay' */
  DW_Subsystem1_Voltage_Control_T Subsystem1;/* '<S60>/Subsystem1' */
  DW_Subsystempi2delay_Voltage__T Subsystempi2delay;/* '<S60>/Subsystem - pi//2 delay' */
  DW_Subsystem1_Voltage_Control_T Subsystem1_h;/* '<S49>/Subsystem1' */
  DW_Subsystempi2delay_Voltage__T Subsystempi2delay_d;/* '<S49>/Subsystem - pi//2 delay' */
} DW_Voltage_Control_T;

/* Continuous states (default storage) */
typedef struct {
  real_T Integrator_CSTATE;            /* '<S4>/Integrator' */
  real_T StateSpace1_CSTATE[3];        /* '<S4>/State-Space1' */
  real_T Integrator_CSTATE_f;          /* '<S37>/Integrator' */
  real_T Integrator_CSTATE_fe;         /* '<S39>/Integrator' */
  real_T VariableTransportDelay_CSTATE;/* '<S59>/Variable Transport Delay' */
  real_T integrator_CSTATE;            /* '<S59>/integrator' */
  real_T TransferFcn_CSTATE;           /* '<S39>/Transfer Fcn' */
  real_T Integrator_x1_CSTATE;         /* '<S55>/Integrator_x1' */
  real_T Integrator_x2_CSTATE;         /* '<S55>/Integrator_x2' */
  real_T VariableTransportDelay_CSTATE_n;/* '<S47>/Variable Transport Delay' */
  real_T integrator_CSTATE_m;          /* '<S47>/integrator' */
  real_T VariableTransportDelay_CSTAT_nz;/* '<S48>/Variable Transport Delay' */
  real_T integrator_CSTATE_l;          /* '<S48>/integrator' */
} X_Voltage_Control_T;

/* Periodic continuous state vector (global) */
typedef int_T PeriodicIndX_Voltage_Control_T[1];
typedef real_T PeriodicRngX_Voltage_Control_T[2];

/* State derivatives (default storage) */
typedef struct {
  real_T Integrator_CSTATE;            /* '<S4>/Integrator' */
  real_T StateSpace1_CSTATE[3];        /* '<S4>/State-Space1' */
  real_T Integrator_CSTATE_f;          /* '<S37>/Integrator' */
  real_T Integrator_CSTATE_fe;         /* '<S39>/Integrator' */
  real_T VariableTransportDelay_CSTATE;/* '<S59>/Variable Transport Delay' */
  real_T integrator_CSTATE;            /* '<S59>/integrator' */
  real_T TransferFcn_CSTATE;           /* '<S39>/Transfer Fcn' */
  real_T Integrator_x1_CSTATE;         /* '<S55>/Integrator_x1' */
  real_T Integrator_x2_CSTATE;         /* '<S55>/Integrator_x2' */
  real_T VariableTransportDelay_CSTATE_n;/* '<S47>/Variable Transport Delay' */
  real_T integrator_CSTATE_m;          /* '<S47>/integrator' */
  real_T VariableTransportDelay_CSTAT_nz;/* '<S48>/Variable Transport Delay' */
  real_T integrator_CSTATE_l;          /* '<S48>/integrator' */
} XDot_Voltage_Control_T;

/* State disabled  */
typedef struct {
  boolean_T Integrator_CSTATE;         /* '<S4>/Integrator' */
  boolean_T StateSpace1_CSTATE[3];     /* '<S4>/State-Space1' */
  boolean_T Integrator_CSTATE_f;       /* '<S37>/Integrator' */
  boolean_T Integrator_CSTATE_fe;      /* '<S39>/Integrator' */
  boolean_T VariableTransportDelay_CSTATE;/* '<S59>/Variable Transport Delay' */
  boolean_T integrator_CSTATE;         /* '<S59>/integrator' */
  boolean_T TransferFcn_CSTATE;        /* '<S39>/Transfer Fcn' */
  boolean_T Integrator_x1_CSTATE;      /* '<S55>/Integrator_x1' */
  boolean_T Integrator_x2_CSTATE;      /* '<S55>/Integrator_x2' */
  boolean_T VariableTransportDelay_CSTATE_n;/* '<S47>/Variable Transport Delay' */
  boolean_T integrator_CSTATE_m;       /* '<S47>/integrator' */
  boolean_T VariableTransportDelay_CSTAT_nz;/* '<S48>/Variable Transport Delay' */
  boolean_T integrator_CSTATE_l;       /* '<S48>/integrator' */
} XDis_Voltage_Control_T;

/* Zero-crossing (trigger) state */
typedef struct {
  ZCSigState Integrator_Reset_ZCE;     /* '<S37>/Integrator' */
} PrevZCX_Voltage_Control_T;

#ifndef ODE4_INTG
#define ODE4_INTG

/* ODE4 Integration Data */
typedef struct {
  real_T *y;                           /* output */
  real_T *f[4];                        /* derivatives */
} ODE4_IntgData;

#endif

/* External inputs (root inport signals with default storage) */
typedef struct {
  real_T action[2];                    /* '<Root>/action' */
  real_T rand_param[3];                /* '<Root>/rand_param' */
  real_T t_step;                       /* '<Root>/t_step' */
  real_T rand_input[4];                /* '<Root>/rand_input' */
} ExtU_Voltage_Control_T;

/* External outputs (root outports fed by signals with default storage) */
typedef struct {
  real_T obs[175];                     /* '<Root>/obs' */
  real_T y_meas[4];                    /* '<Root>/y_meas' */
  real_T y_true[4];                    /* '<Root>/y_true' */
} ExtY_Voltage_Control_T;

/* Parameters for system: '<S49>/Subsystem - pi//2 delay' */
struct P_Subsystempi2delay_Voltage_C_T_ {
  real_T dq_Y0[2];                     /* Expression: [0,0]
                                        * Referenced by: '<S53>/dq'
                                        */
};

/* Parameters for system: '<S49>/Subsystem1' */
struct P_Subsystem1_Voltage_Control_T_ {
  real_T dq_Y0[2];                     /* Expression: [0,0]
                                        * Referenced by: '<S54>/dq'
                                        */
};

/* Parameters (default storage) */
struct P_Voltage_Control_T_ {
  real_T AlphaBetaZerotodq0_Alignment;
                                 /* Mask Parameter: AlphaBetaZerotodq0_Alignment
                                  * Referenced by: '<S49>/Constant'
                                  */
  real_T dq0toAlphaBetaZero_Alignment;
                                 /* Mask Parameter: dq0toAlphaBetaZero_Alignment
                                  * Referenced by: '<S134>/Constant'
                                  */
  real_T AlphaBetaZerotodq0_Alignment_g;
                               /* Mask Parameter: AlphaBetaZerotodq0_Alignment_g
                                * Referenced by: '<S121>/Constant'
                                */
  real_T AlphaBetaZerotodq0_Alignment_i;
                               /* Mask Parameter: AlphaBetaZerotodq0_Alignment_i
                                * Referenced by: '<S127>/Constant'
                                */
  real_T AlphaBetaZerotodq0_Alignment_c;
                               /* Mask Parameter: AlphaBetaZerotodq0_Alignment_c
                                * Referenced by: '<S60>/Constant'
                                */
  real_T Continuous_Init;              /* Mask Parameter: Continuous_Init
                                        * Referenced by: '<S39>/Integrator'
                                        */
  real_T Continuous_Kd;                /* Mask Parameter: Continuous_Kd
                                        * Referenced by: '<S39>/Kp6'
                                        */
  real_T Continuous_Ki;                /* Mask Parameter: Continuous_Ki
                                        * Referenced by: '<S39>/Kp5'
                                        */
  real_T Continuous_Kp;                /* Mask Parameter: Continuous_Kp
                                        * Referenced by: '<S39>/Kp4'
                                        */
  real_T CompareToConstant1_const;   /* Mask Parameter: CompareToConstant1_const
                                      * Referenced by: '<S52>/Constant'
                                      */
  real_T CompareToConstant_const;     /* Mask Parameter: CompareToConstant_const
                                       * Referenced by: '<S51>/Constant'
                                       */
  real_T CompareToConstant_const_e; /* Mask Parameter: CompareToConstant_const_e
                                     * Referenced by: '<S135>/Constant'
                                     */
  real_T CompareToConstant1_const_m;
                                   /* Mask Parameter: CompareToConstant1_const_m
                                    * Referenced by: '<S136>/Constant'
                                    */
  real_T CompareToConstant_const_p; /* Mask Parameter: CompareToConstant_const_p
                                     * Referenced by: '<S123>/Constant'
                                     */
  real_T CompareToConstant1_const_l;
                                   /* Mask Parameter: CompareToConstant1_const_l
                                    * Referenced by: '<S124>/Constant'
                                    */
  real_T CompareToConstant_const_m; /* Mask Parameter: CompareToConstant_const_m
                                     * Referenced by: '<S129>/Constant'
                                     */
  real_T CompareToConstant1_const_h;
                                   /* Mask Parameter: CompareToConstant1_const_h
                                    * Referenced by: '<S130>/Constant'
                                    */
  real_T CompareToConstant_const_d; /* Mask Parameter: CompareToConstant_const_d
                                     * Referenced by: '<S62>/Constant'
                                     */
  real_T CompareToConstant1_const_i;
                                   /* Mask Parameter: CompareToConstant1_const_i
                                    * Referenced by: '<S63>/Constant'
                                    */
  real_T Gain_Y0;                      /* Expression: [1]
                                        * Referenced by: '<S38>/Gain'
                                        */
  real_T Gain3_Gain[9];
  /* Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
   * Referenced by: '<S50>/Gain3'
   */
  real_T Gain1_Gain;                   /* Expression: 2/3
                                        * Referenced by: '<S50>/Gain1'
                                        */
  real_T VariableTransportDelay_MaxDelay;/* Expression: 1/sps.Fmin+eps
                                          * Referenced by: '<S47>/Variable Transport Delay'
                                          */
  real_T VariableTransportDelay_InitOutp;/* Expression: 0
                                          * Referenced by: '<S47>/Variable Transport Delay'
                                          */
  real_T integrator_IC;                /* Expression: 0
                                        * Referenced by: '<S47>/integrator'
                                        */
  real_T Constant_Value;               /* Expression: 1/sps.Finit
                                        * Referenced by: '<S47>/Constant'
                                        */
  real_T Memory_InitialCondition;      /* Expression: sps.Vinit
                                        * Referenced by: '<S47>/Memory'
                                        */
  real_T VariableTransportDelay_MaxDel_e;/* Expression: 1/sps.Fmin+eps
                                          * Referenced by: '<S48>/Variable Transport Delay'
                                          */
  real_T VariableTransportDelay_InitOu_n;/* Expression: 0
                                          * Referenced by: '<S48>/Variable Transport Delay'
                                          */
  real_T integrator_IC_l;              /* Expression: 0
                                        * Referenced by: '<S48>/integrator'
                                        */
  real_T Constant_Value_n;             /* Expression: 1/sps.Finit
                                        * Referenced by: '<S48>/Constant'
                                        */
  real_T Memory_InitialCondition_i;    /* Expression: sps.Vinit
                                        * Referenced by: '<S48>/Memory'
                                        */
  real_T Saturation_UpperSat;          /* Expression: inf
                                        * Referenced by: '<S38>/Saturation'
                                        */
  real_T Saturation_LowerSat;          /* Expression: eps
                                        * Referenced by: '<S38>/Saturation'
                                        */
  real_T Toavoiddivisionbyzero_UpperSat;/* Expression: 1e6
                                         * Referenced by: '<S47>/To avoid division  by zero'
                                         */
  real_T Toavoiddivisionbyzero_LowerSat;/* Expression: eps
                                         * Referenced by: '<S47>/To avoid division  by zero'
                                         */
  real_T Toavoiddivisionbyzero_UpperSa_h;/* Expression: 1e6
                                          * Referenced by: '<S48>/To avoid division  by zero'
                                          */
  real_T Toavoiddivisionbyzero_LowerSa_h;/* Expression: eps
                                          * Referenced by: '<S48>/To avoid division  by zero'
                                          */
  real_T alpha_beta_Y0[2];             /* Expression: [0,0]
                                        * Referenced by: '<S137>/alpha_beta'
                                        */
  real_T alpha_beta_Y0_e[2];           /* Expression: [0,0]
                                        * Referenced by: '<S138>/alpha_beta'
                                        */
  real_T Integrator_IC;                /* Expression: 0
                                        * Referenced by: '<S4>/Integrator'
                                        */
  real_T Integrator_WrappedStateUpperVal;/* Expression: 2*pi
                                          * Referenced by: '<S4>/Integrator'
                                          */
  real_T Integrator_WrappedStateLowerVal;/* Expression: 0
                                          * Referenced by: '<S4>/Integrator'
                                          */
  real_T Constant2_Value;              /* Expression: 2*pi/3
                                        * Referenced by: '<S4>/Constant2'
                                        */
  real_T Gain2_Gain;                   /* Expression: 120*sqrt(2)
                                        * Referenced by: '<S4>/Gain2'
                                        */
  real_T Constant_Value_a;             /* Expression: 60
                                        * Referenced by: '<S4>/Constant'
                                        */
  real_T StateSpace1_A[5];             /* Computed Parameter: StateSpace1_A
                                        * Referenced by: '<S4>/State-Space1'
                                        */
  real_T StateSpace1_B;                /* Computed Parameter: StateSpace1_B
                                        * Referenced by: '<S4>/State-Space1'
                                        */
  real_T StateSpace1_C[3];             /* Computed Parameter: StateSpace1_C
                                        * Referenced by: '<S4>/State-Space1'
                                        */
  real_T StateSpace1_InitialCondition[3];/* Expression: [0;0;0]
                                          * Referenced by: '<S4>/State-Space1'
                                          */
  real_T Gain_Gain;                    /* Expression: 60
                                        * Referenced by: '<S4>/Gain'
                                        */
  real_T Gain4_Gain;                   /* Expression: 1/60
                                        * Referenced by: '<S4>/Gain4'
                                        */
  real_T Gain4_Gain_n;                 /* Expression: 392
                                        * Referenced by: '<Root>/Gain4'
                                        */
  real_T Constant2_Value_h;            /* Expression: 2*pi
                                        * Referenced by: '<S37>/Constant2'
                                        */
  real_T Initial_Value;                /* Expression: sps.Phase_Init*pi/180
                                        * Referenced by: '<S37>/Initial'
                                        */
  real_T Integrator_UpperSat;          /* Expression: inf
                                        * Referenced by: '<S37>/Integrator'
                                        */
  real_T Integrator_LowerSat;          /* Expression: -inf
                                        * Referenced by: '<S37>/Integrator'
                                        */
  real_T Gain3_Gain_l[9];
          /* Expression: [ 1   0   1; -1/2  sqrt(3)/2   1; -1/2  -sqrt(3)/2  1 ]
           * Referenced by: '<S133>/Gain3'
           */
  real_T UnitDelay_InitialCondition;   /* Expression: 0
                                        * Referenced by: '<S67>/Unit Delay'
                                        */
  real_T UnitDelay_InitialCondition_e; /* Expression: 0
                                        * Referenced by: '<S68>/Unit Delay'
                                        */
  real_T UnitDelay_InitialCondition_j; /* Expression: 0
                                        * Referenced by: '<S69>/Unit Delay'
                                        */
  real_T StateSpace_AS_param[36];      /* Expression: S.A
                                        * Referenced by: '<S139>/State-Space'
                                        */
  real_T StateSpace_BS_param[54];      /* Expression: S.B
                                        * Referenced by: '<S139>/State-Space'
                                        */
  real_T StateSpace_CS_param[144];     /* Expression: S.C
                                        * Referenced by: '<S139>/State-Space'
                                        */
  real_T StateSpace_DS_param[216];     /* Expression: S.D
                                        * Referenced by: '<S139>/State-Space'
                                        */
  real_T StateSpace_X0_param[6];       /* Expression: S.x0
                                        * Referenced by: '<S139>/State-Space'
                                        */
  real_T UnitDelay_InitialCondition_l; /* Expression: 0
                                        * Referenced by: '<Root>/Unit Delay'
                                        */
  real_T donotdeletethisgain_Gain;     /* Expression: 1
                                        * Referenced by: '<S91>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_n;   /* Expression: 1
                                        * Referenced by: '<S92>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_f;   /* Expression: 1
                                        * Referenced by: '<S93>/do not delete this gain'
                                        */
  real_T Kv_Gain;                      /* Expression: Ki
                                        * Referenced by: '<S10>/Kv'
                                        */
  real_T Gain3_Gain_o[9];
  /* Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
   * Referenced by: '<S122>/Gain3'
   */
  real_T Gain1_Gain_f;                 /* Expression: 2/3
                                        * Referenced by: '<S122>/Gain1'
                                        */
  real_T Gain_Gain_e;                  /* Expression: 1/392
                                        * Referenced by: '<Root>/Gain'
                                        */
  real_T Gain1_Gain_o;                 /* Expression: 1/392
                                        * Referenced by: '<Root>/Gain1'
                                        */
  real_T donotdeletethisgain_Gain_b;   /* Expression: 1
                                        * Referenced by: '<S94>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_a;   /* Expression: 1
                                        * Referenced by: '<S95>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_h;   /* Expression: 1
                                        * Referenced by: '<S96>/do not delete this gain'
                                        */
  real_T Kv1_Gain;                     /* Expression: Kv
                                        * Referenced by: '<S10>/Kv1'
                                        */
  real_T Gain3_Gain_d[9];
  /* Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
   * Referenced by: '<S128>/Gain3'
   */
  real_T Gain1_Gain_h;                 /* Expression: 2/3
                                        * Referenced by: '<S128>/Gain1'
                                        */
  real_T Gain2_Gain_d;                 /* Expression: 1/170
                                        * Referenced by: '<Root>/Gain2'
                                        */
  real_T Gain3_Gain_dd;                /* Expression: 1/170
                                        * Referenced by: '<Root>/Gain3'
                                        */
  real_T donotdeletethisgain_Gain_d;   /* Expression: 1
                                        * Referenced by: '<S71>/do not delete this gain'
                                        */
  real_T Constant_Value_h;             /* Expression: 208
                                        * Referenced by: '<Root>/Constant'
                                        */
  real_T Gain9_Gain;                   /* Expression: 100e3
                                        * Referenced by: '<Root>/Gain9'
                                        */
  real_T Saturation_UpperSat_e;        /* Expression: inf
                                        * Referenced by: '<S67>/Saturation'
                                        */
  real_T Saturation_LowerSat_c;        /* Expression: 1e-3
                                        * Referenced by: '<S67>/Saturation'
                                        */
  real_T donotdeletethisgain_Gain_nm;  /* Expression: 1
                                        * Referenced by: '<S74>/do not delete this gain'
                                        */
  real_T Saturation_UpperSat_p;        /* Expression: inf
                                        * Referenced by: '<S68>/Saturation'
                                        */
  real_T Saturation_LowerSat_e;        /* Expression: 1e-3
                                        * Referenced by: '<S68>/Saturation'
                                        */
  real_T donotdeletethisgain_Gain_ft;  /* Expression: 1
                                        * Referenced by: '<S77>/do not delete this gain'
                                        */
  real_T Saturation_UpperSat_f;        /* Expression: inf
                                        * Referenced by: '<S69>/Saturation'
                                        */
  real_T Saturation_LowerSat_a;        /* Expression: 1e-3
                                        * Referenced by: '<S69>/Saturation'
                                        */
  real_T HitCrossing_Offset;           /* Expression: pi
                                        * Referenced by: '<S37>/Hit  Crossing'
                                        */
  real_T Memory_InitialCondition_o;    /* Expression: sps.Finit
                                        * Referenced by: '<S37>/Memory'
                                        */
  real_T Constant1_Value;              /* Expression: sps.AGC
                                        * Referenced by: '<S37>/Constant1'
                                        */
  real_T Integrator_UpperSat_a;        /* Expression: Par_Limits(1)
                                        * Referenced by: '<S39>/Integrator'
                                        */
  real_T Integrator_LowerSat_n;        /* Expression: Par_Limits(2)
                                        * Referenced by: '<S39>/Integrator'
                                        */
  real_T VariableTransportDelay_MaxDel_n;/* Expression: 1/sps.Fmin+eps
                                          * Referenced by: '<S59>/Variable Transport Delay'
                                          */
  real_T VariableTransportDelay_InitOu_h;/* Expression: 0
                                          * Referenced by: '<S59>/Variable Transport Delay'
                                          */
  real_T integrator_IC_lo;             /* Expression: 0
                                        * Referenced by: '<S59>/integrator'
                                        */
  real_T Constant_Value_k;             /* Expression: 1/sps.Finit
                                        * Referenced by: '<S59>/Constant'
                                        */
  real_T Memory_InitialCondition_a;    /* Expression: sps.Vinit
                                        * Referenced by: '<S59>/Memory'
                                        */
  real_T TransferFcn_A;                /* Computed Parameter: TransferFcn_A
                                        * Referenced by: '<S39>/Transfer Fcn'
                                        */
  real_T TransferFcn_C;                /* Computed Parameter: TransferFcn_C
                                        * Referenced by: '<S39>/Transfer Fcn'
                                        */
  real_T TransferFcn_D;                /* Computed Parameter: TransferFcn_D
                                        * Referenced by: '<S39>/Transfer Fcn'
                                        */
  real_T Saturation2_UpperSat;         /* Expression: Par_Limits(1)
                                        * Referenced by: '<S39>/Saturation2'
                                        */
  real_T Saturation2_LowerSat;         /* Expression: Par_Limits(2)
                                        * Referenced by: '<S39>/Saturation2'
                                        */
  real_T Gain10_Gain;                  /* Expression: 1/2/pi
                                        * Referenced by: '<S37>/Gain10'
                                        */
  real_T RateLimiter_RisingLim;        /* Expression: sps.MaxRateChangeFreq
                                        * Referenced by: '<S37>/Rate Limiter'
                                        */
  real_T RateLimiter_FallingLim;       /* Expression: -sps.MaxRateChangeFreq
                                        * Referenced by: '<S37>/Rate Limiter'
                                        */
  real_T Integrator_x1_IC;             /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S55>/Integrator_x1'
                                        */
  real_T A11_Gain;                     /* Expression: sps.A11
                                        * Referenced by: '<S56>/A11'
                                        */
  real_T Integrator_x2_IC;             /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S55>/Integrator_x2'
                                        */
  real_T A12_Gain;                     /* Expression: sps.A12
                                        * Referenced by: '<S56>/A12'
                                        */
  real_T A21_Gain;                     /* Expression: sps.A21
                                        * Referenced by: '<S56>/A21'
                                        */
  real_T A22_Gain;                     /* Expression: sps.A22
                                        * Referenced by: '<S56>/A22'
                                        */
  real_T B11_Gain;                     /* Expression: sps.B11
                                        * Referenced by: '<S57>/B11'
                                        */
  real_T B21_Gain;                     /* Expression: sps.B21
                                        * Referenced by: '<S57>/B21'
                                        */
  real_T C11_Gain;                     /* Expression: sps.C11
                                        * Referenced by: '<S58>/C11'
                                        */
  real_T C12_Gain;                     /* Expression: sps.C12
                                        * Referenced by: '<S58>/C12'
                                        */
  real_T Du_Gain;                      /* Expression: sps.D
                                        * Referenced by: '<S55>/D*u'
                                        */
  real_T Toavoiddivisionbyzero_UpperSa_k;/* Expression: 1e6
                                          * Referenced by: '<S59>/To avoid division  by zero'
                                          */
  real_T Toavoiddivisionbyzero_LowerSa_c;/* Expression: eps
                                          * Referenced by: '<S59>/To avoid division  by zero'
                                          */
  real_T Gain3_Gain_m[9];
  /* Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
   * Referenced by: '<S61>/Gain3'
   */
  real_T Gain1_Gain_k;                 /* Expression: 2/3
                                        * Referenced by: '<S61>/Gain1'
                                        */
  real_T Gain1_Gain_i;                 /* Expression: 2*pi
                                        * Referenced by: '<S4>/Gain1'
                                        */
  real_T donotdeletethisgain_Gain_dg;  /* Expression: 1
                                        * Referenced by: '<S28>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_m;   /* Expression: 1
                                        * Referenced by: '<S29>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_j;   /* Expression: 1
                                        * Referenced by: '<S30>/do not delete this gain'
                                        */
  real_T Kv1_Gain_d;                   /* Expression: Kv
                                        * Referenced by: '<S21>/Kv1'
                                        */
  real_T donotdeletethisgain_Gain_p;   /* Expression: 1
                                        * Referenced by: '<S25>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_nw;  /* Expression: 1
                                        * Referenced by: '<S26>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_bf;  /* Expression: 1
                                        * Referenced by: '<S27>/do not delete this gain'
                                        */
  real_T Kv_Gain_j;                    /* Expression: Ki
                                        * Referenced by: '<S21>/Kv'
                                        */
  real_T Gain3_Gain_lb;                /* Expression: 1/100e3
                                        * Referenced by: '<S4>/Gain3'
                                        */
  P_Subsystem1_Voltage_Control_T Subsystem1_ex;/* '<S127>/Subsystem1' */
  P_Subsystempi2delay_Voltage_C_T Subsystempi2delay_n;/* '<S127>/Subsystem - pi//2 delay' */
  P_Subsystem1_Voltage_Control_T Subsystem1_e;/* '<S121>/Subsystem1' */
  P_Subsystempi2delay_Voltage_C_T Subsystempi2delay_b;/* '<S121>/Subsystem - pi//2 delay' */
  P_Subsystem1_Voltage_Control_T Subsystem1;/* '<S60>/Subsystem1' */
  P_Subsystempi2delay_Voltage_C_T Subsystempi2delay;/* '<S60>/Subsystem - pi//2 delay' */
  P_Subsystem1_Voltage_Control_T Subsystem1_h;/* '<S49>/Subsystem1' */
  P_Subsystempi2delay_Voltage_C_T Subsystempi2delay_d;/* '<S49>/Subsystem - pi//2 delay' */
};

/* Real-time Model Data Structure */
struct tag_RTM_Voltage_Control_T {
  const char_T *errorStatus;
  RTWSolverInfo solverInfo;
  X_Voltage_Control_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  boolean_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeY[15];
  real_T odeF[4][15];
  ODE4_IntgData intgData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numSampTimes;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    boolean_T firstInitCondFlag;
    struct {
      uint8_T TID[3];
    } TaskCounters;

    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[3];
  } Timing;
};

/* Block parameters (default storage) */
extern P_Voltage_Control_T Voltage_Control_P;

/* Block signals (default storage) */
extern B_Voltage_Control_T Voltage_Control_B;

/* Continuous states (default storage) */
extern X_Voltage_Control_T Voltage_Control_X;

/* Block states (default storage) */
extern DW_Voltage_Control_T Voltage_Control_DW;

/* Zero-crossing (trigger) state */
extern PrevZCX_Voltage_Control_T Voltage_Control_PrevZCX;

/* External inputs (root inport signals with default storage) */
extern ExtU_Voltage_Control_T Voltage_Control_U;

/* External outputs (root outports fed by signals with default storage) */
extern ExtY_Voltage_Control_T Voltage_Control_Y;

/* Model entry point functions */
extern void Voltage_Control_initialize(void);
extern void Voltage_Control_output(void);
extern void Voltage_Control_update(void);
extern void Voltage_Control_terminate(void);

/* Real-time Model object */
extern RT_MODEL_Voltage_Control_T *const Voltage_Control_M;

/*-
 * These blocks were eliminated from the model due to optimizations:
 *
 * Block '<S20>/Product3' : Unused code path elimination
 * Block '<S20>/Sum' : Unused code path elimination
 * Block '<S20>/Sum1' : Unused code path elimination
 * Block '<S20>/Sum5' : Unused code path elimination
 * Block '<S20>/Sum6' : Unused code path elimination
 * Block '<S20>/pu->V' : Unused code path elimination
 * Block '<S4>/Scope1' : Unused code path elimination
 * Block '<S4>/Scope2' : Unused code path elimination
 * Block '<S4>/Scope3' : Unused code path elimination
 * Block '<S4>/Scope4' : Unused code path elimination
 * Block '<S43>/Rad->Deg.' : Unused code path elimination
 * Block '<Root>/Scope4' : Unused code path elimination
 * Block '<S9>/Kv' : Unused code path elimination
 * Block '<S9>/Kv1' : Unused code path elimination
 * Block '<S82>/do not delete this gain' : Unused code path elimination
 * Block '<S83>/do not delete this gain' : Unused code path elimination
 * Block '<S84>/do not delete this gain' : Unused code path elimination
 * Block '<S11>/Kv' : Unused code path elimination
 * Block '<S11>/Kv1' : Unused code path elimination
 * Block '<S106>/do not delete this gain' : Unused code path elimination
 * Block '<S107>/do not delete this gain' : Unused code path elimination
 * Block '<S108>/do not delete this gain' : Unused code path elimination
 * Block '<S12>/Kv' : Unused code path elimination
 * Block '<S12>/Kv1' : Unused code path elimination
 * Block '<S115>/do not delete this gain' : Unused code path elimination
 * Block '<S116>/do not delete this gain' : Unused code path elimination
 * Block '<S117>/do not delete this gain' : Unused code path elimination
 */

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'Voltage_Control'
 * '<S1>'   : 'Voltage_Control/Controlled Current Source1'
 * '<S2>'   : 'Voltage_Control/Controlled Current Source2'
 * '<S3>'   : 'Voltage_Control/Controlled Current Source6'
 * '<S4>'   : 'Voltage_Control/Generator3'
 * '<S5>'   : 'Voltage_Control/MATLAB Function'
 * '<S6>'   : 'Voltage_Control/MATLAB Function1'
 * '<S7>'   : 'Voltage_Control/PLL (3ph)'
 * '<S8>'   : 'Voltage_Control/Subsystem1'
 * '<S9>'   : 'Voltage_Control/Three-Phase V-I Measurement1'
 * '<S10>'  : 'Voltage_Control/Three-Phase V-I Measurement2'
 * '<S11>'  : 'Voltage_Control/Three-Phase V-I Measurement6'
 * '<S12>'  : 'Voltage_Control/Three-Phase V-I Measurement7'
 * '<S13>'  : 'Voltage_Control/abc to dq0'
 * '<S14>'  : 'Voltage_Control/abc to dq1'
 * '<S15>'  : 'Voltage_Control/dq0 to abc'
 * '<S16>'  : 'Voltage_Control/powergui'
 * '<S17>'  : 'Voltage_Control/Generator3/Controlled Voltage Source'
 * '<S18>'  : 'Voltage_Control/Generator3/Controlled Voltage Source1'
 * '<S19>'  : 'Voltage_Control/Generator3/Controlled Voltage Source2'
 * '<S20>'  : 'Voltage_Control/Generator3/Power (3ph, Instantaneous)'
 * '<S21>'  : 'Voltage_Control/Generator3/Three-Phase V-I Measurement3'
 * '<S22>'  : 'Voltage_Control/Generator3/Three-Phase V-I Measurement3/Mode I'
 * '<S23>'  : 'Voltage_Control/Generator3/Three-Phase V-I Measurement3/Mode V'
 * '<S24>'  : 'Voltage_Control/Generator3/Three-Phase V-I Measurement3/Model'
 * '<S25>'  : 'Voltage_Control/Generator3/Three-Phase V-I Measurement3/Model/I A:'
 * '<S26>'  : 'Voltage_Control/Generator3/Three-Phase V-I Measurement3/Model/I B:'
 * '<S27>'  : 'Voltage_Control/Generator3/Three-Phase V-I Measurement3/Model/I C:'
 * '<S28>'  : 'Voltage_Control/Generator3/Three-Phase V-I Measurement3/Model/U A:'
 * '<S29>'  : 'Voltage_Control/Generator3/Three-Phase V-I Measurement3/Model/U B:'
 * '<S30>'  : 'Voltage_Control/Generator3/Three-Phase V-I Measurement3/Model/U C:'
 * '<S31>'  : 'Voltage_Control/Generator3/Three-Phase V-I Measurement3/Model/I A:/Model'
 * '<S32>'  : 'Voltage_Control/Generator3/Three-Phase V-I Measurement3/Model/I B:/Model'
 * '<S33>'  : 'Voltage_Control/Generator3/Three-Phase V-I Measurement3/Model/I C:/Model'
 * '<S34>'  : 'Voltage_Control/Generator3/Three-Phase V-I Measurement3/Model/U A:/Model'
 * '<S35>'  : 'Voltage_Control/Generator3/Three-Phase V-I Measurement3/Model/U B:/Model'
 * '<S36>'  : 'Voltage_Control/Generator3/Three-Phase V-I Measurement3/Model/U C:/Model'
 * '<S37>'  : 'Voltage_Control/PLL (3ph)/Model'
 * '<S38>'  : 'Voltage_Control/PLL (3ph)/Model/Automatic Gain Control'
 * '<S39>'  : 'Voltage_Control/PLL (3ph)/Model/Continuous'
 * '<S40>'  : 'Voltage_Control/PLL (3ph)/Model/Second-Order Filter'
 * '<S41>'  : 'Voltage_Control/PLL (3ph)/Model/Variable Frequency Mean value'
 * '<S42>'  : 'Voltage_Control/PLL (3ph)/Model/abc to dq0'
 * '<S43>'  : 'Voltage_Control/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)'
 * '<S44>'  : 'Voltage_Control/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/Mean (Variable Frequency)1'
 * '<S45>'  : 'Voltage_Control/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/Mean (Variable Frequency)2'
 * '<S46>'  : 'Voltage_Control/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/abc ->dq0'
 * '<S47>'  : 'Voltage_Control/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/Mean (Variable Frequency)1/Model'
 * '<S48>'  : 'Voltage_Control/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/Mean (Variable Frequency)2/Model'
 * '<S49>'  : 'Voltage_Control/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/abc ->dq0/Alpha-Beta-Zero to dq0'
 * '<S50>'  : 'Voltage_Control/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/abc ->dq0/abc to Alpha-Beta-Zero'
 * '<S51>'  : 'Voltage_Control/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/abc ->dq0/Alpha-Beta-Zero to dq0/Compare To Constant'
 * '<S52>'  : 'Voltage_Control/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/abc ->dq0/Alpha-Beta-Zero to dq0/Compare To Constant1'
 * '<S53>'  : 'Voltage_Control/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/abc ->dq0/Alpha-Beta-Zero to dq0/Subsystem - pi//2 delay'
 * '<S54>'  : 'Voltage_Control/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/abc ->dq0/Alpha-Beta-Zero to dq0/Subsystem1'
 * '<S55>'  : 'Voltage_Control/PLL (3ph)/Model/Second-Order Filter/Model'
 * '<S56>'  : 'Voltage_Control/PLL (3ph)/Model/Second-Order Filter/Model/A*x'
 * '<S57>'  : 'Voltage_Control/PLL (3ph)/Model/Second-Order Filter/Model/B*u'
 * '<S58>'  : 'Voltage_Control/PLL (3ph)/Model/Second-Order Filter/Model/C*x'
 * '<S59>'  : 'Voltage_Control/PLL (3ph)/Model/Variable Frequency Mean value/Model'
 * '<S60>'  : 'Voltage_Control/PLL (3ph)/Model/abc to dq0/Alpha-Beta-Zero to dq0'
 * '<S61>'  : 'Voltage_Control/PLL (3ph)/Model/abc to dq0/abc to Alpha-Beta-Zero'
 * '<S62>'  : 'Voltage_Control/PLL (3ph)/Model/abc to dq0/Alpha-Beta-Zero to dq0/Compare To Constant'
 * '<S63>'  : 'Voltage_Control/PLL (3ph)/Model/abc to dq0/Alpha-Beta-Zero to dq0/Compare To Constant1'
 * '<S64>'  : 'Voltage_Control/PLL (3ph)/Model/abc to dq0/Alpha-Beta-Zero to dq0/Subsystem - pi//2 delay'
 * '<S65>'  : 'Voltage_Control/PLL (3ph)/Model/abc to dq0/Alpha-Beta-Zero to dq0/Subsystem1'
 * '<S66>'  : 'Voltage_Control/Subsystem1/MATLAB Function'
 * '<S67>'  : 'Voltage_Control/Subsystem1/Subsystem'
 * '<S68>'  : 'Voltage_Control/Subsystem1/Subsystem1'
 * '<S69>'  : 'Voltage_Control/Subsystem1/Subsystem2'
 * '<S70>'  : 'Voltage_Control/Subsystem1/Subsystem/Controlled Current Source'
 * '<S71>'  : 'Voltage_Control/Subsystem1/Subsystem/Voltage Measurement'
 * '<S72>'  : 'Voltage_Control/Subsystem1/Subsystem/Voltage Measurement/Model'
 * '<S73>'  : 'Voltage_Control/Subsystem1/Subsystem1/Controlled Current Source'
 * '<S74>'  : 'Voltage_Control/Subsystem1/Subsystem1/Voltage Measurement'
 * '<S75>'  : 'Voltage_Control/Subsystem1/Subsystem1/Voltage Measurement/Model'
 * '<S76>'  : 'Voltage_Control/Subsystem1/Subsystem2/Controlled Current Source'
 * '<S77>'  : 'Voltage_Control/Subsystem1/Subsystem2/Voltage Measurement'
 * '<S78>'  : 'Voltage_Control/Subsystem1/Subsystem2/Voltage Measurement/Model'
 * '<S79>'  : 'Voltage_Control/Three-Phase V-I Measurement1/Mode I'
 * '<S80>'  : 'Voltage_Control/Three-Phase V-I Measurement1/Mode V'
 * '<S81>'  : 'Voltage_Control/Three-Phase V-I Measurement1/Model'
 * '<S82>'  : 'Voltage_Control/Three-Phase V-I Measurement1/Model/I A:'
 * '<S83>'  : 'Voltage_Control/Three-Phase V-I Measurement1/Model/I B:'
 * '<S84>'  : 'Voltage_Control/Three-Phase V-I Measurement1/Model/I C:'
 * '<S85>'  : 'Voltage_Control/Three-Phase V-I Measurement1/Model/I A:/Model'
 * '<S86>'  : 'Voltage_Control/Three-Phase V-I Measurement1/Model/I B:/Model'
 * '<S87>'  : 'Voltage_Control/Three-Phase V-I Measurement1/Model/I C:/Model'
 * '<S88>'  : 'Voltage_Control/Three-Phase V-I Measurement2/Mode I'
 * '<S89>'  : 'Voltage_Control/Three-Phase V-I Measurement2/Mode V'
 * '<S90>'  : 'Voltage_Control/Three-Phase V-I Measurement2/Model'
 * '<S91>'  : 'Voltage_Control/Three-Phase V-I Measurement2/Model/I A:'
 * '<S92>'  : 'Voltage_Control/Three-Phase V-I Measurement2/Model/I B:'
 * '<S93>'  : 'Voltage_Control/Three-Phase V-I Measurement2/Model/I C:'
 * '<S94>'  : 'Voltage_Control/Three-Phase V-I Measurement2/Model/U A:'
 * '<S95>'  : 'Voltage_Control/Three-Phase V-I Measurement2/Model/U B:'
 * '<S96>'  : 'Voltage_Control/Three-Phase V-I Measurement2/Model/U C:'
 * '<S97>'  : 'Voltage_Control/Three-Phase V-I Measurement2/Model/I A:/Model'
 * '<S98>'  : 'Voltage_Control/Three-Phase V-I Measurement2/Model/I B:/Model'
 * '<S99>'  : 'Voltage_Control/Three-Phase V-I Measurement2/Model/I C:/Model'
 * '<S100>' : 'Voltage_Control/Three-Phase V-I Measurement2/Model/U A:/Model'
 * '<S101>' : 'Voltage_Control/Three-Phase V-I Measurement2/Model/U B:/Model'
 * '<S102>' : 'Voltage_Control/Three-Phase V-I Measurement2/Model/U C:/Model'
 * '<S103>' : 'Voltage_Control/Three-Phase V-I Measurement6/Mode I'
 * '<S104>' : 'Voltage_Control/Three-Phase V-I Measurement6/Mode V'
 * '<S105>' : 'Voltage_Control/Three-Phase V-I Measurement6/Model'
 * '<S106>' : 'Voltage_Control/Three-Phase V-I Measurement6/Model/I A:'
 * '<S107>' : 'Voltage_Control/Three-Phase V-I Measurement6/Model/I B:'
 * '<S108>' : 'Voltage_Control/Three-Phase V-I Measurement6/Model/I C:'
 * '<S109>' : 'Voltage_Control/Three-Phase V-I Measurement6/Model/I A:/Model'
 * '<S110>' : 'Voltage_Control/Three-Phase V-I Measurement6/Model/I B:/Model'
 * '<S111>' : 'Voltage_Control/Three-Phase V-I Measurement6/Model/I C:/Model'
 * '<S112>' : 'Voltage_Control/Three-Phase V-I Measurement7/Mode I'
 * '<S113>' : 'Voltage_Control/Three-Phase V-I Measurement7/Mode V'
 * '<S114>' : 'Voltage_Control/Three-Phase V-I Measurement7/Model'
 * '<S115>' : 'Voltage_Control/Three-Phase V-I Measurement7/Model/U A:'
 * '<S116>' : 'Voltage_Control/Three-Phase V-I Measurement7/Model/U B:'
 * '<S117>' : 'Voltage_Control/Three-Phase V-I Measurement7/Model/U C:'
 * '<S118>' : 'Voltage_Control/Three-Phase V-I Measurement7/Model/U A:/Model'
 * '<S119>' : 'Voltage_Control/Three-Phase V-I Measurement7/Model/U B:/Model'
 * '<S120>' : 'Voltage_Control/Three-Phase V-I Measurement7/Model/U C:/Model'
 * '<S121>' : 'Voltage_Control/abc to dq0/Alpha-Beta-Zero to dq0'
 * '<S122>' : 'Voltage_Control/abc to dq0/abc to Alpha-Beta-Zero'
 * '<S123>' : 'Voltage_Control/abc to dq0/Alpha-Beta-Zero to dq0/Compare To Constant'
 * '<S124>' : 'Voltage_Control/abc to dq0/Alpha-Beta-Zero to dq0/Compare To Constant1'
 * '<S125>' : 'Voltage_Control/abc to dq0/Alpha-Beta-Zero to dq0/Subsystem - pi//2 delay'
 * '<S126>' : 'Voltage_Control/abc to dq0/Alpha-Beta-Zero to dq0/Subsystem1'
 * '<S127>' : 'Voltage_Control/abc to dq1/Alpha-Beta-Zero to dq0'
 * '<S128>' : 'Voltage_Control/abc to dq1/abc to Alpha-Beta-Zero'
 * '<S129>' : 'Voltage_Control/abc to dq1/Alpha-Beta-Zero to dq0/Compare To Constant'
 * '<S130>' : 'Voltage_Control/abc to dq1/Alpha-Beta-Zero to dq0/Compare To Constant1'
 * '<S131>' : 'Voltage_Control/abc to dq1/Alpha-Beta-Zero to dq0/Subsystem - pi//2 delay'
 * '<S132>' : 'Voltage_Control/abc to dq1/Alpha-Beta-Zero to dq0/Subsystem1'
 * '<S133>' : 'Voltage_Control/dq0 to abc/Alpha-Beta-Zero to abc'
 * '<S134>' : 'Voltage_Control/dq0 to abc/dq0 to Alpha-Beta-Zero'
 * '<S135>' : 'Voltage_Control/dq0 to abc/dq0 to Alpha-Beta-Zero/Compare To Constant'
 * '<S136>' : 'Voltage_Control/dq0 to abc/dq0 to Alpha-Beta-Zero/Compare To Constant1'
 * '<S137>' : 'Voltage_Control/dq0 to abc/dq0 to Alpha-Beta-Zero/Subsystem - pi//2 delay'
 * '<S138>' : 'Voltage_Control/dq0 to abc/dq0 to Alpha-Beta-Zero/Subsystem1'
 * '<S139>' : 'Voltage_Control/powergui/EquivalentModel1'
 * '<S140>' : 'Voltage_Control/powergui/EquivalentModel1/Sources'
 * '<S141>' : 'Voltage_Control/powergui/EquivalentModel1/Yout'
 */
#endif                                 /* RTW_HEADER_Voltage_Control_h_ */
