/*
 * Voltage_Control_data.c
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

#include "Voltage_Control.h"
#include "Voltage_Control_private.h"

/* Block parameters (default storage) */
P_Voltage_Control_T Voltage_Control_P = {
  /* Mask Parameter: AlphaBetaZerotodq0_Alignment
   * Referenced by: '<S49>/Constant'
   */
  2.0,

  /* Mask Parameter: dq0toAlphaBetaZero_Alignment
   * Referenced by: '<S134>/Constant'
   */
  2.0,

  /* Mask Parameter: AlphaBetaZerotodq0_Alignment_g
   * Referenced by: '<S121>/Constant'
   */
  2.0,

  /* Mask Parameter: AlphaBetaZerotodq0_Alignment_i
   * Referenced by: '<S127>/Constant'
   */
  2.0,

  /* Mask Parameter: AlphaBetaZerotodq0_Alignment_c
   * Referenced by: '<S60>/Constant'
   */
  2.0,

  /* Mask Parameter: Continuous_Init
   * Referenced by: '<S39>/Integrator'
   */
  376.99111843077515,

  /* Mask Parameter: Continuous_Kd
   * Referenced by: '<S39>/Kp6'
   */
  1.0,

  /* Mask Parameter: Continuous_Ki
   * Referenced by: '<S39>/Kp5'
   */
  3200.0,

  /* Mask Parameter: Continuous_Kp
   * Referenced by: '<S39>/Kp4'
   */
  180.0,

  /* Mask Parameter: CompareToConstant1_const
   * Referenced by: '<S52>/Constant'
   */
  2.0,

  /* Mask Parameter: CompareToConstant_const
   * Referenced by: '<S51>/Constant'
   */
  1.0,

  /* Mask Parameter: CompareToConstant_const_e
   * Referenced by: '<S135>/Constant'
   */
  1.0,

  /* Mask Parameter: CompareToConstant1_const_m
   * Referenced by: '<S136>/Constant'
   */
  2.0,

  /* Mask Parameter: CompareToConstant_const_p
   * Referenced by: '<S123>/Constant'
   */
  1.0,

  /* Mask Parameter: CompareToConstant1_const_l
   * Referenced by: '<S124>/Constant'
   */
  2.0,

  /* Mask Parameter: CompareToConstant_const_m
   * Referenced by: '<S129>/Constant'
   */
  1.0,

  /* Mask Parameter: CompareToConstant1_const_h
   * Referenced by: '<S130>/Constant'
   */
  2.0,

  /* Mask Parameter: CompareToConstant_const_d
   * Referenced by: '<S62>/Constant'
   */
  1.0,

  /* Mask Parameter: CompareToConstant1_const_i
   * Referenced by: '<S63>/Constant'
   */
  2.0,

  /* Expression: [1]
   * Referenced by: '<S38>/Gain'
   */
  1.0,

  /* Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
   * Referenced by: '<S50>/Gain3'
   */
  { 1.0, 0.0, 0.5, -0.5, 0.8660254037844386, 0.5, -0.5, -0.8660254037844386, 0.5
  },

  /* Expression: 2/3
   * Referenced by: '<S50>/Gain1'
   */
  0.66666666666666663,

  /* Expression: 1/sps.Fmin+eps
   * Referenced by: '<S47>/Variable Transport Delay'
   */
  0.022222222222222445,

  /* Expression: 0
   * Referenced by: '<S47>/Variable Transport Delay'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S47>/integrator'
   */
  0.0,

  /* Expression: 1/sps.Finit
   * Referenced by: '<S47>/Constant'
   */
  0.016666666666666666,

  /* Expression: sps.Vinit
   * Referenced by: '<S47>/Memory'
   */
  1.0,

  /* Expression: 1/sps.Fmin+eps
   * Referenced by: '<S48>/Variable Transport Delay'
   */
  0.022222222222222445,

  /* Expression: 0
   * Referenced by: '<S48>/Variable Transport Delay'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S48>/integrator'
   */
  0.0,

  /* Expression: 1/sps.Finit
   * Referenced by: '<S48>/Constant'
   */
  0.016666666666666666,

  /* Expression: sps.Vinit
   * Referenced by: '<S48>/Memory'
   */
  0.0,

  /* Expression: inf
   * Referenced by: '<S38>/Saturation'
   */
  0.0,

  /* Expression: eps
   * Referenced by: '<S38>/Saturation'
   */
  2.2204460492503131E-16,

  /* Expression: 1e6
   * Referenced by: '<S47>/To avoid division  by zero'
   */
  1.0E+6,

  /* Expression: eps
   * Referenced by: '<S47>/To avoid division  by zero'
   */
  2.2204460492503131E-16,

  /* Expression: 1e6
   * Referenced by: '<S48>/To avoid division  by zero'
   */
  1.0E+6,

  /* Expression: eps
   * Referenced by: '<S48>/To avoid division  by zero'
   */
  2.2204460492503131E-16,

  /* Expression: [0,0]
   * Referenced by: '<S137>/alpha_beta'
   */
  { 0.0, 0.0 },

  /* Expression: [0,0]
   * Referenced by: '<S138>/alpha_beta'
   */
  { 0.0, 0.0 },

  /* Expression: 0
   * Referenced by: '<S4>/Integrator'
   */
  0.0,

  /* Expression: 2*pi
   * Referenced by: '<S4>/Integrator'
   */
  6.2831853071795862,

  /* Expression: 0
   * Referenced by: '<S4>/Integrator'
   */
  0.0,

  /* Expression: 2*pi/3
   * Referenced by: '<S4>/Constant2'
   */
  2.0943951023931953,

  /* Expression: 120*sqrt(2)
   * Referenced by: '<S4>/Gain2'
   */
  169.70562748477141,

  /* Expression: 60
   * Referenced by: '<S4>/Constant'
   */
  60.0,

  /* Computed Parameter: StateSpace1_A
   * Referenced by: '<S4>/State-Space1'
   */
  { 1.0, 1.0, -2.5, -26.874999999999996, -5.375 },

  /* Computed Parameter: StateSpace1_B
   * Referenced by: '<S4>/State-Space1'
   */
  -1.25,

  /* Computed Parameter: StateSpace1_C
   * Referenced by: '<S4>/State-Space1'
   */
  { 1.0, 1.0, 1.0 },

  /* Expression: [0;0;0]
   * Referenced by: '<S4>/State-Space1'
   */
  { 0.0, 0.0, 0.0 },

  /* Expression: 60
   * Referenced by: '<S4>/Gain'
   */
  60.0,

  /* Expression: 1/60
   * Referenced by: '<S4>/Gain4'
   */
  0.016666666666666666,

  /* Expression: 392
   * Referenced by: '<Root>/Gain4'
   */
  392.0,

  /* Expression: 2*pi
   * Referenced by: '<S37>/Constant2'
   */
  6.2831853071795862,

  /* Expression: sps.Phase_Init*pi/180
   * Referenced by: '<S37>/Initial'
   */
  0.0,

  /* Expression: inf
   * Referenced by: '<S37>/Integrator'
   */
  0.0,

  /* Expression: -inf
   * Referenced by: '<S37>/Integrator'
   */
  0.0,

  /* Expression: [ 1   0   1; -1/2  sqrt(3)/2   1; -1/2  -sqrt(3)/2  1 ]
   * Referenced by: '<S133>/Gain3'
   */
  { 1.0, -0.5, -0.5, 0.0, 0.8660254037844386, -0.8660254037844386, 1.0, 1.0, 1.0
  },

  /* Expression: 0
   * Referenced by: '<S67>/Unit Delay'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S68>/Unit Delay'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S69>/Unit Delay'
   */
  0.0,

  /* Expression: S.A
   * Referenced by: '<S139>/State-Space'
   */
  { 0.99896934842207252, 0.0, 0.0, 0.0, 0.045348669428812896, 0.0, 0.0,
    0.99896934842207252, 0.0, 0.045348669428812896, 0.0, 0.0, 0.0, 0.0,
    0.99896934842207252, 0.0, 0.0, 0.045348669428812896, 0.0,
    -0.045348669428812896, 0.0, 0.99534145486776748, 0.0, 0.0,
    -0.045348669428812896, 0.0, 0.0, 0.0, 0.99534145486776748, 0.0, 0.0, 0.0,
    -0.045348669428812896, 0.0, 0.0, 0.99534145486776748 },

  /* Expression: S.B
   * Referenced by: '<S139>/State-Space'
   */
  { 0.0, 103.06515779275658, 0.0, -4534.8669428812891, 0.0, 0.0,
    103.06515779275658, 0.0, 0.0, 0.0, -4534.8669428812891, 0.0, 0.0, 0.0,
    103.06515779275658, 0.0, 0.0, -4534.8669428812891, 0.0, 4543.11215550471,
    0.0, 103.06515779275658, 0.0, 0.0, 4543.11215550471, 0.0, 0.0, 0.0,
    103.06515779275658, 0.0, 0.0, 0.0, 4543.11215550471, 0.0, 0.0,
    103.06515779275658, 0.0, -4543.11215550471, 0.0, -103.06515779275658, 0.0,
    0.0, -4543.11215550471, 0.0, 0.0, 0.0, -103.06515779275658, 0.0, 0.0, 0.0,
    -4543.11215550471, 0.0, 0.0, -103.06515779275658 },

  /* Expression: S.C
   * Referenced by: '<S139>/State-Space'
   */
  { 0.0, 0.0, 0.0, 0.0, 9.9948467421103633E-6, 0.0, 0.0, 9.9948467421103633E-6,
    0.0, 0.0, 0.0, 0.0, 0.0, -2.267433471440645E-7, 0.0, 0.0, 0.0, 0.0, 0.0,
    2.267433471440645E-7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    9.9948467421103633E-6, 0.0, 0.0, 9.9948467421103633E-6, 0.0, 0.0, 0.0, 0.0,
    0.0, -2.267433471440645E-7, 0.0, 0.0, 0.0, 0.0, 0.0, 2.267433471440645E-7,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.9948467421103633E-6, 0.0,
    0.0, 9.9948467421103633E-6, 0.0, 0.0, 0.0, 0.0, 0.0, -2.267433471440645E-7,
    0.0, 0.0, 0.0, 0.0, 0.0, 2.267433471440645E-7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -2.267433471440645E-7, 0.0, 0.0, -2.267433471440645E-7, 0.0, 0.0, 0.0, 0.0,
    0.0, -9.9767072743388382E-6, 0.0, 0.0, 0.0, 0.0, 0.0, 9.9767072743388382E-6,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.267433471440645E-7, 0.0, 0.0,
    -2.267433471440645E-7, 0.0, 0.0, 0.0, 0.0, 0.0, -9.9767072743388382E-6, 0.0,
    0.0, 0.0, 0.0, 0.0, 9.9767072743388382E-6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -2.267433471440645E-7, 0.0, 0.0, -2.267433471440645E-7, 0.0, 0.0,
    0.0, 0.0, 0.0, -9.9767072743388382E-6, 0.0, 0.0, 0.0, 0.0, 0.0,
    9.9767072743388382E-6, 0.0, 0.0, 0.0 },

  /* Expression: S.D
   * Referenced by: '<S139>/State-Space'
   */
  { 1.0, 0.0, 0.0, 0.00051532578896378294, 0.0, 0.0, 0.00051532578896378294, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.022674334714406448, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.022674334714406448, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.00051532578896378294, 0.0, 0.0, 0.00051532578896378294, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.022674334714406448, 0.0, 0.0, 0.0, 0.0, 0.0, -0.022674334714406448,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.00051532578896378294, 0.0,
    0.0, 0.00051532578896378294, 0.0, 0.0, 1.0, 0.0, 0.0, 0.022674334714406448,
    0.0, 0.0, 0.0, 0.0, 0.0, -0.022674334714406448, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.022715560777523552, 0.0, 0.0, 0.022715560777523552, 0.0, 0.0, 0.0, 0.0,
    0.0, -0.00051532578896378294, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.00051532578896378294, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.022715560777523552, 0.0, 0.0, 0.022715560777523552, 0.0, 0.0, 0.0, 0.0,
    0.0, -0.00051532578896378294, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.00051532578896378294, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.022715560777523552, 0.0, 0.0, 0.022715560777523552, 0.0, 0.0, 0.0, 0.0,
    0.0, -0.00051532578896378294, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.00051532578896378294, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.022715560777523552,
    0.0, 0.0, -0.022715560777523552, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.00051532578896378294, 0.0, 0.0, 0.0, 0.0, 0.0, 0.99948467421103626, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.022715560777523552, 0.0, 0.0,
    -0.022715560777523552, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00051532578896378294, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.99948467421103626, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.022715560777523552, 0.0, 0.0, -0.022715560777523552, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.00051532578896378294, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.99948467421103626, 0.0, 0.0, 1.0 },

  /* Expression: S.x0
   * Referenced by: '<S139>/State-Space'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  /* Expression: 0
   * Referenced by: '<Root>/Unit Delay'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S91>/do not delete this gain'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S92>/do not delete this gain'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S93>/do not delete this gain'
   */
  1.0,

  /* Expression: Ki
   * Referenced by: '<S10>/Kv'
   */
  1.0,

  /* Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
   * Referenced by: '<S122>/Gain3'
   */
  { 1.0, 0.0, 0.5, -0.5, 0.8660254037844386, 0.5, -0.5, -0.8660254037844386, 0.5
  },

  /* Expression: 2/3
   * Referenced by: '<S122>/Gain1'
   */
  0.66666666666666663,

  /* Expression: 1/392
   * Referenced by: '<Root>/Gain'
   */
  0.0025510204081632651,

  /* Expression: 1/392
   * Referenced by: '<Root>/Gain1'
   */
  0.0025510204081632651,

  /* Expression: 1
   * Referenced by: '<S94>/do not delete this gain'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S95>/do not delete this gain'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S96>/do not delete this gain'
   */
  1.0,

  /* Expression: Kv
   * Referenced by: '<S10>/Kv1'
   */
  1.0,

  /* Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
   * Referenced by: '<S128>/Gain3'
   */
  { 1.0, 0.0, 0.5, -0.5, 0.8660254037844386, 0.5, -0.5, -0.8660254037844386, 0.5
  },

  /* Expression: 2/3
   * Referenced by: '<S128>/Gain1'
   */
  0.66666666666666663,

  /* Expression: 1/170
   * Referenced by: '<Root>/Gain2'
   */
  0.0058823529411764705,

  /* Expression: 1/170
   * Referenced by: '<Root>/Gain3'
   */
  0.0058823529411764705,

  /* Expression: 1
   * Referenced by: '<S71>/do not delete this gain'
   */
  1.0,

  /* Expression: 208
   * Referenced by: '<Root>/Constant'
   */
  208.0,

  /* Expression: 100e3
   * Referenced by: '<Root>/Gain9'
   */
  100000.0,

  /* Expression: inf
   * Referenced by: '<S67>/Saturation'
   */
  0.0,

  /* Expression: 1e-3
   * Referenced by: '<S67>/Saturation'
   */
  0.001,

  /* Expression: 1
   * Referenced by: '<S74>/do not delete this gain'
   */
  1.0,

  /* Expression: inf
   * Referenced by: '<S68>/Saturation'
   */
  0.0,

  /* Expression: 1e-3
   * Referenced by: '<S68>/Saturation'
   */
  0.001,

  /* Expression: 1
   * Referenced by: '<S77>/do not delete this gain'
   */
  1.0,

  /* Expression: inf
   * Referenced by: '<S69>/Saturation'
   */
  0.0,

  /* Expression: 1e-3
   * Referenced by: '<S69>/Saturation'
   */
  0.001,

  /* Expression: pi
   * Referenced by: '<S37>/Hit  Crossing'
   */
  3.1415926535897931,

  /* Expression: sps.Finit
   * Referenced by: '<S37>/Memory'
   */
  60.0,

  /* Expression: sps.AGC
   * Referenced by: '<S37>/Constant1'
   */
  1.0,

  /* Expression: Par_Limits(1)
   * Referenced by: '<S39>/Integrator'
   */
  0.0,

  /* Expression: Par_Limits(2)
   * Referenced by: '<S39>/Integrator'
   */
  0.0,

  /* Expression: 1/sps.Fmin+eps
   * Referenced by: '<S59>/Variable Transport Delay'
   */
  0.022222222222222445,

  /* Expression: 0
   * Referenced by: '<S59>/Variable Transport Delay'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S59>/integrator'
   */
  0.0,

  /* Expression: 1/sps.Finit
   * Referenced by: '<S59>/Constant'
   */
  0.016666666666666666,

  /* Expression: sps.Vinit
   * Referenced by: '<S59>/Memory'
   */
  0.0,

  /* Computed Parameter: TransferFcn_A
   * Referenced by: '<S39>/Transfer Fcn'
   */
  -10000.0,

  /* Computed Parameter: TransferFcn_C
   * Referenced by: '<S39>/Transfer Fcn'
   */
  -1.0E+8,

  /* Computed Parameter: TransferFcn_D
   * Referenced by: '<S39>/Transfer Fcn'
   */
  10000.0,

  /* Expression: Par_Limits(1)
   * Referenced by: '<S39>/Saturation2'
   */
  0.0,

  /* Expression: Par_Limits(2)
   * Referenced by: '<S39>/Saturation2'
   */
  0.0,

  /* Expression: 1/2/pi
   * Referenced by: '<S37>/Gain10'
   */
  0.15915494309189535,

  /* Expression: sps.MaxRateChangeFreq
   * Referenced by: '<S37>/Rate Limiter'
   */
  12.0,

  /* Expression: -sps.MaxRateChangeFreq
   * Referenced by: '<S37>/Rate Limiter'
   */
  -12.0,

  /* Expression: sps.x0(1,:)
   * Referenced by: '<S55>/Integrator_x1'
   */
  0.0024317084074161068,

  /* Expression: sps.A11
   * Referenced by: '<S56>/A11'
   */
  0.0,

  /* Expression: sps.x0(2,:)
   * Referenced by: '<S55>/Integrator_x2'
   */
  0.0,

  /* Expression: sps.A12
   * Referenced by: '<S56>/A12'
   */
  1.0,

  /* Expression: sps.A21
   * Referenced by: '<S56>/A21'
   */
  -24674.011002723397,

  /* Expression: sps.A22
   * Referenced by: '<S56>/A22'
   */
  -222.11060060879836,

  /* Expression: sps.B11
   * Referenced by: '<S57>/B11'
   */
  0.0,

  /* Expression: sps.B21
   * Referenced by: '<S57>/B21'
   */
  1.0,

  /* Expression: sps.C11
   * Referenced by: '<S58>/C11'
   */
  24674.011002723397,

  /* Expression: sps.C12
   * Referenced by: '<S58>/C12'
   */
  0.0,

  /* Expression: sps.D
   * Referenced by: '<S55>/D*u'
   */
  0.0,

  /* Expression: 1e6
   * Referenced by: '<S59>/To avoid division  by zero'
   */
  1.0E+6,

  /* Expression: eps
   * Referenced by: '<S59>/To avoid division  by zero'
   */
  2.2204460492503131E-16,

  /* Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
   * Referenced by: '<S61>/Gain3'
   */
  { 1.0, 0.0, 0.5, -0.5, 0.8660254037844386, 0.5, -0.5, -0.8660254037844386, 0.5
  },

  /* Expression: 2/3
   * Referenced by: '<S61>/Gain1'
   */
  0.66666666666666663,

  /* Expression: 2*pi
   * Referenced by: '<S4>/Gain1'
   */
  6.2831853071795862,

  /* Expression: 1
   * Referenced by: '<S28>/do not delete this gain'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S29>/do not delete this gain'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S30>/do not delete this gain'
   */
  1.0,

  /* Expression: Kv
   * Referenced by: '<S21>/Kv1'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S25>/do not delete this gain'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S26>/do not delete this gain'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S27>/do not delete this gain'
   */
  1.0,

  /* Expression: Ki
   * Referenced by: '<S21>/Kv'
   */
  1.0,

  /* Expression: 1/100e3
   * Referenced by: '<S4>/Gain3'
   */
  1.0E-5,

  /* Start of '<S127>/Subsystem1' */
  {
    /* Expression: [0,0]
     * Referenced by: '<S132>/dq'
     */
    { 0.0, 0.0 }
  }
  ,

  /* End of '<S127>/Subsystem1' */

  /* Start of '<S127>/Subsystem - pi//2 delay' */
  {
    /* Expression: [0,0]
     * Referenced by: '<S131>/dq'
     */
    { 0.0, 0.0 }
  }
  ,

  /* End of '<S127>/Subsystem - pi//2 delay' */

  /* Start of '<S121>/Subsystem1' */
  {
    /* Expression: [0,0]
     * Referenced by: '<S126>/dq'
     */
    { 0.0, 0.0 }
  }
  ,

  /* End of '<S121>/Subsystem1' */

  /* Start of '<S121>/Subsystem - pi//2 delay' */
  {
    /* Expression: [0,0]
     * Referenced by: '<S125>/dq'
     */
    { 0.0, 0.0 }
  }
  ,

  /* End of '<S121>/Subsystem - pi//2 delay' */

  /* Start of '<S60>/Subsystem1' */
  {
    /* Expression: [0,0]
     * Referenced by: '<S65>/dq'
     */
    { 0.0, 0.0 }
  }
  ,

  /* End of '<S60>/Subsystem1' */

  /* Start of '<S60>/Subsystem - pi//2 delay' */
  {
    /* Expression: [0,0]
     * Referenced by: '<S64>/dq'
     */
    { 0.0, 0.0 }
  }
  ,

  /* End of '<S60>/Subsystem - pi//2 delay' */

  /* Start of '<S49>/Subsystem1' */
  {
    /* Expression: [0,0]
     * Referenced by: '<S54>/dq'
     */
    { 0.0, 0.0 }
  }
  ,

  /* End of '<S49>/Subsystem1' */

  /* Start of '<S49>/Subsystem - pi//2 delay' */
  {
    /* Expression: [0,0]
     * Referenced by: '<S53>/dq'
     */
    { 0.0, 0.0 }
  }
  /* End of '<S49>/Subsystem - pi//2 delay' */
};
