/*
 * Voltage_Control.c
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

/* Block signals (default storage) */
B_Voltage_Control_T Voltage_Control_B;

/* Continuous states */
X_Voltage_Control_T Voltage_Control_X;

/* Periodic continuous states */
PeriodicIndX_Voltage_Control_T Voltage_Control_PeriodicIndX;
PeriodicRngX_Voltage_Control_T Voltage_Control_PeriodicRngX;

/* Block states (default storage) */
DW_Voltage_Control_T Voltage_Control_DW;

/* Previous zero-crossings (trigger) states */
PrevZCX_Voltage_Control_T Voltage_Control_PrevZCX;

/* External inputs (root inport signals with default storage) */
ExtU_Voltage_Control_T Voltage_Control_U;

/* External outputs (root outports fed by signals with default storage) */
ExtY_Voltage_Control_T Voltage_Control_Y;

/* Real-time model */
RT_MODEL_Voltage_Control_T Voltage_Control_M_;
RT_MODEL_Voltage_Control_T *const Voltage_Control_M = &Voltage_Control_M_;
static void rate_scheduler(void);

/* For variable transport delay block, find the real delay time */
real_T rt_VTDelayfindtDInterpolate(
  real_T x,real_T* tBuf,real_T* uBuf,real_T* xBuf,int_T bufSz,int_T head,int_T
  tail,int_T* pLast,real_T t,real_T tStart,boolean_T discrete,boolean_T
  minorStepAndTAtLastMajorOutput,real_T initOutput,real_T* appliedDelay)
{
  int_T n, k;
  real_T f;
  int_T kp1;
  real_T tminustD, tL, tR, uD, uL, uR, fU;
  if (minorStepAndTAtLastMajorOutput) {
    /* pretend that the entry at head has not been added */
    if (*pLast == head) {
      *pLast = (*pLast == 0) ? bufSz-1 : *pLast-1;
    }

    head = (head == 0) ? bufSz-1 : head-1;
  }

  /*
   * The loop below finds k such that:
   *      x(t)-x(tminustD) =1 or
   *      x - xBuf[k+1] <= 1.0 < x - xBuf[k]
   *
   * Note that we have:
   *
   * tStart = tBuf[0] < tBuf[1] < ... < tBuf[tail] < ... tBuf[head] <= t
   *      0 = xBuf[0] < xBuf[1] < ... < xBuf[tail] < ... xBuf[head] <  x
   *
   * This is true if we assume the direction of transport is always positive
   * such as a flow goes through a pipe from one end to another. However, for
   * model such as convey belt, the transportation can change direction. For
   * this case, there will be more than one solution to x(t)-x(tminustD) = 1,
   * should found the minimum tminustD and tminustD > 0. The search will not
   * be as efficient as the following code.
   */

  /*
   * when x<=1, physically it means the flow didn't reach the output yet,
   * t-tD will be less then zero, so force output to be the initial output
   */
  if (x <= 1) {
    return initOutput;
  }

  /*
   * if the x is monoton increase, only one solution. use k=pLast for now
   */
  k= *pLast;
  n = 0;
  for (;;) {
    n++;
    if (n>bufSz)
      break;
    if (x - xBuf[k] > 1.0) {
      /* move k forward, unless k = head */
      if (k == head) {
        /* xxx this situation means tD= appliedDelay = 0
         *
         * linearly interpolate using (tBuf[head], xBuf[head])
         * and (t,x) to find (tD,xD) such that: x - xD = 1.0
         */
        int_T km1;
        f = (x - 1.0 - xBuf[k]) / (x - xBuf[k]);
        tminustD = (1.0-f)*tBuf[k] + f*t;
        km1 = k-1;
        if (km1 < 0)
          km1 = bufSz-1;
        tL = tBuf[km1];
        tR = tBuf[k];
        uL = uBuf[km1];
        uR = uBuf[k];
        break;
      }

      kp1 = k+1;
      if (kp1 == bufSz)
        kp1 = 0;
      if (x - xBuf[kp1] <= 1.0) {
        /*
         * linearly interpolate using (tBuf[k], xBuf[k])
         * and  (tBuf[k+1], xBuf[k+1]) to find (tminustD,xD)
         * such that: x - xD = 1.0
         */
        f = (x - 1.0 - xBuf[k]) / (xBuf[kp1] - xBuf[k]);
        tL = tBuf[k];
        tR = tBuf[kp1];
        uL = uBuf[k];
        uR = uBuf[kp1];
        tminustD = (1.0-f)*tL + f*tR;
        break;
      }

      k = kp1;
    } else {
      /* moved k backward, unless k = tail */
      if (k == tail) {
        /* This situation means tminustD <= Ttail*/
        f = (x - 1.0)/xBuf[k];
        if (discrete) {
          return(uBuf[tail]);
        }

        kp1 = k+1;
        if (kp1 == bufSz)
          kp1 = 0;

        /* * linearly interpolate using (tStart, 0)
         * and  (tBuf[tail], xBuf[tail]) to find (tminustD,xD)
         * such that: x - xD = 1.0
         */

        /* Here it is better to use Tstart because since x>1, tminustD
         * must > 0. Since x is monotone increase, its linearity is
         * better.
         */
        tminustD = (1-f)*tStart + f*tBuf[k];

        /* linearly interpolate using (t[tail], x[tail])
         * and  (tBuf[tail+1], xBuf[tail+1]) to find (tminustD,xD)
         * such that: x - xD = 1.0.
         * For time delay block, use t[tail] and t[tail+1], not good
         * for transport delay block since it may give tminstD < 0
         */

        /*  f = (tBuf[kp1]-tBuf[k])/(xBuf[kp1]-xBuf[k]);
         *  tminustD = tBuf[kp1]-f*(1+xBuf[kp1]-x);
         */
        tL = tBuf[k];
        tR = tBuf[kp1];
        uL = uBuf[k];
        uR = uBuf[kp1];
        break;
      }

      k = k - 1;
      if (k < 0)
        k = bufSz-1;
    }
  }

  *pLast = k;
  if (tR == tL) {
    fU = 1.0;
  } else {
    fU = (tminustD-tL)/(tR-tL);
  }

  /* for discrete signal, no interpolation, use either uL or uR
   * depend on wehre tminustD is.
   */
  if (discrete) {
    uD= (fU > (1.0-fU))? uR: uL;
  } else {
    uD = (1.0-fU)*uL + fU*uR;
  }

  /* we want return tD= t-(t-tD);*/
  *appliedDelay = t-tminustD;
  return uD;
}

/*
 *   This function updates active task flag for each subrate.
 * The function is called at model base rate, hence the
 * generated code self-manages all its subrates.
 */
static void rate_scheduler(void)
{
  /* Compute which subrates run during the next base time step.  Subrates
   * are an integer multiple of the base rate counter.  Therefore, the subtask
   * counter is reset when it reaches its limit (zero means run).
   */
  (Voltage_Control_M->Timing.TaskCounters.TID[2])++;
  if ((Voltage_Control_M->Timing.TaskCounters.TID[2]) > 9) {/* Sample time: [0.0001s, 0.0s] */
    Voltage_Control_M->Timing.TaskCounters.TID[2] = 0;
  }
}

/* State reduction function */
void local_stateReduction(real_T* x, int_T* p, int_T n, real_T* r)
{
  int_T i, j;
  for (i = 0, j = 0; i < n; ++i, ++j) {
    int_T k = p[i];
    real_T lb = r[j++];
    real_T xk = x[k]-lb;
    real_T rk = r[j]-lb;
    int_T q = (int_T) floor(xk/rk);
    if (q) {
      x[k] = xk-q*rk+lb;
    }
  }
}

/*
 * This function updates continuous states using the ODE4 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE4_IntgData *id = (ODE4_IntgData *)rtsiGetSolverData(si);
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T *f3 = id->f[3];
  real_T temp;
  int_T i;
  int_T nXc = 15;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) memcpy(y, x,
                (uint_T)nXc*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  Voltage_Control_derivatives();

  /* f1 = f(t + (h/2), y + (h/2)*f0) */
  temp = 0.5 * h;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f0[i]);
  }

  rtsiSetT(si, t + temp);
  rtsiSetdX(si, f1);
  Voltage_Control_output();
  Voltage_Control_derivatives();

  /* f2 = f(t + (h/2), y + (h/2)*f1) */
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f1[i]);
  }

  rtsiSetdX(si, f2);
  Voltage_Control_output();
  Voltage_Control_derivatives();

  /* f3 = f(t + h, y + h*f2) */
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (h*f2[i]);
  }

  rtsiSetT(si, tnew);
  rtsiSetdX(si, f3);
  Voltage_Control_output();
  Voltage_Control_derivatives();

  /* tnew = t + h
     ynew = y + (h/6)*(f0 + 2*f1 + 2*f2 + 2*f3) */
  temp = h / 6.0;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + temp*(f0[i] + 2.0*f1[i] + 2.0*f2[i] + f3[i]);
  }

  local_stateReduction(x, rtsiGetPeriodicContStateIndices(si), 1,
                       rtsiGetPeriodicContStateRanges(si));
  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/*
 * System initialize for enable system:
 *    '<S49>/Subsystem - pi//2 delay'
 *    '<S60>/Subsystem - pi//2 delay'
 *    '<S121>/Subsystem - pi//2 delay'
 *    '<S127>/Subsystem - pi//2 delay'
 */
void Voltage__Subsystempi2delay_Init(B_Subsystempi2delay_Voltage_C_T *localB,
  P_Subsystempi2delay_Voltage_C_T *localP)
{
  /* SystemInitialize for Outport: '<S53>/dq' */
  localB->Fcn = localP->dq_Y0[0];
  localB->Fcn1 = localP->dq_Y0[1];
}

/*
 * Disable for enable system:
 *    '<S49>/Subsystem - pi//2 delay'
 *    '<S60>/Subsystem - pi//2 delay'
 *    '<S121>/Subsystem - pi//2 delay'
 *    '<S127>/Subsystem - pi//2 delay'
 */
void Volta_Subsystempi2delay_Disable(DW_Subsystempi2delay_Voltage__T *localDW)
{
  localDW->Subsystempi2delay_MODE = false;
}

/*
 * Output and update for enable system:
 *    '<S49>/Subsystem - pi//2 delay'
 *    '<S60>/Subsystem - pi//2 delay'
 *    '<S121>/Subsystem - pi//2 delay'
 *    '<S127>/Subsystem - pi//2 delay'
 */
void Voltage_Contr_Subsystempi2delay(RT_MODEL_Voltage_Control_T * const
  Voltage_Control_M, uint8_T rtu_Enable, const real_T rtu_alpha_beta[2], real_T
  rtu_wt, B_Subsystempi2delay_Voltage_C_T *localB,
  DW_Subsystempi2delay_Voltage__T *localDW)
{
  real_T Fcn_tmp;
  real_T Fcn_tmp_0;

  /* Outputs for Enabled SubSystem: '<S49>/Subsystem - pi//2 delay' incorporates:
   *  EnablePort: '<S53>/Enable'
   */
  if ((rtmIsMajorTimeStep(Voltage_Control_M) &&
       Voltage_Control_M->Timing.TaskCounters.TID[1] == 0) && rtmIsMajorTimeStep
      (Voltage_Control_M)) {
    if (rtu_Enable > 0) {
      localDW->Subsystempi2delay_MODE = true;
    } else {
      if (localDW->Subsystempi2delay_MODE) {
        Volta_Subsystempi2delay_Disable(localDW);
      }
    }
  }

  if (localDW->Subsystempi2delay_MODE) {
    /* Fcn: '<S53>/Fcn' incorporates:
     *  Fcn: '<S53>/Fcn1'
     */
    Fcn_tmp = cos(rtu_wt);
    Fcn_tmp_0 = sin(rtu_wt);
    localB->Fcn = rtu_alpha_beta[0] * Fcn_tmp_0 - rtu_alpha_beta[1] * Fcn_tmp;

    /* Fcn: '<S53>/Fcn1' */
    localB->Fcn1 = rtu_alpha_beta[0] * Fcn_tmp + rtu_alpha_beta[1] * Fcn_tmp_0;
  }

  /* End of Outputs for SubSystem: '<S49>/Subsystem - pi//2 delay' */
}

/*
 * System initialize for enable system:
 *    '<S49>/Subsystem1'
 *    '<S60>/Subsystem1'
 *    '<S121>/Subsystem1'
 *    '<S127>/Subsystem1'
 */
void Voltage_Control_Subsystem1_Init(B_Subsystem1_Voltage_Control_T *localB,
  P_Subsystem1_Voltage_Control_T *localP)
{
  /* SystemInitialize for Outport: '<S54>/dq' */
  localB->Fcn = localP->dq_Y0[0];
  localB->Fcn1 = localP->dq_Y0[1];
}

/*
 * Disable for enable system:
 *    '<S49>/Subsystem1'
 *    '<S60>/Subsystem1'
 *    '<S121>/Subsystem1'
 *    '<S127>/Subsystem1'
 */
void Voltage_Cont_Subsystem1_Disable(DW_Subsystem1_Voltage_Control_T *localDW)
{
  localDW->Subsystem1_MODE = false;
}

/*
 * Output and update for enable system:
 *    '<S49>/Subsystem1'
 *    '<S60>/Subsystem1'
 *    '<S121>/Subsystem1'
 *    '<S127>/Subsystem1'
 */
void Voltage_Control_Subsystem1(RT_MODEL_Voltage_Control_T * const
  Voltage_Control_M, uint8_T rtu_Enable, const real_T rtu_alpha_beta[2], real_T
  rtu_wt, B_Subsystem1_Voltage_Control_T *localB,
  DW_Subsystem1_Voltage_Control_T *localDW)
{
  real_T Fcn_tmp;
  real_T Fcn_tmp_0;

  /* Outputs for Enabled SubSystem: '<S49>/Subsystem1' incorporates:
   *  EnablePort: '<S54>/Enable'
   */
  if ((rtmIsMajorTimeStep(Voltage_Control_M) &&
       Voltage_Control_M->Timing.TaskCounters.TID[1] == 0) && rtmIsMajorTimeStep
      (Voltage_Control_M)) {
    if (rtu_Enable > 0) {
      localDW->Subsystem1_MODE = true;
    } else {
      if (localDW->Subsystem1_MODE) {
        Voltage_Cont_Subsystem1_Disable(localDW);
      }
    }
  }

  if (localDW->Subsystem1_MODE) {
    /* Fcn: '<S54>/Fcn' incorporates:
     *  Fcn: '<S54>/Fcn1'
     */
    Fcn_tmp = sin(rtu_wt);
    Fcn_tmp_0 = cos(rtu_wt);
    localB->Fcn = rtu_alpha_beta[0] * Fcn_tmp_0 + rtu_alpha_beta[1] * Fcn_tmp;

    /* Fcn: '<S54>/Fcn1' */
    localB->Fcn1 = -rtu_alpha_beta[0] * Fcn_tmp + rtu_alpha_beta[1] * Fcn_tmp_0;
  }

  /* End of Outputs for SubSystem: '<S49>/Subsystem1' */
}

real_T rt_hypotd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T a;
  a = fabs(u0);
  y = fabs(u1);
  if (a < y) {
    a /= y;
    y *= sqrt(a * a + 1.0);
  } else if (a > y) {
    y /= a;
    y = sqrt(y * y + 1.0) * a;
  } else {
    if (!rtIsNaN(y)) {
      y = a * 1.4142135623730951;
    }
  }

  return y;
}

/* Model output function */
void Voltage_Control_output(void)
{
  /* local block i/o variables */
  real_T rtb_VariableTransportDelay;
  real_T rtb_Gain1[3];
  real_T rtb_Gain1_b[3];
  real_T rtb_VariableTransportDelay_h;
  real_T rtb_VariableTransportDelay_e;
  boolean_T didZcEventOccur;
  int32_T k;
  real_T rtb_Gain4;
  real_T rtb_Gain1_l;
  real_T rtb_Add;
  real_T rtb_Gain3_idx_2;
  real_T rtb_Gain3_idx_1;
  real_T rtb_Gain3_idx_0;
  real_T rtb_UnitDelay_idx_0;
  real_T rtb_UnitDelay_idx_1;
  real_T Fcn_tmp;
  int32_T uData_tmp;
  int32_T uData_tmp_0;
  if (rtmIsMajorTimeStep(Voltage_Control_M)) {
    /* set solver stop time */
    if (!(Voltage_Control_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&Voltage_Control_M->solverInfo,
                            ((Voltage_Control_M->Timing.clockTickH0 + 1) *
        Voltage_Control_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&Voltage_Control_M->solverInfo,
                            ((Voltage_Control_M->Timing.clockTick0 + 1) *
        Voltage_Control_M->Timing.stepSize0 +
        Voltage_Control_M->Timing.clockTickH0 *
        Voltage_Control_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(Voltage_Control_M)) {
    Voltage_Control_M->Timing.t[0] = rtsiGetT(&Voltage_Control_M->solverInfo);
  }

  /* Sum: '<S4>/Add' incorporates:
   *  Constant: '<S4>/Constant'
   *  Gain: '<S4>/Gain'
   *  StateSpace: '<S4>/State-Space1'
   */
  rtb_Add = Voltage_Control_P.StateSpace1_C[1] *
    Voltage_Control_X.StateSpace1_CSTATE[1] * Voltage_Control_P.Gain_Gain +
    Voltage_Control_P.Constant_Value_a;

  /* Gain: '<S4>/Gain4' */
  rtb_Gain4 = Voltage_Control_P.Gain4_Gain * rtb_Add;

  /* Product: '<S4>/Product' incorporates:
   *  Constant: '<S4>/Constant2'
   *  Gain: '<S4>/Gain2'
   *  Integrator: '<S4>/Integrator'
   *  Sum: '<S4>/Add1'
   *  Sum: '<S4>/Add2'
   *  Trigonometry: '<S4>/Trigonometric Function'
   *  Trigonometry: '<S4>/Trigonometric Function1'
   *  Trigonometry: '<S4>/Trigonometric Function2'
   */
  Voltage_Control_B.Product[0] = Voltage_Control_P.Gain2_Gain * sin
    (Voltage_Control_X.Integrator_CSTATE) * rtb_Gain4;
  Voltage_Control_B.Product[1] = sin(Voltage_Control_X.Integrator_CSTATE -
    Voltage_Control_P.Constant2_Value) * Voltage_Control_P.Gain2_Gain *
    rtb_Gain4;
  Voltage_Control_B.Product[2] = sin(Voltage_Control_X.Integrator_CSTATE +
    Voltage_Control_P.Constant2_Value) * Voltage_Control_P.Gain2_Gain *
    rtb_Gain4;

  /* Gain: '<Root>/Gain4' incorporates:
   *  Inport: '<Root>/action'
   */
  rtb_Gain3_idx_0 = Voltage_Control_P.Gain4_Gain_n * Voltage_Control_U.action[0];
  rtb_Gain3_idx_1 = Voltage_Control_P.Gain4_Gain_n * Voltage_Control_U.action[1];
  rtb_Gain3_idx_2 = Voltage_Control_P.Gain4_Gain_n * 0.0;

  /* RelationalOperator: '<S37>/Relational Operator' incorporates:
   *  Constant: '<S37>/Constant2'
   */
  Voltage_Control_B.RelationalOperator = (Voltage_Control_X.Integrator_CSTATE_f >
    Voltage_Control_P.Constant2_Value_h);

  /* InitialCondition: '<S37>/Initial' incorporates:
   *  Constant: '<S37>/Constant2'
   *  Sum: '<S37>/Subtract'
   */
  rtb_Gain4 = Voltage_Control_M->Timing.t[0];
  if ((Voltage_Control_DW.Initial_FirstOutputTime == (rtMinusInf)) ||
      (Voltage_Control_DW.Initial_FirstOutputTime == rtb_Gain4)) {
    Voltage_Control_DW.Initial_FirstOutputTime = rtb_Gain4;
    Voltage_Control_B.Initial = Voltage_Control_P.Initial_Value;
  } else {
    Voltage_Control_B.Initial = Voltage_Control_X.Integrator_CSTATE_f -
      Voltage_Control_P.Constant2_Value_h;
  }

  /* Integrator: '<S37>/Integrator' */
  /* Limited  Integrator  */
  if (rtmIsMajorTimeStep(Voltage_Control_M)) {
    didZcEventOccur = (Voltage_Control_B.RelationalOperator &&
                       (Voltage_Control_PrevZCX.Integrator_Reset_ZCE != 1));
    Voltage_Control_PrevZCX.Integrator_Reset_ZCE =
      Voltage_Control_B.RelationalOperator;

    /* evaluate zero-crossings */
    if (didZcEventOccur || (Voltage_Control_DW.Integrator_IWORK != 0)) {
      Voltage_Control_X.Integrator_CSTATE_f = Voltage_Control_B.Initial;
      rtsiSetBlockStateForSolverChangedAtMajorStep
        (&Voltage_Control_M->solverInfo, true);
      rtsiSetContTimeOutputInconsistentWithStateAtMajorStep
        (&Voltage_Control_M->solverInfo, true);
    }
  }

  if (Voltage_Control_X.Integrator_CSTATE_f >=
      Voltage_Control_P.Integrator_UpperSat) {
    if (Voltage_Control_X.Integrator_CSTATE_f !=
        Voltage_Control_P.Integrator_UpperSat) {
      Voltage_Control_X.Integrator_CSTATE_f =
        Voltage_Control_P.Integrator_UpperSat;
      rtsiSetBlockStateForSolverChangedAtMajorStep
        (&Voltage_Control_M->solverInfo, true);
      rtsiSetContTimeOutputInconsistentWithStateAtMajorStep
        (&Voltage_Control_M->solverInfo, true);
    }
  } else {
    if ((Voltage_Control_X.Integrator_CSTATE_f <=
         Voltage_Control_P.Integrator_LowerSat) &&
        (Voltage_Control_X.Integrator_CSTATE_f !=
         Voltage_Control_P.Integrator_LowerSat)) {
      Voltage_Control_X.Integrator_CSTATE_f =
        Voltage_Control_P.Integrator_LowerSat;
      rtsiSetBlockStateForSolverChangedAtMajorStep
        (&Voltage_Control_M->solverInfo, true);
      rtsiSetContTimeOutputInconsistentWithStateAtMajorStep
        (&Voltage_Control_M->solverInfo, true);
    }
  }

  Voltage_Control_B.Integrator = Voltage_Control_X.Integrator_CSTATE_f;

  /* End of Integrator: '<S37>/Integrator' */
  if (rtmIsMajorTimeStep(Voltage_Control_M) &&
      Voltage_Control_M->Timing.TaskCounters.TID[1] == 0) {
    /* RelationalOperator: '<S135>/Compare' incorporates:
     *  Constant: '<S134>/Constant'
     *  Constant: '<S135>/Constant'
     */
    Voltage_Control_B.Compare = (uint8_T)
      (Voltage_Control_P.dq0toAlphaBetaZero_Alignment ==
       Voltage_Control_P.CompareToConstant_const_e);

    /* Outputs for Enabled SubSystem: '<S134>/Subsystem1' incorporates:
     *  EnablePort: '<S138>/Enable'
     */
    if (rtmIsMajorTimeStep(Voltage_Control_M)) {
      Voltage_Control_DW.Subsystem1_MODE = (Voltage_Control_B.Compare > 0);
    }

    /* End of Outputs for SubSystem: '<S134>/Subsystem1' */

    /* RelationalOperator: '<S136>/Compare' incorporates:
     *  Constant: '<S134>/Constant'
     *  Constant: '<S136>/Constant'
     */
    Voltage_Control_B.Compare_p = (uint8_T)
      (Voltage_Control_P.dq0toAlphaBetaZero_Alignment ==
       Voltage_Control_P.CompareToConstant1_const_m);

    /* Outputs for Enabled SubSystem: '<S134>/Subsystem - pi//2 delay' incorporates:
     *  EnablePort: '<S137>/Enable'
     */
    if (rtmIsMajorTimeStep(Voltage_Control_M)) {
      Voltage_Control_DW.Subsystempi2delay_MODE = (Voltage_Control_B.Compare_p >
        0);
    }

    /* End of Outputs for SubSystem: '<S134>/Subsystem - pi//2 delay' */
  }

  /* Outputs for Enabled SubSystem: '<S134>/Subsystem1' incorporates:
   *  EnablePort: '<S138>/Enable'
   */
  if (Voltage_Control_DW.Subsystem1_MODE) {
    /* Fcn: '<S138>/Fcn' incorporates:
     *  Fcn: '<S138>/Fcn1'
     */
    rtb_Gain1_l = sin(Voltage_Control_B.Integrator);
    Fcn_tmp = cos(Voltage_Control_B.Integrator);
    Voltage_Control_B.Fcn = rtb_Gain3_idx_0 * Fcn_tmp - rtb_Gain3_idx_1 *
      rtb_Gain1_l;

    /* Fcn: '<S138>/Fcn1' */
    Voltage_Control_B.Fcn1 = rtb_Gain3_idx_0 * rtb_Gain1_l + rtb_Gain3_idx_1 *
      Fcn_tmp;
  }

  /* End of Outputs for SubSystem: '<S134>/Subsystem1' */

  /* Outputs for Enabled SubSystem: '<S134>/Subsystem - pi//2 delay' incorporates:
   *  EnablePort: '<S137>/Enable'
   */
  if (Voltage_Control_DW.Subsystempi2delay_MODE) {
    /* Fcn: '<S137>/Fcn' */
    Voltage_Control_B.Fcn_b = rtb_Gain3_idx_0 * sin(Voltage_Control_B.Integrator)
      + rtb_Gain3_idx_1 * cos(Voltage_Control_B.Integrator);

    /* Fcn: '<S137>/Fcn1' */
    Voltage_Control_B.Fcn1_b = -rtb_Gain3_idx_0 * cos
      (Voltage_Control_B.Integrator) + rtb_Gain3_idx_1 * sin
      (Voltage_Control_B.Integrator);
  }

  /* End of Outputs for SubSystem: '<S134>/Subsystem - pi//2 delay' */

  /* Switch: '<S134>/Switch' */
  if (Voltage_Control_B.Compare != 0) {
    /* SignalConversion generated from: '<S133>/Gain3' */
    rtb_Gain3_idx_0 = Voltage_Control_B.Fcn;
    rtb_Gain3_idx_1 = Voltage_Control_B.Fcn1;
  } else {
    /* SignalConversion generated from: '<S133>/Gain3' */
    rtb_Gain3_idx_0 = Voltage_Control_B.Fcn_b;
    rtb_Gain3_idx_1 = Voltage_Control_B.Fcn1_b;
  }

  /* End of Switch: '<S134>/Switch' */

  /* Gain: '<S133>/Gain3' incorporates:
   *  SignalConversion generated from: '<S133>/Gain3'
   */
  for (uData_tmp = 0; uData_tmp < 3; uData_tmp++) {
    Voltage_Control_B.Gain3[uData_tmp] = 0.0;
    Voltage_Control_B.Gain3[uData_tmp] +=
      Voltage_Control_P.Gain3_Gain_l[uData_tmp] * rtb_Gain3_idx_0;
    Voltage_Control_B.Gain3[uData_tmp] +=
      Voltage_Control_P.Gain3_Gain_l[uData_tmp + 3] * rtb_Gain3_idx_1;
    Voltage_Control_B.Gain3[uData_tmp] +=
      Voltage_Control_P.Gain3_Gain_l[uData_tmp + 6] * rtb_Gain3_idx_2;
  }

  /* End of Gain: '<S133>/Gain3' */
  if (rtmIsMajorTimeStep(Voltage_Control_M) &&
      Voltage_Control_M->Timing.TaskCounters.TID[1] == 0) {
    /* UnitDelay: '<S67>/Unit Delay' */
    Voltage_Control_B.UnitDelay = Voltage_Control_DW.UnitDelay_DSTATE;

    /* UnitDelay: '<S68>/Unit Delay' */
    Voltage_Control_B.UnitDelay_b = Voltage_Control_DW.UnitDelay_DSTATE_h;

    /* UnitDelay: '<S69>/Unit Delay' */
    Voltage_Control_B.UnitDelay_i = Voltage_Control_DW.UnitDelay_DSTATE_i;

    /* S-Function (sfun_spssw_discc): '<S139>/State-Space' */

    /* S-Function block: <S139>/State-Space */
    {
      real_T accum;

      /*
       * Compute outputs:
       * ---------------
       */
      real_T *Cs = (real_T*)Voltage_Control_DW.StateSpace_PWORK.CS;
      real_T *Ds = (real_T*)Voltage_Control_DW.StateSpace_PWORK.DS;

      {
        int_T i1;
        real_T *y0 = &Voltage_Control_B.StateSpace[0];
        for (i1=0; i1 < 24; i1++) {
          accum = 0.0;

          {
            int_T i2;
            real_T *xd = &Voltage_Control_DW.StateSpace_DSTATE[0];
            for (i2=0; i2 < 6; i2++) {
              accum += *(Cs++) * xd[i2];
            }
          }

          accum += *(Ds++) * Voltage_Control_B.Product[0];
          accum += *(Ds++) * Voltage_Control_B.Product[1];
          accum += *(Ds++) * Voltage_Control_B.Product[2];
          accum += *(Ds++) * Voltage_Control_B.Gain3[0];
          accum += *(Ds++) * Voltage_Control_B.Gain3[1];
          accum += *(Ds++) * Voltage_Control_B.Gain3[2];
          accum += *(Ds++) * Voltage_Control_B.UnitDelay;
          accum += *(Ds++) * Voltage_Control_B.UnitDelay_b;
          accum += *(Ds++) * Voltage_Control_B.UnitDelay_i;
          y0[i1] = accum;
        }
      }
    }

    /* Gain: '<S10>/Kv' incorporates:
     *  Gain: '<S91>/do not delete this gain'
     *  Gain: '<S92>/do not delete this gain'
     *  Gain: '<S93>/do not delete this gain'
     */
    rtb_Gain3_idx_0 = Voltage_Control_P.donotdeletethisgain_Gain *
      Voltage_Control_B.StateSpace[18] * Voltage_Control_P.Kv_Gain;
    rtb_Gain3_idx_1 = Voltage_Control_P.donotdeletethisgain_Gain_n *
      Voltage_Control_B.StateSpace[19] * Voltage_Control_P.Kv_Gain;
    rtb_Gain3_idx_2 = Voltage_Control_P.donotdeletethisgain_Gain_f *
      Voltage_Control_B.StateSpace[20] * Voltage_Control_P.Kv_Gain;
    for (uData_tmp = 0; uData_tmp < 3; uData_tmp++) {
      /* Gain: '<S122>/Gain1' incorporates:
       *  Gain: '<S122>/Gain3'
       */
      Voltage_Control_B.Gain1[uData_tmp] = Voltage_Control_P.Gain1_Gain_f *
        (Voltage_Control_P.Gain3_Gain_o[uData_tmp + 6] * rtb_Gain3_idx_2 +
         (Voltage_Control_P.Gain3_Gain_o[uData_tmp + 3] * rtb_Gain3_idx_1 +
          Voltage_Control_P.Gain3_Gain_o[uData_tmp] * rtb_Gain3_idx_0));
    }

    /* RelationalOperator: '<S123>/Compare' incorporates:
     *  Constant: '<S121>/Constant'
     *  Constant: '<S123>/Constant'
     */
    Voltage_Control_B.Compare_m = (uint8_T)
      (Voltage_Control_P.AlphaBetaZerotodq0_Alignment_g ==
       Voltage_Control_P.CompareToConstant_const_p);
  }

  if (rtmIsMajorTimeStep(Voltage_Control_M) &&
      Voltage_Control_M->Timing.TaskCounters.TID[2] == 0) {
    /* UnitDelay: '<Root>/Unit Delay' */
    rtb_UnitDelay_idx_0 = Voltage_Control_DW.UnitDelay_DSTATE_n[0];
    rtb_UnitDelay_idx_1 = Voltage_Control_DW.UnitDelay_DSTATE_n[1];
  }

  /* Outputs for Enabled SubSystem: '<S121>/Subsystem1' */
  Voltage_Control_Subsystem1(Voltage_Control_M, Voltage_Control_B.Compare_m,
    &Voltage_Control_B.Gain1[0], Voltage_Control_B.Integrator,
    &Voltage_Control_B.Subsystem1_e, &Voltage_Control_DW.Subsystem1_e);

  /* End of Outputs for SubSystem: '<S121>/Subsystem1' */
  if (rtmIsMajorTimeStep(Voltage_Control_M) &&
      Voltage_Control_M->Timing.TaskCounters.TID[1] == 0) {
    /* RelationalOperator: '<S124>/Compare' incorporates:
     *  Constant: '<S121>/Constant'
     *  Constant: '<S124>/Constant'
     */
    Voltage_Control_B.Compare_a = (uint8_T)
      (Voltage_Control_P.AlphaBetaZerotodq0_Alignment_g ==
       Voltage_Control_P.CompareToConstant1_const_l);
  }

  /* Outputs for Enabled SubSystem: '<S121>/Subsystem - pi//2 delay' */
  Voltage_Contr_Subsystempi2delay(Voltage_Control_M, Voltage_Control_B.Compare_a,
    &Voltage_Control_B.Gain1[0], Voltage_Control_B.Integrator,
    &Voltage_Control_B.Subsystempi2delay_b,
    &Voltage_Control_DW.Subsystempi2delay_b);

  /* End of Outputs for SubSystem: '<S121>/Subsystem - pi//2 delay' */

  /* Switch: '<S121>/Switch' */
  if (Voltage_Control_B.Compare_m != 0) {
    rtb_Gain3_idx_2 = Voltage_Control_B.Subsystem1_e.Fcn;
    rtb_Gain3_idx_0 = Voltage_Control_B.Subsystem1_e.Fcn1;
  } else {
    rtb_Gain3_idx_2 = Voltage_Control_B.Subsystempi2delay_b.Fcn;
    rtb_Gain3_idx_0 = Voltage_Control_B.Subsystempi2delay_b.Fcn1;
  }

  /* End of Switch: '<S121>/Switch' */

  /* Gain: '<Root>/Gain' */
  rtb_Gain3_idx_1 = Voltage_Control_P.Gain_Gain_e * rtb_Gain3_idx_2;

  /* Gain: '<Root>/Gain1' */
  rtb_Gain1_l = Voltage_Control_P.Gain1_Gain_o * rtb_Gain3_idx_0;
  if (rtmIsMajorTimeStep(Voltage_Control_M) &&
      Voltage_Control_M->Timing.TaskCounters.TID[1] == 0) {
    /* Gain: '<S10>/Kv1' incorporates:
     *  Gain: '<S94>/do not delete this gain'
     *  Gain: '<S95>/do not delete this gain'
     *  Gain: '<S96>/do not delete this gain'
     */
    Voltage_Control_B.Kv1[0] = Voltage_Control_P.donotdeletethisgain_Gain_b *
      Voltage_Control_B.StateSpace[6] * Voltage_Control_P.Kv1_Gain;
    Voltage_Control_B.Kv1[1] = Voltage_Control_P.donotdeletethisgain_Gain_a *
      Voltage_Control_B.StateSpace[7] * Voltage_Control_P.Kv1_Gain;
    Voltage_Control_B.Kv1[2] = Voltage_Control_P.donotdeletethisgain_Gain_h *
      Voltage_Control_B.StateSpace[8] * Voltage_Control_P.Kv1_Gain;
    for (uData_tmp = 0; uData_tmp < 3; uData_tmp++) {
      /* Gain: '<S128>/Gain1' incorporates:
       *  Gain: '<S128>/Gain3'
       */
      Voltage_Control_B.Gain1_l[uData_tmp] = Voltage_Control_P.Gain1_Gain_h *
        (Voltage_Control_P.Gain3_Gain_d[uData_tmp + 6] * Voltage_Control_B.Kv1[2]
         + (Voltage_Control_P.Gain3_Gain_d[uData_tmp + 3] *
            Voltage_Control_B.Kv1[1] + Voltage_Control_P.Gain3_Gain_d[uData_tmp]
            * Voltage_Control_B.Kv1[0]));
    }

    /* RelationalOperator: '<S129>/Compare' incorporates:
     *  Constant: '<S127>/Constant'
     *  Constant: '<S129>/Constant'
     */
    Voltage_Control_B.Compare_f = (uint8_T)
      (Voltage_Control_P.AlphaBetaZerotodq0_Alignment_i ==
       Voltage_Control_P.CompareToConstant_const_m);
  }

  /* Outputs for Enabled SubSystem: '<S127>/Subsystem1' */
  Voltage_Control_Subsystem1(Voltage_Control_M, Voltage_Control_B.Compare_f,
    &Voltage_Control_B.Gain1_l[0], Voltage_Control_B.Integrator,
    &Voltage_Control_B.Subsystem1_ex, &Voltage_Control_DW.Subsystem1_ex);

  /* End of Outputs for SubSystem: '<S127>/Subsystem1' */
  if (rtmIsMajorTimeStep(Voltage_Control_M) &&
      Voltage_Control_M->Timing.TaskCounters.TID[1] == 0) {
    /* RelationalOperator: '<S130>/Compare' incorporates:
     *  Constant: '<S127>/Constant'
     *  Constant: '<S130>/Constant'
     */
    Voltage_Control_B.Compare_a5 = (uint8_T)
      (Voltage_Control_P.AlphaBetaZerotodq0_Alignment_i ==
       Voltage_Control_P.CompareToConstant1_const_h);
  }

  /* Outputs for Enabled SubSystem: '<S127>/Subsystem - pi//2 delay' */
  Voltage_Contr_Subsystempi2delay(Voltage_Control_M,
    Voltage_Control_B.Compare_a5, &Voltage_Control_B.Gain1_l[0],
    Voltage_Control_B.Integrator, &Voltage_Control_B.Subsystempi2delay_n,
    &Voltage_Control_DW.Subsystempi2delay_n);

  /* End of Outputs for SubSystem: '<S127>/Subsystem - pi//2 delay' */

  /* Switch: '<S127>/Switch' */
  if (Voltage_Control_B.Compare_f != 0) {
    rtb_Gain3_idx_2 = Voltage_Control_B.Subsystem1_ex.Fcn;
    rtb_Gain3_idx_0 = Voltage_Control_B.Subsystem1_ex.Fcn1;
  } else {
    rtb_Gain3_idx_2 = Voltage_Control_B.Subsystempi2delay_n.Fcn;
    rtb_Gain3_idx_0 = Voltage_Control_B.Subsystempi2delay_n.Fcn1;
  }

  /* End of Switch: '<S127>/Switch' */

  /* Gain: '<Root>/Gain2' */
  rtb_Gain3_idx_2 *= Voltage_Control_P.Gain2_Gain_d;

  /* Gain: '<Root>/Gain3' */
  rtb_Gain3_idx_0 *= Voltage_Control_P.Gain3_Gain_dd;

  /* Sum: '<Root>/Sum' incorporates:
   *  Inport: '<Root>/rand_input'
   */
  Voltage_Control_B.Sum[0] = Voltage_Control_U.rand_input[0] + rtb_Gain3_idx_1;
  Voltage_Control_B.Sum[1] = Voltage_Control_U.rand_input[1] + rtb_Gain1_l;
  Voltage_Control_B.Sum[2] = Voltage_Control_U.rand_input[2] + rtb_Gain3_idx_2;
  Voltage_Control_B.Sum[3] = Voltage_Control_U.rand_input[3] + rtb_Gain3_idx_0;
  if (rtmIsMajorTimeStep(Voltage_Control_M) &&
      Voltage_Control_M->Timing.TaskCounters.TID[2] == 0) {
    /* MATLAB Function: '<Root>/MATLAB Function1' incorporates:
     *  Inport: '<Root>/rand_param'
     *  SignalConversion generated from: '<S6>/ SFunction '
     */
    for (k = 0; k < 24; k++) {
      for (uData_tmp = 0; uData_tmp < 5; uData_tmp++) {
        Voltage_Control_DW.yData[uData_tmp + 5 * k] = Voltage_Control_DW.yData
          [(k + 1) * 5 + uData_tmp];
      }

      uData_tmp = (k + 1) << 1;
      uData_tmp_0 = k << 1;
      Voltage_Control_DW.uData[uData_tmp_0] = Voltage_Control_DW.uData[uData_tmp];
      Voltage_Control_DW.uData[uData_tmp_0 + 1] =
        Voltage_Control_DW.uData[uData_tmp + 1];
    }

    Voltage_Control_DW.yData[120] = Voltage_Control_B.Sum[0];
    Voltage_Control_DW.yData[121] = Voltage_Control_B.Sum[1];
    Voltage_Control_DW.yData[122] = Voltage_Control_B.Sum[2];
    Voltage_Control_DW.yData[123] = Voltage_Control_B.Sum[3];
    Voltage_Control_DW.yData[124] = Voltage_Control_U.rand_param[2];
    Voltage_Control_DW.uData[48] = rtb_UnitDelay_idx_0;
    Voltage_Control_DW.uData[49] = rtb_UnitDelay_idx_1;

    /* Outport: '<Root>/obs' incorporates:
     *  MATLAB Function: '<Root>/MATLAB Function1'
     */
    for (uData_tmp = 0; uData_tmp < 25; uData_tmp++) {
      Voltage_Control_Y.obs[uData_tmp] = Voltage_Control_DW.yData[5 * uData_tmp];
      Voltage_Control_Y.obs[uData_tmp + 25] = Voltage_Control_DW.yData[5 *
        uData_tmp + 1];
      Voltage_Control_Y.obs[uData_tmp + 50] = Voltage_Control_DW.yData[5 *
        uData_tmp + 2];
      Voltage_Control_Y.obs[uData_tmp + 75] = Voltage_Control_DW.yData[5 *
        uData_tmp + 3];
      Voltage_Control_Y.obs[uData_tmp + 100] = Voltage_Control_DW.yData[5 *
        uData_tmp + 4];

      /* MATLAB Function: '<Root>/MATLAB Function1' */
      k = uData_tmp << 1;
      Voltage_Control_Y.obs[uData_tmp + 125] = Voltage_Control_DW.uData[k];
      Voltage_Control_Y.obs[uData_tmp + 150] = Voltage_Control_DW.uData[k + 1];
    }

    /* End of Outport: '<Root>/obs' */
  }

  /* Outport: '<Root>/y_meas' */
  Voltage_Control_Y.y_meas[0] = Voltage_Control_B.Sum[0];
  Voltage_Control_Y.y_meas[1] = Voltage_Control_B.Sum[1];
  Voltage_Control_Y.y_meas[2] = Voltage_Control_B.Sum[2];
  Voltage_Control_Y.y_meas[3] = Voltage_Control_B.Sum[3];

  /* Outport: '<Root>/y_true' */
  Voltage_Control_Y.y_true[0] = rtb_Gain3_idx_1;
  Voltage_Control_Y.y_true[1] = rtb_Gain1_l;
  Voltage_Control_Y.y_true[2] = rtb_Gain3_idx_2;
  Voltage_Control_Y.y_true[3] = rtb_Gain3_idx_0;
  if (rtmIsMajorTimeStep(Voltage_Control_M) &&
      Voltage_Control_M->Timing.TaskCounters.TID[1] == 0) {
    /* Gain: '<S71>/do not delete this gain' */
    Voltage_Control_B.donotdeletethisgain =
      Voltage_Control_P.donotdeletethisgain_Gain_d *
      Voltage_Control_B.StateSpace[3];
  }

  /* MATLAB Function: '<S8>/MATLAB Function' incorporates:
   *  Constant: '<Root>/Constant'
   */
  rtb_Gain3_idx_2 = Voltage_Control_P.Constant_Value_h / 1.7320508075688772;

  /* Outputs for Enabled SubSystem: '<S37>/Automatic Gain Control' incorporates:
   *  EnablePort: '<S38>/Enable'
   */
  /* Clock: '<Root>/Clock' incorporates:
   *  Clock: '<S47>/Clock'
   *  Clock: '<S59>/Clock'
   */
  rtb_UnitDelay_idx_0 = Voltage_Control_M->Timing.t[0];

  /* End of Outputs for SubSystem: '<S37>/Automatic Gain Control' */

  /* MATLAB Function: '<Root>/MATLAB Function' incorporates:
   *  Clock: '<Root>/Clock'
   *  Inport: '<Root>/rand_param'
   *  Inport: '<Root>/t_step'
   */
  if (rtb_UnitDelay_idx_0 < Voltage_Control_U.t_step - 1.0E-5) {
    rtb_UnitDelay_idx_1 = Voltage_Control_U.rand_param[0];
  } else {
    rtb_UnitDelay_idx_1 = Voltage_Control_U.rand_param[0] +
      Voltage_Control_U.rand_param[1];
  }

  /* End of MATLAB Function: '<Root>/MATLAB Function' */

  /* MATLAB Function: '<S8>/MATLAB Function' incorporates:
   *  Gain: '<Root>/Gain9'
   */
  rtb_UnitDelay_idx_1 = rtb_Gain3_idx_2 * rtb_Gain3_idx_2 /
    (Voltage_Control_P.Gain9_Gain * rtb_UnitDelay_idx_1 / 3.0);

  /* Saturate: '<S67>/Saturation' */
  if (rtb_UnitDelay_idx_1 > Voltage_Control_P.Saturation_UpperSat_e) {
    rtb_Gain3_idx_2 = Voltage_Control_P.Saturation_UpperSat_e;
  } else if (rtb_UnitDelay_idx_1 < Voltage_Control_P.Saturation_LowerSat_c) {
    rtb_Gain3_idx_2 = Voltage_Control_P.Saturation_LowerSat_c;
  } else {
    rtb_Gain3_idx_2 = rtb_UnitDelay_idx_1;
  }

  /* End of Saturate: '<S67>/Saturation' */

  /* Product: '<S67>/Divide' */
  Voltage_Control_B.Divide = Voltage_Control_B.donotdeletethisgain /
    rtb_Gain3_idx_2;
  if (rtmIsMajorTimeStep(Voltage_Control_M) &&
      Voltage_Control_M->Timing.TaskCounters.TID[1] == 0) {
    /* Gain: '<S74>/do not delete this gain' */
    Voltage_Control_B.donotdeletethisgain_a =
      Voltage_Control_P.donotdeletethisgain_Gain_nm *
      Voltage_Control_B.StateSpace[4];
  }

  /* Saturate: '<S68>/Saturation' */
  if (rtb_UnitDelay_idx_1 > Voltage_Control_P.Saturation_UpperSat_p) {
    rtb_Gain3_idx_2 = Voltage_Control_P.Saturation_UpperSat_p;
  } else if (rtb_UnitDelay_idx_1 < Voltage_Control_P.Saturation_LowerSat_e) {
    rtb_Gain3_idx_2 = Voltage_Control_P.Saturation_LowerSat_e;
  } else {
    rtb_Gain3_idx_2 = rtb_UnitDelay_idx_1;
  }

  /* End of Saturate: '<S68>/Saturation' */

  /* Product: '<S68>/Divide' */
  Voltage_Control_B.Divide_b = Voltage_Control_B.donotdeletethisgain_a /
    rtb_Gain3_idx_2;
  if (rtmIsMajorTimeStep(Voltage_Control_M) &&
      Voltage_Control_M->Timing.TaskCounters.TID[1] == 0) {
    /* Gain: '<S77>/do not delete this gain' */
    Voltage_Control_B.donotdeletethisgain_p =
      Voltage_Control_P.donotdeletethisgain_Gain_ft *
      Voltage_Control_B.StateSpace[5];
  }

  /* Saturate: '<S69>/Saturation' */
  if (rtb_UnitDelay_idx_1 > Voltage_Control_P.Saturation_UpperSat_f) {
    rtb_UnitDelay_idx_1 = Voltage_Control_P.Saturation_UpperSat_f;
  } else {
    if (rtb_UnitDelay_idx_1 < Voltage_Control_P.Saturation_LowerSat_a) {
      rtb_UnitDelay_idx_1 = Voltage_Control_P.Saturation_LowerSat_a;
    }
  }

  /* End of Saturate: '<S69>/Saturation' */

  /* Product: '<S69>/Divide' */
  Voltage_Control_B.Divide_n = Voltage_Control_B.donotdeletethisgain_p /
    rtb_UnitDelay_idx_1;
  if (rtmIsMajorTimeStep(Voltage_Control_M) &&
      Voltage_Control_M->Timing.TaskCounters.TID[1] == 0) {
    /* Memory: '<S37>/Memory' */
    Voltage_Control_B.Memory = Voltage_Control_DW.Memory_PreviousInput;

    /* Outputs for Enabled SubSystem: '<S37>/Automatic Gain Control' incorporates:
     *  EnablePort: '<S38>/Enable'
     */
    if (rtmIsMajorTimeStep(Voltage_Control_M)) {
      /* Constant: '<S37>/Constant1' */
      if (Voltage_Control_P.Constant1_Value > 0.0) {
        Voltage_Control_DW.AutomaticGainControl_MODE = true;
      } else {
        if (Voltage_Control_DW.AutomaticGainControl_MODE) {
          /* Disable for Enabled SubSystem: '<S49>/Subsystem - pi//2 delay' */
          if (Voltage_Control_DW.Subsystempi2delay_d.Subsystempi2delay_MODE) {
            Volta_Subsystempi2delay_Disable
              (&Voltage_Control_DW.Subsystempi2delay_d);
          }

          /* End of Disable for SubSystem: '<S49>/Subsystem - pi//2 delay' */

          /* Disable for Enabled SubSystem: '<S49>/Subsystem1' */
          if (Voltage_Control_DW.Subsystem1_h.Subsystem1_MODE) {
            Voltage_Cont_Subsystem1_Disable(&Voltage_Control_DW.Subsystem1_h);
          }

          /* End of Disable for SubSystem: '<S49>/Subsystem1' */
          Voltage_Control_DW.AutomaticGainControl_MODE = false;
        }
      }

      /* End of Constant: '<S37>/Constant1' */
    }

    /* End of Outputs for SubSystem: '<S37>/Automatic Gain Control' */
  }

  /* Outputs for Enabled SubSystem: '<S37>/Automatic Gain Control' incorporates:
   *  EnablePort: '<S38>/Enable'
   */
  if (Voltage_Control_DW.AutomaticGainControl_MODE) {
    for (uData_tmp = 0; uData_tmp < 3; uData_tmp++) {
      /* Gain: '<S50>/Gain1' incorporates:
       *  Gain: '<S50>/Gain3'
       */
      rtb_Gain1_b[uData_tmp] = Voltage_Control_P.Gain1_Gain *
        (Voltage_Control_P.Gain3_Gain[uData_tmp + 6] * Voltage_Control_B.Kv1[2]
         + (Voltage_Control_P.Gain3_Gain[uData_tmp + 3] * Voltage_Control_B.Kv1
            [1] + Voltage_Control_P.Gain3_Gain[uData_tmp] *
            Voltage_Control_B.Kv1[0]));
    }

    if (rtmIsMajorTimeStep(Voltage_Control_M) &&
        Voltage_Control_M->Timing.TaskCounters.TID[1] == 0) {
      /* RelationalOperator: '<S52>/Compare' incorporates:
       *  Constant: '<S49>/Constant'
       *  Constant: '<S52>/Constant'
       */
      Voltage_Control_B.Compare_l = (uint8_T)
        (Voltage_Control_P.AlphaBetaZerotodq0_Alignment ==
         Voltage_Control_P.CompareToConstant1_const);
    }

    /* Outputs for Enabled SubSystem: '<S49>/Subsystem - pi//2 delay' */
    Voltage_Contr_Subsystempi2delay(Voltage_Control_M,
      Voltage_Control_B.Compare_l, &rtb_Gain1_b[0], Voltage_Control_B.Integrator,
      &Voltage_Control_B.Subsystempi2delay_d,
      &Voltage_Control_DW.Subsystempi2delay_d);

    /* End of Outputs for SubSystem: '<S49>/Subsystem - pi//2 delay' */
    if (rtmIsMajorTimeStep(Voltage_Control_M) &&
        Voltage_Control_M->Timing.TaskCounters.TID[1] == 0) {
      /* RelationalOperator: '<S51>/Compare' incorporates:
       *  Constant: '<S49>/Constant'
       *  Constant: '<S51>/Constant'
       */
      Voltage_Control_B.Compare_pa = (uint8_T)
        (Voltage_Control_P.AlphaBetaZerotodq0_Alignment ==
         Voltage_Control_P.CompareToConstant_const);
    }

    /* Outputs for Enabled SubSystem: '<S49>/Subsystem1' */
    Voltage_Control_Subsystem1(Voltage_Control_M, Voltage_Control_B.Compare_pa,
      &rtb_Gain1_b[0], Voltage_Control_B.Integrator,
      &Voltage_Control_B.Subsystem1_h, &Voltage_Control_DW.Subsystem1_h);

    /* End of Outputs for SubSystem: '<S49>/Subsystem1' */

    /* Switch: '<S49>/Switch' */
    if (Voltage_Control_B.Compare_pa != 0) {
      Voltage_Control_B.Switch_ow[0] = Voltage_Control_B.Subsystem1_h.Fcn;
      Voltage_Control_B.Switch_ow[1] = Voltage_Control_B.Subsystem1_h.Fcn1;
    } else {
      Voltage_Control_B.Switch_ow[0] = Voltage_Control_B.Subsystempi2delay_d.Fcn;
      Voltage_Control_B.Switch_ow[1] =
        Voltage_Control_B.Subsystempi2delay_d.Fcn1;
    }

    /* End of Switch: '<S49>/Switch' */

    /* VariableTransportDelay: '<S47>/Variable Transport Delay' */
    {
      real_T **uBuffer = (real_T**)
        &Voltage_Control_DW.VariableTransportDelay_PWORK_p.TUbufferPtrs[0];
      real_T **tBuffer = (real_T**)
        &Voltage_Control_DW.VariableTransportDelay_PWORK_p.TUbufferPtrs[1];
      real_T **xBuffer = (real_T**)
        &Voltage_Control_DW.VariableTransportDelay_PWORK_p.TUbufferPtrs[2];
      real_T simTime = Voltage_Control_M->Timing.t[0];
      real_T appliedDelay;

      /* For variable transport dealy, find the real applied dealy
       * here and then output
       */
      rtb_VariableTransportDelay_h= rt_VTDelayfindtDInterpolate
        (Voltage_Control_X.VariableTransportDelay_CSTATE_n,*tBuffer,*uBuffer,
         *xBuffer,
         Voltage_Control_DW.VariableTransportDelay_IWORK_l.CircularBufSize,
         Voltage_Control_DW.VariableTransportDelay_IWORK_l.Head,
         Voltage_Control_DW.VariableTransportDelay_IWORK_l.Tail,
         &Voltage_Control_DW.VariableTransportDelay_IWORK_l.Last, simTime, 0.0,0,
         0, Voltage_Control_P.VariableTransportDelay_InitOutp,
         &appliedDelay);
    }

    /* Integrator: '<S47>/integrator' */
    Voltage_Control_B.integrator_p = Voltage_Control_X.integrator_CSTATE_m;
    if (rtmIsMajorTimeStep(Voltage_Control_M) &&
        Voltage_Control_M->Timing.TaskCounters.TID[1] == 0) {
      /* Memory: '<S47>/Memory' */
      Voltage_Control_B.Memory_p = Voltage_Control_DW.Memory_PreviousInput_d;
    }

    /* Switch: '<S47>/Switch' incorporates:
     *  Constant: '<S47>/Constant'
     *  Product: '<S47>/Product'
     *  RelationalOperator: '<S47>/Relational Operator'
     *  Sum: '<S47>/Sum7'
     */
    if (rtb_UnitDelay_idx_0 >= Voltage_Control_P.Constant_Value) {
      Voltage_Control_B.Switch_j = (Voltage_Control_B.integrator_p -
        rtb_VariableTransportDelay_h) * Voltage_Control_B.Memory;
    } else {
      Voltage_Control_B.Switch_j = Voltage_Control_B.Memory_p;
    }

    /* End of Switch: '<S47>/Switch' */

    /* VariableTransportDelay: '<S48>/Variable Transport Delay' */
    {
      real_T **uBuffer = (real_T**)
        &Voltage_Control_DW.VariableTransportDelay_PWORK_f.TUbufferPtrs[0];
      real_T **tBuffer = (real_T**)
        &Voltage_Control_DW.VariableTransportDelay_PWORK_f.TUbufferPtrs[1];
      real_T **xBuffer = (real_T**)
        &Voltage_Control_DW.VariableTransportDelay_PWORK_f.TUbufferPtrs[2];
      real_T simTime = Voltage_Control_M->Timing.t[0];
      real_T appliedDelay;

      /* For variable transport dealy, find the real applied dealy
       * here and then output
       */
      rtb_VariableTransportDelay_e= rt_VTDelayfindtDInterpolate
        (Voltage_Control_X.VariableTransportDelay_CSTAT_nz,*tBuffer,*uBuffer,
         *xBuffer,
         Voltage_Control_DW.VariableTransportDelay_IWORK_g.CircularBufSize,
         Voltage_Control_DW.VariableTransportDelay_IWORK_g.Head,
         Voltage_Control_DW.VariableTransportDelay_IWORK_g.Tail,
         &Voltage_Control_DW.VariableTransportDelay_IWORK_g.Last, simTime, 0.0,0,
         0, Voltage_Control_P.VariableTransportDelay_InitOu_n,
         &appliedDelay);
    }

    /* Integrator: '<S48>/integrator' */
    Voltage_Control_B.integrator_m = Voltage_Control_X.integrator_CSTATE_l;
    if (rtmIsMajorTimeStep(Voltage_Control_M) &&
        Voltage_Control_M->Timing.TaskCounters.TID[1] == 0) {
      /* Memory: '<S48>/Memory' */
      Voltage_Control_B.Memory_j = Voltage_Control_DW.Memory_PreviousInput_o;
    }

    /* Switch: '<S48>/Switch' incorporates:
     *  Constant: '<S48>/Constant'
     *  Product: '<S48>/Product'
     *  RelationalOperator: '<S48>/Relational Operator'
     *  Sum: '<S48>/Sum7'
     */
    if (rtb_UnitDelay_idx_0 >= Voltage_Control_P.Constant_Value_n) {
      Voltage_Control_B.Switch_i = (Voltage_Control_B.integrator_m -
        rtb_VariableTransportDelay_e) * Voltage_Control_B.Memory;
    } else {
      Voltage_Control_B.Switch_i = Voltage_Control_B.Memory_j;
    }

    /* End of Switch: '<S48>/Switch' */

    /* ComplexToMagnitudeAngle: '<S43>/Complex to Magnitude-Angle' incorporates:
     *  RealImagToComplex: '<S43>/Real-Imag to Complex'
     */
    rtb_UnitDelay_idx_1 = rt_hypotd_snf(Voltage_Control_B.Switch_j,
      Voltage_Control_B.Switch_i);

    /* Saturate: '<S38>/Saturation' */
    if (rtb_UnitDelay_idx_1 > Voltage_Control_P.Saturation_UpperSat) {
      rtb_UnitDelay_idx_1 = Voltage_Control_P.Saturation_UpperSat;
    } else {
      if (rtb_UnitDelay_idx_1 < Voltage_Control_P.Saturation_LowerSat) {
        rtb_UnitDelay_idx_1 = Voltage_Control_P.Saturation_LowerSat;
      }
    }

    /* End of Saturate: '<S38>/Saturation' */

    /* Math: '<S38>/Math Function'
     *
     * About '<S38>/Math Function':
     *  Operator: reciprocal
     */
    Voltage_Control_B.MathFunction = 1.0 / rtb_UnitDelay_idx_1;

    /* Saturate: '<S47>/To avoid division  by zero' */
    if (Voltage_Control_B.Memory >
        Voltage_Control_P.Toavoiddivisionbyzero_UpperSat) {
      rtb_UnitDelay_idx_1 = Voltage_Control_P.Toavoiddivisionbyzero_UpperSat;
    } else if (Voltage_Control_B.Memory <
               Voltage_Control_P.Toavoiddivisionbyzero_LowerSat) {
      rtb_UnitDelay_idx_1 = Voltage_Control_P.Toavoiddivisionbyzero_LowerSat;
    } else {
      rtb_UnitDelay_idx_1 = Voltage_Control_B.Memory;
    }

    /* End of Saturate: '<S47>/To avoid division  by zero' */

    /* Fcn: '<S47>/period' */
    Voltage_Control_B.period_p = 1.0 / rtb_UnitDelay_idx_1;

    /* Saturate: '<S48>/To avoid division  by zero' */
    if (Voltage_Control_B.Memory >
        Voltage_Control_P.Toavoiddivisionbyzero_UpperSa_h) {
      rtb_UnitDelay_idx_1 = Voltage_Control_P.Toavoiddivisionbyzero_UpperSa_h;
    } else if (Voltage_Control_B.Memory <
               Voltage_Control_P.Toavoiddivisionbyzero_LowerSa_h) {
      rtb_UnitDelay_idx_1 = Voltage_Control_P.Toavoiddivisionbyzero_LowerSa_h;
    } else {
      rtb_UnitDelay_idx_1 = Voltage_Control_B.Memory;
    }

    /* End of Saturate: '<S48>/To avoid division  by zero' */

    /* Fcn: '<S48>/period' */
    Voltage_Control_B.period_a = 1.0 / rtb_UnitDelay_idx_1;
  }

  /* End of Outputs for SubSystem: '<S37>/Automatic Gain Control' */

  /* Integrator: '<S39>/Integrator' */
  /* Limited  Integrator  */
  if (Voltage_Control_X.Integrator_CSTATE_fe >=
      Voltage_Control_P.Integrator_UpperSat_a) {
    if (Voltage_Control_X.Integrator_CSTATE_fe >
        Voltage_Control_P.Integrator_UpperSat_a) {
      rtsiSetBlockStateForSolverChangedAtMajorStep
        (&Voltage_Control_M->solverInfo, true);
    }

    Voltage_Control_X.Integrator_CSTATE_fe =
      Voltage_Control_P.Integrator_UpperSat_a;
  } else {
    if (Voltage_Control_X.Integrator_CSTATE_fe <=
        Voltage_Control_P.Integrator_LowerSat_n) {
      if (Voltage_Control_X.Integrator_CSTATE_fe <
          Voltage_Control_P.Integrator_LowerSat_n) {
        rtsiSetBlockStateForSolverChangedAtMajorStep
          (&Voltage_Control_M->solverInfo, true);
      }

      Voltage_Control_X.Integrator_CSTATE_fe =
        Voltage_Control_P.Integrator_LowerSat_n;
    }
  }

  rtb_Gain3_idx_0 = Voltage_Control_X.Integrator_CSTATE_fe;

  /* End of Integrator: '<S39>/Integrator' */

  /* VariableTransportDelay: '<S59>/Variable Transport Delay' */
  {
    real_T **uBuffer = (real_T**)
      &Voltage_Control_DW.VariableTransportDelay_PWORK.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)
      &Voltage_Control_DW.VariableTransportDelay_PWORK.TUbufferPtrs[1];
    real_T **xBuffer = (real_T**)
      &Voltage_Control_DW.VariableTransportDelay_PWORK.TUbufferPtrs[2];
    real_T simTime = Voltage_Control_M->Timing.t[0];
    real_T appliedDelay;

    /* For variable transport dealy, find the real applied dealy
     * here and then output
     */
    rtb_VariableTransportDelay= rt_VTDelayfindtDInterpolate
      (Voltage_Control_X.VariableTransportDelay_CSTATE,*tBuffer,*uBuffer,
       *xBuffer,
       Voltage_Control_DW.VariableTransportDelay_IWORK.CircularBufSize,
       Voltage_Control_DW.VariableTransportDelay_IWORK.Head,
       Voltage_Control_DW.VariableTransportDelay_IWORK.Tail,
       &Voltage_Control_DW.VariableTransportDelay_IWORK.Last, simTime, 0.0,0,
       0, Voltage_Control_P.VariableTransportDelay_InitOu_h,
       &appliedDelay);
  }

  /* Integrator: '<S59>/integrator' */
  Voltage_Control_B.integrator = Voltage_Control_X.integrator_CSTATE;
  if (rtmIsMajorTimeStep(Voltage_Control_M) &&
      Voltage_Control_M->Timing.TaskCounters.TID[1] == 0) {
    /* Memory: '<S59>/Memory' */
    Voltage_Control_B.Memory_h = Voltage_Control_DW.Memory_PreviousInput_f;
  }

  /* Switch: '<S59>/Switch' incorporates:
   *  Constant: '<S59>/Constant'
   *  Product: '<S59>/Product'
   *  RelationalOperator: '<S59>/Relational Operator'
   *  Sum: '<S59>/Sum7'
   */
  if (rtb_UnitDelay_idx_0 >= Voltage_Control_P.Constant_Value_k) {
    Voltage_Control_B.Switch = (Voltage_Control_B.integrator -
      rtb_VariableTransportDelay) * Voltage_Control_B.Memory;
  } else {
    Voltage_Control_B.Switch = Voltage_Control_B.Memory_h;
  }

  /* End of Switch: '<S59>/Switch' */

  /* Product: '<S37>/Divide' */
  rtb_Gain3_idx_2 = Voltage_Control_B.Switch * Voltage_Control_B.MathFunction;

  /* Gain: '<S39>/Kp5' */
  Voltage_Control_B.Kp5 = Voltage_Control_P.Continuous_Ki * rtb_Gain3_idx_2;

  /* Gain: '<S39>/Kp6' */
  Voltage_Control_B.Kp6 = Voltage_Control_P.Continuous_Kd * rtb_Gain3_idx_2;

  /* Sum: '<S39>/Sum6' incorporates:
   *  Gain: '<S39>/Kp4'
   *  TransferFcn: '<S39>/Transfer Fcn'
   */
  rtb_Gain3_idx_0 = (Voltage_Control_P.Continuous_Kp * rtb_Gain3_idx_2 +
                     rtb_Gain3_idx_0) + (Voltage_Control_P.TransferFcn_C *
    Voltage_Control_X.TransferFcn_CSTATE + Voltage_Control_P.TransferFcn_D *
    Voltage_Control_B.Kp6);

  /* Saturate: '<S39>/Saturation2' */
  if (rtb_Gain3_idx_0 > Voltage_Control_P.Saturation2_UpperSat) {
    Voltage_Control_B.Saturation2 = Voltage_Control_P.Saturation2_UpperSat;
  } else if (rtb_Gain3_idx_0 < Voltage_Control_P.Saturation2_LowerSat) {
    Voltage_Control_B.Saturation2 = Voltage_Control_P.Saturation2_LowerSat;
  } else {
    Voltage_Control_B.Saturation2 = rtb_Gain3_idx_0;
  }

  /* End of Saturate: '<S39>/Saturation2' */

  /* Gain: '<S37>/Gain10' */
  Voltage_Control_B.RateLimiter = Voltage_Control_P.Gain10_Gain *
    Voltage_Control_B.Saturation2;

  /* RateLimiter: '<S37>/Rate Limiter' incorporates:
   *  InitialCondition: '<S37>/Initial'
   */
  if (!(Voltage_Control_DW.LastMajorTime == (rtInf))) {
    rtb_UnitDelay_idx_0 = rtb_Gain4 - Voltage_Control_DW.LastMajorTime;
    rtb_UnitDelay_idx_1 = rtb_UnitDelay_idx_0 *
      Voltage_Control_P.RateLimiter_RisingLim;
    rtb_Gain4 = Voltage_Control_B.RateLimiter - Voltage_Control_DW.PrevY;
    if (rtb_Gain4 > rtb_UnitDelay_idx_1) {
      Voltage_Control_B.RateLimiter = Voltage_Control_DW.PrevY +
        rtb_UnitDelay_idx_1;
    } else {
      rtb_UnitDelay_idx_0 *= Voltage_Control_P.RateLimiter_FallingLim;
      if (rtb_Gain4 < rtb_UnitDelay_idx_0) {
        Voltage_Control_B.RateLimiter = Voltage_Control_DW.PrevY +
          rtb_UnitDelay_idx_0;
      }
    }
  }

  /* End of RateLimiter: '<S37>/Rate Limiter' */

  /* Sum: '<S55>/A*x1+ B* u1' incorporates:
   *  Gain: '<S56>/A11'
   *  Gain: '<S56>/A12'
   *  Gain: '<S57>/B11'
   *  Integrator: '<S55>/Integrator_x1'
   *  Integrator: '<S55>/Integrator_x2'
   *  Sum: '<S56>/sum2'
   */
  Voltage_Control_B.x1 = (Voltage_Control_P.A11_Gain *
    Voltage_Control_X.Integrator_x1_CSTATE + Voltage_Control_P.A12_Gain *
    Voltage_Control_X.Integrator_x2_CSTATE) + Voltage_Control_P.B11_Gain *
    Voltage_Control_B.RateLimiter;

  /* Sum: '<S55>/A*x2+ B*u2' incorporates:
   *  Gain: '<S56>/A21'
   *  Gain: '<S56>/A22'
   *  Gain: '<S57>/B21'
   *  Integrator: '<S55>/Integrator_x1'
   *  Integrator: '<S55>/Integrator_x2'
   *  Sum: '<S56>/sum3'
   */
  Voltage_Control_B.x2 = (Voltage_Control_P.A21_Gain *
    Voltage_Control_X.Integrator_x1_CSTATE + Voltage_Control_P.A22_Gain *
    Voltage_Control_X.Integrator_x2_CSTATE) + Voltage_Control_P.B21_Gain *
    Voltage_Control_B.RateLimiter;

  /* Sum: '<S55>/C*x + D*u' incorporates:
   *  Gain: '<S55>/D*u'
   *  Gain: '<S58>/C11'
   *  Gain: '<S58>/C12'
   *  Integrator: '<S55>/Integrator_x1'
   *  Integrator: '<S55>/Integrator_x2'
   *  Sum: '<S58>/sum2'
   */
  Voltage_Control_B.y = (Voltage_Control_P.C11_Gain *
    Voltage_Control_X.Integrator_x1_CSTATE + Voltage_Control_P.C12_Gain *
    Voltage_Control_X.Integrator_x2_CSTATE) + Voltage_Control_P.Du_Gain *
    Voltage_Control_B.RateLimiter;

  /* Saturate: '<S59>/To avoid division  by zero' */
  if (Voltage_Control_B.Memory >
      Voltage_Control_P.Toavoiddivisionbyzero_UpperSa_k) {
    rtb_Gain4 = Voltage_Control_P.Toavoiddivisionbyzero_UpperSa_k;
  } else if (Voltage_Control_B.Memory <
             Voltage_Control_P.Toavoiddivisionbyzero_LowerSa_c) {
    rtb_Gain4 = Voltage_Control_P.Toavoiddivisionbyzero_LowerSa_c;
  } else {
    rtb_Gain4 = Voltage_Control_B.Memory;
  }

  /* End of Saturate: '<S59>/To avoid division  by zero' */

  /* Fcn: '<S59>/period' */
  Voltage_Control_B.period = 1.0 / rtb_Gain4;
  if (rtmIsMajorTimeStep(Voltage_Control_M) &&
      Voltage_Control_M->Timing.TaskCounters.TID[1] == 0) {
    /* RelationalOperator: '<S62>/Compare' incorporates:
     *  Constant: '<S60>/Constant'
     *  Constant: '<S62>/Constant'
     */
    Voltage_Control_B.Compare_o = (uint8_T)
      (Voltage_Control_P.AlphaBetaZerotodq0_Alignment_c ==
       Voltage_Control_P.CompareToConstant_const_d);

    /* RelationalOperator: '<S63>/Compare' incorporates:
     *  Constant: '<S60>/Constant'
     *  Constant: '<S63>/Constant'
     */
    Voltage_Control_B.Compare_n = (uint8_T)
      (Voltage_Control_P.AlphaBetaZerotodq0_Alignment_c ==
       Voltage_Control_P.CompareToConstant1_const_i);
  }

  for (uData_tmp = 0; uData_tmp < 3; uData_tmp++) {
    /* Gain: '<S61>/Gain1' incorporates:
     *  Gain: '<S61>/Gain3'
     */
    rtb_Gain1[uData_tmp] = Voltage_Control_P.Gain1_Gain_k *
      (Voltage_Control_P.Gain3_Gain_m[uData_tmp + 6] * Voltage_Control_B.Kv1[2]
       + (Voltage_Control_P.Gain3_Gain_m[uData_tmp + 3] * Voltage_Control_B.Kv1
          [1] + Voltage_Control_P.Gain3_Gain_m[uData_tmp] *
          Voltage_Control_B.Kv1[0]));
  }

  /* Outputs for Enabled SubSystem: '<S60>/Subsystem - pi//2 delay' */
  Voltage_Contr_Subsystempi2delay(Voltage_Control_M, Voltage_Control_B.Compare_n,
    &rtb_Gain1[0], Voltage_Control_B.Integrator,
    &Voltage_Control_B.Subsystempi2delay, &Voltage_Control_DW.Subsystempi2delay);

  /* End of Outputs for SubSystem: '<S60>/Subsystem - pi//2 delay' */

  /* Outputs for Enabled SubSystem: '<S60>/Subsystem1' */
  Voltage_Control_Subsystem1(Voltage_Control_M, Voltage_Control_B.Compare_o,
    &rtb_Gain1[0], Voltage_Control_B.Integrator, &Voltage_Control_B.Subsystem1,
    &Voltage_Control_DW.Subsystem1);

  /* End of Outputs for SubSystem: '<S60>/Subsystem1' */

  /* Switch: '<S60>/Switch' */
  if (Voltage_Control_B.Compare_o != 0) {
    Voltage_Control_B.Switch_o[0] = Voltage_Control_B.Subsystem1.Fcn;
    Voltage_Control_B.Switch_o[1] = Voltage_Control_B.Subsystem1.Fcn1;
  } else {
    Voltage_Control_B.Switch_o[0] = Voltage_Control_B.Subsystempi2delay.Fcn;
    Voltage_Control_B.Switch_o[1] = Voltage_Control_B.Subsystempi2delay.Fcn1;
  }

  /* End of Switch: '<S60>/Switch' */

  /* Gain: '<S4>/Gain1' */
  Voltage_Control_B.Gain1_a = Voltage_Control_P.Gain1_Gain_i * rtb_Add;
  if (rtmIsMajorTimeStep(Voltage_Control_M) &&
      Voltage_Control_M->Timing.TaskCounters.TID[1] == 0) {
    /* Gain: '<S4>/Gain3' incorporates:
     *  Gain: '<S21>/Kv'
     *  Gain: '<S21>/Kv1'
     *  Gain: '<S25>/do not delete this gain'
     *  Gain: '<S26>/do not delete this gain'
     *  Gain: '<S27>/do not delete this gain'
     *  Gain: '<S28>/do not delete this gain'
     *  Gain: '<S29>/do not delete this gain'
     *  Gain: '<S30>/do not delete this gain'
     *  Product: '<S20>/Product1'
     *  Sum: '<S20>/Sum4'
     */
    Voltage_Control_B.Gain3_a = ((Voltage_Control_P.donotdeletethisgain_Gain_dg *
      Voltage_Control_B.StateSpace[0] * Voltage_Control_P.Kv1_Gain_d *
      (Voltage_Control_P.donotdeletethisgain_Gain_p *
       Voltage_Control_B.StateSpace[12] * Voltage_Control_P.Kv_Gain_j) +
      Voltage_Control_P.donotdeletethisgain_Gain_m *
      Voltage_Control_B.StateSpace[1] * Voltage_Control_P.Kv1_Gain_d *
      (Voltage_Control_P.donotdeletethisgain_Gain_nw *
       Voltage_Control_B.StateSpace[13] * Voltage_Control_P.Kv_Gain_j)) +
      Voltage_Control_P.donotdeletethisgain_Gain_j *
      Voltage_Control_B.StateSpace[2] * Voltage_Control_P.Kv1_Gain_d *
      (Voltage_Control_P.donotdeletethisgain_Gain_bf *
       Voltage_Control_B.StateSpace[14] * Voltage_Control_P.Kv_Gain_j)) *
      Voltage_Control_P.Gain3_Gain_lb;
  }
}

/* Model update function */
void Voltage_Control_update(void)
{
  /* Update for Integrator: '<S37>/Integrator' */
  Voltage_Control_DW.Integrator_IWORK = 0;
  if (rtmIsMajorTimeStep(Voltage_Control_M) &&
      Voltage_Control_M->Timing.TaskCounters.TID[1] == 0) {
    /* Update for UnitDelay: '<S67>/Unit Delay' */
    Voltage_Control_DW.UnitDelay_DSTATE = Voltage_Control_B.Divide;

    /* Update for UnitDelay: '<S68>/Unit Delay' */
    Voltage_Control_DW.UnitDelay_DSTATE_h = Voltage_Control_B.Divide_b;

    /* Update for UnitDelay: '<S69>/Unit Delay' */
    Voltage_Control_DW.UnitDelay_DSTATE_i = Voltage_Control_B.Divide_n;

    /* Update for S-Function (sfun_spssw_discc): '<S139>/State-Space' */

    /* S-Function block: <S139>/State-Space */
    {
      const real_T *As = (real_T*)Voltage_Control_DW.StateSpace_PWORK.AS;
      const real_T *Bs = (real_T*)Voltage_Control_DW.StateSpace_PWORK.BS;
      real_T *xtmp = (real_T*)Voltage_Control_DW.StateSpace_PWORK.XTMP;
      real_T accum;

      /* Calculate new states... */
      {
        int_T i1;
        real_T *xd = &Voltage_Control_DW.StateSpace_DSTATE[0];
        for (i1=0; i1 < 6; i1++) {
          accum = 0.0;

          {
            int_T i2;
            real_T *xd = &Voltage_Control_DW.StateSpace_DSTATE[0];
            for (i2=0; i2 < 6; i2++) {
              accum += *(As++) * xd[i2];
            }
          }

          accum += *(Bs++) * Voltage_Control_B.Product[0];
          accum += *(Bs++) * Voltage_Control_B.Product[1];
          accum += *(Bs++) * Voltage_Control_B.Product[2];
          accum += *(Bs++) * Voltage_Control_B.Gain3[0];
          accum += *(Bs++) * Voltage_Control_B.Gain3[1];
          accum += *(Bs++) * Voltage_Control_B.Gain3[2];
          accum += *(Bs++) * Voltage_Control_B.UnitDelay;
          accum += *(Bs++) * Voltage_Control_B.UnitDelay_b;
          accum += *(Bs++) * Voltage_Control_B.UnitDelay_i;
          xtmp[i1] = accum;
        }
      }

      {
        int_T i1;
        real_T *xd = &Voltage_Control_DW.StateSpace_DSTATE[0];
        for (i1=0; i1 < 6; i1++) {
          xd[i1] = xtmp[i1];
        }
      }
    }

    /* Update for Memory: '<S37>/Memory' */
    Voltage_Control_DW.Memory_PreviousInput = Voltage_Control_B.y;
  }

  if (rtmIsMajorTimeStep(Voltage_Control_M) &&
      Voltage_Control_M->Timing.TaskCounters.TID[2] == 0) {
    /* Update for UnitDelay: '<Root>/Unit Delay' incorporates:
     *  Inport: '<Root>/action'
     */
    Voltage_Control_DW.UnitDelay_DSTATE_n[0] = Voltage_Control_U.action[0];
    Voltage_Control_DW.UnitDelay_DSTATE_n[1] = Voltage_Control_U.action[1];
  }

  /* Update for Enabled SubSystem: '<S37>/Automatic Gain Control' incorporates:
   *  EnablePort: '<S38>/Enable'
   */
  if (Voltage_Control_DW.AutomaticGainControl_MODE) {
    /* Update for VariableTransportDelay: '<S47>/Variable Transport Delay' */
    {
      real_T **uBuffer = (real_T**)
        &Voltage_Control_DW.VariableTransportDelay_PWORK_p.TUbufferPtrs[0];
      real_T **tBuffer = (real_T**)
        &Voltage_Control_DW.VariableTransportDelay_PWORK_p.TUbufferPtrs[1];
      real_T **xBuffer = (real_T**)
        &Voltage_Control_DW.VariableTransportDelay_PWORK_p.TUbufferPtrs[2];
      real_T simTime = Voltage_Control_M->Timing.t[0];
      Voltage_Control_DW.VariableTransportDelay_IWORK_l.Head =
        ((Voltage_Control_DW.VariableTransportDelay_IWORK_l.Head <
          (Voltage_Control_DW.VariableTransportDelay_IWORK_l.CircularBufSize-1))
         ? (Voltage_Control_DW.VariableTransportDelay_IWORK_l.Head+1) : 0);
      if (Voltage_Control_DW.VariableTransportDelay_IWORK_l.Head ==
          Voltage_Control_DW.VariableTransportDelay_IWORK_l.Tail) {
        Voltage_Control_DW.VariableTransportDelay_IWORK_l.Tail =
          ((Voltage_Control_DW.VariableTransportDelay_IWORK_l.Tail <
            (Voltage_Control_DW.VariableTransportDelay_IWORK_l.CircularBufSize-1))
           ? (Voltage_Control_DW.VariableTransportDelay_IWORK_l.Tail+1) : 0);
      }

      (*tBuffer)[Voltage_Control_DW.VariableTransportDelay_IWORK_l.Head] =
        simTime;
      (*uBuffer)[Voltage_Control_DW.VariableTransportDelay_IWORK_l.Head] =
        Voltage_Control_B.integrator_p;
      (*xBuffer)[Voltage_Control_DW.VariableTransportDelay_IWORK_l.Head] =
        Voltage_Control_X.VariableTransportDelay_CSTATE_n;

      /* when use fixed buffer, reset solver at when buffer is updated
       * to avoid output consistency fail.
       */
    }

    if (rtmIsMajorTimeStep(Voltage_Control_M) &&
        Voltage_Control_M->Timing.TaskCounters.TID[1] == 0) {
      /* Update for Memory: '<S47>/Memory' */
      Voltage_Control_DW.Memory_PreviousInput_d = Voltage_Control_B.Switch_j;
    }

    /* Update for VariableTransportDelay: '<S48>/Variable Transport Delay' */
    {
      real_T **uBuffer = (real_T**)
        &Voltage_Control_DW.VariableTransportDelay_PWORK_f.TUbufferPtrs[0];
      real_T **tBuffer = (real_T**)
        &Voltage_Control_DW.VariableTransportDelay_PWORK_f.TUbufferPtrs[1];
      real_T **xBuffer = (real_T**)
        &Voltage_Control_DW.VariableTransportDelay_PWORK_f.TUbufferPtrs[2];
      real_T simTime = Voltage_Control_M->Timing.t[0];
      Voltage_Control_DW.VariableTransportDelay_IWORK_g.Head =
        ((Voltage_Control_DW.VariableTransportDelay_IWORK_g.Head <
          (Voltage_Control_DW.VariableTransportDelay_IWORK_g.CircularBufSize-1))
         ? (Voltage_Control_DW.VariableTransportDelay_IWORK_g.Head+1) : 0);
      if (Voltage_Control_DW.VariableTransportDelay_IWORK_g.Head ==
          Voltage_Control_DW.VariableTransportDelay_IWORK_g.Tail) {
        Voltage_Control_DW.VariableTransportDelay_IWORK_g.Tail =
          ((Voltage_Control_DW.VariableTransportDelay_IWORK_g.Tail <
            (Voltage_Control_DW.VariableTransportDelay_IWORK_g.CircularBufSize-1))
           ? (Voltage_Control_DW.VariableTransportDelay_IWORK_g.Tail+1) : 0);
      }

      (*tBuffer)[Voltage_Control_DW.VariableTransportDelay_IWORK_g.Head] =
        simTime;
      (*uBuffer)[Voltage_Control_DW.VariableTransportDelay_IWORK_g.Head] =
        Voltage_Control_B.integrator_m;
      (*xBuffer)[Voltage_Control_DW.VariableTransportDelay_IWORK_g.Head] =
        Voltage_Control_X.VariableTransportDelay_CSTAT_nz;

      /* when use fixed buffer, reset solver at when buffer is updated
       * to avoid output consistency fail.
       */
    }

    if (rtmIsMajorTimeStep(Voltage_Control_M) &&
        Voltage_Control_M->Timing.TaskCounters.TID[1] == 0) {
      /* Update for Memory: '<S48>/Memory' */
      Voltage_Control_DW.Memory_PreviousInput_o = Voltage_Control_B.Switch_i;
    }
  }

  /* End of Update for SubSystem: '<S37>/Automatic Gain Control' */

  /* Update for VariableTransportDelay: '<S59>/Variable Transport Delay' */
  {
    real_T **uBuffer = (real_T**)
      &Voltage_Control_DW.VariableTransportDelay_PWORK.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)
      &Voltage_Control_DW.VariableTransportDelay_PWORK.TUbufferPtrs[1];
    real_T **xBuffer = (real_T**)
      &Voltage_Control_DW.VariableTransportDelay_PWORK.TUbufferPtrs[2];
    real_T simTime = Voltage_Control_M->Timing.t[0];
    Voltage_Control_DW.VariableTransportDelay_IWORK.Head =
      ((Voltage_Control_DW.VariableTransportDelay_IWORK.Head <
        (Voltage_Control_DW.VariableTransportDelay_IWORK.CircularBufSize-1)) ?
       (Voltage_Control_DW.VariableTransportDelay_IWORK.Head+1) : 0);
    if (Voltage_Control_DW.VariableTransportDelay_IWORK.Head ==
        Voltage_Control_DW.VariableTransportDelay_IWORK.Tail) {
      Voltage_Control_DW.VariableTransportDelay_IWORK.Tail =
        ((Voltage_Control_DW.VariableTransportDelay_IWORK.Tail <
          (Voltage_Control_DW.VariableTransportDelay_IWORK.CircularBufSize-1)) ?
         (Voltage_Control_DW.VariableTransportDelay_IWORK.Tail+1) : 0);
    }

    (*tBuffer)[Voltage_Control_DW.VariableTransportDelay_IWORK.Head] = simTime;
    (*uBuffer)[Voltage_Control_DW.VariableTransportDelay_IWORK.Head] =
      Voltage_Control_B.integrator;
    (*xBuffer)[Voltage_Control_DW.VariableTransportDelay_IWORK.Head] =
      Voltage_Control_X.VariableTransportDelay_CSTATE;

    /* when use fixed buffer, reset solver at when buffer is updated
     * to avoid output consistency fail.
     */
  }

  if (rtmIsMajorTimeStep(Voltage_Control_M) &&
      Voltage_Control_M->Timing.TaskCounters.TID[1] == 0) {
    /* Update for Memory: '<S59>/Memory' */
    Voltage_Control_DW.Memory_PreviousInput_f = Voltage_Control_B.Switch;
  }

  /* Update for RateLimiter: '<S37>/Rate Limiter' */
  Voltage_Control_DW.PrevY = Voltage_Control_B.RateLimiter;
  Voltage_Control_DW.LastMajorTime = Voltage_Control_M->Timing.t[0];

  /* ContTimeOutputInconsistentWithStateAtMajorOutputFlag is set, need to run a minor output */
  if (rtmIsMajorTimeStep(Voltage_Control_M)) {
    if (rtsiGetContTimeOutputInconsistentWithStateAtMajorStep
        (&Voltage_Control_M->solverInfo)) {
      rtsiSetSimTimeStep(&Voltage_Control_M->solverInfo,MINOR_TIME_STEP);
      rtsiSetContTimeOutputInconsistentWithStateAtMajorStep
        (&Voltage_Control_M->solverInfo, false);
      Voltage_Control_output();
      rtsiSetSimTimeStep(&Voltage_Control_M->solverInfo, MAJOR_TIME_STEP);
    }
  }

  if (rtmIsMajorTimeStep(Voltage_Control_M)) {
    rt_ertODEUpdateContinuousStates(&Voltage_Control_M->solverInfo);
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++Voltage_Control_M->Timing.clockTick0)) {
    ++Voltage_Control_M->Timing.clockTickH0;
  }

  Voltage_Control_M->Timing.t[0] = rtsiGetSolverStopTime
    (&Voltage_Control_M->solverInfo);

  {
    /* Update absolute timer for sample time: [1.0E-5s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The resolution of this integer timer is 1.0E-5, which is the step size
     * of the task. Size of "clockTick1" ensures timer will not overflow during the
     * application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    Voltage_Control_M->Timing.clockTick1++;
    if (!Voltage_Control_M->Timing.clockTick1) {
      Voltage_Control_M->Timing.clockTickH1++;
    }
  }

  rate_scheduler();
}

/* Derivatives for root system: '<Root>' */
void Voltage_Control_derivatives(void)
{
  boolean_T lsat;
  boolean_T usat;
  XDot_Voltage_Control_T *_rtXdot;
  _rtXdot = ((XDot_Voltage_Control_T *) Voltage_Control_M->derivs);

  /* Derivatives for Integrator: '<S4>/Integrator' */
  _rtXdot->Integrator_CSTATE = Voltage_Control_B.Gain1_a;

  /* Derivatives for StateSpace: '<S4>/State-Space1' */
  _rtXdot->StateSpace1_CSTATE[0] = 0.0;
  _rtXdot->StateSpace1_CSTATE[1] = 0.0;
  _rtXdot->StateSpace1_CSTATE[2] = 0.0;
  _rtXdot->StateSpace1_CSTATE[0] += Voltage_Control_P.StateSpace1_A[0] *
    Voltage_Control_X.StateSpace1_CSTATE[1];
  _rtXdot->StateSpace1_CSTATE[1] += Voltage_Control_P.StateSpace1_A[1] *
    Voltage_Control_X.StateSpace1_CSTATE[2];
  _rtXdot->StateSpace1_CSTATE[2] += Voltage_Control_P.StateSpace1_A[2] *
    Voltage_Control_X.StateSpace1_CSTATE[0];
  _rtXdot->StateSpace1_CSTATE[2] += Voltage_Control_P.StateSpace1_A[3] *
    Voltage_Control_X.StateSpace1_CSTATE[1];
  _rtXdot->StateSpace1_CSTATE[2] += Voltage_Control_P.StateSpace1_A[4] *
    Voltage_Control_X.StateSpace1_CSTATE[2];
  _rtXdot->StateSpace1_CSTATE[2] += Voltage_Control_P.StateSpace1_B *
    Voltage_Control_B.Gain3_a;

  /* Derivatives for Integrator: '<S37>/Integrator' */
  lsat = (Voltage_Control_X.Integrator_CSTATE_f <=
          Voltage_Control_P.Integrator_LowerSat);
  usat = (Voltage_Control_X.Integrator_CSTATE_f >=
          Voltage_Control_P.Integrator_UpperSat);
  if (((!lsat) && (!usat)) || (lsat && (Voltage_Control_B.Saturation2 > 0.0)) ||
      (usat && (Voltage_Control_B.Saturation2 < 0.0))) {
    _rtXdot->Integrator_CSTATE_f = Voltage_Control_B.Saturation2;
  } else {
    /* in saturation */
    _rtXdot->Integrator_CSTATE_f = 0.0;
  }

  /* End of Derivatives for Integrator: '<S37>/Integrator' */

  /* Derivatives for Enabled SubSystem: '<S37>/Automatic Gain Control' */
  if (Voltage_Control_DW.AutomaticGainControl_MODE) {
    /* Derivatives for VariableTransportDelay: '<S47>/Variable Transport Delay' */
    {
      real_T instantDelay;
      instantDelay = Voltage_Control_B.period_p;
      if (instantDelay > (Voltage_Control_P.VariableTransportDelay_MaxDelay)) {
        instantDelay = (Voltage_Control_P.VariableTransportDelay_MaxDelay);
      }

      if (instantDelay < 0.0) {
        ((XDot_Voltage_Control_T *) Voltage_Control_M->derivs)
          ->VariableTransportDelay_CSTATE_n = 0;
      } else {
        ((XDot_Voltage_Control_T *) Voltage_Control_M->derivs)
          ->VariableTransportDelay_CSTATE_n = 1.0/instantDelay;
      }
    }

    /* Derivatives for Integrator: '<S47>/integrator' */
    _rtXdot->integrator_CSTATE_m = Voltage_Control_B.Switch_ow[0];

    /* Derivatives for VariableTransportDelay: '<S48>/Variable Transport Delay' */
    {
      real_T instantDelay;
      instantDelay = Voltage_Control_B.period_a;
      if (instantDelay > (Voltage_Control_P.VariableTransportDelay_MaxDel_e)) {
        instantDelay = (Voltage_Control_P.VariableTransportDelay_MaxDel_e);
      }

      if (instantDelay < 0.0) {
        ((XDot_Voltage_Control_T *) Voltage_Control_M->derivs)
          ->VariableTransportDelay_CSTAT_nz = 0;
      } else {
        ((XDot_Voltage_Control_T *) Voltage_Control_M->derivs)
          ->VariableTransportDelay_CSTAT_nz = 1.0/instantDelay;
      }
    }

    /* Derivatives for Integrator: '<S48>/integrator' */
    _rtXdot->integrator_CSTATE_l = Voltage_Control_B.Switch_ow[1];
  } else {
    {
      real_T *dx;
      int_T i;
      dx = &(((XDot_Voltage_Control_T *) Voltage_Control_M->derivs)
             ->VariableTransportDelay_CSTATE_n);
      for (i=0; i < 4; i++) {
        dx[i] = 0.0;
      }
    }
  }

  /* End of Derivatives for SubSystem: '<S37>/Automatic Gain Control' */

  /* Derivatives for Integrator: '<S39>/Integrator' */
  lsat = (Voltage_Control_X.Integrator_CSTATE_fe <=
          Voltage_Control_P.Integrator_LowerSat_n);
  usat = (Voltage_Control_X.Integrator_CSTATE_fe >=
          Voltage_Control_P.Integrator_UpperSat_a);
  if (((!lsat) && (!usat)) || (lsat && (Voltage_Control_B.Kp5 > 0.0)) || (usat &&
       (Voltage_Control_B.Kp5 < 0.0))) {
    _rtXdot->Integrator_CSTATE_fe = Voltage_Control_B.Kp5;
  } else {
    /* in saturation */
    _rtXdot->Integrator_CSTATE_fe = 0.0;
  }

  /* End of Derivatives for Integrator: '<S39>/Integrator' */

  /* Derivatives for VariableTransportDelay: '<S59>/Variable Transport Delay' */
  {
    real_T instantDelay;
    instantDelay = Voltage_Control_B.period;
    if (instantDelay > (Voltage_Control_P.VariableTransportDelay_MaxDel_n)) {
      instantDelay = (Voltage_Control_P.VariableTransportDelay_MaxDel_n);
    }

    if (instantDelay < 0.0) {
      ((XDot_Voltage_Control_T *) Voltage_Control_M->derivs)
        ->VariableTransportDelay_CSTATE = 0;
    } else {
      ((XDot_Voltage_Control_T *) Voltage_Control_M->derivs)
        ->VariableTransportDelay_CSTATE = 1.0/instantDelay;
    }
  }

  /* Derivatives for Integrator: '<S59>/integrator' */
  _rtXdot->integrator_CSTATE = Voltage_Control_B.Switch_o[1];

  /* Derivatives for TransferFcn: '<S39>/Transfer Fcn' */
  _rtXdot->TransferFcn_CSTATE = 0.0;
  _rtXdot->TransferFcn_CSTATE += Voltage_Control_P.TransferFcn_A *
    Voltage_Control_X.TransferFcn_CSTATE;
  _rtXdot->TransferFcn_CSTATE += Voltage_Control_B.Kp6;

  /* Derivatives for Integrator: '<S55>/Integrator_x1' */
  _rtXdot->Integrator_x1_CSTATE = Voltage_Control_B.x1;

  /* Derivatives for Integrator: '<S55>/Integrator_x2' */
  _rtXdot->Integrator_x2_CSTATE = Voltage_Control_B.x2;
}

/* Model initialize function */
void Voltage_Control_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  Voltage_Control_P.Saturation_UpperSat = rtInf;
  Voltage_Control_P.Integrator_UpperSat = rtInf;
  Voltage_Control_P.Integrator_LowerSat = rtMinusInf;
  Voltage_Control_P.Saturation_UpperSat_e = rtInf;
  Voltage_Control_P.Saturation_UpperSat_p = rtInf;
  Voltage_Control_P.Saturation_UpperSat_f = rtInf;
  Voltage_Control_P.Integrator_UpperSat_a = rtInf;
  Voltage_Control_P.Saturation2_UpperSat = rtInf;

  /* initialize real-time model */
  (void) memset((void *)Voltage_Control_M, 0,
                sizeof(RT_MODEL_Voltage_Control_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&Voltage_Control_M->solverInfo,
                          &Voltage_Control_M->Timing.simTimeStep);
    rtsiSetTPtr(&Voltage_Control_M->solverInfo, &rtmGetTPtr(Voltage_Control_M));
    rtsiSetStepSizePtr(&Voltage_Control_M->solverInfo,
                       &Voltage_Control_M->Timing.stepSize0);
    rtsiSetdXPtr(&Voltage_Control_M->solverInfo, &Voltage_Control_M->derivs);
    rtsiSetContStatesPtr(&Voltage_Control_M->solverInfo, (real_T **)
                         &Voltage_Control_M->contStates);
    rtsiSetNumContStatesPtr(&Voltage_Control_M->solverInfo,
      &Voltage_Control_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&Voltage_Control_M->solverInfo,
      &Voltage_Control_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&Voltage_Control_M->solverInfo,
      &Voltage_Control_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&Voltage_Control_M->solverInfo,
      &Voltage_Control_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&Voltage_Control_M->solverInfo, (&rtmGetErrorStatus
      (Voltage_Control_M)));
    rtsiSetRTModelPtr(&Voltage_Control_M->solverInfo, Voltage_Control_M);
  }

  rtsiSetSimTimeStep(&Voltage_Control_M->solverInfo, MAJOR_TIME_STEP);
  Voltage_Control_M->intgData.y = Voltage_Control_M->odeY;
  Voltage_Control_M->intgData.f[0] = Voltage_Control_M->odeF[0];
  Voltage_Control_M->intgData.f[1] = Voltage_Control_M->odeF[1];
  Voltage_Control_M->intgData.f[2] = Voltage_Control_M->odeF[2];
  Voltage_Control_M->intgData.f[3] = Voltage_Control_M->odeF[3];
  Voltage_Control_M->contStates = ((X_Voltage_Control_T *) &Voltage_Control_X);
  Voltage_Control_M->periodicContStateIndices = ((int_T*)
    Voltage_Control_PeriodicIndX);
  Voltage_Control_M->periodicContStateRanges = ((real_T*)
    Voltage_Control_PeriodicRngX);
  rtsiSetSolverData(&Voltage_Control_M->solverInfo, (void *)
                    &Voltage_Control_M->intgData);
  rtsiSetSolverName(&Voltage_Control_M->solverInfo,"ode4");
  rtmSetTPtr(Voltage_Control_M, &Voltage_Control_M->Timing.tArray[0]);
  Voltage_Control_M->Timing.stepSize0 = 1.0E-5;
  rtmSetFirstInitCond(Voltage_Control_M, 1);

  /* block I/O */
  (void) memset(((void *) &Voltage_Control_B), 0,
                sizeof(B_Voltage_Control_T));

  /* states (continuous) */
  {
    (void) memset((void *)&Voltage_Control_X, 0,
                  sizeof(X_Voltage_Control_T));
  }

  /* Periodic continuous states */
  {
    (void) memset((void*) Voltage_Control_PeriodicIndX, 0,
                  1*sizeof(int_T));
    (void) memset((void*) Voltage_Control_PeriodicRngX, 0,
                  2*sizeof(real_T));
  }

  /* states (dwork) */
  (void) memset((void *)&Voltage_Control_DW, 0,
                sizeof(DW_Voltage_Control_T));

  /* external inputs */
  (void)memset(&Voltage_Control_U, 0, sizeof(ExtU_Voltage_Control_T));

  /* external outputs */
  (void) memset((void *)&Voltage_Control_Y, 0,
                sizeof(ExtY_Voltage_Control_T));

  /* Start for InitialCondition: '<S37>/Initial' */
  Voltage_Control_B.Initial = Voltage_Control_P.Initial_Value;
  Voltage_Control_DW.Initial_FirstOutputTime = (rtMinusInf);

  /* Start for S-Function (sfun_spssw_discc): '<S139>/State-Space' */

  /* S-Function block: <S139>/State-Space */
  {
    Voltage_Control_DW.StateSpace_PWORK.AS = (real_T*)calloc(6 * 6, sizeof
      (real_T));
    Voltage_Control_DW.StateSpace_PWORK.BS = (real_T*)calloc(6 * 9, sizeof
      (real_T));
    Voltage_Control_DW.StateSpace_PWORK.CS = (real_T*)calloc(24 * 6, sizeof
      (real_T));
    Voltage_Control_DW.StateSpace_PWORK.DS = (real_T*)calloc(24 * 9, sizeof
      (real_T));
    Voltage_Control_DW.StateSpace_PWORK.DX_COL = (real_T*)calloc(24, sizeof
      (real_T));
    Voltage_Control_DW.StateSpace_PWORK.TMP2 = (real_T*)calloc(9, sizeof(real_T));
    Voltage_Control_DW.StateSpace_PWORK.BD_COL = (real_T*)calloc(6, sizeof
      (real_T));
    Voltage_Control_DW.StateSpace_PWORK.TMP1 = (real_T*)calloc(6, sizeof(real_T));
    Voltage_Control_DW.StateSpace_PWORK.XTMP = (real_T*)calloc(6, sizeof(real_T));
  }

  /* Start for VariableTransportDelay: '<S59>/Variable Transport Delay' */
  {
    real_T *pBuffer =
      &Voltage_Control_DW.VariableTransportDelay_RWORK.TUbufferArea[0];
    int_T j;
    Voltage_Control_DW.VariableTransportDelay_IWORK.Tail = 0;
    Voltage_Control_DW.VariableTransportDelay_IWORK.Head = 0;
    Voltage_Control_DW.VariableTransportDelay_IWORK.Last = 0;
    Voltage_Control_DW.VariableTransportDelay_IWORK.CircularBufSize = 8192;
    for (j=0; j < 8192; j++) {
      pBuffer[j] = Voltage_Control_P.VariableTransportDelay_InitOu_h;
      pBuffer[8192 + j] = Voltage_Control_M->Timing.t[0];
    }

    pBuffer[2*8192] = 0.0;
    Voltage_Control_DW.VariableTransportDelay_PWORK.TUbufferPtrs[0] = (void *)
      &pBuffer[0];
    Voltage_Control_DW.VariableTransportDelay_PWORK.TUbufferPtrs[1] = (void *)
      &pBuffer[8192];
    Voltage_Control_DW.VariableTransportDelay_PWORK.TUbufferPtrs[2] = (void *)
      &pBuffer[2*8192];
  }

  Voltage_Control_PrevZCX.Integrator_Reset_ZCE = UNINITIALIZED_ZCSIG;

  /* InitializeConditions for Integrator: '<S4>/Integrator' */
  Voltage_Control_X.Integrator_CSTATE = Voltage_Control_P.Integrator_IC;

  /* InitializeConditions for StateSpace: '<S4>/State-Space1' */
  Voltage_Control_X.StateSpace1_CSTATE[0] =
    Voltage_Control_P.StateSpace1_InitialCondition[0];
  Voltage_Control_X.StateSpace1_CSTATE[1] =
    Voltage_Control_P.StateSpace1_InitialCondition[1];
  Voltage_Control_X.StateSpace1_CSTATE[2] =
    Voltage_Control_P.StateSpace1_InitialCondition[2];

  /* InitializeConditions for Integrator: '<S37>/Integrator' */
  if (rtmIsFirstInitCond(Voltage_Control_M)) {
    Voltage_Control_X.Integrator_CSTATE_f = 0.0;
  }

  Voltage_Control_DW.Integrator_IWORK = 1;

  /* End of InitializeConditions for Integrator: '<S37>/Integrator' */

  /* InitializeConditions for UnitDelay: '<S67>/Unit Delay' */
  Voltage_Control_DW.UnitDelay_DSTATE =
    Voltage_Control_P.UnitDelay_InitialCondition;

  /* InitializeConditions for UnitDelay: '<S68>/Unit Delay' */
  Voltage_Control_DW.UnitDelay_DSTATE_h =
    Voltage_Control_P.UnitDelay_InitialCondition_e;

  /* InitializeConditions for UnitDelay: '<S69>/Unit Delay' */
  Voltage_Control_DW.UnitDelay_DSTATE_i =
    Voltage_Control_P.UnitDelay_InitialCondition_j;

  /* InitializeConditions for S-Function (sfun_spssw_discc): '<S139>/State-Space' */
  {
    int32_T i, j;
    real_T *As = (real_T*)Voltage_Control_DW.StateSpace_PWORK.AS;
    real_T *Bs = (real_T*)Voltage_Control_DW.StateSpace_PWORK.BS;
    real_T *Cs = (real_T*)Voltage_Control_DW.StateSpace_PWORK.CS;
    real_T *Ds = (real_T*)Voltage_Control_DW.StateSpace_PWORK.DS;
    real_T *X0 = (real_T*)&Voltage_Control_DW.StateSpace_DSTATE[0];
    for (i = 0; i < 6; i++) {
      X0[i] = (Voltage_Control_P.StateSpace_X0_param[i]);
    }

    /* Copy and transpose A and B */
    for (i=0; i<6; i++) {
      for (j=0; j<6; j++)
        As[i*6 + j] = (Voltage_Control_P.StateSpace_AS_param[i + j*6]);
      for (j=0; j<9; j++)
        Bs[i*9 + j] = (Voltage_Control_P.StateSpace_BS_param[i + j*6]);
    }

    /* Copy and transpose C */
    for (i=0; i<24; i++) {
      for (j=0; j<6; j++)
        Cs[i*6 + j] = (Voltage_Control_P.StateSpace_CS_param[i + j*24]);
    }

    /* Copy and transpose D */
    for (i=0; i<24; i++) {
      for (j=0; j<9; j++)
        Ds[i*9 + j] = (Voltage_Control_P.StateSpace_DS_param[i + j*24]);
    }
  }

  /* InitializeConditions for UnitDelay: '<Root>/Unit Delay' */
  Voltage_Control_DW.UnitDelay_DSTATE_n[0] =
    Voltage_Control_P.UnitDelay_InitialCondition_l;
  Voltage_Control_DW.UnitDelay_DSTATE_n[1] =
    Voltage_Control_P.UnitDelay_InitialCondition_l;

  /* InitializeConditions for Memory: '<S37>/Memory' */
  Voltage_Control_DW.Memory_PreviousInput =
    Voltage_Control_P.Memory_InitialCondition_o;

  /* InitializeConditions for Integrator: '<S39>/Integrator' */
  Voltage_Control_X.Integrator_CSTATE_fe = Voltage_Control_P.Continuous_Init;

  /* InitializeConditions for VariableTransportDelay: '<S59>/Variable Transport Delay' */
  Voltage_Control_X.VariableTransportDelay_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S59>/integrator' */
  Voltage_Control_X.integrator_CSTATE = Voltage_Control_P.integrator_IC_lo;

  /* InitializeConditions for Memory: '<S59>/Memory' */
  Voltage_Control_DW.Memory_PreviousInput_f =
    Voltage_Control_P.Memory_InitialCondition_a;

  /* InitializeConditions for TransferFcn: '<S39>/Transfer Fcn' */
  Voltage_Control_X.TransferFcn_CSTATE = 0.0;

  /* InitializeConditions for RateLimiter: '<S37>/Rate Limiter' */
  Voltage_Control_DW.LastMajorTime = (rtInf);

  /* InitializeConditions for Integrator: '<S55>/Integrator_x1' */
  Voltage_Control_X.Integrator_x1_CSTATE = Voltage_Control_P.Integrator_x1_IC;

  /* InitializeConditions for Integrator: '<S55>/Integrator_x2' */
  Voltage_Control_X.Integrator_x2_CSTATE = Voltage_Control_P.Integrator_x2_IC;

  /* SystemInitialize for Enabled SubSystem: '<S134>/Subsystem1' */
  /* SystemInitialize for Outport: '<S138>/alpha_beta' */
  Voltage_Control_B.Fcn = Voltage_Control_P.alpha_beta_Y0_e[0];
  Voltage_Control_B.Fcn1 = Voltage_Control_P.alpha_beta_Y0_e[1];

  /* End of SystemInitialize for SubSystem: '<S134>/Subsystem1' */

  /* SystemInitialize for Enabled SubSystem: '<S134>/Subsystem - pi//2 delay' */
  /* SystemInitialize for Outport: '<S137>/alpha_beta' */
  Voltage_Control_B.Fcn_b = Voltage_Control_P.alpha_beta_Y0[0];
  Voltage_Control_B.Fcn1_b = Voltage_Control_P.alpha_beta_Y0[1];

  /* End of SystemInitialize for SubSystem: '<S134>/Subsystem - pi//2 delay' */

  /* SystemInitialize for Enabled SubSystem: '<S121>/Subsystem1' */
  Voltage_Control_Subsystem1_Init(&Voltage_Control_B.Subsystem1_e,
    &Voltage_Control_P.Subsystem1_e);

  /* End of SystemInitialize for SubSystem: '<S121>/Subsystem1' */

  /* SystemInitialize for Enabled SubSystem: '<S121>/Subsystem - pi//2 delay' */
  Voltage__Subsystempi2delay_Init(&Voltage_Control_B.Subsystempi2delay_b,
    &Voltage_Control_P.Subsystempi2delay_b);

  /* End of SystemInitialize for SubSystem: '<S121>/Subsystem - pi//2 delay' */

  /* SystemInitialize for Enabled SubSystem: '<S127>/Subsystem1' */
  Voltage_Control_Subsystem1_Init(&Voltage_Control_B.Subsystem1_ex,
    &Voltage_Control_P.Subsystem1_ex);

  /* End of SystemInitialize for SubSystem: '<S127>/Subsystem1' */

  /* SystemInitialize for Enabled SubSystem: '<S127>/Subsystem - pi//2 delay' */
  Voltage__Subsystempi2delay_Init(&Voltage_Control_B.Subsystempi2delay_n,
    &Voltage_Control_P.Subsystempi2delay_n);

  /* End of SystemInitialize for SubSystem: '<S127>/Subsystem - pi//2 delay' */

  /* SystemInitialize for MATLAB Function: '<Root>/MATLAB Function1' */
  memset(&Voltage_Control_DW.yData[0], 0, 125U * sizeof(real_T));
  memset(&Voltage_Control_DW.uData[0], 0, 50U * sizeof(real_T));

  /* SystemInitialize for Enabled SubSystem: '<S37>/Automatic Gain Control' */
  /* Start for VariableTransportDelay: '<S47>/Variable Transport Delay' */
  {
    real_T *pBuffer =
      &Voltage_Control_DW.VariableTransportDelay_RWORK_k.TUbufferArea[0];
    int_T j;
    Voltage_Control_DW.VariableTransportDelay_IWORK_l.Tail = 0;
    Voltage_Control_DW.VariableTransportDelay_IWORK_l.Head = 0;
    Voltage_Control_DW.VariableTransportDelay_IWORK_l.Last = 0;
    Voltage_Control_DW.VariableTransportDelay_IWORK_l.CircularBufSize = 8192;
    for (j=0; j < 8192; j++) {
      pBuffer[j] = Voltage_Control_P.VariableTransportDelay_InitOutp;
      pBuffer[8192 + j] = Voltage_Control_M->Timing.t[0];
    }

    pBuffer[2*8192] = 0.0;
    Voltage_Control_DW.VariableTransportDelay_PWORK_p.TUbufferPtrs[0] = (void *)
      &pBuffer[0];
    Voltage_Control_DW.VariableTransportDelay_PWORK_p.TUbufferPtrs[1] = (void *)
      &pBuffer[8192];
    Voltage_Control_DW.VariableTransportDelay_PWORK_p.TUbufferPtrs[2] = (void *)
      &pBuffer[2*8192];
  }

  /* Start for VariableTransportDelay: '<S48>/Variable Transport Delay' */
  {
    real_T *pBuffer =
      &Voltage_Control_DW.VariableTransportDelay_RWORK_c.TUbufferArea[0];
    int_T j;
    Voltage_Control_DW.VariableTransportDelay_IWORK_g.Tail = 0;
    Voltage_Control_DW.VariableTransportDelay_IWORK_g.Head = 0;
    Voltage_Control_DW.VariableTransportDelay_IWORK_g.Last = 0;
    Voltage_Control_DW.VariableTransportDelay_IWORK_g.CircularBufSize = 8192;
    for (j=0; j < 8192; j++) {
      pBuffer[j] = Voltage_Control_P.VariableTransportDelay_InitOu_n;
      pBuffer[8192 + j] = Voltage_Control_M->Timing.t[0];
    }

    pBuffer[2*8192] = 0.0;
    Voltage_Control_DW.VariableTransportDelay_PWORK_f.TUbufferPtrs[0] = (void *)
      &pBuffer[0];
    Voltage_Control_DW.VariableTransportDelay_PWORK_f.TUbufferPtrs[1] = (void *)
      &pBuffer[8192];
    Voltage_Control_DW.VariableTransportDelay_PWORK_f.TUbufferPtrs[2] = (void *)
      &pBuffer[2*8192];
  }

  /* InitializeConditions for VariableTransportDelay: '<S47>/Variable Transport Delay' */
  Voltage_Control_X.VariableTransportDelay_CSTATE_n = 0.0;

  /* InitializeConditions for Integrator: '<S47>/integrator' */
  Voltage_Control_X.integrator_CSTATE_m = Voltage_Control_P.integrator_IC;

  /* InitializeConditions for Memory: '<S47>/Memory' */
  Voltage_Control_DW.Memory_PreviousInput_d =
    Voltage_Control_P.Memory_InitialCondition;

  /* InitializeConditions for VariableTransportDelay: '<S48>/Variable Transport Delay' */
  Voltage_Control_X.VariableTransportDelay_CSTAT_nz = 0.0;

  /* InitializeConditions for Integrator: '<S48>/integrator' */
  Voltage_Control_X.integrator_CSTATE_l = Voltage_Control_P.integrator_IC_l;

  /* InitializeConditions for Memory: '<S48>/Memory' */
  Voltage_Control_DW.Memory_PreviousInput_o =
    Voltage_Control_P.Memory_InitialCondition_i;

  /* SystemInitialize for Enabled SubSystem: '<S49>/Subsystem - pi//2 delay' */
  Voltage__Subsystempi2delay_Init(&Voltage_Control_B.Subsystempi2delay_d,
    &Voltage_Control_P.Subsystempi2delay_d);

  /* End of SystemInitialize for SubSystem: '<S49>/Subsystem - pi//2 delay' */

  /* SystemInitialize for Enabled SubSystem: '<S49>/Subsystem1' */
  Voltage_Control_Subsystem1_Init(&Voltage_Control_B.Subsystem1_h,
    &Voltage_Control_P.Subsystem1_h);

  /* End of SystemInitialize for SubSystem: '<S49>/Subsystem1' */

  /* SystemInitialize for Outport: '<S38>/Gain' */
  Voltage_Control_B.MathFunction = Voltage_Control_P.Gain_Y0;

  /* End of SystemInitialize for SubSystem: '<S37>/Automatic Gain Control' */

  /* SystemInitialize for Enabled SubSystem: '<S60>/Subsystem - pi//2 delay' */
  Voltage__Subsystempi2delay_Init(&Voltage_Control_B.Subsystempi2delay,
    &Voltage_Control_P.Subsystempi2delay);

  /* End of SystemInitialize for SubSystem: '<S60>/Subsystem - pi//2 delay' */

  /* SystemInitialize for Enabled SubSystem: '<S60>/Subsystem1' */
  Voltage_Control_Subsystem1_Init(&Voltage_Control_B.Subsystem1,
    &Voltage_Control_P.Subsystem1);

  /* End of SystemInitialize for SubSystem: '<S60>/Subsystem1' */

  /* InitializeConditions for root-level periodic continuous states */
  {
    int_T rootPeriodicContStateIndices[1] = { 0 };

    real_T rootPeriodicContStateRanges[2] = { 0.0, 6.2831853071795862 };

    (void) memcpy((void*)Voltage_Control_PeriodicIndX,
                  rootPeriodicContStateIndices,
                  1*sizeof(int_T));
    (void) memcpy((void*)Voltage_Control_PeriodicRngX,
                  rootPeriodicContStateRanges,
                  2*sizeof(real_T));
  }

  /* set "at time zero" to false */
  if (rtmIsFirstInitCond(Voltage_Control_M)) {
    rtmSetFirstInitCond(Voltage_Control_M, 0);
  }
}

/* Model terminate function */
void Voltage_Control_terminate(void)
{
  /* Terminate for S-Function (sfun_spssw_discc): '<S139>/State-Space' */

  /* S-Function block: <S139>/State-Space */
  {
    /* Free memory */
    free(Voltage_Control_DW.StateSpace_PWORK.AS);
    free(Voltage_Control_DW.StateSpace_PWORK.BS);
    free(Voltage_Control_DW.StateSpace_PWORK.CS);
    free(Voltage_Control_DW.StateSpace_PWORK.DS);
    free(Voltage_Control_DW.StateSpace_PWORK.DX_COL);
    free(Voltage_Control_DW.StateSpace_PWORK.TMP2);
    free(Voltage_Control_DW.StateSpace_PWORK.BD_COL);
    free(Voltage_Control_DW.StateSpace_PWORK.TMP1);
    free(Voltage_Control_DW.StateSpace_PWORK.XTMP);
  }
}
