"""
Filename:       Voltage_Control_gyminterface.py
Written by:     Niranjan Bhujel
Date:           2022-11-30
"""


import os
import numpy as np
import gym
from gym.spaces import Box
try:
    from .Voltage_Control_pyinterface import Voltage_Control
except:
    from Voltage_Control_pyinterface import Voltage_Control


class Voltage_Control_gym(gym.Env):
    def __init__(self, obs_space_lower: np.ndarray, obs_space_upper: np.ndarray, action_space_lower: np.ndarray, action_space_upper: np.ndarray, rand_param,  rand_input , t_step: np.ndarray, reset_time: float, stop_time: float, compile: bool=True):
        """
        OpenAI gym interface to `Voltage_Control`.
        
        Parameters
        ----------
        obs_space_lower: np.ndarray
            Lower bound of observation space of shape 175
        obs_space_upper: np.ndarray
            Upper bound of observation space of shape 175
        action_space_lower: np.ndarray
            Lower bound of action space of shape 2
        action_space_upper: np.ndarray
            Upper bound of action space of shape 2
        rand_param: np.ndarray or callable
            Value of the parameter or callable (with no arguments) that generate random parameter (called at the start of episode). Shape is 3
        rand_input: scipy.stats.continuous_rv
            Value of the input or callable (with no arguments) that generate random input (called every time steps). Shape is 4
        reset_time: float
            `reset()` function reset and run the simulation upto `reset_time`
        stop_time: float
            Time upto which simulation to run
        t_step: np.ndarray
            Parameter `t_step`
        compile: bool
            Whether to compile while creating object.
        """
        super().__init__()
        
        self._rand_param = rand_param
        
        self._rand_input = rand_input

        self.reset_time = reset_time
        self.stop_time = stop_time

        
        self._t_step = t_step.copy()

        self._compile = compile

        self.observation_space = Box(
            low=obs_space_lower,
            high=obs_space_upper,
            dtype=np.float64
        )

        self.action_space = Box(
            low=action_space_lower,
            high=action_space_upper,
            dtype=np.float64
        )

        self.model = Voltage_Control(
            compile=self._compile
        )

        self._reset_state = False

    def reset(self):
        action = self.action_space.sample()*0.0
        
        if hasattr(self._rand_param, "__call__"):
            self.rand_param = self._rand_param()
        else:
            self.rand_param = self._rand_param
        rand_param = self.rand_param.copy()
        
        if hasattr(self._rand_input, "__call__"):
            rand_input = self._rand_input()
        else:
            rand_input = self._rand_input
        
        t_step = self._t_step.copy()
    
        self.model.reset(action, rand_param, t_step, rand_input)
        while self.model.current_time < self.reset_time:
            self.model.one_step(action.copy(), rand_param, t_step, rand_input)

        self._reset_state = True
        
        return self.model.obs.copy().reshape(self.observation_space.shape, order="F")


    def step(self, action, steps=10):
        action_ = action.astype(dtype=np.float64)
        if self._reset_state is False:
            raise Exception("Reset function not called. Please call reset() function first.")
        
        rand_param = self.rand_param.copy()
        
        if hasattr(self._rand_input, "__call__"):
            rand_input = self._rand_input()
        else:
            rand_input = self._rand_input
        
        t_step = self._t_step.copy()
        
        # Define reward here
        reward = -(1*(self.model.y_meas[2] - self.rand_param[2])
                   ** 2 + 0.05*action_[0]**2 + 0.005*action_[1]**2)

        if reward is None:
            raise Exception("Reward not defined!!!")
            
        for k in range(steps):
            self.model.one_step(action_.copy(), rand_param, t_step, rand_input)
            # if k==0:
            #     reward = self.model.reward[0]
        
        done = True if self.model.current_time >= self.stop_time else False
        if done:
            self._reset_state = False

        return self.model.obs.copy().reshape(self.observation_space.shape, order="F"), reward, done, {}
