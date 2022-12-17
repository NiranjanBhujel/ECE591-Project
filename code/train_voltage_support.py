import tensorflow as tf
from tensorflow import keras
from sac_algo.sac import SACAgent
import gym
import os
import sys
from grid_model.Voltage_Control_gyminterface import Voltage_Control_gym
import numpy as np

os.environ["CUDA_VISIBLE_DEVICES"] = "-1"


def get_pi_network(env: gym.Env):
    """
    Function to create actor network.
    """
    obs_dim = env.observation_space.shape
    action_dim = env.action_space.shape[0]
    lower_bound = env.action_space.low
    upper_bound = env.action_space.high

    # Input and shared
    input = keras.layers.Input(shape=obs_dim)
    output = keras.layers.LSTM(10) (input)
    
    pi_out = keras.layers.Dense(
        64,
        activation="relu")(output)
    pi_out = keras.layers.Dense(
        64,
        activation="relu")(pi_out)

    action_out = keras.layers.Dense(
        action_dim)(pi_out)
    
    return keras.Model(inputs=input, outputs=action_out)

def get_q_network(env: gym.Env):
    obs_dim = env.observation_space.shape
    action_dim = env.action_space.shape[0]
    lower_bound = env.action_space.low
    upper_bound = env.action_space.high

    obs_input = keras.layers.Input(shape=obs_dim)
    action_input = keras.layers.Input(shape=action_dim)

    obs_feature = keras.layers.LSTM(10)(obs_input)

    combined_feature = keras.layers.Concatenate()([obs_feature, action_input])
    # Value only
    q_out = keras.layers.Dense(
        64,
        activation="relu")(combined_feature)
    
    q_out = keras.layers.Dense(
        64,
        activation="relu")(q_out)

    q_out = keras.layers.Dense(
        1)(q_out)


    model = keras.Model(inputs=[obs_input, action_input], outputs=q_out)
    return model


env = Voltage_Control_gym(
    obs_space_lower=-np.ones((25, 7), dtype=np.float64),
    obs_space_upper=np.ones((25, 7), dtype=np.float64),
    action_space_lower=np.array([-0.5, -0.5]),
    action_space_upper=np.array([0.5, 0.5]),
    rand_param=lambda : np.random.uniform(
        np.array([0.2, -0.2, 0.8]), 
        np.array([0.8, 0.2, 1.2])),
    rand_input=lambda : np.random.normal(
        np.array([0.0, 0.0, 0.0, 0.0]),
        np.array([1e-5, 1e-5, 1e-4, 1e-4])**0.5 / 2,
    ),
    t_step=np.array([0.15]),
    reset_time=0.140,
    stop_time=0.170,
    compile=True,
    )

pi_net = get_pi_network(env)
q1_net = get_q_network(env)
q2_net = get_q_network(env)

print("pi net:\n")
print(pi_net.summary())
print("q net:\n")
print(q1_net.summary())

agent = SACAgent(
    env=env,
    pi_network=pi_net,
    q1_network=q1_net,
    q2_network=q2_net,
)


agent.train(
    total_timesteps=int(1e6),
    tensorboard_log="train_log",
    tensorboard_filename="SAC_test_noise",
)

agent.pi_network.save("saved_network/pi_network.h5")

