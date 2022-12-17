import tensorflow as tf
from tensorflow import keras
from sac_algo.sac import SACAgent
import gym
import os
import sys

os.environ["CUDA_VISIBLE_DEVICES"] = "-1"


def get_q_network(obs_space, action_space):
    state_input = keras.layers.Input(shape=(obs_space.shape[0]))
    action_input = keras.layers.Input(shape=(action_space.shape[0]))

    q_concat = keras.layers.Concatenate()([state_input, action_input])
    q_out = keras.layers.Dense(256, activation="relu")(q_concat)
    q_out = keras.layers.Dense(256, activation="relu")(q_out)
    q_out = keras.layers.Dense(1)(q_out)

    model = keras.Model(inputs=[state_input, action_input], outputs=q_out)
    return model


def get_pi_network(obs_space, action_space):
    input = keras.layers.Input(shape=obs_space.shape[0])
    output = keras.layers.Dense(
        256,
        activation="relu")(input)
    output = keras.layers.Dense(
        256,
        activation="relu")(output)

    mu_out = keras.layers.Dense(
        action_space.shape[0],
    )(output)
    model = keras.Model(inputs=input, outputs=mu_out)
    return model


env = gym.make("Pendulum-v1")

if sys.argv[1] == "train":
    agent = SACAgent(
        env=env,
        pi_network=get_pi_network(env.observation_space, env.action_space),
        q1_network=get_q_network(env.observation_space, env.action_space),
        q2_network=get_q_network(env.observation_space, env.action_space)
    )

    agent.train(
        total_timesteps=int(200*200),
        tensorboard_log="train_test"
    )
    agent.pi_network.save("pi_network.h5")

elif sys.argv[1] == "test":
    agent = SACAgent(
        env=env,
        pi_network=keras.models.load_model("pi_network.h5"),
        q1_network=get_q_network(env.observation_space, env.action_space),
        q2_network=get_q_network(env.observation_space, env.action_space)
    )

    obs = env.reset()
    done = False

    while not done:
        env.render(mode="human")
        action = agent.policy.predict(obs)
        obs, reward, done, info = env.step(action)
