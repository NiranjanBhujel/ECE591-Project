import tensorflow as tf
from tensorflow import keras
from sac_algo.sac import SACAgent
from sac_algo.utils import get_parsed_metrics
from grid_model.Voltage_Control_gyminterface import Voltage_Control_gym
import numpy as np
from grid_model.Voltage_Control_gyminterface import Voltage_Control_gym
import matplotlib.pyplot as plt
import matplotlib as mpl
import gym

# os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

plt.rcParams['mathtext.fontset'] = 'stix'
mpl.rc('font', family='times new roman')
mpl.rcParams['savefig.dpi'] = 600
mpl.rcParams['axes.unicode_minus'] = False
mpl.rcParams['agg.path.chunksize'] = 10000
plt.rcParams['font.size'] = 10

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
    # output = keras.layers.Dense(
    #     64,
    #     activation="tanh",
    #     kernel_initializer=keras.initializers.orthogonal(gain=np.sqrt(2)))(input)
    # output = keras.layers.Flatten()(input)
    output = keras.layers.LSTM(10) (input)
    
    # Policy only
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

    # model = keras.Model(inputs=input, outputs=[action_out, output])
    model = keras.Model(inputs=[obs_input, action_input], outputs=q_out)
    return model

t0, t1 = 0.1495, 0.154
t_step = 0.15

env = Voltage_Control_gym(
    obs_space_lower=-np.ones((25, 7), dtype=np.float64),
    obs_space_upper=np.ones((25, 7), dtype=np.float64),
    action_space_lower=np.array([-0.5, -0.5]),
    action_space_upper=np.array([0.5, 0.5]),
    rand_param=np.array([0.3, 0.2, 0.95]),
    rand_input=lambda : np.random.normal(
        np.array([0.0, 0.0, 0.0, 0.0]),
        np.array([0.2e-5, 0.2e-5, 0.2e-4, 0.2e-4])**0.5,
    ),
    t_step=np.array([t_step]),
    reset_time=0.0,
    stop_time=t1,
    compile=False,
)

pi_net = keras.models.load_model("saved_network/pi_network.h5")
q1_net = get_q_network(env)
q2_net = get_q_network(env)

agent = SACAgent(
    env=env,
    pi_network=pi_net,
    q1_network=q1_net,
    q2_network=q2_net,
)

done = False
data_sac = {
    "t": [],
    "igd": [],
    "vcd": [],
    "iinvd": [],
    "iinvq": [],
}

# With control action
obs = env.reset()
count = 0
while not done:
    if count % 10==0:
        action = agent.policy.predict(obs, deterministic=True).numpy()
        # action = action.numpy()
        # pi_out = agent.policy.pi_network(tf.expand_dims(obs, 0), training=False)
        # mu = pi_out
        # action = np.tanh(mu[0].numpy().astype(dtype=np.float64)) * 0.5

    data_sac["t"].append(env.model.current_time)
    data_sac["igd"].append(env.model.y_true[0])
    data_sac["vcd"].append(env.model.y_true[2])
    data_sac["iinvd"].append(action[0])
    data_sac["iinvq"].append(action[1])

    obs, reward, done, _ = env.step(action, steps=1)
    count += 1


done = False
data_no_control = {
    "t": [],
    "igd": [],
    "vcd": [],
    "iinvd": [],
    "iinvq": [],
}


obs = env.reset()
count = 0
while not done:
    if count % 10==0:
        # action = agent.predict(obs, deterministic=True)
        action = np.array([0.0, 0.0])

    data_no_control["t"].append(env.model.current_time)
    data_no_control["igd"].append(env.model.y_true[0])
    data_no_control["vcd"].append(env.model.y_true[2])
    data_no_control["iinvd"].append(action[0])
    data_no_control["iinvq"].append(action[1])

    obs, reward, done, _ = env.step(action, steps=1)
    count += 1


fig, ax = plt.subplots(2, 3, figsize=(7.1, 1.95*2), constrained_layout=True)

data_ = get_parsed_metrics("train_log/SAC_test_8")
ax[0,0].plot(
    data_["rollout/ep_rew_mean"][:,0],
    data_["rollout/ep_rew_mean"][:,1]
)

ax[0,0].set_xlabel("timesteps\n(a)")
ax[0,0].set_ylabel("mean episodic reward")

ax[0,0].grid(alpha=0.3, linestyle="-.")

ax[0,0].ticklabel_format(axis="x", style="scientific", scilimits=(6, 6), useMathText=True)


ax[0,1].plot(
    data_["train/entropy"][:,0],
    data_["train/entropy"][:,1]
)

ax[0,1].set_xlabel("timesteps\n(b)")
ax[0,1].set_ylabel("entropy")

ax[0,1].grid(alpha=0.3, linestyle="-.")

ax[0,1].ticklabel_format(axis="x", style="scientific", scilimits=(6, 6), useMathText=True)


ax[0,2].plot(
    data_["train/std1"][:,0],
    data_["train/std1"][:,1],
    label="$i_{invd}$"
)

ax[0,2].plot(
    data_["train/std2"][:,0],
    data_["train/std2"][:,1],
    label="$i_{invq}$"
)

ax[0,2].set_xlabel("timesteps\n(c)")
ax[0,2].set_ylabel("$\sigma$")

ax[0,2].grid(alpha=0.3, linestyle="-.")
ax[0,2].legend(loc=0, labelspacing=0.15, borderpad=0.15)
ax[0,2].ticklabel_format(axis="x", style="scientific", scilimits=(6, 6), useMathText=True)

tk = (t0 <= np.array(data_sac["t"])) & (np.array(data_sac["t"]) <= t1)
ax[1,0].plot(
    np.array(data_no_control["t"])[tk],
    np.array(data_no_control["igd"])[tk],
    label="no control"
)
ax[1,0].plot(
    np.array(data_sac["t"])[tk],
    np.array(data_sac["igd"])[tk],
    label="with SAC"
)

ax[1,0].set_xlabel("time [s]\n(d)")
ax[1,0].set_ylabel("$i_{gd}$ [p.u.]")
ax[1,0].legend(loc=0, labelspacing=0.15, borderpad=0.15)
ax[1,0].grid(True, linestyle="--", alpha=0.3)
ax[1,0].set_xlim(left=t0, right=t1)

ax[1,1].plot(
    np.array(data_no_control["t"])[tk],
    np.array(data_no_control["vcd"])[tk],
    label="no control"
)
ax[1,1].plot(
    np.array(data_sac["t"])[tk],
    np.array(data_sac["vcd"])[tk],
    label="with SAC"
)

ax[1,1].set_xlabel("time [s]\n(e)")
ax[1,1].set_ylabel("$v_{cd}$ [p.u.]")
ax[1,1].legend(loc=0, labelspacing=0.15, borderpad=0.15)
ax[1,1].grid(True, linestyle="--", alpha=0.3)
ax[1,1].set_xlim(left=t0, right=t1)


tk = (t0 <= np.array(data_sac["t"])) & (np.array(data_sac["t"]) <= t1)
ax[1,2].plot(
    np.array(data_sac["t"])[tk],
    np.array(data_sac["iinvd"])[tk],
    label="$i_{invd}$"
)
ax[1,2].plot(
    np.array(data_sac["t"])[tk],
    np.array(data_sac["iinvq"])[tk],
    label="$i_{invq}$"
)

ax[1,2].set_xlabel("time [s]\n(f)")
ax[1,2].set_ylabel("$i_{inv}$ [p.u.]")
ax[1,2].legend(loc=0, labelspacing=0.15, borderpad=0.15)
ax[1,2].grid(True, linestyle="--", alpha=0.3)
ax[1,2].set_xlim(left=t0, right=t1)

plt.savefig("Results/perf.svg", dpi=600)
plt.show()