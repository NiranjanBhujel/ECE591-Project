"""
Author:         Niranjan Bhujel
Description:    Implementation of soft actor critic (SAC) algorithm.
"""


import tensorflow as tf
from tensorflow import keras
import numpy as np
import gym
from .policies import ActorOnly
from .buffer import Buffer
from .utils import explained_variance, MetricsBuffer, SummaryLogger, to_batchtensor


class SACAgent:
    def __init__(
            self,
            env: gym.Env,
            pi_network,
            q1_network,
            q2_network,
            pi_optimizer=None,
            q_optimizer=None,
            ent_optimizer=None,
            ent_coeff=0.1,
            target_ent=None,
            auto_ent=True,
            initial_std=1.0,
            gamma=0.99,
            tau=0.005,
            buffer_capacity=1000000,
            batch_size=256,
            buffer_seed=None,
            train_interval=1,
            grad_steps=None) -> None:

        self.algo_name = "SAC"

        (
            self.train_env,
            self.pi_network,
            self.q1_network,
            self.q2_network,
            self.pi_optimizer,
            self.q_optimizer,
            self.ent_optimizer,
            self.ent_coeff,
            self.target_ent,
            self.auto_ent,
            self.initial_std,
            self.gamma,
            self.tau,
            self.buffer_capacity,
            self.batch_size,
            self.buffer_seed,
            self.train_interval,
            self.grad_steps,
        ) = (
            env,
            pi_network,
            q1_network,
            q2_network,
            pi_optimizer if pi_optimizer is not None else keras.optimizers.Adam(
                learning_rate=3e-4),
            q_optimizer if q_optimizer is not None else keras.optimizers.Adam(
                learning_rate=3e-4),
            ent_optimizer if ent_optimizer is not None else keras.optimizers.Adam(
                learning_rate=3e-4),
            ent_coeff,
            target_ent if target_ent is not None else -env.action_space.shape[0],
            auto_ent,
            initial_std,
            gamma,
            tau,
            buffer_capacity,
            batch_size,
            buffer_seed,
            train_interval,
            grad_steps,
        )

        self.obs_shape = self.train_env.observation_space.shape
        self.action_shape = self.train_env.action_space.shape

        self.setup_model()

    def setup_model(self):
        self.q1_target = keras.models.clone_model(self.q1_network)
        self.q2_target = keras.models.clone_model(self.q2_network)

        self.policy = ActorOnly(
            self.train_env.observation_space,
            self.train_env.action_space,
            pi_network=self.pi_network,
            initial_std=self.initial_std)

        self.buffer = Buffer(
            buffer_size=self.buffer_capacity,
            var_name=[
                "obs",
                "action",
                "reward",
                "next_obs",
                "done"],
            var_shape=[
                self.obs_shape,
                self.action_shape,
                (1,),
                self.obs_shape,
                (1,)],
            var_dtype=[np.float32] * 5
        )

        self.log_ent_coeff = tf.Variable(
            initial_value=float(tf.math.log(self.ent_coeff)) *
            tf.ones(1, dtype=tf.float32),
            trainable=True,
            dtype=tf.float32
        )

        self.q_trainable_params = self.q1_network.trainable_variables + \
            self.q2_network.trainable_variables


    @tf.function(experimental_follow_type_hints=True)
    def compute_q_loss(
            self,
            obs_batch: tf.Tensor,
            action_batch: tf.Tensor,
            reward_batch: tf.Tensor,
            next_value_batch: tf.Tensor,
            done_batch: tf.Tensor):

        q1_value = self.q1_network([obs_batch, self.policy.unscale_action(
            action_batch)], training=self.policy._training)
        q2_value = self.q2_network([obs_batch, self.policy.unscale_action(
            action_batch)], training=self.policy._training)

        q1_loss = 1/2*tf.reduce_mean(tf.square(q1_value - (
            reward_batch + self.gamma * (1-done_batch) * next_value_batch)))
        q2_loss = 1/2*tf.reduce_mean(tf.square(q2_value - (
            reward_batch + self.gamma * (1-done_batch) * next_value_batch)))

        explained_var1 = explained_variance(
            (reward_batch + self.gamma * (1-done_batch) * next_value_batch), q1_value)
        explained_var2 = explained_variance(
            (reward_batch + self.gamma * (1-done_batch) * next_value_batch), q2_value)
        return {
            "train/q_loss": q1_loss + q2_loss,
            "train/explained_variance": 0.5*explained_var1 + 0.5*explained_var2,
            "train/q_learning_rate": self.q_optimizer.learning_rate
        }

    @tf.function(experimental_follow_type_hints=True)
    def update_q(
            self,
            obs_batch: tf.Tensor,
            action_batch: tf.Tensor,
            reward_batch: tf.Tensor,
            next_obs_batch: tf.Tensor,
            done_batch: tf.Tensor):

        self.policy._training = True
        sampled_action, log_prob = self.policy.sample_action(
            next_obs_batch, batched=True)
        q1_value = self.q1_target(
            [next_obs_batch, self.policy.unscale_action(sampled_action)], training=False)
        q2_value = self.q2_target(
            [next_obs_batch, self.policy.unscale_action(sampled_action)], training=False)
        ent_coeff = tf.stop_gradient(tf.exp(self.log_ent_coeff))
        next_value_batch = tf.minimum(
            q1_value, q2_value) - ent_coeff * log_prob

        with tf.GradientTape(watch_accessed_variables=False) as tape:
            tape.watch(self.q_trainable_params)
            losses_metrics = self.compute_q_loss(
                obs_batch, action_batch, reward_batch, next_value_batch, done_batch)
            q_loss = losses_metrics["train/q_loss"]

        grads = tape.gradient(q_loss, self.q_trainable_params)

        self.q_optimizer.apply_gradients(
            zip(
                grads,
                self.q_trainable_params
            )
        )
        self.policy._training = False

        return losses_metrics

    @tf.function(experimental_follow_type_hints=True)
    def compute_pi_loss(
            self,
            obs_batch: tf.Tensor):

        sampled_action, log_prob = self.policy.sample_action(
            obs_batch, batched=True)
        q1_value = self.q1_target(
            [obs_batch, self.policy.unscale_action(sampled_action)], training=False)
        q2_value = self.q2_target(
            [obs_batch, self.policy.unscale_action(sampled_action)], training=False)
        q_value = tf.minimum(q1_value, q2_value)

        ent_coeff = tf.stop_gradient(tf.exp(self.log_ent_coeff))
        pi_loss = tf.reduce_mean(ent_coeff * log_prob - q_value)
        ent_coeff_loss = - tf.reduce_mean(self.log_ent_coeff * tf.stop_gradient(log_prob + self.target_ent))

        return {
            "train/pi_loss": pi_loss,
            "train/ent_coeff_loss": ent_coeff_loss,
            "train/ent_coeff": tf.exp(self.log_ent_coeff),
            "train/entropy": -tf.reduce_mean(log_prob),
        }

    @tf.function(experimental_follow_type_hints=True)
    def update_pi(
        self,
        obs_batch: tf.Tensor
    ):

        self.policy._training = True
        with tf.GradientTape(watch_accessed_variables=False, persistent=True) as tape:
            tape.watch(self.policy.trainable_params + [self.log_ent_coeff])
            losses_metrics = self.compute_pi_loss(obs_batch)
            pi_loss = losses_metrics["train/pi_loss"]
            ent_coeff_loss = losses_metrics["train/ent_coeff_loss"]

        grads_pi = tape.gradient(pi_loss, self.policy.trainable_params)

        self.pi_optimizer.apply_gradients(
            zip(
                grads_pi,
                self.policy.trainable_params
            )
        )

        if self.auto_ent:
            grads_ent = tape.gradient(ent_coeff_loss, [self.log_ent_coeff])

            self.ent_optimizer.apply_gradients(
                zip(
                    grads_ent,
                    [self.log_ent_coeff]
                )
            )

        self.policy._training = False

        if self.train_env.action_space.shape[0] == 1:
            losses_metrics['train/std'] = tf.exp(self.policy.log_std)
        else:
            for k in range(self.train_env.action_space.shape[0]):
                losses_metrics[f'train/std{k+1}'] = tf.exp(
                    self.policy.log_std)[k]

        return losses_metrics

    def update_network(self):
        grad_steps = self.grad_steps if self.grad_steps is not None else self.train_interval

        metrics_step = MetricsBuffer(action=np.mean)
        for _ in range(grad_steps):
            batch_data = self.buffer.sample_batch(
                self.batch_size, output='tensorflow')

            q_out_step = self.update_q(
                batch_data['obs'],
                batch_data['action'],
                batch_data['reward'],
                batch_data['next_obs'],
                batch_data['done'],
            )

            pi_out_step = self.update_pi(batch_data['obs'])

            self.update_target(self.q1_network, self.q1_target, tau=self.tau)
            self.update_target(self.q2_network, self.q2_target, tau=self.tau)

            metrics = {**pi_out_step, **q_out_step}
            metrics_step.add_data(metrics)

        metrics_step.add_data(
            {
                "train/pi_learning_rate": self.pi_optimizer.learning_rate,
                "train/q_learning_rate": self.q_optimizer.learning_rate,
            }
        )
        self.initial_std = tf.exp(self.policy.log_std).numpy()
        self.ent_coeff = tf.exp(self.log_ent_coeff).numpy()
        metrics_step.add_data({"train/target_entropy": self.target_ent})

        return metrics_step

    @tf.function
    def update_target(self, main_network, target_network, tau):
        """
        Update the provided target_network with main network with softupdate.
        """
        for (a, b) in zip(target_network.variables, main_network.variables):
            a.assign(b * tau + a * (1 - tau))

    def run_rollout(self, rollout_steps):
        n_steps = 0

        steps_to_take = rollout_steps

        for _ in range(steps_to_take):
            action, log_prob = self.policy.sample_action(self.obs)
            action, log_prob = action.numpy(), log_prob.numpy()

            next_obs, reward, done, _ = self.train_env.step(action)

            self.t += 1
            self.ep_rew += reward
            self.ep_len += 1
            n_steps += 1

            self.buffer.record(
                {
                    "obs": self.obs,
                    "action": action,
                    "reward": reward,
                    "next_obs": next_obs,
                    "done": done
                }
            )
            self.obs = next_obs

            if done or self.t >= self.total_timesteps:
                if done:
                    self.obs = self.train_env.reset()
                    self.running_reward.record(
                        {
                            "ep_rew": self.ep_rew,
                            "ep_len": self.ep_len,
                        }
                    )
                    self.ep_rew = 0.0
                    self.ep_len = 0
                    self.ep_count += 1
                break

        return n_steps, done

    def train(
            self,
            total_timesteps,
            smooth_horizon=40,
            tensorboard_log=None,
            tensorboard_filename="",
            log_interval=1):
        """
        Train the agents.

        Parameters
        ----------
        total_timesteps : int
            Total timesteps to run the simulation. If `done` is returned `True` from environment, the environment is reset and trained again.
        smooth_horizon : int, optional
            Number of episode to consider to calculate mean episodic length and reward.
        tensorboard_log : str, optional
            Folder name where log is to be written, by default None. If set `None`, log is not written. 
        tensorboard_filename : str, optional
            File name for tensorboard file, by default algorithm name is used
        verbose : bool, optional
            Whether to print training info, by default False
        log_interval: int, optional
            Log every `log_interval` episodes, by default 1
        """

        self.running_reward = Buffer(
            buffer_size=smooth_horizon,
            var_name=["ep_rew", "ep_len"],
            var_shape=[(1,), (1,)],
            var_dtype=[np.float32, np.float32]
        )
        self.ep_rew = 0.0
        self.ep_len = 0
        self.ep_count = 0

        self.total_timesteps = total_timesteps
        tensorboard_filename = self.algo_name if tensorboard_filename == "" else tensorboard_filename
        if tensorboard_log is not None:
            summary_writer = SummaryLogger(
                tensorboard_log, tensorboard_filename)

        self.t = 0

        metrics_episode = MetricsBuffer(action=np.mean)
        self.obs = self.train_env.reset()

        print("Training started")
        while True:
            rollout_steps = min(self.total_timesteps -
                                self.t, self.train_interval)
            done = True
            if rollout_steps > 0:
                _, done = self.run_rollout(rollout_steps)

            if self.t >= self.batch_size:
                metrics_step = self.update_network()
                metrics_episode.add_data(metrics_step.get_data())

            if done:
                self.obs = self.train_env.reset()
                if self.running_reward.num_data > 0:
                    data = self.running_reward.get_data()
                    ep_rew_mean = np.mean(data["ep_rew"])
                    ep_len_mean = np.mean(data["ep_len"])
                else:
                    ep_rew_mean = 0.0
                    ep_len_mean = 0

                if tensorboard_log is not None and self.ep_count % log_interval == 0:
                    summary_writer.add_scalar(
                        "rollout/ep_rew_mean", ep_rew_mean, self.t)
                    summary_writer.add_scalar(
                        "rollout/ep_len_mean", ep_len_mean, self.t)
                    summary_writer.add_scalars(
                        metrics_episode.get_data(), self.t)

                metrics_episode.clear()

            if self.t >= self.total_timesteps:
                break
            
            if self.t % int(0.05 * self.total_timesteps) == 0:
                print(f"{int((self.t + 1) / self.total_timesteps * 100)}% complete!!!")

        if tensorboard_log is not None:
            summary_writer.close()

        print("Training ended")
