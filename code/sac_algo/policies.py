import tensorflow as tf
from tensorflow import keras
import tensorflow_probability as tfp
from gym import Space, spaces
from .utils import to_batchtensor


class SquashedDiagGaussian:
    def __init__(self, mu, std, epsilon=1e-6):
        """
        Guassian distribution with tanh squashing

        Parameters
        ----------
        mu : tf.Tensor
            Mean
        std : tf.Tensor
            Standard deviation
        epsilon : float, optional
            Small value to find stable log_prob, by default 1e-6
        """        
        self._mu = mu
        self._std = std
        self._epsilon = epsilon
        self._distrib = tfp.distributions.Normal(mu, std)

    def mean(self):
        return tf.math.tanh(self._mu)

    def sample(self, log_prob=True):
        gaussian_x = self._distrib.sample()
        x = tf.math.tanh(gaussian_x)
        if log_prob:
            log_prob = self._distrib.log_prob(
                gaussian_x) - tf.math.log(1-tf.square(x) + self._epsilon)
            return x, tf.reduce_sum(log_prob, axis=1, keepdims=True)
        else:
            return x

    def log_prob(self, x, gaussian_x=None):
        x_safe = tf.clip_by_value(x, -1+self._epsilon, 1-self._epsilon)
        if gaussian_x is None:
            gaussian_x = 1/2 * (tf.math.log1p(x_safe) - tf.math.log1p(-x_safe))
        log_prob = self._distrib.log_prob(
            gaussian_x) - tf.math.log(1-tf.square(x) + self._epsilon)
        return tf.reduce_sum(log_prob, axis=1, keepdims=True)


class ActorOnly:
    def __init__(self, obs_space, action_space, pi_network, initial_std=1):
        self.obs_space = obs_space
        self.action_space = action_space
        self.action_shape = action_space.shape
        self.pi_network = pi_network
        self.log_std = tf.Variable(
            initial_value=tf.cast(tf.math.log(
                initial_std), dtype=tf.float32) * tf.ones(self.action_shape, dtype=tf.float32),
            trainable=True,
            dtype=tf.float32
        )

        self.trainable_params = self.pi_network.trainable_variables + \
            [self.log_std]
        self._training = False
        self._scale = (action_space.high - action_space.low) / 2
        self._offset = (action_space.high + action_space.low) / 2

    @tf.function
    def _forward(self, obs):
        """
        Function to find shape parameter pf distribution from actor network.

        Parameters
        ----------
        obs : tf.Tensor
            Observation

        Returns
        -------
        tuple
            Mean and standard deviation
        """
        out = self.pi_network(obs, training=self._training)
        mu = out
        distrib_param = (mu, tf.ones_like(mu) * tf.exp(self.log_std))

        return distrib_param

    def forward(self, obs):
        """
        Function to create distribution

        Parameters
        ----------
        obs : tf.Tensor
            Observation

        Returns
        -------
        Squashed Gaussian distributiomn
        """
        distrib_param = self._forward(obs)
        distrib = SquashedDiagGaussian(*distrib_param)
        return distrib


    @tf.function
    def sample_action(self, obs, batched=False):
        obs_tf = obs if batched else to_batchtensor(obs)
        distrib = self.forward(obs_tf)
        action, log_prob = distrib.sample(log_prob=True)
        action_scaled = self.scale_action(action)

        if batched:
            return action_scaled, log_prob
        else:
            return action_scaled[0], log_prob[0]

    @tf.function
    def evaluate_action(self, obs, action):
        distrib = self.forward(obs)
        log_prob = distrib.log_prob(self.unscale_action(action))

        return log_prob


    @tf.function
    def scale_action(self, action):
        return action * self._scale + self._offset

    @tf.function
    def unscale_action(self, action):
        return (action - self._offset) / self._scale

    @tf.function
    def predict(self, obs, deterministic=True):
        obs_tf = to_batchtensor(obs)
        distrib = self.forward(obs_tf)
        if deterministic:
            action = distrib.mean()
        else:
            action = distrib.sample(log_prob=False)

        return self.scale_action(action)[0]