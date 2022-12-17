"""
Relay buffer to store experience.

Author:     Niranjan Bhujel

"""

import tensorflow as tf
import numpy as np


class Buffer:
    def __init__(
        self, 
        buffer_size: int,
        var_name: list, 
        var_shape: list,
        var_dtype: list,
        seed=None) -> None:

        self.buffer_size = buffer_size
        self.var_name = var_name
        self._data = {}
        for k in range(len(var_name)):
            self._data[var_name[k]] = np.zeros(shape=(self.buffer_size, *var_shape[k]), dtype=var_dtype[k])

        self.pointer = 0
        self.num_data = 0

        self.rng = np.random.default_rng(seed=seed)
        
    def record(self, data_dict):   
        for k, v in data_dict.items():
            self._data[k][self.pointer] = v

        self.pointer += 1
        if self.num_data < self.buffer_size:
            self.num_data += 1

        if self.pointer == self.buffer_size:
            self.pointer = 0

    def get_data(self, output="numpy"):
        whole_slice = slice(0, self.num_data)
        data_out = {}

        for k, v in self._data.items():
            if output=='tensorflow':
                data_out[k] = tf.convert_to_tensor(v[whole_slice])
            else:
                 data_out[k] = v[whole_slice]

        return data_out
        
    def sample_batch(self, batch_size, output='numpy', probs=None):
        """
        Sample data of size `batch_size` from the buffer

        Parameters
        ----------
        batch_size : int
            Batch size of data to be sampled
        output : str, optional
            Whether output should be numpy or tensor, by default 'numpy'

        Returns
        -------
        dict
            Dictionary of data with following keywords: `obs`, `action`, `reward`, `next_obs`, `done`
        """        

        batch_indices = self.rng.choice(self.num_data, size=batch_size, p=probs)
        batch_out = {}
        for k, v in self._data.items():
            if output=='tensorflow':
                batch_out[k] = tf.convert_to_tensor(v[batch_indices])
            else:
                batch_out[k] = v[batch_indices]
        return batch_out

    def sample_mini_batch(self, batch_size):
        indices = np.arange(self.num_data)
        self.rng.shuffle(indices)
        
        split_indices = []
        point = batch_size
        while point < self.num_data:
            split_indices.append(point)
            point += batch_size

        temp_data = {}

        for k, v in self._data.items():
            temp_data[k] = np.split(v[indices], split_indices)

        n = len(temp_data[self.var_name[0]])
        data_out = []
        for i in range(n):
            tmp = {}
            for k, v in temp_data.items():
                tmp[k] = v[i]
            
            data_out.append(tmp)
        
        return data_out

    def clear(self):
        self.pointer = 0
        self.num_data = 0