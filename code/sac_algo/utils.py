"""
Collection of utility function/classes.

Author:     Niranjan Bhujel
Date:       July 29, 2022
"""

import numpy as np
import tensorflow as tf
import os
import tbparse


class MetricsBuffer:
    def __init__(self, action=np.sum):
        """
        Create buffer for metrics. 

        It stores calculated metrics in the buffer which can be accessed by performing desired reduction.

        Parameters
        ----------
        action : callable, optional
            Callable that return reduced metric, by default np.sum
        """        
        self.metrics_data = {}
        self.action = action
        
    def add_data(self, metrics):
        """
        Add new metrics

        Parameters
        ----------
        metrics : dict
            Dictionary of metrics
        """        
        for k in metrics:
            if k in self.metrics_data:
                self.metrics_data[k].append(float(metrics[k]))
            else:
                self.metrics_data[k] = [float(metrics[k])]
    
    def get_data(self, action_dict={}) -> dict:
        """
        Get reduced metrics.

        Parameters
        ----------
        action_dict : dict, optional
            Reduction to apply for specific metrics. Default reduction is applied for metrics not specified here., by default {}

        Returns
        -------
        dict
            Dictionary of reduced metrics
        """        
        metrics_out = {}
        for k in self.metrics_data:
            if k not in action_dict.keys():
                metrics_out[k] = self.action(np.array(self.metrics_data[k])).copy()
            else:
                metrics_out[k] = action_dict[k](np.array(self.metrics_data[k])).copy()
        
        return metrics_out

    def clear(self):
        self.metrics_data = {}
    

def to_batchtensor(obs):
    return tf.expand_dims(tf.cast(obs, dtype=tf.float32), 0)


def explained_variance(y_true, y_pred):
    """
    Calculate explained variance

    Parameters
    ----------
    y_true : array like
        True value or target value
    y_pred : array like
        Predicted value

    Returns
    -------
    float
        Explained variance
    """    
    var_y = tf.math.reduce_variance(y_true)
    var_d = tf.math.reduce_variance(y_true - y_pred)
    expl_var = 1-var_d / var_y
    return np.nan if var_y==0 else tf.clip_by_value(expl_var, 0.0, 1.0)


class SummaryLogger:
    def __init__(self, tensorboard_log, tensorboard_filename):
        if tensorboard_log not in os.listdir('./') and tensorboard_log != "./":
            os.mkdir(tensorboard_log)
        file_counter = 1
        while f"{tensorboard_filename}_{file_counter}" in os.listdir(tensorboard_log):
            file_counter += 1
        self.summary_writer = tf.summary.create_file_writer(
            os.path.join(tensorboard_log, f"{tensorboard_filename}_{file_counter}")
        )

        print(f"Logging to " + os.path.join(tensorboard_log, f"{tensorboard_filename}_{file_counter}"))

    def add_scalar(self, key, val, step):
        with self.summary_writer.as_default():
            tf.summary.scalar(key, val, step)
    
    def add_scalars(self, value_dict, step):
        with self.summary_writer.as_default():
            for k, val in value_dict.items():
                    tf.summary.scalar(k, val, step)

    def close(self):
        self.summary_writer.close()

def get_parsed_metrics(filenames):
    """
    Read tensorboard log file and returns steps and metrics.

    Parameters
    ----------
    filenames : str or List(str)
        Filename or list of filename

    Returns
    -------
    dict
        Dictionary of numpy array. Key to the dictionary is same as name of metrics. Value is numpy array whose first column is step index and second column is corresponding metrics value
    """  
    def get_parsed_single(f):
        reader = tbparse.SummaryReader(f)
        df = reader.tensors
        colnames = df['tag'].unique()

        df_ = {}

        for col in colnames:
            steps = np.array(df.loc[df['tag']==col]['step']).reshape((-1, 1))
            values = np.array(df.loc[df['tag']==col]['value']).reshape((-1, 1))
            df_[col] = np.hstack((steps, values))
        
        return df_
    
    if isinstance(filenames, str):
        return get_parsed_single(filenames)
    elif isinstance(filenames, list):
        return [get_parsed_single(f) for f in filenames]
    else:
        raise Exception("Not supported data types for filenames argument. Must be either str or list of str.")