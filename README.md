# ECE591-Project
Code for the project ECE591. The implementation of soft-actor algorithm for dynamic voltage support in microgrid.

## Required Packages
- Tensorflow 2.x
- OpenAI Gym 0.21.0 (newer version might not work)
- tbparse

All the codes are inside `code` folder.
Implementation of SAC algorithm is in `sac_algo` folder.
Model of grid (exported from simulink) is in `grid_model` folder.
Tensorboard log for some of the trainings are in `train_log` folder.

Example to perform training on general environment is shown in `test_pendulum.py`. 

To compile the model of environment, go to `grid_model` folder and enter following command on terminal
```{shell}
make shared_lib
```
Note: GNU make and gcc compiler should be installed to compile.



To perform training on voltage support, enter following command on terminal
```{shell}
python train_voltage_support.py
```


To check the performance of voltage support, enter following command on terminal
```{shell}
python test_voltage_support.py
```

To visualize training log, go to `train_log` and enter following command on terminal
```{shell}
tensorboard --logdir=./
```