"""
Filename:       Voltage_Control_pyinterface.py
Written by:     Niranjan Bhujel
Date:           2022-11-30
"""


import os
import ctypes
import numpy as np
import h5py
import time
import argparse


class Voltage_Control:
    def __init__(self, out_lib_name: str="Voltage_Control_lib.so", compile: bool=True):
        """
        Class to simulate the Voltage_Control system. The model input and output are shown below:
        
        Input:
            action [2]
            rand_param [3]
            t_step [1]
            rand_input [4]
            
        Output:
            obs [175]
            y_meas [4]
            y_true [4]
            
        """
        filepath = __file__
        filepath = os.path.split(filepath)[0]
        if filepath=="":
            filepath = "./"

        if compile:
            compile_cmd = "gcc -fPIC -shared -O3 -lm " + '"' +\
                os.path.join(filepath, 'Voltage_Control_interface.c') + '"' + ' "' +\
                os.path.join(filepath, 'sim_code"/*.c') + ' -o ' + '"' +\
                os.path.join(filepath, out_lib_name) + '"'
            # print(compile_cmd)
            os.system(compile_cmd)

        self.start_time = 0.0
        self.step_time = 1.0E-5
        self.current_time = self.start_time
        self.next_time = self.current_time + self.step_time
        
        self.Voltage_Control_system = ctypes.cdll.LoadLibrary(os.path.join(filepath, out_lib_name))

        self.Voltage_Control_system.initialize.restype = None
        self.Voltage_Control_system.initialize.argtypes = [
            np.ctypeslib.ndpointer(dtype=np.float64),
            np.ctypeslib.ndpointer(dtype=np.float64),
            np.ctypeslib.ndpointer(dtype=np.float64),
            np.ctypeslib.ndpointer(dtype=np.float64),
            np.ctypeslib.ndpointer(dtype=np.float64),
            np.ctypeslib.ndpointer(dtype=np.float64),
            np.ctypeslib.ndpointer(dtype=np.float64),
            
        ]

        self.Voltage_Control_system.one_step.restype = None
        self.Voltage_Control_system.one_step.argtypes = [
            np.ctypeslib.ndpointer(dtype=np.float64),
            np.ctypeslib.ndpointer(dtype=np.float64),
            np.ctypeslib.ndpointer(dtype=np.float64),
            np.ctypeslib.ndpointer(dtype=np.float64),
            np.ctypeslib.ndpointer(dtype=np.float64),
            np.ctypeslib.ndpointer(dtype=np.float64),
            np.ctypeslib.ndpointer(dtype=np.float64),
            
        ]

        self.obs = np.zeros((175,), dtype=np.float64)
        self.y_meas = np.zeros((4,), dtype=np.float64)
        self.y_true = np.zeros((4,), dtype=np.float64)
        
    
    def reset(self, action: np.ndarray, rand_param: np.ndarray, t_step: np.ndarray, rand_input: np.ndarray):
        """
        Reset the simulation back to initial state.
        """
        self.current_time = self.start_time
        self.next_time = self.current_time + self.step_time
        self.Voltage_Control_system.initialize(action, rand_param, t_step, rand_input,  self.obs, self.y_meas, self.y_true)
        
        # return {
        #     "obs": self.obs,
        #     "y_meas": self.y_meas,
        #     "y_true": self.y_true,
        # 
        # }

    def one_step(self, action: np.ndarray, rand_param: np.ndarray, t_step: np.ndarray, rand_input: np.ndarray):
        """
        Run one step of simulation
        """
        self.Voltage_Control_system.one_step(action, rand_param, t_step, rand_input,  self.obs, self.y_meas, self.y_true)
        self.current_time += self.step_time
        self.next_time = self.current_time + self.step_time

        # return {
        #     "obs": self.obs,
        #     "y_meas": self.y_meas,
        #     "y_true": self.y_true,
        # 
        # }


# if __name__=="__main__":
#     parser = argparse.ArgumentParser(
#         prog="Voltage_Control_pyinterface",
#         description="Python interface to Voltage_Control.",
#     )
#     parser.add_argument("--stop_time", help="Simulation stop time. Total time of simulation.", type=float)
#     parser.add_argument("--decimation", help="Decimation factor. Data is saved once every specified decimation factor time steps.", default=1, type=int)
#     args = parser.parse_args()
#     if args.stop_time is None:
#         raise Exception("Stop time not specified. Type `python Voltage_Control_pyinterface.py --help` for more info!!!")

#     STEP_TIME = 1.0E-5
#     STOP_TIME = args.stop_time
#     DECIMATION = args.decimation

#     if STEP_TIME==-1:
#         raise Exception("STEP_TIME not specified.")
#     if STOP_TIME==-1:
#         raise Exception("STOP_TIME not specified.")
    
#     def to_numpy(*args):
#         return np.array(args, dtype=np.float64)

#     m = Voltage_Control("Voltage_Control_lib.so", step_time=STEP_TIME, compile=True)
    
#     DATA_OUT = np.zeros((1+int(STOP_TIME/(STEP_TIME*DECIMATION)), 1+175 + 4 + 4), dtype=np.float64)

#     # Specify value of input here
#     # action = 
#     # rand_param = 
#     # t_step = 
#     # rand_input = 
#     
#     m.reset(action, rand_param, t_step, rand_input)
#     DATA_OUT[0,:] = to_numpy(m.current_time, m.obs[0], m.obs[1], m.obs[2], m.obs[3], m.obs[4], m.obs[5], m.obs[6], m.obs[7], m.obs[8], m.obs[9], m.obs[10], m.obs[11], m.obs[12], m.obs[13], m.obs[14], m.obs[15], m.obs[16], m.obs[17], m.obs[18], m.obs[19], m.obs[20], m.obs[21], m.obs[22], m.obs[23], m.obs[24], m.obs[25], m.obs[26], m.obs[27], m.obs[28], m.obs[29], m.obs[30], m.obs[31], m.obs[32], m.obs[33], m.obs[34], m.obs[35], m.obs[36], m.obs[37], m.obs[38], m.obs[39], m.obs[40], m.obs[41], m.obs[42], m.obs[43], m.obs[44], m.obs[45], m.obs[46], m.obs[47], m.obs[48], m.obs[49], m.obs[50], m.obs[51], m.obs[52], m.obs[53], m.obs[54], m.obs[55], m.obs[56], m.obs[57], m.obs[58], m.obs[59], m.obs[60], m.obs[61], m.obs[62], m.obs[63], m.obs[64], m.obs[65], m.obs[66], m.obs[67], m.obs[68], m.obs[69], m.obs[70], m.obs[71], m.obs[72], m.obs[73], m.obs[74], m.obs[75], m.obs[76], m.obs[77], m.obs[78], m.obs[79], m.obs[80], m.obs[81], m.obs[82], m.obs[83], m.obs[84], m.obs[85], m.obs[86], m.obs[87], m.obs[88], m.obs[89], m.obs[90], m.obs[91], m.obs[92], m.obs[93], m.obs[94], m.obs[95], m.obs[96], m.obs[97], m.obs[98], m.obs[99], m.obs[100], m.obs[101], m.obs[102], m.obs[103], m.obs[104], m.obs[105], m.obs[106], m.obs[107], m.obs[108], m.obs[109], m.obs[110], m.obs[111], m.obs[112], m.obs[113], m.obs[114], m.obs[115], m.obs[116], m.obs[117], m.obs[118], m.obs[119], m.obs[120], m.obs[121], m.obs[122], m.obs[123], m.obs[124], m.obs[125], m.obs[126], m.obs[127], m.obs[128], m.obs[129], m.obs[130], m.obs[131], m.obs[132], m.obs[133], m.obs[134], m.obs[135], m.obs[136], m.obs[137], m.obs[138], m.obs[139], m.obs[140], m.obs[141], m.obs[142], m.obs[143], m.obs[144], m.obs[145], m.obs[146], m.obs[147], m.obs[148], m.obs[149], m.obs[150], m.obs[151], m.obs[152], m.obs[153], m.obs[154], m.obs[155], m.obs[156], m.obs[157], m.obs[158], m.obs[159], m.obs[160], m.obs[161], m.obs[162], m.obs[163], m.obs[164], m.obs[165], m.obs[166], m.obs[167], m.obs[168], m.obs[169], m.obs[170], m.obs[171], m.obs[172], m.obs[173], m.obs[174], m.y_meas[0], m.y_meas[1], m.y_meas[2], m.y_meas[3], m.y_true[0], m.y_true[1], m.y_true[2], m.y_true[3])
    
#     LAST_CURRENT_TIME = m.current_time
#     WALL_TIME_ORIGIN = time.time_ns()
#     __counter__ = 0
#     while m.current_time < STOP_TIME:
#         # Specify value of input here
#         # action = 
#         # rand_param = 
#         # t_step = 
#         # rand_input = 
#         

#         m.one_step(action, rand_param, t_step, rand_input)
#         __counter__ += 1

#         if __counter__ % DECIMATION == 0:
#             try:
#                 DATA_OUT[int(__counter__ / DECIMATION),:] = to_numpy(m.current_time, m.obs[0], m.obs[1], m.obs[2], m.obs[3], m.obs[4], m.obs[5], m.obs[6], m.obs[7], m.obs[8], m.obs[9], m.obs[10], m.obs[11], m.obs[12], m.obs[13], m.obs[14], m.obs[15], m.obs[16], m.obs[17], m.obs[18], m.obs[19], m.obs[20], m.obs[21], m.obs[22], m.obs[23], m.obs[24], m.obs[25], m.obs[26], m.obs[27], m.obs[28], m.obs[29], m.obs[30], m.obs[31], m.obs[32], m.obs[33], m.obs[34], m.obs[35], m.obs[36], m.obs[37], m.obs[38], m.obs[39], m.obs[40], m.obs[41], m.obs[42], m.obs[43], m.obs[44], m.obs[45], m.obs[46], m.obs[47], m.obs[48], m.obs[49], m.obs[50], m.obs[51], m.obs[52], m.obs[53], m.obs[54], m.obs[55], m.obs[56], m.obs[57], m.obs[58], m.obs[59], m.obs[60], m.obs[61], m.obs[62], m.obs[63], m.obs[64], m.obs[65], m.obs[66], m.obs[67], m.obs[68], m.obs[69], m.obs[70], m.obs[71], m.obs[72], m.obs[73], m.obs[74], m.obs[75], m.obs[76], m.obs[77], m.obs[78], m.obs[79], m.obs[80], m.obs[81], m.obs[82], m.obs[83], m.obs[84], m.obs[85], m.obs[86], m.obs[87], m.obs[88], m.obs[89], m.obs[90], m.obs[91], m.obs[92], m.obs[93], m.obs[94], m.obs[95], m.obs[96], m.obs[97], m.obs[98], m.obs[99], m.obs[100], m.obs[101], m.obs[102], m.obs[103], m.obs[104], m.obs[105], m.obs[106], m.obs[107], m.obs[108], m.obs[109], m.obs[110], m.obs[111], m.obs[112], m.obs[113], m.obs[114], m.obs[115], m.obs[116], m.obs[117], m.obs[118], m.obs[119], m.obs[120], m.obs[121], m.obs[122], m.obs[123], m.obs[124], m.obs[125], m.obs[126], m.obs[127], m.obs[128], m.obs[129], m.obs[130], m.obs[131], m.obs[132], m.obs[133], m.obs[134], m.obs[135], m.obs[136], m.obs[137], m.obs[138], m.obs[139], m.obs[140], m.obs[141], m.obs[142], m.obs[143], m.obs[144], m.obs[145], m.obs[146], m.obs[147], m.obs[148], m.obs[149], m.obs[150], m.obs[151], m.obs[152], m.obs[153], m.obs[154], m.obs[155], m.obs[156], m.obs[157], m.obs[158], m.obs[159], m.obs[160], m.obs[161], m.obs[162], m.obs[163], m.obs[164], m.obs[165], m.obs[166], m.obs[167], m.obs[168], m.obs[169], m.obs[170], m.obs[171], m.obs[172], m.obs[173], m.obs[174], m.y_meas[0], m.y_meas[1], m.y_meas[2], m.y_meas[3], m.y_true[0], m.y_true[1], m.y_true[2], m.y_true[3])
#             except:
#                 pass

#         if (m.current_time - LAST_CURRENT_TIME) / STOP_TIME * 100 > 5:
#             print(f"{round(m.current_time / STOP_TIME * 100, 0)}% complete. Total time elapsed: {(time.time_ns() - WALL_TIME_ORIGIN)/1e9} !!!")
#             LAST_CURRENT_TIME = m.current_time

#     if round(LAST_CURRENT_TIME / STOP_TIME * 100, 0)!=5:
#         print(f"{round(100.0, 0)}% complete. Total time elapsed: {(time.time_ns() - WALL_TIME_ORIGIN)/1e9}!!!")
#     hf = h5py.File("Voltage_Control_out.h5", "w")
#     hf.create_dataset("sim_out", data=DATA_OUT)
#     hf.close()