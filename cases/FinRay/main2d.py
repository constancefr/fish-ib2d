import numpy as np
import os
import sys

# IB2d solver comes from the external/IB2d submodule (repo root is two levels up)
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(_REPO_ROOT, 'external', 'IB2d', 'pyIB2d', 'IBM_Blackbox'))
import IBM_Driver as Driver
from please_Initialize_Simulation import please_Initialize_Simulation

def give_Me_input2d_Parameters():
    ''' Function to read in input files from "input2d"
    
    Returns:
        params: ndarray of parameter values
        struct_name: name of structure'''
        
    filename = 'input2d' #Name of file to read in
    
    #This is more sophisticated than what was here before the port.
    #If desired, names could be double checked for consistancy.
    names,params = np.loadtxt(filename,dtype={'names': ('param','value'),\
            'formats': ('S25','f8')},comments=('%','string_name'),\
            delimiter='=',unpack=True)
            
    #now get string_name
    with open('input2d','r') as f:
        for line in f:
            #look for 'string_name'. make sure it comes before a comment char
            if ('string_name' in line) and \
                line.find('string_name') < line.find('%'):
                #split the line by whitespace
                words = line.split()
                #look for the equals sign
                for n,word in enumerate(words):
                    if word == '=':
                        struct_name = words[n+1]
                        break
                break

    return (params,struct_name)

def main2d():
    
    '''This is the "main" function, which gets called to run the 
    Immersed Boundary Simulation. It reads in all the parameters from 
    "input2d", and sends them off to the "IBM_Driver" function to actually 
    perform the simulation.'''
    
    # OLD FORMAT #
    # READ-IN INPUT2d PARAMTERS #
    # params,struct_name = give_Me_input2d_Parameters()
    
    # NEW FORMAT #
    # READ-IN INPUT2d PARAMETERS #
    Fluid_Params, Grid_Params, Time_Params, Lag_Struct_Params, Output_Params, Lag_Name_Params = please_Initialize_Simulation()    


    #-#-#-# DO THE IMMERSED BOUNDARY SOLVE!!!!!!!! #-#-#-#
    #[X, Y, U, V, xLags, yLags] = Driver.main(struct_name, mu, rho, grid_Info, dt, T_final, model_Info)
    
    #For debugging only!
    #Driver.main(struct_name, mu, rho, grid_Info, dt, T_final, model_Info)
    Driver.main(Fluid_Params,Grid_Params,Time_Params,Lag_Struct_Params,Output_Params,Lag_Name_Params)
    
if __name__ == "__main__":
    main2d()