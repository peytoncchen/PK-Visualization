###
# DO NOT RUN... I mean you theoretically could but there are things that are pyodide specific
# i.e. interfacing with Javascript objects and functions. This is just a place to document the functions
# and provide comments as to what they do that isn't in the middle of function call in a 1000-line HTML file.
###

import numpy as np
from scipy.optimize import brute
import math
import js

# Importing all the variables needed to make the calculations we need on the python side.
in_csv_data = js.in_csv_data
js_table = js.k_table2
dydt2 = js.dydt2
represents = js.represents

# Interfacing between Python and Javascript is weird. So to avoid any Typeerrors "Jsproxy" is not compatible with
# <insert random Python type> I just copy the whole k_table array to make it a Python array.
def copy_array(array):
    copy = [[],[],[]]
    copy[0] = [[lst[0], lst[1]] if lst else [] for lst in array[0]]
    copy[1] = [num for num in array[1]]
    copy[2] = [num for num in array[2]]
    return copy

k_table2 = copy_array(js_table)

# Turns the CSV data array into a dictionary that we can manipulate.
def arr_to_dic(arr):
    result = {}
    for i in range(len(arr)):
        if arr[i][0] not in result and len(arr[i]) != 1:
            result[arr[i][0]] = [arr[i][1]]
        elif len(arr[i]) != 1:
            result[arr[i][0]].append(arr[i][1])
    return result

dic = arr_to_dic(in_csv_data)
print(dic)

# Sum of squared errors, least squares function that is called from scipy.optimize
def SSE(fit_dic, actual_dic):
    """
    input:
    fit_dic is a dictionary {time: value, time2: value2}
    actual_dic is a dictionary for the actual values i.e. {0:[0,0,0,0,0],... 90:[0.2,0.4,0.6,0.8,0.3]}
    --------
    output: 
    Sum of squared errors
    """
    SSE = 0
    for time in fit_dic: # iterates through the keys for time values
        fit_conc = fit_dic[time]
        for actual_value in actual_dic[time]: # get the actual value
            diff = actual_value - fit_conc
            SSE += math.sqrt(diff**2) # add to SSE
    return SSE

total_time = max(list(dic))
timepoints = 100 # Arbitrary value

# Constructs time_measured which is a list of times measured from the CSV input and 
# constructs time_index which is a list of all indexes generated from np.linspace that match up to the index 
# corresponding to the times from the actual dictionary.
def construct_time(dic):
    time_range = np.linspace(0,total_time,total_time*timepoints+1)
    time_measured = list(dic)

    time_range_list = list(time_range)
    time_index = []  # index of the correct time points
    for timepoint in time_measured:
        for index in range(len(time_range_list)):
            if time_range_list[index] >= timepoint: # accounts for floating point weirdness
                time_index.append(index)
                break
    return time_measured, time_index


# Makes dictionary to fit to in SSE.
def make_fit_dic(solve, actual_dic):
    """
    input:
    solve: output from ordinary differential equation (ode) solver for given ks
    actual_dic: actual dictionary constructed from CSV data inputted
    --------
    output:
    a normalized dictionary of the fit values at the correct timepoints based on the input data
    """
    returnDic = {}
    timemeasured,timeindex = construct_time(actual_dic)
    plas_solve = []
    for i in range(len(solve)):
        plas_solve.append(solve[i][represents - 1]) # Only data for the compartment CSV data represents
    non_normal_conc = [] 
    for index in timeindex:
        non_normal_conc.append(plas_solve[index])
    max_conc = max(non_normal_conc)
    normal_conc = [x/max_conc for x in non_normal_conc] # Normalize
    for i in range(len(timeindex)):
        returnDic[timemeasured[i]] = normal_conc[i] # Populate the return dictionary
    return returnDic

### Global variables for our Integrator
y0 = k_table2[1][:]
k = []
dt = total_time/(total_time*timepoints+1)
###

# From the given k, integrates with the same 4th order Runge-Kutta method and returns the SSE value.
def SSEfromK(ks):
    global k # Accesses global variable k
    k = ks
    if k_table2[2]:
        for i in range(len(k_table2[2])):
            entry = k_table2[2][i]
            k[entry[1] - 1] = entry[0]
    
    integrator = js.Integrator.new(y0, dydt2, 0, dt)
    integrator.steps(total_time*timepoints+1)
    full_dataset2 = integrator.y_collection

    fit = make_fit_dic(full_dataset2, dic)
    return SSE(fit, dic)

# Given the k_table array processes the array and returns a tuple of the ranges we want to brute for
# each k-value.
def process_bounds(array):
    result = []
    bounds = array[0]
    k_try = array[2]
    index_k_try = 0
    for i in range(len(bounds)):
        if bounds[i]:
            result.append((bounds[i][0], bounds[i][1]))
        else:
            result.append((k_try[index_k_try][0], k_try[index_k_try][0]))
            index_k_try += 1
    return tuple(result)


bounds_tup = tuple(process_bounds(k_table2))
final_result = list(brute(SSEfromK, bounds_tup, Ns=2))
print(final_result)