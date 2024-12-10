import pandas as pd
import math
from scipy.integrate import quad
import numpy as np
import matplotlib.pyplot as plt

adjustment_dict = {
    "1": [-4.18, -0.37, -9.31],
    "2": [2.16, -1.20, 2.66],
    "3": [2.02, 1.57, 6.65],
    "2006": [-0.37, 0.72, 0.83],
    "2007": [-0.59, 1.00, 1.23],
    "2008": [-0.31, 0.47, 0.98],
    "2009": [-0.24, 0.24, 0.60],
    "2010": [-0.11, -0.14, 0.31],
    "2011": [0.10, -0.58, -0.56],
    "2012": [0.33, -0.71, -0.83],
    "2013": [0.27, -0.51, -0.73],
    "2014": [0.12, 0.069, -0.37],
    "2015": [0.28, -0.12, -0.68],
    "2016": [0.52, -0.44, -0.78],
    "January": [-0.46, 1.81, 3.13],
    "February": [0.18, 0.76, 2.43],
    "March": [1.05, -0.77, 1.04],
    "April": [1.58, -2.03, -0.56],
    "May": [1.49, -2.47, -1.95],
    "June": [0.74, -2.01, -2.75],
    "July": [-0.41, -0.81, -2.68],
    "August": [-0.96, 0.11, -2.06],
    "September": [-1.08, 0.78, -1.08],
    "October": [-0.85, 1.20, 0.27],
    "November": [-0.63, 1.45, 1.51],
    "December": [-0.65, 1.98, 2.70],
    "Appalachian": [-0.22, -0.042, -0.89],
    "Corn Belt": [0.55, -0.58, -1.12],
    "Delta": [-2.56, 0.59, 1.47],
    "Lake": [0.61, -0.40, -0.64],
    "Mountain": [-0.96, 3.13, 1.50],
    "Northeast": [1.04, -1.99, -1.13],
    "Northern Plains": [-0.26, 0.19, -0.79],
    "New York": [0.67, -1.21, -0.45],
    "Pennsylvania": [1.15, -0.96, 0.06],
    "Southeast": [-2.00, 2.59, 2.60],
    "Southern Plains": [-0.51, -1.02, -0.93],
    "West Coast": [1.09, 0.53, 0.52],
    "Wisconsin": [1.4, -0.83, -0.2],
    "2x/d": [-0.74, 0.090, 0.15],
    "3x/d": [0.74, -0.090, -0.15]
}

def calculate_lactation_group_yield(AHMP, num_milking_cows, parity_structure_pct = (0.386, 0.281, 0.333)):
    herdM305 = AHMP * 305 / (365 * num_milking_cows) 
    P1_pct, P2_pct, P3_pct = parity_structure_pct

    # 3 equations solve 3 unknown variables. Constant 1.18 and 1.25 from Manfei lact curve paper.
    P1_M305  = herdM305 / (P1_pct + 1.18 * P2_pct + 1.25 * P3_pct)
    P2_M305 = P1_M305 * 1.18
    P3_M305 = P1_M305 * 1.25

    # Returning the results as a tuple
    return (P1_M305, P2_M305, P3_M305)

def get_t_values(start=1, end=305, step=1):
    return np.arange(start, end, step)

def get_y_values_wood_curve(t, parameter_a, parameter_b, parameter_c):
    return parameter_a * np.power(t, parameter_b) * np.exp(-1 * parameter_c * t)

def calc_integral_wood_curve(parameter_a, parameter_b, parameter_c):
    result, _ = quad(get_y_values_wood_curve, 1, 305, args=(parameter_a, parameter_b, parameter_c))
    return result

def plot_graph(t, y):
    plt.plot(t, y)
    plt.xlabel("DIM (d)")
    plt.ylabel("Milk yield (kg/d)")
    plt.grid(linestyle=':')
    plt.show()

def plot_from_parameters(parameter_a, parameter_b, parameter_c):
    t = get_t_values()
    y = get_y_values_wood_curve(t, parameter_a, parameter_b, parameter_c)
    return t, y

def get_wood_parameters(lactation_group = None, year = None, month = None, region = None, milking_frequency = None, MY_305d = None):
    parameter_a = 19.9
    parameter_b = 24.7 * 1e-2
    parameter_c = 33.76 * 1e-4
    t = get_t_values()

    for category in [lactation_group, year, month, region, milking_frequency]:
        if category:
            parameter_a += adjustment_dict[category][0]
            parameter_b += adjustment_dict[category][1] * 1e-2
            parameter_c += adjustment_dict[category][2] * 1e-4

    if MY_305d == None:
        MY_305d = calc_integral_wood_curve(parameter_a, parameter_b, parameter_c)
        # y = get_y_values_wood_curve(t, parameter_a, parameter_b, parameter_c)
        # return parameter_a, parameter_b, parameter_c, MY_305d
        print("parameter_a:", parameter_a,
              "\nparameter_b:", parameter_b,
              "\nparameter_c:", parameter_c,
              "\nMY_305d:", MY_305d)
        return parameter_a, parameter_b, parameter_c
    else:
        min_diff = float('inf')
        parameter_a_best_etimate = 0
        MY_305d_best_etimate = 0
        for parameter_a_error in np.arange(-20, 20, 0.01):
            parameter_a_vary = parameter_a + parameter_a_error
            MY_305d_vary = calc_integral_wood_curve(parameter_a_vary, parameter_b, parameter_c)
            if abs(MY_305d_vary - MY_305d) < min_diff:
                min_diff =  abs(MY_305d_vary - MY_305d)
                parameter_a_best_etimate = parameter_a_vary
                MY_305d_best_etimate = MY_305d_vary
        print("parameter_a:", parameter_a_best_etimate,
              "\nparameter_b:", parameter_b,
              "\nparameter_c:", parameter_c,
              "\nMY_305d:", MY_305d_best_etimate)
        return parameter_a_best_etimate, parameter_b, parameter_c


def method1(year, region, milking_frequency):
    print("====METHOD1====")
    print(">> Parity 1:")
    p1_parameter_a, p1_parameter_b, p1_parameter_c = get_wood_parameters(lactation_group = '1', year = year, region = region, milking_frequency = milking_frequency)
    p1_t, p1_y = plot_from_parameters(p1_parameter_a, p1_parameter_b, p1_parameter_c)
    print(">> Parity 2:")
    p2_parameter_a, p2_parameter_b, p2_parameter_c = get_wood_parameters(lactation_group = '2', year = year, region = region, milking_frequency = milking_frequency)
    p2_t, p2_y = plot_from_parameters(p2_parameter_a, p2_parameter_b, p2_parameter_c)
    print(">> Parity 3+:")
    p3_parameter_a, p3_parameter_b, p3_parameter_c = get_wood_parameters(lactation_group = '3', year = year, region = region, milking_frequency = milking_frequency)
    p3_t, p3_y = plot_from_parameters(p3_parameter_a, p3_parameter_b, p3_parameter_c)
    plt.plot(p1_t, p1_y, label = "Parity 1") 
    plt.plot(p2_t, p2_y, label = "Parity 2") 
    plt.plot(p3_t, p3_y, label = "Parity 3") 
    plt.legend() 
    plt.show()

def method2(year, region, milking_frequency, AHMP, parity_structure_pct, num_milking_cows):
    P1_M305, P2_M305, P3_M305 = calculate_lactation_group_yield(AHMP, num_milking_cows, parity_structure_pct)
    print("====MODEL2====")
    print(">> Parity 1:")
    p1_parameter_a_best_etimate, p1_parameter_b, p1_parameter_c = get_wood_parameters(lactation_group = '1', year = year, region = region, milking_frequency = milking_frequency, MY_305d = P1_M305)
    p1_t, p1_y = plot_from_parameters(p1_parameter_a_best_etimate, p1_parameter_b, p1_parameter_c)
    print(">> Parity 2:")
    p2_parameter_a_best_etimate, p2_parameter_b, p2_parameter_c = get_wood_parameters(lactation_group = '2', year = year, region = region, milking_frequency = milking_frequency, MY_305d = P2_M305)
    p2_t, p2_y = plot_from_parameters(p2_parameter_a_best_etimate, p2_parameter_b, p2_parameter_c)
    print(">> Parity 3+:")
    p3_parameter_a_best_etimate, p3_parameter_b, p3_parameter_c = get_wood_parameters(lactation_group = '3', year = year, region = region, milking_frequency = milking_frequency, MY_305d = P3_M305)
    p3_t, p3_y = plot_from_parameters(p3_parameter_a_best_etimate, p3_parameter_b, p3_parameter_c)
    plt.plot(p1_t, p1_y, label = "Parity 1") 
    plt.plot(p2_t, p2_y, label = "Parity 2") 
    plt.plot(p3_t, p3_y, label = "Parity 3") 
    plt.legend() 
    plt.show()

### METHOD1 ###
milking_frequency = '3x/d' # Either '2x/d' or '3x/d'
region = 'New York' 
# method1('2016', region, milking_frequency)

### METHOD2 ###
parity_structure_pct = (0.359525, 0.263485, 0.3769925)
num_milking_cows = 1284
AHMP = 18405332   # Unit: kg
method2('2016', region, milking_frequency, AHMP, parity_structure_pct, num_milking_cows)
