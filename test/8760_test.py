import numpy as np
from sympy import *
from matrix_cal import *
from scipy.optimize import minimize
from scipy.optimize import SR1

def cal_12(x):
    a = np.zeros(6)
    for i in range(6):
        a[i] = x[i] - x[i+1]*x[i+1]
    return a

para = [1,2,3,4,5,6,1]
b = cal_12(para)
print(b)