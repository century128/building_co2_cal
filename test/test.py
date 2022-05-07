import numpy as np
from sympy import *
from matrix_cal import *
from scipy.optimize import minimize

def objective(x):
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]
    x4 = x[3]
    ob,cb,db = 0,0,0
    for i in range(2):
        ob += x[i]
    for i in range(2):
        cb += x[2+i]
    db = cb+ob
    return db

def constraint1(x):
    return  x[0] * x[1] * x[2] * x[3] -25

def constraint2(x):
    sum_sq = 40
    for i in range(4):
        sum_sq = sum_sq - x[i] ** 2
    return sum_sq

def cal_12(x):
    a = np.zeros(4)
    for i in range(3):
        a[i] = x[i] -2* x[i+1]
    return a

# x0 = np.zeros(4)
# for i in range(4):
#     x0[i] += np.random.randn()
x0 = np.zeros(4)
for i in range(4):
    x0[i] += np.random.randn()
b = (1.0 , 5.0)
bnds = (b,b,b,b)
cons1 = {'type' : 'ineq' , 'fun' : constraint1}
cons2 = {'type' : 'eq' , 'fun' : constraint2}
cons3 = {'type' : 'ineq' , 'fun' :cal_12}
cons = [cons1,cons2,cons3]

sol = minimize(objective, x0 = x0,constraints=cons)
print(sol)
