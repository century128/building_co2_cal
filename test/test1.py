from scipy.optimize import Bounds
from scipy.optimize import LinearConstraint
from scipy.optimize import NonlinearConstraint
import numpy as np
from scipy.optimize import minimize
import matrix_cal

def cons_j(x):
    f =rosen(x)
    return

def rosen(x):
    return 100*(x[1]-x[0]**2)**2 + (1-x[0])**2
def cons_f(x):
    return [x[0]**2 + x[1], x[0]**2 - x[1]]
def cons_J(x):
    return [[2*x[0], 1], [2*x[0], -1]]
def cons_H(x, v):
    return v[0]*np.array([[2, 0], [0, 0]]) + v[1]*np.array([[2, 0], [0, 0]])
def jac(x):
    return [-400*x[0]*x[1] +2*x[0] + 400*x[0]**3-2+2*x[0],]

bounds = Bounds([0, -0.5], [1.0, 2.0])
linear_constraint = LinearConstraint([[1, 2], [2, 1]], [-np.inf, 1], [1, 1])
nonlinear_constraint = NonlinearConstraint(cons_f, -np.inf, 1, jac=cons_J, hess=cons_H)
x0 = np.array([0.5, 0])
res = minimize(rosen, x0, method='trust-constr', jac=cons_J, hess=cons_H,
               constraints=[linear_constraint, nonlinear_constraint],
               options={'verbose': 1}, bounds=bounds)