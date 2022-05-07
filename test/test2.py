from scipy.optimize import Bounds
from scipy.optimize import LinearConstraint
from scipy.optimize import NonlinearConstraint
import numpy as np
from scipy.optimize import minimize
from matrix_cal import *
from scipy.optimize import BFGS
from scipy.optimize import SR1

def cons_jj(x):
    f =rosen(x)
    return jacobian_cal(x,f)

def rosen(x):
    return 100*(x[1]-x[0]**2)**2 + (1-x[0])**2
def cons_f(x):
    return [x[0]**2 + x[1], x[0]**2 - x[1]]
def cons_J(x):
    return [[2*x[0], 1], [2*x[0], -1]]
def cons_H(x):
    return np.array([[2, 0], [0, 0]])
ineq_cons = {'type': 'ineq',
             'fun' : lambda x: np.array([1 - x[0] - 2*x[1],
                                         1 - x[0]**2 - x[1],
                                        ]),
             }
ineq_cons2 = {'type': 'ineq',
             'fun' : lambda x: np.array([ 1 - x[0]**2 + x[1]]),
             }
eq_cons = {'type': 'eq',
           'fun' : lambda x: np.array([2*x[0] + x[1] - 1]),
           }

bounds = Bounds([0, -0.5], [1.0, 2.0])
# linear_constraint = LinearConstraint([[1, 2], [2, 1]], [-np.inf, 1], [1, 1])
# nonlinear_constraint = NonlinearConstraint(cons_f, -np.inf, 1,jac=cons_J, hess=BFGS())
# x0 = np.array([0.5, 0])
# res = minimize(rosen, x0, method='trust-constr', jac="2-point", hess=SR1(),
#                constraints=[linear_constraint, nonlinear_constraint],
#                options={'verbose': 1}, bounds=bounds)
# print(res.x)
x0 = np.array([-200, 100])
res = minimize(rosen, x0, method='SLSQP', jac="2-point",
               constraints=[eq_cons, ineq_cons,ineq_cons2], options={'ftol': 1e-9, 'disp': True},
               )
print(res.x)