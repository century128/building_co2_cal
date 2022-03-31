from sympy import *
import numpy as np

def para_to_symbol(x):
    for i in range(len(x)):
        x[i] = Symbol(x[i])
    return x

def jacobian_cal(vars,func ):
    #vars:输入['x','y','z',...]形式的自变量列表
    #func:输入"func"形式函数

    f = Matrix(func)
    v = Matrix(vars)
    return np.array(f.jacobian(v)).reshape(len(vars))

def hessian_cal(vars,func):
    # vars:输入"var var ……"形式的自变量字符串
    # func:输入"func"形式函数

    f = sympify(func)
    H = zeros(len(vars),len(vars))

    for i, fi in enumerate(f):
        for j, r in enumerate(vars):
            for k, s in enumerate(vars):
                H[j, k] = diff(diff(fi, r), s)
    return np.array(H)
