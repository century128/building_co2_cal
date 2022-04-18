from math import *
import numpy as np
# for i in range(4):
#     print(i)

# global a
# a = [1,2,3,4]
# def p():
#     for i in range(len(a)) :
#         print(a[i])
#         a[i] = i+1
# p()
# print("\n")
# print(a)
a = [-0.7, -0.5, 3.3, 0, -0.7]
sum = 0
for i in range(len(a)):
    sum += exp(a[i])
b = np.zeros(len(a))
for j in range(len(a)):
    b[j] = exp(a[j])/sum

print(b)
