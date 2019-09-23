import numpy as np
import math
#from scipy.optimize import fsolve
#import scipy
from scipy import optimize
from scipy import integrate

def w(x, my):
    return np.cos(my*x)
def d_w (x,my):
    return -my*np.sin(my*x)

def p(x):
    return 2- np.sin(math.pi*x)

def q(x):
    return 5
h1, h2 = 0, 1
H1, H2 = 1, 4

def a(x):
    return np.sin(math.pi*x)

def f(x):
    return 2*x**2 + np.sin(2*x)

def  get_my(x):
    return np.tan(x) - 1/4/x

def create_n_my (n):
    res = np.zeros(n)
    for i in range(n):
        res[i] = optimize.fsolve(get_my,0.5+math.pi*i)
    return res

def right_function(x, i,j):
    return p(x)*d_w(x,i)*d_w(x,j) + a(x)*d_w(x,i)*w(x,j) + q(x)*w(x,i)*w(x,j)
def left_function(x,j):
    return f(x)*w(x,j)

def A_ij (i,j):
    res = integrate.quad(right_function,0,1,args=(i,j))[0] + p(1)*w(1,i)*w(1,j)*H1/H2 + p(0)*w(0,i)*w(0,j)*h1/h2
    return res

def create_matrix_A_F (n):
    A = np.zeros((n,n))
    F = np.zeros(n)
    my = create_n_my(n)
    for i in range(n):
        for j in range(n):
            A[i][j] = A_ij(my[i],my[j])
        F[i] = integrate.quad(left_function,0,1,args=(my[i]))[0]
    return A,F

def result_SLAR(n):
    A, F = create_matrix_A_F(n)
    c = np.zeros(n)
    c = np.linalg.solve(A,F)
    return c

def result_function (x, n):
    c = result_SLAR(n)
    my = create_n_my(n)
    res = 0
    for i in range(n):
        res += c[i]*w(x,i)
    return res


if __name__ == '__main__':
    #print(result_SLAR(10))

    x = np.arange(0.1,1,0.1)
    print(result_function(x, 3)  - result_function(x, 4))




