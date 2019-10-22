import numpy as np
import math
#from scipy.optimize import fsolve
import scipy.misc
from scipy import optimize
from scipy import integrate
import matplotlib.pyplot as plt

def w_for_normalise(x,my):
    return np.cos(my*x)**2

def w(x, my):
    return np.cos(my*x)/ np.sqrt(integrate.quad(w_for_normalise,0,1,args=(my))[0])
def d_w (x,my):
    return -my*np.sin(my*x) / np.sqrt(integrate.quad(w_for_normalise,0,1,args=(my))[0])

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
            A[j][i] = A_ij(my[i],my[j])
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

def create_result_coefficients(n):
    return result_SLAR(n), create_n_my(n)

def y(x,n,c,my):
    res = 0
    for i in range(n):
        res = res + c[i]*w(x,my[i])
    return res

def first_dim (x,n,c,my):
    return p(x)*scipy.misc.derivative(y,x,n=1,args=(n,c,my), dx = 1e-8)


def diffyr(x,n,c,my):
    return -1*scipy.misc.derivative(first_dim,x,n=1,args=(n,c,my), dx = 1e-8)+ a(x)*scipy.misc.derivative(y,x,n=1,args=(n,c,my), dx = 1e-8) + q(x)*y(x,n,c,my) - f(x)
def first_condition(n,c,my):
    return h1*y(0,n,c,my) - h2*scipy.misc.derivative(y,0,n=1,args=(n,c,my), dx = 1e-8)
def second_condition(n,c,my):
    return H1*y(1,n,c,my) + H2*scipy.misc.derivative(y,1,n=1,args=(n,c,my), dx = 1e-8)


def test_y (n):
    c, my = create_result_coefficients(n)
    print(c,my)
    # print(create_n_my(n))
    # print("first condition: ", first_condition(n,c,my))
    # print("second condition: ", second_condition(n,c,my))

    # mass = np.linspace(0, 1)

    # for i in range(mass.size):
    #     print(mass[i],'  ',diffyr(mass[i],n,c,my))

    xx = np.linspace(0,1,100)
    yy = np.zeros(100)
    for i in range(100):
        yy[i] = y(xx[i],n,c,my)
    fig,ax = plt.subplots()
    
    ax.plot(xx,yy)
    plt.show()

def create_data_for_plot( N: int, x: np.ndarray):
    c, my = create_result_coefficients(N)
    yy = np.zeros(x.shape[0])
    for i in range(x.shape[0]):
        yy[i] = y(x[i],N,c,my)
    return yy


if __name__ == '__main__':
    #print(result_SLAR(10))
    xx = np.linspace(0,1,100)
    yy = create_data_for_plot(10,xx)
    
    fig,ax = plt.subplots()
    #plt.axis('scaled')
    
    ax.plot(xx,yy)
    plt.show()






