import numpy as np
import math
#from scipy.optimize import fsolve
import scipy.misc
from scipy import optimize
from scipy import integrate
import matplotlib.pyplot as plt


def phi(x,i, N, h):
    if 0<i<N :
        if h*(i-1) <=x<=h*i :
            return (x-(i-1)*h)/h
        elif i*h <= x<=(i+1)*h:
            return ((i+1)*h - x)/h
        else:
            return 0
    elif i==0:
        if 0<=x <h:
            return (h-x)/h
        else:
            return 0
    else:
        if (N-1)*h < x <= N:
            return (x - (N-1)*h)/h
        else:
            return 0


def d_phi(x,i, N, h):
    if 0<i<N :
        if h*(i-1) <=x<=h*i :
            return 1/h
        elif i*h <= x<=(i+1)*h:
            return -1/h
        else:
            return 0
    elif i==0:
        if 0<=x <h:
            return -1/h
        else:
            return 0
    else:
        if (N-1)*h < x <= N:
            return 1/h
        else:
            return 0


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


def right_function(x, i,j,N):
    return p(x)*d_phi(x,i,N, 1/N)*d_phi(x,j,N,1/N) + a(x)*d_phi(x,i, N, 1/N)*phi(x,j, N, 1/N) + q(x)*phi(x,i,N, 1/N)*phi(x,j, N, 1/N)
def left_function(x,j,N):
    return f(x)*phi(x,j, N, 1/N)

def A_ij (i,j,N):
    res = integrate.quad(right_function,max((i-1)/N,0),min(1,(i+1)/N),args=(i,j, N),limit = 200)[0] + p(1)*phi(1,i, N, 1/N)*phi(1,j, N,1/N)*H1/H2\
         + p(0)*phi(0,i, N, 1/N)*phi(0,j, N, 1/N)*h1/h2
    print(integrate.quad(right_function,max((i-1)/N,0),min(1,(i+1)/N),args=(i,j, N), limit = 100)[0],'  ', i, j)
    #print(p(1)*phi(1,i, N, 1/N)*phi(1,j, N,1/N)*H1/H2 + p(0)*phi(0,i, N, 1/N)*phi(0,j, N, 1/N)*h1/h2,'  ', i, j)
    return res

#max((i-1)/N,0),min(1,(i+1)/N)

def create_matrix_A_F (N):
    A = np.zeros((N+1,N+1))
    F = np.zeros(N+1)

    for i in range(N+1):
        for j in range(1,-2,-1):
            if 0<=i-j <=N:
                A[i-j][i] = A_ij(i,i-j,N)
        # for j in range(N+1):
        #     A[j][i] = A_ij(i,j,N)

        F[i] = integrate.quad(left_function,max((i-1)/N,0),min(1,(i+1)/N),args=(i,N))[0]
    return A,F

def result_SLAR(n):
    A, F = create_matrix_A_F(n)
    c = np.zeros(n+1)
    c = np.linalg.solve(A,F)
    return c

def result_function (x, n):
    c = result_SLAR(n)
    res = 0
    for i in range(n+1):
        res += c[i]*phi(x,i,n,1/n)
    return res

def y(x,n,c):
    res = 0
    for i in range(n+1):
        res = res + c[i]*phi(x,i,n,1/n)
        #res = res + phi(x,i,n,1/n)
    return res


def create_data_for_plot( N: int, x: np.ndarray):
    c = result_SLAR(N)
    yy = np.zeros(x.shape[0])
    print(y(0,N,c))
    for i in range(x.shape[0]):
        yy[i] = y(x[i],N,c)
        #yy[i] = phi(x[i],N,N,1/N)
    return yy


def test_2_plots(N1, N2):
    xx = np.linspace(0,1,50)
    y1 = create_data_for_plot(N1,xx)
    y2 = create_data_for_plot(N2,xx)

    fig,ax = plt.subplots()
    ax.scatter(xx,y1, marker = 'X')
    ax.plot(xx,y2, color = 'red')
    #plt.legend("n1="+str(N1), "n2="+str(N2))
    plt.show()

if __name__ == '__main__':
    #print(result_SLAR(10))
    #print(right_function(0.1,2,0,10))
    # xx = np.linspace(0,1,100)
    # yy = create_data_for_plot(100,xx)
    
    # fig,ax = plt.subplots()
    # #plt.axis('scaled')
    
    # ax.plot(xx,yy)
    # plt.show()

    test_2_plots(50,100)

