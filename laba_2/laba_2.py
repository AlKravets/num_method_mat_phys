import numpy as np
import math
#from scipy.optimize import fsolve
import scipy.misc
from scipy import optimize
from scipy import integrate
import matplotlib.pyplot as plt

import laba_1

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
    # return 0

def f(x):
    return 2*x**2 + np.sin(2*x)


## Изменим вид диф. Уравнения для работы 2 метода,для 1 метода форма не важна

def mu(x):
    # return np.exp(x + (4 * np.arctan((1-2*np.tan(math.pi*x/2))/3**0.5))/(3**0.5 * math.pi)\
    #     +np.log(2-np.sin(math.pi *x)))
    return np.exp(x + (4 * np.arctan((1-2*np.tan(math.pi*x/2))/3**0.5))/(3**0.5 * math.pi))

def p_1 (x):
    # return -1*(np.sin(math.pi*x)-2)*mu(x)
    if x< 1/2**0.5:
        return -1*(np.sin(math.pi*x)-2)*mu(x)
    else: return -1*(np.sin(math.pi*x)-2)*mu(x)+3
    
    
def p_1_for_integrate(x):
    return 1/p_1(x)

def q_1(x):
    return q(x)*mu(x)
    #  return q(x)

def f_1(x):
    return f(x)*mu(x)

def a_1 (x):
    return 0

def right_function(x, i,j,N):
    return p_1(x)*d_phi(x,i,N, 1/N)*d_phi(x,j,N,1/N) + a_1(x)*d_phi(x,i, N, 1/N)*phi(x,j, N, 1/N) + q_1(x)*phi(x,i,N, 1/N)*phi(x,j, N, 1/N)
def left_function(x,j,N):
    return f_1(x)*phi(x,j, N, 1/N)

def A_ij (i,j,N):
    res = integrate.quad(right_function,max((i-1)/N,0),min(1,(i+1)/N),args=(i,j, N),limit = 200)[0] + p_1(1)*phi(1,i, N, 1/N)*phi(1,j, N,1/N)*H1/H2\
         + p_1(0)*phi(0,i, N, 1/N)*phi(0,j, N, 1/N)*h1/h2
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
    xx = np.linspace(0,1,100)
    y1 = create_data_for_plot(N1,xx)
    y2 = create_data_for_plot(N2,xx)
    c1 = result_SLAR(N1)
    c2 = result_SLAR(N2)

    fig,ax = plt.subplots()
    ax.scatter(xx,y1, marker = 'X', s = 5, label = '50 точек')
    ax.plot(xx,y2, color = 'red', label = '100 точек')
    #ax.plot(xx,y1, color = 'blue')
    plt.legend()
    plt.show()
    fig.savefig('photo3.jpg')

    x1 = np.linspace(0,1,11)
    for i in range(x1.shape[0]):
        print('{0:.2f} & {1:.6f} \\'.format(x1[i], abs(y(x1[i],N1,c1)-y(x1[i],N2,c2))))




# Метод 2





def _p_ (x_1,x_2,h):
    return (integrate.quad(p_1_for_integrate,x_1,x_2)[0]/h)**-1

def _q_ (x_1,x_2,h):
    return integrate.quad(q_1,x_1,x_2)[0]/h

def _f_ (x_1,x_2,h):
    return integrate.quad(f_1,x_1,x_2)[0]/h

def _p_0_N (x_1,x_2,h):
    return 2*_p_ (x_1,x_2,h)

def _q_0_N (x_1,x_2,h):
    return 2*_q_ (x_1,x_2,h) 


def _f_0_N (x_1,x_2,h):
    return 2*_f_ (x_1,x_2,h)

alfa_1 = p_1(0)*h1/h2
alfa_2 = p_1(1)*H1/H2

def create_SLAR_for_2_metod(x:np.ndarray):
    N = x.shape[0]-1
    h = 1/N
    u =np.zeros(N+1)
    F =np.zeros(N+1)
    A = np.zeros((N+1,N+1))


    for i in range(1, N):
        A[i][i-1] = -1* _p_(h*(i-1), h*(i), h)/h**2 # p_i
        A[i][i] = (_p_(h*(i), h*(i+1), h)+_p_(h*(i-1), h*(i), h))/h**2 +_q_(h*(i-0.5),h*(i+0.5),h)  # p_1+1 + p_i
        A[i][i+1] = -1*_p_(h*(i), h*(i+1), h)/h**2 #p_i+1
        F[i] = _f_(h*(i-0.5),h*(i+0.5),h)


    ##
    A[0][1] = -_p_(0,h,h )/h
    A[0][0] = _p_(0,h,h )/h +alfa_1 + h/2*_q_0_N(0,0.5*h,h)
    F[0] = h/2*_f_0_N(0,0.5*h,h)


    ##
    A[N][N-1] =  -_p_(h*(N-1), h*N, h)/h
    A[N][N] = _p_(h*(N-1), h*N, h)/h + alfa_2 + h/2*_q_0_N(h*(N-0.5), h*N, h)
    F[N] = h/2*_f_0_N(h*(N-0.5), h*N, h)

    u = np.linalg.solve(A,F)
    #print(A)
    return u

def test_metod_2(x):
    u = create_SLAR_for_2_metod(x)

    #y1 = laba_1.create_data_for_plot(10,x)

    y1 = create_data_for_plot(x.shape[0],x)
    fig,ax = plt.subplots()
    ax.scatter(x,u, color = 'red', label = 'Интегро-интер. метод')
    ax.plot(x,y1, color = 'blue', label = 'Метод конечных разниц')
    plt.legend()
    plt.show()
    #fig.savefig('photo2_1.jpg')

def delta_3_metods(x):
    y_0 = laba_1.create_data_for_plot(10, x)
    y_1 = create_data_for_plot(x.shape[0],x)
    y_2 = create_SLAR_for_2_metod(x)

    for i in range(0,x.shape[0],10):
        delta_0_1 = abs(y_0[i] - y_1[i])
        delta_0_2 = abs(y_0[i] - y_2[i])
        delta_1_2 = abs(y_2[i] - y_1[i])
        print('{0:.2f} & {1:.6f} & {2:.6f} & {3:.6f} & {4:.6f} & {5:.6f} & {6:.6f}'.format(x[i], y_0[i], y_1[i], y_2[i], delta_0_1, delta_0_2, delta_1_2))



def test_2_plots_2_metods():
    x1 = np.linspace(0,1,51)
    x2 = np.linspace(0,1,101)
    y1 = create_SLAR_for_2_metod(x1)
    y2 = create_SLAR_for_2_metod(x2)
    

    fig,ax = plt.subplots()
    ax.scatter(x1,y1, marker = 'X', s = 5, label = '50 точек')
    ax.plot(x2,y2, color = 'red', label = '100 точек')
    #ax.plot(xx,y1, color = 'blue')
    plt.legend()
    plt.show()
    fig.savefig('photo3_1.jpg')

    x1 = np.linspace(0,1,11)
    for i in range(x1.shape[0]):
        print('{0:.2f} & {1:.6f} \\'.format(x1[i], abs(y1[i*5] - y2[i*10])))




if __name__ == '__main__':
    #print(result_SLAR(10))
    #print(right_function(0.1,2,0,10))
    # xx = np.linspace(0,1,100)
    # yy = create_data_for_plot(100,xx)
    
    # fig,ax = plt.subplots()
    # #plt.axis('scaled')
    
    # ax.plot(xx,yy)
    # plt.show()

    #test_2_plots(50,100)

    x = np.linspace(0,1,101)
    test_metod_2(x)
    #print(p_1(0.3))
    #delta_3_metods(x)

    #test_2_plots_2_metods()



