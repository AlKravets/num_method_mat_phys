import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate


# R = 0.05
# чит метод 
R = 0.1
#
u_kp = 50

def return_to_norm(x,y):
    return x*R, (y*u_kp + u_kp)

lam =398
c = 380
ro =8900

def print_time(t, tay):
    return t*tay*R**2/(lam/c/ro)

def show_res(h, y):
    x = np.arange(0,1+h,h)
    x1, y1 = return_to_norm(x,y)
    fig,ax = plt.subplots()
    ax.scatter(x1,y1, marker = 'X', s = 1)
    
    plt.legend()
    plt.show()
    fig.savefig('photo3_1.jpg')

def show_all_res(h, tay, y_res, t_end):
    x = np.arange(0,1+h,h)
    y = np.array(y_res)
    x1,y1 = return_to_norm(x,y)
    
    fig,ax = plt.subplots()

    for i in range(y1.shape[0]):
        if i == y1.shape[0]-1:
            time = t_end
        else:
            time = i*40
        # ax.scatter(x1,y1[i], marker = 'X', s = 1,label = 't = {0:.4f}'.format(print_time(time,tay)) )
        ax.plot(x1,y1[i], label = 't = {0:.4f}'.format(print_time(time,tay)) )
    plt.legend()    
    plt.show()
    fig.savefig('photo3_1.jpg')




def x_i_plus_1_2 (i, h):
    return integrate.quad(lambda x : x**2, i*h , (i+1)*h)[0] / h
    # return ((i-0.5)*h)**2

def x_i (i,h):
    # return integrate.quad(lambda x : x**2, i*h , (i+1)*h)[0] / h
    return h*i


def create_matrix (h, tay):
    N = int(1/h)
    A = np.zeros((N+1,N+1))
    for i in range(1, N):
        A[i][i-1] = -0.5/h**2 * (x_i_plus_1_2(i-1, h))**2
        A[i][i] = x_i(i,h)**2/tay + 0.5/h**2 *((x_i_plus_1_2(i-1, h))**2 + (x_i_plus_1_2(i, h))**2)
        A[i][i+1] = -0.5/h**2 * (x_i_plus_1_2(i, h))**2
    
    A[N,N] = 1

    # A[0,0] = h/tay/2 + 1/h
    # A[0,1] = h/tay/2 -1/h

    # A[0][0] = 1/tay

    # A[0,0] = 1/tay + 0.5/h
    # A[0,1] = -0.5/h

    # Чит метод
    A[0,0]  = 1
    #

    return A

def result (h, tay):
    N = int(1/h)
    y_j = np.zeros(N+1)
    #y_j_1 = np.zeros(N+1)
    F = np.zeros(N+1)
    A = create_matrix(h, tay)
    y_res = [y_j.tolist(),]

    t = 0



    while t<10000 and y_j[int(0.5/h)]> -4/5:

        for i in range(1, N):
            F[i] = y_j[i-1]*(0.5/h**2 * (x_i_plus_1_2(i-1, h))**2) \
                + y_j[i]*(x_i(i,h)**2/tay - 0.5/h**2 *((x_i_plus_1_2(i-1, h))**2 + (x_i_plus_1_2(i, h))**2))\
                    + y_j[i+1]*(0.5/h**2 * (x_i_plus_1_2(i, h))**2)
        
        F[N] = -1
        # F[0] = y_j[0]*(h/tay/2 - 1/h) + y_j[1]*(h/2/tay + 1/h)
        # F[0] = 1/h*(y_j[1] - y_j[0]) + y_j[0]/tay


        # F[0] = y_j[0]*(1/tay - 0.5/h) + y_j[1]*(0.5/h)

        # Чит метод
        F[0] = -1
        # 

        y_j = np.linalg.solve(A,F)

        t+=1
        print(t, '  ', y_j[0])
        if t% 40 ==0:
            y_res.append(y_j.tolist())
    #show_res(h,y_j)
    y_res.append(y_j.tolist())
    show_all_res(h, tay,y_res, t)



if __name__ == '__main__':
    A = create_matrix(0.1,0.1)
    print(A[-1,:])
    result(0.01,0.01)
    print(print_time(100,0.01))


