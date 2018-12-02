import matplotlib.pyplot as plt
import numpy as np
from cvxopt import matrix, solvers
import warnings
import time
from sklearn.preprocessing import normalize
warnings.filterwarnings("ignore")


def IHT(A, b, K, epislon = 1e-6):
    m, N = np.shape(A)
    A_ = normalize(A, axis = 0)
    A_norm = 1 / np.linalg.norm(A, axis = 0)
    W = np.zeros(shape = (N, N))
    np.fill_diagonal(W, A_norm)
    xk = np.zeros(N)
    
    for _ in range(100):
        x = xk.copy()
        r = b - np.dot(A_, x)
        xk = x + np.dot(A_.T, r)
        y = abs(xk)
        I = np.argsort(y)[::-1][K:]
        xk[I] = np.zeros(len(I))
        if np.linalg.norm(r) < epislon:
            return np.dot(W, xk)
        
    return np.dot(W, xk)   
## --------- Function of OMP --------- ##

def OMP(A, b, epislon = 1e-6):
    
    m, N = np.shape(A)
    
    ''' Initialize '''
    index_set = []
    rk = b.copy()
    ''''''''''''''''''
    
    while True:
        r = rk.copy()
        corr = abs(np.dot(A.T, r))
        idx = np.argsort(corr)[-1]
        
        if idx not in index_set:
            index_set.append(idx)
        
        xk = np.linalg.lstsq(A[:, index_set], b)[0]
        rk = b - np.dot(A[:, index_set], xk)
        
        if np.linalg.norm(rk) < epislon: 
            break
    
    x_ = np.zeros(N)
    x_[index_set] = xk
    
    return x_    

## ------------------------------------ ##

## -------- Fucntion of CoSaMP -------- ##

def CoSaMP(A, y, K, alpha):
    
    M, N = np.shape(A)
    theta = np.zeros(N) #
    Pos_theta = [] 
    r_n = y 
    theta_ls = None
    for kk in range(K): 
        # (1) Identification
        product = np.dot(A.T, r_n) 
        pos = np.argsort(abs(product))[::-1]
        Js = pos[:alpha*K] 
        # (2) Support Merger
        Is = np.union1d(Pos_theta, Js).astype('int') 
        # (3) Estimation
        
        if len(Is) <= M:
            At = A[:,Is] 
        else: 
            if kk == 1:
                theta_ls = 0
            break
        
        
        theta_ls = np.linalg.lstsq(At, y)[0] 
        
        # (4) Pruning
        pos = np.argsort(abs(theta_ls))[::-1]
        
        # (5) Sample Update
        
        Pos_theta = Is[pos[:K]];
        theta_ls = theta_ls[pos[:K]];
        r_n = y - np.dot(At[:,pos[:K]], theta_ls)
        if np.linalg.norm(r_n) < 1e-6: 
            break

    
    theta[Pos_theta] = theta_ls 
    
    return theta

## ------------------------------------ ##

    
## --------- Function of BP --------- ##

def Basis_pursuit(A, b):
    
    weights, N = np.shape(A)
    
    P = matrix(np.zeros(shape = (2*N, 2*N)))
    q = matrix(np.concatenate((np.zeros(shape = (N,1)), np.ones(shape = (N,1))), axis = 0))
    A_ = matrix(np.concatenate((A, np.zeros(shape = (weights, N))), axis = 1)) 
    I1 = np.concatenate((np.eye(N), -np.eye(N)), axis = 1)
    I2 = np.concatenate((-np.eye(N), -np.eye(N)), axis = 1)
    G = matrix(np.concatenate((I1, I2), axis = 0))
    h = matrix(np.zeros(2*N))
    b = matrix(b)
    
    x = solvers.qp(P, q, G, h, A_,b)['x'][:N]
    
    return x

## ---------------------------------- ##


## --------- Extra Function --------- ## 
def Generate_x(N, k): # k = nonzero numbers
    x = np.zeros(N)
    I = np.random.permutation(np.arange(0, N))[:k]
    
    for index in I:
        x[index] = np.random.randn()
    
    return x 
    

def check(x, x_, epislon = 1e-6):
    
    l = len(x)

    for i in range(l):
        if abs(x[i] - x_[i]) > epislon:
            return False
    return True    
    
## ---------------------------------- ##  
    


if __name__ == '__main__':
        
    accuracy_matrix_bp = np.zeros(71)
    accuracy_matrix_omp = np.zeros(71)
    accuracy_matrix_cosamp1 = np.zeros(71)
    accuracy_matrix_cosamp2 = np.zeros(71)
    accuracy_matrix_cosamp3 = np.zeros(71)
    accuracy_matrix_iht = np.zeros(71)
    N = 1000
    m = 500
    A = np.random.randn(m, N)
    bound = m/(2*np.log(N))
    
    bp_time_mtx = np.zeros(71)
    omp_time_mtx = np.zeros(71)
    cosamp1_time_mtx = np.zeros(71)
    cosamp2_time_mtx = np.zeros(71)
    cosamp3_time_mtx = np.zeros(71)
    iht_time_mtx = np.zeros(71)

    
    for cf_num in range(1, 71):
        
        count_bp = 0
        count_omp = 0
        count_cosamp1 = 0
        count_cosamp2 = 0
        count_cosamp3 = 0
        count_iht = 0
        
        bp_time = 0
        omp_time = 0
        cosamp1_time = 0
        cosamp2_time = 0
        cosamp3_time = 0
        iht_time = 0
        
        for _ in range(100):
            
            x = Generate_x(N, cf_num)
            b = np.dot(A, x)
            
#            '************ BP ***********'
#            bp_start = time.time()
#            x_bp = Basis_pursuit(A, b)
#            bp_end = time.time()
#            
#            if check(x, x_bp) == True:
#                count_bp += 1
#                
#            bp_time += (bp_end - bp_start)
#            '***************************'
#            
#            '*********** OMP ***********'
#            
#            omp_start = time.time()
#            x_omp = OMP(A, b)
#            omp_end = time.time()
#
#            if check(x, x_omp) == True:
#                count_omp += 1
#            
#            omp_time += (omp_end - omp_start)
#            '***************************'
#            
#            '********** CoSaMP2 *********'
#
#            cosamp2_start = time.time()
#            x_cosamp2 = CoSaMP(A, b, cf_num, alpha = 2)
#            cosamp2_end = time.time()
#
#            if check(x, x_cosamp2) == True:
#                count_cosamp2 += 1
#                
#            cosamp2_time += (cosamp2_end - cosamp2_start)
#            
#            '***************************'
#            
#            '********** CoSaMP 1 & 3 *********'
#            cosamp1_start = time.time()
#            x_cosamp1 = CoSaMP(A, b, cf_num, alpha = 1)
#            cosamp1_end = time.time()
#            if check(x, x_cosamp1) == True:
#                count_cosamp1 += 1
#            cosamp1_time += (cosamp1_end - cosamp1_start)
#    
#            cosamp3_start = time.time()
#            x_cosamp3 = CoSaMP(A, b, cf_num, alpha = 3)
#            cosamp3_end = time.time()
#            if check(x, x_cosamp3) == True:
#                count_cosamp3 += 1
#            cosamp3_time += (cosamp3_end - cosamp3_start)
#            '*********************************'
            
            
            '*********** IHT ***********'
            
            iht_start = time.time()
            x_iht = IHT(A, b, cf_num)
            iht_end = time.time()

            if check(x, x_iht) == True:
                count_iht += 1
            
            iht_time += (iht_end - iht_start)
            '***************************'
            
        accuracy_matrix_bp[cf_num] = count_bp / 100
        accuracy_matrix_omp[cf_num] = count_omp / 100
        accuracy_matrix_cosamp1[cf_num] = count_cosamp1 / 100
        accuracy_matrix_cosamp2[cf_num] = count_cosamp2 / 100
        accuracy_matrix_cosamp3[cf_num] = count_cosamp3 / 100
        accuracy_matrix_iht[cf_num] = count_iht / 100
        
        bp_time_mtx[cf_num] = bp_time
        omp_time_mtx[cf_num] = omp_time
        cosamp1_time_mtx[cf_num] = cosamp1_time
        cosamp2_time_mtx[cf_num] = cosamp2_time
        cosamp3_time_mtx[cf_num] = cosamp3_time
        iht_time_mtx[cf_num] = iht_time

#        
    plt.figure(figsize = (12, 6))
#    plt.plot(range(1, 71), accuracy_matrix_bp[1:71], '-ro', label = 'BP')
#    plt.plot(range(1, 71), accuracy_matrix_omp[1:71], '-bo', label = 'OMP')
#    
#    plt.plot(range(1, 71), accuracy_matrix_cosamp1[1:71], '-ro', label = 'CoSaMP(α = 1)')
#    plt.plot(range(1, 71), accuracy_matrix_cosamp2[1:71], '-go', label = 'CoSaMP(α = 2)')
#    plt.plot(range(1, 71), accuracy_matrix_cosamp3[1:71], '-bo', label = 'CoSaMP(α = 3)')
#
    plt.plot(range(1, 71), accuracy_matrix_iht[1:71], '-co', label = 'IHT')
#    plt.plot([bound, bound], [0, 1], 'k--', linewidth = 5 , label = 'm/(2logN)')
    plt.xlabel('# sparsity')
    plt.ylabel('Accuracy rate')
#    plt.title('CoSaMP(α = 1,2,3)')       
    plt.legend()
    plt.show()
#    
#    method_time = [bp_time, omp_time, cosamp1_time, cosamp2_time, cosamp3_time, iht_time]
#    method = ['BP', 'OMP', 'CoSaMP1', 'CoSaMP2', 'CoSaMP3', 'IHT']
#    for name, time_cost in zip(method, method_time):
#        print('{} Time Cost : {}'.format(name, round(time_cost, 4)))

    