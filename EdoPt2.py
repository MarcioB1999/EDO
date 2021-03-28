import numpy as np
import matplotlib.pyplot as plt

def PVC1(N,yi,yf):
    f_left = lambda dx: 1/(dx**2)
    f_mid = lambda dx: -(2/(dx**2) + 1)
    f_right = lambda dx: 1/(dx**2)
    
    x = np.linspace(yi,yf,N+1)
    dx = x[1] - x[0]

    A = np.zeros((N-1,N-1))
    b = np.zeros(N-1)
    
    
    #Primeira Linha
    A[0,0:2] = np.array([f_mid(dx),f_right(dx)])
    b[0] = -1 * (yi * f_left(dx))
    b[N-2] = -1 * (yf * f_right(dx)) 

    #Linhas do meio
    for i in range(1,N-2):
        A[i,i-1:i+2] = np.array([f_left(dx),f_mid(dx),f_right(dx)])

    #Última Linha
    A[N-2,N-3:N-1] = np.array([f_left(dx),f_mid(dx)])
    
    print("Matriz do PVC1 = \n",A,'\n') 
    
    y = np.linalg.solve(A,b)
    y_ = np.hstack([0,y,1])
    
    return x,y_
       
            
def PVC2(N=8):
    def get_neighbors(i, j, N):
        neighbors = [(i-1, j), (i, j+1), (i+1, j), (i, j-1)]
        neighbors = [(i, j) for (i, j) in neighbors if 0 <= i < N and 0 <= j < N]
    
        return neighbors
    
    A = np.zeros(((N-1)**2, (N-1)**2))
    b = np.ones((N-1)**2) * 4
    dx = 1/N
    idx = lambda i, j: i*(N-1) + j
    
    for i in range(N-1):
        for j in range(N-1):
            neighbors = get_neighbors(i, j, N-1)
            for x, y in neighbors:
                A[idx(i, j)][idx(x, y)] = 1/(dx**2)
            A[idx(i, j)][idx(i, j)] = -2* (1/(dx**2) + 1/(dx**2))
            

    y = np.linalg.solve(A,b)
    
    print("Matriz do PVC2 = \n",A,'\n') 
    return y

def relative_error(y_approx,y):
    return abs(y_approx - y) / y

if __name__ == "__main__":  
    N,yi,yf = 8,0,1
    x,y = PVC1(N,yi,yf)
    
    F = lambda x: (1/(np.exp(-1)-np.exp(1))) * (np.exp(-x) - np.exp(x))        #Sol. exata
    y_exato = [F(i) for i in x]
    errors = [relative_error(y[i],y_exato[i]) for i in range(1,len(x)-1)]
    print('x = ',x,'\n')
    print(f"Meus y = {y} \n\n Y exatos = {y_exato}\n\n")
    print(f"Erros relativos sem o 0 e o 1 = {errors}\n\n")
    
    u = PVC2(N)
    y_pdf = [-0.171875,-0.21875,-0.171875,-0.21875,-0.28125,-0.21875,-0.171875,-0.21875,-0.171875]
    #errors2 = [relative_error(y2,y_pdf[i]) for y in y2]
    print(f"Meus u = {u} \n\n")
    #print(f"Erros relativos = {errors2}\n\n")
        
        
        
        
        
       
    
    
