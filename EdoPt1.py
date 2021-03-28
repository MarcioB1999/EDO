import numpy as np

def runge_kutta4(F,S_0,dt,t = 0):
    curr_S = []
    curr_S.append(S_0)
    prev_S = S_0 - np.ones(S_0.shape[0])

    #Só 3 estados
    for i in range(3):
        prev_S = curr_S[i]
        S_2 = curr_S[i] + dt/2 * F(curr_S[i])
        S_3 = curr_S[i] + dt/2 * F(S_2)
        S_4 = curr_S[i] + dt * F(S_3)        
        curr_S.append(curr_S[i] + dt/6 * (F(curr_S[i]) + 2 * F(S_2) + 2*F(S_3) +F(S_4)))

    
    return curr_S                           #Retorna o S_0 inicial e os 3 gerados
    

def pred_cor_4(F,S_0,dt,t = 0,stop_cond=0):
    s = runge_kutta4(F,S_0,dt)                             # 4 estados iniciais. i-3,i-2,i-1,i

    F_3,F_2,F_1,F_0 = (F(s[0]),F(s[1]),F(s[2]),F(s[3]))    # 4 estados inicias na F
    
    curr_S = S_0
    prev_S = S_0 - np.ones(S_0.shape)
   
    i = 0
    while stop_cond(prev_S,curr_S):
        if (i+1) % 100000 == 0: print("curr_s = ",curr_S) 
        prev_S = curr_S
        
        F_4,F_3,F_2,F_1 = F_3,F_2,F_1,F_0
        curr_S_pred = prev_S + (dt/24) * (55*F_1 -59*F_2 +37*F_3 -9*F_4)                #Predicao
        F_0 = F(curr_S_pred)

        curr_S = prev_S + (dt/24) * (9*F_0 +19*F_1 -5*F_2 + F_3)                        #Correcao
        F_0 = F(curr_S)  
        i+=1

    return i,prev_S,curr_S

if __name__ == "__main__":
    #Dados do exemplo
    k, m, g = (.25, 2,10)
    S_0 = np.array([5, 200])            # S_0 = [v_0, y_0]

    F = lambda S_t: np.array([          # F = F(S(t),t)
        -g - k/m * S_t[0],          
        S_t[0]
    ])


    # Altura máxima da trajetória e tempo decorrido até a altura máxima
    stop_cond = lambda prev_S, curr_S: curr_S[1] > prev_S[1]

    for i in range(2,6):
        dt = 10 ** (-i)
        print('--'*8+f' delta t = {dt} '+'--'*8+'\n')
        iter, prev_S, curr_S = pred_cor_4(F,S_0,dt,stop_cond=stop_cond)
        print(f"{iter} estados")
        print(f"Penúltimo S = {prev_S} \nÚltimo S = {curr_S}")
        print(f'A altura maxima ocorre em {iter * dt:.03f} segundos com altura de {prev_S[1]:.03f} metros \n')

    # # Tempo total até a queda no mar e a velocidade no momento do impacto com o mar
    stop_cond = lambda prev_S, curr_S: curr_S[1] > 0

    for i in range(2,6):
        dt = 10 ** (-i)
        print('--'*8+f' delta t = {dt} '+'--'*8+'\n')
        iter, prev_S, curr_S = pred_cor_4(F,S_0,dt,stop_cond=stop_cond)
        print(f"{iter} estados")
        print(f"Penúltimo S = {prev_S} \nÚltimo S = {curr_S}")
        print(f'Caiu no mar apos {iter * dt: .04f} segundos a uma velocidade de {curr_S[0]: .04f} m/s')
        
