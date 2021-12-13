import numpy as np
import matplotlib.pyplot as plt

def seir_model(current_values):
    alpha = 0.25
    beta = 0.2 #?
    beta_H = 0.2 #?
    gamma = 0.05
    epsilon = 0.00001
    f = 0.3
    p = 0.4
    mu_n = 0
    mu_i = 0.001
    mu_hn = 0.025
    mu_hi = 0.025 
    h_n = 0.0001
    h_i = 0.01
    d_n = 0.2
    d_i = 0.05

    S_prev, I0_prev, Ip_prev, Im_prev, H_s, R_prev, H_im, H_ip, H_R = current_values

    S_next = -(mu_n*S_prev)-beta*(I0_prev + Im_prev) - (S_prev*h_n) + (H_s*d_n)
    I0_next = (S_prev*beta*Im_prev) + I0_prev*(-mu_n + S_prev*beta -(1-f)*(alpha - gamma)*gamma*((alpha*(1-f)+gamma*f)/((gamma**2)*(1-f)+(alpha**2)*f)) + f*gamma*(2*p-1)-h_n)
    Ip_next = - (Ip_prev*mu_i) + (I0_prev*f*gamma*p) - Ip_prev*(h_i + (alpha*gamma)/(alpha-gamma))
    Im_next = - Im_prev*(mu_n + h_n + (alpha*gamma)/(alpha-gamma)) + I0_prev*f*gamma*(1-p)
    R_next = - R_prev*(mu_n+h_n)+Im_prev*((alpha*gamma)/(alpha-gamma))+I0_prev*((1-f)*(alpha-gamma)*gamma*((alpha*(1-f)+gamma*f)/(alpha**2*(1-f)+gamma**2*f))) + (H_im * d_n) + (H_ip*d_i) + (H_R*d_n)
    Hs_next = - H_s*d_n + S_prev*h_n - H_s* mu_hn - H_s * beta_H * H_im
    Him_next = I0_prev*(1-p)*h_n - H_im*h_n*(1-p) - H_im*d_n + H_s*beta_H*H_im - H_im*mu_hn
    Hip_next = - H_ip*d_i + Im_prev*h_n*p + Ip_prev*h_i + I0_prev*p*h_n - H_im*(mu_hi + epsilon)
    HR_next = - (H_R*d_n) + (R_prev*h_n) - (H_R* mu_hn)
    return [S_next, I0_next, Ip_next, Im_next, Hs_next, R_next, Him_next, Hip_next, HR_next]

def euler(current_values, h):
    numerical_values = []
    estimated_values = seir_model(current_values)
    for i in range(0, len(current_values)):
        numerical_values.append(current_values[i] + h*estimated_values[i])
    return numerical_values

# set as initial values first
numerical_solutions = []
approx_values = [[1-0.0010001, 0.0000001, 0, 0, 0.0001, 0, 0.1, 0, 0], 
        [1-0.0010001, 0.0000001, 0, 0.0000001, 0.0001, 0, 0.1, 0.001, 0.00001],
        [1-0.0010001, 0.0001, 0, 0.00001, 0.0001, 0, 0.1, 0.003, 0.0001],
        [1-0.0010001, 0.01, 0, 0.001, 0.0001, 0, 0.1, 0.002, 0.001]]
h = 0.1
for initial_conds in approx_values:
    current_values = initial_conds
    current_numerical = []
    for i in range(0,500):
        updated_values = euler(current_values,h)
        current_numerical.append(updated_values)
        current_values = updated_values
    numerical_solutions.append(current_numerical)


for sol_list in numerical_solutions:
    s_list = []
    i0_list = []
    ip_list = []
    im_list = []
    h_list = []
    r_list = []
    him_list = []
    hip_list = []
    hr_list = []
    for sol in sol_list:
        for i in range(0, len(sol)):
            if i == 0:
                s_list.append(sol[i])
            elif i == 1:
                i0_list.append(sol[i])
            elif i ==2:
                ip_list.append(sol[i])
            elif i == 3:
                im_list.append(sol[i])
            elif i == 4:
                h_list.append(sol[i])
            elif i == 5:
                r_list.append(sol[i])
            elif i == 6:
                him_list.append(sol[i])
            elif i == 7:
                hip_list.append(sol[i])
            else:
                hr_list.append(sol[i])

    #plt.plot(s_list, c = 'b', label = 'Susceptible Population')
    plt.plot(i0_list, c ='c', label = 'Infected Not Tested')
    plt.plot(ip_list, c = 'y', label = 'Tested Positive')
    plt.plot(im_list, c = 'purple', label = 'Tested Negative')
    #plt.plot(h_list, c = 'c', label = 'In Hospital (not related)')
    plt.plot(r_list, c = 'orange', label = 'Recovered')
    #plt.plot(him_list, c = 'y', label = 'In Hospital (not related)')
    plt.plot(hip_list, c = 'r', label = 'In Hospital (COVID-19 Related)')
    #plt.plot(hr_list, c = 'orange', label = 'Recovered (not related)')
    plt.title('Model Results with B=0.2')
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.tight_layout()
    plt.show()
    plt.clf()







