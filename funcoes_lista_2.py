import numpy as np
import math
import plotly.graph_objects as go

def funcao_1(t,y):
    dydt = (t-3.2)*y + 8*t*np.exp(((t-3.2)**2)/2)*math.cos(4*t**2)
    return dydt

def analitica_1(t, y):
    t0 = 0
    y0 = 0.75
    C = y0*np.exp(-((t0-3.2)**2)/2) - math.sin(4*t0**2)
    yt = np.exp(((t-3.2)**2)/2)*math.sin(4*t**2) + C

    return yt

def funcao_2_a(t,y):
    y1t = -(8/5)*t + (3/10)*y[0]
    y2t = (8/5)*t - (4/5)*y[1]

    return y1t, y2t

def funcao_2_b(t,y):
    xt = y[0]
    yt = y[1]
    zt = y[2]

    c1t = -0.5*xt
    c2t = 0.5*xt - (1/4)*yt
    c3t = (1/4)*yt - (1/6)*zt

    return c1t, c2t, c3t

def analitica_2_b(t,y):
    x0 = 20
    y0 = 40
    z0 = 60
    c1t = x0*np.exp(-t/2)
    c2t = -2*x0*np.exp(-t/2) + (y0+2*x0)*np.exp(-t/4)
    c3t = (3/2)*x0*np.exp(-t/2) - 3*(y0+2*x0)*np.exp(-t/4) + (z0-(3/2)*x0+3*(y0+2*x0))*np.exp(-t/6)

    return c1t, c2t, c3t

def funcao_2_c(t,y):
    xt = y[0]
    yt = y[1]
    zt = y[2]

    c1t = -(1/6)*xt + (1/6)*zt
    c2t = (1/6)*xt - (1/3)*yt
    c3t = (1/3)*yt - (1/6)*zt

    return c1t, c2t, c3t

def analitica_2_c(t,y):
    c1 = y[0]
    c2 = y[1]
    c3 = y[2]

    c1t = c1 + (c2+2*c3)*np.exp(-(t/3))*math.cos(t/6) + (2*c2+c3)*np.exp(-(t/3))*math.sin(t/6)
    c2t = 0.5*c1 + (-2*c2-c3)*np.exp(-(t/3))*math.cos(t/6) + (c2-2*c3)*np.exp(-(t/3))*math.sin(t/6)
    c3t = c1 + (c2+3*c3)*np.exp(-(t/3))*math.cos(t/6) + (-3*c2+c3)*np.exp(-(t/3))*math.sin(t/6)

    return c1t, c2t, c3t


def funcao_3_a_menos40(t,v):
    v1sobre2 = 1.5 
    k = 16
    gL = 19 
    eL = -67
    gNa = 74
    c = 10 
    eNa = 60 
    I = -40

    m_inf_v = 1/(1 + np.exp((v1sobre2-v)/k))
    vt = (I - gL*(v-eL) - gNa*m_inf_v*(v-eNa))/c

    return vt

    
def funcao_3_a_menos40(t,v):
    v1sobre2, k, gL, eL, gNa, c, eNa, I = 1.5, 16, 19, -67, 74, 10, 60, -40

    m_inf_v = 1/(1 + np.exp((v1sobre2-v)/k))
    vt = (I - gL*(v-eL) - gNa*m_inf_v*(v-eNa))/c

    return vt

def funcao_3_a_menos30(t,v):
    v1sobre2, k, gL, eL, gNa, c, eNa, I = 1.5, 16, 19, -67, 74, 10, 60, -30

    m_inf_v = 1/(1 + np.exp((v1sobre2-v)/k))
    vt = (I - gL*(v-eL) - gNa*m_inf_v*(v-eNa))/c

    return vt

def funcao_3_a_menos20(t,v):
    v1sobre2, k, gL, eL, gNa, c, eNa, I = 1.5, 16, 19, -67, 74, 10, 60, -20

    m_inf_v = 1/(1 + np.exp((v1sobre2-v)/k))
    vt = (I - gL*(v-eL) - gNa*m_inf_v*(v-eNa))/c

    return vt

def funcao_3_a_menos10(t,v):
    v1sobre2, k, gL, eL, gNa, c, eNa, I = 1.5, 16, 19, -67, 74, 10, 60, -10

    m_inf_v = 1/(1 + np.exp((v1sobre2-v)/k))
    vt = (I - gL*(v-eL) - gNa*m_inf_v*(v-eNa))/c

    return vt

def funcao_3_a_0(t,v):
    v1sobre2, k, gL, eL, gNa, c, eNa, I = 1.5, 16, 19, -67, 74, 10, 60, 0

    m_inf_v = 1/(1 + np.exp((v1sobre2-v)/k))
    vt = (I - gL*(v-eL) - gNa*m_inf_v*(v-eNa))/c

    return vt

def funcao_3_a_mais10(t,v):
    v1sobre2, k, gL, eL, gNa, c, eNa, I = 1.5, 16, 19, -67, 74, 10, 60, 10

    m_inf_v = 1/(1 + np.exp((v1sobre2-v)/k))
    vt = (I - gL*(v-eL) - gNa*m_inf_v*(v-eNa))/c

    return vt

def funcao_3_a_mais20(t,v):
    v1sobre2, k, gL, eL, gNa, c, eNa, I = 1.5, 16, 19, -67, 74, 10, 60, 20

    m_inf_v = 1/(1 + np.exp((v1sobre2-v)/k))
    vt = (I - gL*(v-eL) - gNa*m_inf_v*(v-eNa))/c

    return vt

def funcao_3_a_mais30(t,v):
    v1sobre2, k, gL, eL, gNa, c, eNa, I = 1.5, 16, 19, -67, 74, 10, 60, 30

    m_inf_v = 1/(1 + np.exp((v1sobre2-v)/k))
    vt = (I - gL*(v-eL) - gNa*m_inf_v*(v-eNa))/c

    return vt

def funcao_3_a_mais40(t,v):
    v1sobre2, k, gL, eL, gNa, c, eNa, I = 1.5, 16, 19, -67, 74, 10, 60, 40

    m_inf_v = 1/(1 + np.exp((v1sobre2-v)/k))
    vt = (I - gL*(v-eL) - gNa*m_inf_v*(v-eNa))/c

    return vt

def funcao_3_b(t,v):
    v1sobre2, k, gL, eL, gNa, c, eNa, I = 1.5, 16, 19, -67, 74, 10, 60, -1000

    m_inf_v = 1/(1 + np.exp((v1sobre2-v)/k))
    vt = (I - gL*(v-eL) - gNa*m_inf_v*(v-eNa))/c

    return vt

def funcao_4_a(t,y):
    xt = y[0]
    yt = y[1]
    zt = y[2]

    x_t = 10*yt - 10*xt
    y_t = 28*xt - yt - xt*zt
    z_t = xt*yt - (8/3)*zt

    return x_t, y_t, z_t

def funcao_4_b(t,y):
    xt = y[0]
    yt = y[1]
    zt = y[2]

    x_t = 77.27*(yt + xt*(1-8.375*(10**-6)*xt-yt))
    y_t = (1/77.27)*(zt-(1+xt)*yt)
    z_t = 0.161*(xt - zt)

    return x_t, y_t, z_t

def funcao_5_1(t,y):
    # Copy 1
    k0, k1, k2, k3, k4 =  0.1425, 0.4785, 0.2568, 0.3691, 0.7655
    wt = y[0]
    xt = y[1]
    pt = y[2]

    w_t = k0*wt - k1*wt + k2*xt
    x_t = k1*wt - k2*xt - k3*xt
    p_t = k3*xt - k4*pt

    return w_t, x_t, p_t

def analitica_5_a_1(t,y):
    # Copy 1
    k0, k1, k2, k3, k4 =  0.1425, 0.4785, 0.2568, 0.3691, 0.7655
    A = (-k0 + k1 + k2 + k3)
    delta = (A**2 - 4*(-k0*k2 - k0*k3 + k1*k3))
    lambda_1 = -k4
    lambda_2 = (-A + np.sqrt(delta))/2
    lambda_3 = (-A - np.sqrt(delta))/2

    w0 = 25

    wt = (w0/(lambda_3 - lambda_2))*(np.exp(lambda_3*t)*(lambda_3 + k2 + k3) - np.exp(lambda_2*t)*(lambda_2 + k2 + k3))
    xt = ((w0*k1)/(lambda_3-lambda_2))*(np.exp(lambda_3*t) - np.exp(lambda_2*t))
    pt = ((w0*k1*k3)/(lambda_2*lambda_3*(lambda_3 - lambda_2)))*np.exp(lambda_1*t)*(((np.exp((lambda_3-lambda_1)*t)-1)/(lambda_3-lambda_1))-((np.exp((lambda_2-lambda_1)*t)-1)/(lambda_2-lambda_1)))
    
    return wt, xt, pt

def funcao_5_2(t,y):
    # Copy 2
    k0, k1, k2, k3, k4 =  0.3232, 0.7696, 0.2341, 0.7404, 0.7952
    wt = y[0]
    xt = y[1]
    pt = y[2]

    w_t = k0*wt - k1*wt + k2*xt
    x_t = k1*wt - k2*xt - k3*xt
    p_t = k3*xt - k4*pt

    return w_t, x_t, p_t

def analitica_5_a_2(t,y):
    # Copy 2
    k0, k1, k2, k3, k4 =  0.3232, 0.7696, 0.2341, 0.7404, 0.7952
    A = (-k0 + k1 + k2 + k3)
    delta = (A**2 - 4*(-k0*k2 - k0*k3 + k1*k3))
    lambda_1 = -k4
    lambda_2 = (-A + np.sqrt(delta))/2
    lambda_3 = (-A - np.sqrt(delta))/2

    w0 = 25

    wt = (w0/(lambda_3 - lambda_2))*(np.exp(lambda_3*t)*(lambda_3 + k2 + k3) - np.exp(lambda_2*t)*(lambda_2 + k2 + k3))
    xt = ((w0*k1)/(lambda_3-lambda_2))*(np.exp(lambda_3*t) - np.exp(lambda_2*t))
    pt = ((w0*k1*k3)/(lambda_2*lambda_3*(lambda_3 - lambda_2)))*np.exp(lambda_1*t)*(((np.exp((lambda_3-lambda_1)*t)-1)/(lambda_3-lambda_1))-((np.exp((lambda_2-lambda_1)*t)-1)/(lambda_2-lambda_1)))
    
    return wt, xt, pt

def funcao_5_3(t,y):
    # Copy 3
    k0, k1, k2, k3, k4 =  0.1084, 0.8301, 0.2142, 0.4756, 0.1869
    wt = y[0]
    xt = y[1]
    pt = y[2]

    w_t = k0*wt - k1*wt + k2*xt
    x_t = k1*wt - k2*xt - k3*xt
    p_t = k3*xt - k4*pt

    return w_t, x_t, p_t

def analitica_5_a_3(t,y):
    # Copy 3
    k0, k1, k2, k3, k4 =  0.1084, 0.8301, 0.2142, 0.4756, 0.1869
    A = (-k0 + k1 + k2 + k3)
    delta = (A**2 - 4*(-k0*k2 - k0*k3 + k1*k3))
    lambda_1 = -k4
    lambda_2 = (-A + np.sqrt(delta))/2
    lambda_3 = (-A - np.sqrt(delta))/2

    w0 = 25

    wt = (w0/(lambda_3 - lambda_2))*(np.exp(lambda_3*t)*(lambda_3 + k2 + k3) - np.exp(lambda_2*t)*(lambda_2 + k2 + k3))
    xt = ((w0*k1)/(lambda_3-lambda_2))*(np.exp(lambda_3*t) - np.exp(lambda_2*t))
    pt = ((w0*k1*k3)/(lambda_2*lambda_3*(lambda_3 - lambda_2)))*np.exp(lambda_1*t)*(((np.exp((lambda_3-lambda_1)*t)-1)/(lambda_3-lambda_1))-((np.exp((lambda_2-lambda_1)*t)-1)/(lambda_2-lambda_1)))
    
    return wt, xt, pt

def analitica_5_b_1(t,y):
    # Copy 1
    k0, k1, k2, k3, k4 =  0.0, 0.01, 0.0, 0.01, 0.1

    A = (-k0 + k1 + k2 + k3)
    delta = (A**2 - 4*(-k0*k2 - k0*k3 + k1*k3))
    lambda_1 = -k4
    lambda_2 = (-A + np.sqrt(delta))/2
    lambda_3 = (-A - np.sqrt(delta))/2

    w0 = 25

    w_t = w0*np.exp(lambda_2*t)*(1+(lambda_2 + k2 + k3)*t)
    x_t = w0*k1*t*np.exp(lambda_2*t)
    p_t = ((w0*k1*k3)/(lambda_2-lambda_1))*np.exp(lambda_1*t)*(np.exp((lambda_2-lambda_1)*t)*(t-(1/(lambda_2-lambda_1))+(1/(lambda_2-lambda_1))))

    return w_t, x_t, p_t

def analitica_5_b_2(t,y):
    # Copy 2
    k0, k1, k2, k3, k4 =  0.0, 0.02, 0.0, 0.02, 0.11

    A = (-k0 + k1 + k2 + k3)
    delta = (A**2 - 4*(-k0*k2 - k0*k3 + k1*k3))
    lambda_1 = -k4
    lambda_2 = (-A + np.sqrt(delta))/2
    lambda_3 = (-A - np.sqrt(delta))/2

    w0 = 25

    w_t = w0*np.exp(lambda_2*t)*(1+(lambda_2 + k2 + k3)*t)
    x_t = w0*k1*t*np.exp(lambda_2*t)
    p_t = ((w0*k1*k3)/(lambda_2-lambda_1))*np.exp(lambda_1*t)*(np.exp((lambda_2-lambda_1)*t)*(t-(1/(lambda_2-lambda_1))+(1/(lambda_2-lambda_1))))

    return w_t, x_t, p_t

def analitica_5_b_3(t,y):
    # Copy 1
    k0, k1, k2, k3, k4 =  0.0, 0.03, 0.0, 0.03, 0.12

    A = (-k0 + k1 + k2 + k3)
    delta = (A**2 - 4*(-k0*k2 - k0*k3 + k1*k3))
    lambda_1 = -k4
    lambda_2 = (-A + np.sqrt(delta))/2
    lambda_3 = (-A - np.sqrt(delta))/2

    w0 = 25

    w_t = w0*np.exp(lambda_2*t)*(1+(lambda_2 + k2 + k3)*t)
    x_t = w0*k1*t*np.exp(lambda_2*t)
    p_t = ((w0*k1*k3)/(lambda_2-lambda_1))*np.exp(lambda_1*t)*(np.exp((lambda_2-lambda_1)*t)*(t-(1/(lambda_2-lambda_1))+(1/(lambda_2-lambda_1))))

    return w_t, x_t, p_t

def funcao_6(t,y):
    k1, k2, k3, k4, k5, k6 = 1, 3, 2, 1, 50, 1

    if (t <= 50):
        input = 0.5
    elif (50< t <= 100):
        input = 1
    elif(100 < t <= 150):
        input = 1.5
    elif(150 < t <= 200):
        input = 1
    elif(200 < t <= 250):
        input = 0.5

    xt, yt, zt, wt = y[0], y[1], y[2], y[3]

    x_t = k1*wt - k2*xt
    y_t = k3*input*xt - k4*yt
    z_t = k4*yt - k5*zt*wt
    w_t = k6 - k5*zt*wt

    return x_t, y_t, z_t, w_t

def funcao_7(t,y):
    y0, y1, y2 = y[0], y[1], y[2]

    f, kd, ktc, ktd = 0.6, 10**-5, 10**8, 10**8
    kp, kfm = 165.9, 4.021
    I, M = 0.0083, 1.96

    dy0dt = 2*f*kd*I - (ktc + ktd)*y0**2
    dy1dt = kp*M*y0 + kfm*M*(y0 - y1) - (ktc + ktd)*y0*y1
    dy2dt = kp*M*(2*y1+y0) + kfm*M*(y0-y2) - (ktc+ktd)*y0*y2
    dq0dt = kfm*M*y0 + ((ktc/2)+ktd)*y0**2
    dq1dt = kfm*M*y1 + (ktc+ktd)*y0*y1
    dq2dt = kfm*M*y2 + (ktc+ktd)*y0*y2 + ktc*y1**2
    
    return dy0dt, dy1dt, dy2dt, dq0dt, dq1dt, dq2dt


def special_plot(variable, dataset, title):
    fig = go.Figure()
    for column in dataset.columns:
        if column[-1:] == "x":
            dataset = dataset.rename(columns={column:f"{column[:-2]} Copy 1"})
        elif column[-1:] == "y":
            dataset = dataset.rename(columns={column:f"{column[:-2]} Copy 2"})
        else:
            if column != "Passo":
                dataset = dataset.rename(columns={column:f"{column} Copy 3"})
            
    for column in dataset.columns:
        if column[:1] == variable:
            fig.add_trace(go.Scatter(x=dataset['Passo'], y=dataset[column],
                        mode='lines',
                        name=column))
          
    fig.update_layout(title_text='Resultados por Método', title_x=0.5,\
                xaxis_title='t', yaxis_title=f'{title}',\
                height = 400, width = 600, font={'size':10})

    fig.show()

def special_plot_q3(all_var, all_var_datasets, title, v0=None, I=None):
    title = "Variações de I"
    fig = go.Figure()
    if I==None:
        for i in range(len(all_var)):
            fig.add_trace(go.Scatter(x=all_var_datasets[i]['Passo'], y=all_var_datasets[i]['Y(t) Runge-Kutta 4ª Ordem'],
                                mode='lines',
                                name=f"V0 = {v0} | I = {all_var[i]}"))
    if v0 == None:
        for i in range(len(all_var)):
            fig.add_trace(go.Scatter(x=all_var_datasets[i]['Passo'], y=all_var_datasets[i]['Y(t) Runge-Kutta 4ª Ordem'],
                                mode='lines',
                                name=f"I = {I} | V0 = {all_var[i]}")) 
        
    fig.update_layout(title_text='Resultados por Método', title_x=0.5,\
                xaxis_title='t', yaxis_title=f'{title}',\
                height = 400, width = 600, font={'size':10})

    fig.show()  