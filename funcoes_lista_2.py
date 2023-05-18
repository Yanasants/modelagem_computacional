import numpy as np
import math

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

def funcao_3_I40neg(t,v):
    v1sobre2 = 1.5 #*10**-3
    k = 16 #*10**-3
    gL = 19 #*10**-3
    eL = -67 #*10**-3
    gNa = 74 #*10**-3
    c = 10 #*10**-6
    eNa = 60 #*10**-3
    I = -40
    m_inf_v = 1/(1+ np.exp((v1sobre2)-v)/k)
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

    x_t = 77.27*(yt + xt*(1-8.375*(10)**-6*xt-yt))
    y_t = (1/77.27)*(zt-(1+xt)*yt)
    z_t = 0.161*(xt - zt)

    return x_t, y_t, z_t

def funcao_5_a_1(t,y):
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

def funcao_5_a_2(t,y):
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

def funcao_5_a_3(t,y):
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
    # Copy 2
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

