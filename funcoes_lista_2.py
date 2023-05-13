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

def analitica_3_c(t,y):
    c1 = y[0]
    c2 = y[1]
    c3 = y[2]

    c1t = c1 + (c2+2*c3)*np.exp(-(t/3))*math.cos(t/6) + (2*c2+c3)*np.exp(-(t/3))*math.sin(t/6)
    c2t = 0.5*c1 + (-2*c2-c3)*np.exp(-(t/3))*math.cos(t/6) + (c2-2*c3)*np.exp(-(t/3))*math.sin(t/6)
    c3t = c1 + (c2+3*c3)*np.exp(-(t/3))*math.cos(t/6) + (-3*c2+c3)*np.exp(-(t/3))*math.sin(t/6)

    return c1t, c2t, c3t