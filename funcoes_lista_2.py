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

