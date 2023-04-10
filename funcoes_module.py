import numpy as np
import math

def funcao_1(x, y):
    return -1.2*y +7*np.exp(-0.3*x)

def analitica_1(x):
    return (70/9)*np.exp(-0.3*x)-(43/9)*np.exp(-1.2*x)

def funcao_1_lista(x,y):
    k = 0.5493
    return y*np.exp(k*x)

def funcao_2_lista(x,y):
    k = 0.15
    ym = 0.01
    return -k*(y-ym)

def funcao_2_analitica(x):
    k = 0.15
    y0 = 50
    ym = 0.01
    C = y0 - ym
    return C*np.exp(-k*x) + ym

def funcao_3_lista(x, y):
    T0=50
    Tf=0.01
    return T0 - Tf*y

def funcao_3_analitica(x):
    T0=50
    Tf=0.01
    C = -4900
    return T0/Tf + C*np.exp(-Tf*x)

def funcao_4_lista(x,y):
    ki = 0.01
    return ki*(1-y)*y

def funcao_4_analitica(x):
    b0 =10**-14
    ki = 0.01
    return 1/(1+((1/b0)-1)*np.exp(-ki*x))

def funcao_5_lista(x,y):
    return x*y**2+2*y

def funcao_5_analitica(x):
    return (-20*np.exp(2*x))/(9 - 5*np.exp(2*x)+10*np.exp(2*x)*x)

def funcao_6_lista(x,y):
    num = 100+50*math.cos((2*np.pi*x)/365)
    den = 10000
    return (num/den)*(5*np.exp((-2*x)/1000)-y)  
    