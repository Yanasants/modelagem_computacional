import numpy as np

delta = 0.05
p = 0.85
v = 0.05
gama = 0.5
u = 0.02
I = 6*10**6

theta =  gama + u + v + delta

sigma = p + theta + u
b = p*v + p*delta + p*u + u*theta
a = 2*p + u

n1 = (1/2)*(sigma - np.sqrt(sigma**2 - 4*b))
n2 = (1/2)*(sigma + np.sqrt(sigma**2 - 4*b))

C0 = (((2*p + u)*I)/(v*p + u*theta + p*delta + p*u)) + 500
N0 = ((((2*(p+theta) - u+ delta))*I)/(v*p + u*theta + p*delta + p*u)) + 500


K1 = (b*(p + theta - n2)*C0 + I*(a*n2 - b) - p*b*N0)/(b*(n1-n2))
K2 = (-b*(p + theta - n1)*C0 + I*(b - a*n1) + p*b*N0)/(b*(n1-n2))

def equations_system(x,y):

  t = x
  Ct = y[0]
  Nt = y[1]

  dCdt = I - (p + theta)*Ct + p*Nt
  dNdt = 2*I - (v + delta)*Ct - u*Nt

  return dCdt, dNdt

def analitica(t,y):

  Ct = (K1*np.exp(-n1*t)) + (K2*np.exp(-n2*t)) + ((a/b)*I)
  Nt = Ct + ((theta/p)*K1*np.exp(-n1*t)) + ((theta/p)*K2*np.exp(-n2*t)) + \
      (((theta*a)/(p*b))*I) - (I/p) \
      - ((1/p)*(n1*K1*np.exp(-n1*t)+n2*K2*np.exp(-n2*t)))

  return Ct, Nt