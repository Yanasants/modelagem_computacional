import numpy as np

def equations_system_1(t,y):
  Xt = y[0]
  Yt = y[1] 

  dXdt = 998*Xt + 1998*Yt
  dYdt = -999*Xt - 1999*Yt

  return dXdt, dYdt

def analytic_functions_1(t,y):
  Xt = 2*np.exp(-t) - np.exp(-1000*t)
  Yt = np.exp(-1000*t) - np.exp(-t)

  return Xt, Yt

def equations_system_2(t,y):
  
  Ut = y[0]
  Vt = y[1]

  dUdt = 1195*Ut - 1995*Vt
  dVdt = 1197*Ut - 1997*Vt

  return dUdt, dVdt

def analytic_functions_2(t,y):

  Ut = 10*np.exp(-2*t) - 8*np.exp(-800*t)
  Vt = 6*np.exp(-2*t) - 8*np.exp(-800*t)

  return Ut, Vt

def equations_system_3(t,y):
  Ut = y[0]
  Vt = y[1]
  dUdt = -2000*Ut + 999.75*Vt + 1000.25
  dVdt = Ut - Vt

  return dUdt, dVdt

def analytic_functions_3(t, y):
  Ut = (3999/8000)*np.exp(-(4001/2)*t) - (11999/8000)*np.exp(-0.5*t) + 1
  Vt = -(1/4000)*np.exp(-(4001/2)*t)

  return Ut, Vt
  print("Não apresenta valores conforme a tabela")

def equations_system_4(t,y):
  Ut = y[0]
  Vt = y[1]
  mu = 1000

  dUdt = Vt
  dVdt = mu*(1-Ut**2)*(Vt) - Ut

  return dUdt, dVdt

