import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from funcoes_trabalho_final import *

class Edo:

  def __init__(self, y_inicial, dom, a, b, I_app, I_function, methods=[], output=pd.DataFrame([])):
    self.y = y_inicial
    self.dom = dom
    self.a = a
    self.b = b
    self.I_app = I_app
    self.I_function = I_function 
    self.current_density = None
    self.methods = methods

    self.step = (self.b - self.a)/self.dom
    self.size = self.dom-1
    self.all_x = np.arange(self.a, self.b, self.step)
    self.V, self.n, self.m, self.h = np.zeros(len(self.all_x)), np.zeros(len(self.all_x)), np.zeros(len(self.all_x)), np.zeros(len(self.all_x))
    self.V[0], self.n[0], self.m[0], self.h[0] = 0, 0, 0, 0

    self.output = pd.DataFrame([])
    self.df_output()

  def adbash2(self):
    self.all_y_adbash2 = np.zeros(len(self.all_x))
    self.all_y_adbash2[0] = self.y

    k1 = self.f(self.all_x[0],self.all_y_adbash2[0])
    k2 = self.f(self.all_x[1], self.all_y_adbash2[0]+self.step*k1)

    self.all_y_adbash2[1] = self.all_y_adbash2[0]  + 0.5*self.step*(k1+k2)
    for i in range(len(self.all_x)-1):
      k2 = k1
      k1 = self.f(self.all_x[i], self.all_y_adbash2[i])
      yp = self.all_y_adbash2[i] + self.h*(3*k1 - k2)/2
      self.all_y_adbash2[i+1] = self.all_y_adbash2[i] + 0.5*self.h*(k1 + self.f(self.all_x[i+1], yp))
  

  def rk4_step(self, V, n, m, h, dt, I_app):

    k1_V = dt * dV_dt(V, m, h, n, I_app)
    k1_n = dt * dn_dt(V, n)
    k1_m = dt * dm_dt(V, m)
    k1_h = dt * dh_dt(V, h)
    
    k2_V = dt * dV_dt(V + 0.5 * k1_V, m, h, n, I_app)
    k2_n = dt * dn_dt(V + 0.5 * dt, n + 0.5 * k1_n)
    k2_m = dt * dm_dt(V + 0.5 * dt, m + 0.5 * k1_m)
    k2_h = dt * dh_dt(V + 0.5 * dt, h + 0.5 * k1_h)

    k3_V = dt * dV_dt(V + 0.5 * k2_V, m, h, n, I_app)
    k3_n = dt * dn_dt(V + 0.5 * dt, n + 0.5 * k2_n)
    k3_m = dt * dm_dt(V + 0.5 * dt, m + 0.5 * k2_m)
    k3_h = dt * dh_dt(V + 0.5 * dt, h + 0.5 * k2_h)

    k4_V = dt * dV_dt(V + k3_V, m, h, n, I_app)
    k4_n = dt * dn_dt(V + dt, n + k3_n)
    k4_m = dt * dm_dt(V + dt, m + k3_m)
    k4_h = dt * dh_dt(V + dt, h + k3_h)

    V += (k1_V + 2.0 * k2_V + 2.0 * k3_V + k4_V) / 6.0
    n += (k1_n + 2.0 * k2_n + 2.0 * k3_n + k4_n) / 6.0
    m += (k1_m + 2.0 * k2_m + 2.0 * k3_m + k4_m) / 6.0
    h += (k1_h + 2.0 * k2_h + 2.0 * k3_h + k4_h) / 6.0

    return V, n, m, h
  
  def I_function_to_output(self):
    self.current_density = self.I_function(self.V, self.n, self.m, self.h)
  
  def df_output(self):
    self.output['Passo'] = self.all_x
    if 'rk4_final' in self.methods:
      for i in range(1, len(self.all_x)):
        self.V[i], self.n[i], self.m[i], self.h[i] = self.rk4_step(self.V[i-1], self.n[i-1], self.m[i-1], self.h[i-1], self.step, self.I_app)
      self.output['V(t)'] = self.V
      self.output['n(t)'] = self.n
      self.output['m(t)'] = self.m
      self.output['h(t)'] = self.h
      self.I_function_to_output()
      self.output['I(t)'] = self.current_density

    return self.output



 