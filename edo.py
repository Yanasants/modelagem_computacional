import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from funcoes_trabalho_final import *

class Edo:

  def __init__(self, y_inicial, dom, a, b, I_app, I_function, v0, n0, m0, h0, fixo=False, graph=1, methods=[], output=pd.DataFrame([])):
    self.y = y_inicial
    self.dom = dom
    self.a = a
    self.b = b
    self.I_app = I_app
    self.I_function = I_function 
    self.current_density = None
    self.methods = methods
    self.fixo = fixo
    self.graph = graph

    self.step = (self.b - self.a)/self.dom
    self.size = self.dom-1
    self.all_x = np.arange(self.a, self.b, self.step)
    self.all_V, self.all_n, self.all_m, self.all_h = None, None, None, None
    self.all_V, self.all_n, self.all_m, self.all_h = np.zeros(len(self.all_x)), np.zeros(len(self.all_x)), np.zeros(len(self.all_x)), np.zeros(len(self.all_x))
    self.V, self.n, self.m, self.h = np.zeros(len(self.all_x)), np.zeros(len(self.all_x)), np.zeros(len(self.all_x)), np.zeros(len(self.all_x))
    self.v0, self.n0, self.m0, self.h0 = v0, n0, m0, h0
    self.all_V[0], self.all_n[0], self.all_m[0], self.all_h[0] = self.v0, self.n0, self.m0, self.h0

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
  
    
  def pred_corretor(self, V, n, m, h, dt, I_app):
    
    k1_V = dt * dV_dt(V, n, m, h, I_app, self.graph)
    k1_n = dt * dn_dt(V, n)
    k1_m = dt * dm_dt(V, m)
    k1_h = dt * dh_dt(V, h)
    
    k2_V = dt * dV_dt(V + 0.5 * k1_V, n, m, h, I_app, self.graph)
    k2_n = dt * dn_dt(V + 0.5 * dt, n + 0.5 * k1_n)
    k2_m = dt * dm_dt(V + 0.5 * dt, m + 0.5 * k1_m)
    k2_h = dt * dh_dt(V + 0.5 * dt, h + 0.5 * k1_h)

    self.all_V[1] = self.all_V[0]  + 0.5*dt*(k1_V+k2_V)
    self.all_n[1] = self.all_n[0]  + 0.5*dt*(k1_n+k2_n)
    self.all_m[1] = self.all_m[0]  + 0.5*dt*(k1_m+k2_m)
    self.all_h[1] = self.all_h[0]  + 0.5*dt*(k1_h+k2_h)
    

    for i in range(1, len(self.all_x)-1):
      k2_V = k1_V
      k1_V = dt * dV_dt(self.all_V[i], self.all_n[i], self.all_m[i], self.all_h[i], I_app, self.graph)
      yp_V = self.all_V[i] + dt*(3*k1_V - k2_V)/2
      self.all_V[i+1] = self.all_V[i] + 0.5*dt*(k1_V + dV_dt(yp_V, self.all_n[i], self.all_m[i], self.all_h[i], I_app, self.graph))

      k2_n = k1_n
      k1_n = dt * dn_dt(self.all_V[i], self.all_n[i])
      yp_n = self.all_n[i] + dt*(3*k1_n - k2_n)/2
      self.all_n[i+1] = self.all_n[i] + 0.5*dt*(k1_n + dn_dt(self.all_V[i+1], yp_n))

      k2_m = k1_m
      k1_m = dt * dm_dt(self.all_V[i], self.all_m[i])
      yp_m = self.all_m[i] + dt*(3*k1_m - k2_m)/2
      self.all_m[i+1] = self.all_m[i] + 0.5*dt*(k1_m + dm_dt(self.all_V[i+1], yp_m))
    
      k2_h = k1_h
      k1_h = dt * dh_dt(self.all_V[i], self.all_h[i])
      yp_h = self.all_h[i] + dt*(3*k1_h - k2_h)/2
      self.all_h[i+1] = self.all_h[i] + 0.5*dt*(k1_h + dh_dt(self.all_V[i+1], yp_h)) 
     



  # Alpha e Beta variando
  def rk4_step(self, V, n, m, h, dt, I_app):

    k1_V = dt * dV_dt(V, n, m, h, I_app, self.graph)
    k1_n = dt * dn_dt(V, n)
    k1_m = dt * dm_dt(V, m)
    k1_h = dt * dh_dt(V, h)
    
    k2_V = dt * dV_dt(V + 0.5 * k1_V, n, m, h, I_app, self.graph)
    k2_n = dt * dn_dt(V + 0.5 * dt, n + 0.5 * k1_n)
    k2_m = dt * dm_dt(V + 0.5 * dt, m + 0.5 * k1_m)
    k2_h = dt * dh_dt(V + 0.5 * dt, h + 0.5 * k1_h)

    k3_V = dt * dV_dt(V + 0.5 * k2_V, n, m, h, I_app, self.graph)
    k3_n = dt * dn_dt(V + 0.5 * dt, n + 0.5 * k2_n)
    k3_m = dt * dm_dt(V + 0.5 * dt, m + 0.5 * k2_m)
    k3_h = dt * dh_dt(V + 0.5 * dt, h + 0.5 * k2_h)

    k4_V = dt * dV_dt(V + k3_V, n, m, h, I_app, self.graph)
    k4_n = dt * dn_dt(V + dt, n + k3_n)
    k4_m = dt * dm_dt(V + dt, m + k3_m)
    k4_h = dt * dh_dt(V + dt, h + k3_h)

    V += (k1_V + 2.0 * k2_V + 2.0 * k3_V + k4_V) / 6.0
    n += (k1_n + 2.0 * k2_n + 2.0 * k3_n + k4_n) / 6.0
    m += (k1_m + 2.0 * k2_m + 2.0 * k3_m + k4_m) / 6.0
    h += (k1_h + 2.0 * k2_h + 2.0 * k3_h + k4_h) / 6.0

    return V, n, m, h
  
  # Alpha e Beta fixos 
  def rk4_step_fixo(self, V, n, m, h, dt, I_app):

    k1_V = dt * dV_dt(V, n, m, h,I_app, self.graph)
    k1_n = dt * dn_dt(V, n)
    k1_m = dt * dm_params_fixo(V, m)
    k1_h = dt * dh_params_fixo(V, h)
    
    k2_V = dt * dV_dt(V + 0.5 * k1_V,n, m, h, I_app, self.graph)
    k2_n = dt * dn_dt(V + 0.5 * dt, n + 0.5 * k1_n)
    k2_m = dt * dm_params_fixo(V + 0.5 * dt, m + 0.5 * k1_m)
    k2_h = dt * dh_params_fixo(V + 0.5 * dt, h + 0.5 * k1_h)

    k3_V = dt * dV_dt(V + 0.5 * k2_V, n, m, h, I_app, self.graph)
    k3_n = dt * dn_dt(V + 0.5 * dt, n + 0.5 * k2_n)
    k3_m = dt * dm_params_fixo(V + 0.5 * dt, m + 0.5 * k2_m)
    k3_h = dt * dh_params_fixo(V + 0.5 * dt, h + 0.5 * k2_h)

    k4_V = dt * dV_dt(V + k3_V, n, m, h, I_app, self.graph)
    k4_n = dt * dn_dt(V + dt, n + k3_n)
    k4_m = dt * dm_params_fixo(V + dt, m + k3_m)
    k4_h = dt * dh_params_fixo(V + dt, h + k3_h)

    V += (k1_V + 2.0 * k2_V + 2.0 * k3_V + k4_V) / 6.0
    n += (k1_n + 2.0 * k2_n + 2.0 * k3_n + k4_n) / 6.0
    m += (k1_m + 2.0 * k2_m + 2.0 * k3_m + k4_m) / 6.0
    h += (k1_h + 2.0 * k2_h + 2.0 * k3_h + k4_h) / 6.0

    return V, n, m, h
  
  def I_function_to_output(self):
    self.current_density = self.I_function(self.V, self.n, self.m, self.h)

  def I_function_to_output_2(self):
    self.current_density = self.I_function(self.all_V, self.all_n, self.all_m, self.all_h)

  def df_output(self):
    self.output['Passo'] = self.all_x
    if 'rk4_final' in self.methods:
      if self.fixo:
        for i in range(1, len(self.all_x)):
          self.V[i], self.n[i], self.m[i], self.h[i] = self.rk4_step_fixo(self.V[i-1], self.n[i-1], self.m[i-1], self.h[i-1], self.step, self.I_app)
      else:
        for i in range(1, len(self.all_x)):
          self.V[i], self.n[i], self.m[i], self.h[i] = self.rk4_step(self.V[i-1], self.n[i-1], self.m[i-1], self.h[i-1], self.step, self.I_app)    
      self.I_function_to_output()
      self.output['V(t)'] = self.V
      self.output['n(t)'] = self.n
      self.output['m(t)'] = self.m
      self.output['h(t)'] = self.h
      self.output['I(t)'] = self.current_density

    if 'pred_corretor' in self.methods:
      if self.fixo:
          self.pred_corretor(V=self.v0, n=self.n0, m=self.m0, h=self.h0, dt=self.step, I_app=self.I_app)  
      else:
          self.pred_corretor(V=self.v0, n=self.n0, m=self.m0, h=self.h0, dt=self.step, I_app=self.I_app)   
      self.I_function_to_output_2()
      
      self.output['V(t)'] = self.all_V
      self.output['n(t)'] = self.all_n
      self.output['m(t)'] = self.all_m
      self.output['h(t)'] = self.all_h
      self.output['I(t)'] = self.current_density

    return self.output



 