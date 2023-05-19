import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

class Edo:

  def __init__(self, y_inicial, dom, a,b, f_function, methods=[], \
               analytic_function=None, all_error=[], output=pd.DataFrame([]), \
              analytic_solution=False, bidimensional=False, tridimensional=False, \
              quadridimensional=False):
    self.y = y_inicial
    self.dom = dom
    self.a = a
    self.b = b
    self.methods = methods
    self.f_function = f_function
    self.analytic_solution = analytic_solution
    self.analytic_function = analytic_function
    self.bidimensional = bidimensional
    self.tridimensional = tridimensional
    self.quadridimensional = quadridimensional

    self.h = (self.b - self.a)/self.dom
    self.size = self.dom-1
    self.all_x = np.arange(self.a, self.b, self.h)
    self.all_y = []
    self.all_y_euller =[self.y]
    self.all_y_analitica = []
    self.all_y_euller_implicito = [self.y]
    self.all_y_ponto_medio = [self.y]
    self.y_newton = [self.y]
    self.all_y_rk4 = [self.y]
    self.all_y_heuen = [self.y]

    self.all_y_adbash2 = None
    self.all_error = all_error
    self.output = pd.DataFrame([])

  def f(self, x, y):
    return self.f_function(x, y)

  def report_euller_1d(self):
    for i in range(len(self.all_x)-1):
        self.all_y_euller.append(self.all_y_euller[i]+self.h*self.f(self.all_x[i],self.all_y_euller[i]))

  def report_euller_implicito_1d(self):
    for i in np.arange(0, len(self.all_x)):
      if len(self.all_y_euller_implicito)<=self.size:
          self.k1 = self.h*self.f(self.all_x[i],self.all_y_euller[i])
          self.k2 = self.h*self.f(self.all_x[i+1],self.all_y_euller[i]+self.k1)
          self.all_y_euller_implicito.append(self.all_y_euller_implicito[i] + ((self.k1+self.k2)/2))

  def report_ponto_medio_1d(self):
    for i in np.arange(0, len(self.all_x)):
      if len(self.all_y_ponto_medio)<=self.size:
          self.k1_pm = self.f(self.all_x[i],self.all_y_euller[i])
          self.k2_pm = self.f(self.all_x[i]+(self.h)/2,self.all_y_euller[i]+self.k1_pm*(self.h/2))
          self.all_y_ponto_medio.append(self.all_y_ponto_medio[i] + self.k2_pm*self.h)
          
  def report_newton_1d(self):
    for i in np.arange(0, len(self.all_x)):
      if len(self.y_newton)<=self.size:
          z = self.y_newton[i]+self.f(self.all_x[i], self.y_newton[i])*self.h
          count = 0
          while count <= 5:
            F = self.y_newton[i]+self.h*self.f(self.all_x[i+1], z) - z
            dF = self.h*self.df(self.all_x[i+1], z) - 1
            z = z - (F/dF)
            count+=1
          self.y_newton.append(z)

  def report_heuen_1d(self):
    self.all_y_heuen = np.zeros(len(self.all_x))
    self.all_y_heuen[0] = self.y
    for i in range(1, len(self.all_x)):
      k1 = self.h*self.f(self.all_x[i-1], self.all_y_heuen[i-1])
      k2 = self.h*self.f(self.all_x[i], self.all_y_heuen[i-1]+k1)
      self.all_y_heuen[i] = self.all_y_heuen[i-1] + (k1+k2)/2

  def report_rk4_1d(self):
    for i in range(len(self.all_x)-1):
      k1 = self.h*np.array((self.f(self.all_x[i], self.all_y_rk4[i])))
      k2 = self.h*np.array(self.f(self.all_x[i]+self.h/2, self.all_y_rk4[i]+k1/2))
      k3 = self.h*np.array(self.f(self.all_x[i]+self.h/2, self.all_y_rk4[i]+k2/2))
      k4 = self.h*np.array(self.f(self.all_x[i]+self.h, self.all_y_rk4[i]+k3))
      self.all_y_rk4.append(self.all_y_rk4[i] + (k1 + 2*k2 +2*k3 + k4)/6)

  def analitica_1d(self):
    self.all_y_analitica.append(self.y)
    for i in range(len(self.all_x)-1):
        self.all_y_analitica.append(self.f_analitica(self.all_x[i], self.all_y_analitica[i]))

  def report_euller(self):
    self.all_y_euller = np.zeros((len(self.all_x), len(self.y)))
    k = np.zeros((len(self.all_x), len(self.y)))
    self.all_y_euller[0] = self.y

    for i in range(len(self.all_x)-1):
      k[i] = self.f(self.all_x[i], self.all_y_euller[i])
      self.all_y_euller[i+1] = self.h*k[i]+self.all_y_euller[i]

  def report_heuen(self):
    self.all_y_heuen = np.zeros((len(self.all_x), len(self.y)))
    self.all_y_heuen[0] = self.y
    for i in range(1, len(self.all_x)):
      k1 = self.h*np.array(self.f(self.all_x[i-1], self.all_y_heuen[i-1]))
      k2 = self.h*np.array(self.f(self.all_x[i], self.all_y_heuen[i-1]+k1))
      self.all_y_heuen[i] = self.all_y_heuen[i-1] + (k1+k2)/2

  def report_rk4(self):
    self.all_y_rk4 = np.zeros((len(self.all_x), len(self.y)))
    self.all_y_rk4[0] = self.y
    for i in range(len(self.all_x)-1):
      k1 = self.h*np.array((self.f(self.all_x[i], self.all_y_rk4[i])))
      k2 = self.h*np.array(self.f(self.all_x[i]+self.h/2, self.all_y_rk4[i]+k1/2))
      k3 = self.h*np.array(self.f(self.all_x[i]+self.h/2, self.all_y_rk4[i]+k2/2))
      k4 = self.h*np.array(self.f(self.all_x[i]+self.h, self.all_y_rk4[i]+k3))
      self.all_y_rk4[i+1] = self.all_y_rk4[i] + (k1 + 2*k2 +2*k3 + k4)/6

  def adbash2(self):
    self.all_y_adbash2 = np.zeros(len(self.all_x))
    self.all_y_adbash2[0] = self.y

    k1 = self.f(self.all_x[0],self.all_y_adbash2[0])
    k2 = self.f(self.all_x[1], self.all_y_adbash2[0]+self.h*k1)

    self.all_y_adbash2[1] = self.all_y_adbash2[0]  + 0.5*self.h*(k1+k2)
    for i in range(len(self.all_x)-1):
      k2 = k1
      k1 = self.f(self.all_x[i], self.all_y_adbash2[i])
      yp = self.all_y_adbash2[i] + self.h*(3*k1 - k2)/2
      self.all_y_adbash2[i+1] = self.all_y_adbash2[i] + 0.5*self.h*(k1 + self.f(self.all_x[i+1], yp))

  def f_analitica(self, x, y):
    return self.analytic_function(x, y)

  def analitica(self):
    self.all_y_analitica = np.zeros((len(self.all_x), len(self.y)))
    self.all_y_analitica[0] = self.y
    for i in range(len(self.all_x)-1):
      self.all_y_analitica[i+1] = self.f_analitica(self.all_x[i], self.all_y_analitica[i])
      

  @property
  def df_output(self):
    self.output['Passo'] = self.all_x
    if self.bidimensional == False and self.tridimensional == False and self.quadridimensional==False:
      if 'euller' in self.methods:
        self.report_euller_1d()
        self.output['Y(t) Euller'] = self.all_y_euller
      if 'euller_implicito' in self.methods:
        self.report_euller_implicito_1d()
        self.output['Y(t) Euller Implícito'] = self.all_y_euller_implicito
      if 'ponto_medio' in self.methods:
        self.report_euller_1d()
        self.report_ponto_medio_1d()
        self.output['Y(t) Ponto médio'] = self.all_y_ponto_medio
      if 'heuen' in self.methods:
        self.report_heuen_1d()
        self.output['Y(t) Heuen'] = self.all_y_heuen
      if "rk4" in self.methods:
        self.report_rk4_1d()
        self.output['Y(t) Runge-Kutta 4ª Ordem'] = self.all_y_rk4
      if self.analytic_solution:
        self.analitica_1d()
        self.output['Y(t) Analítica'] = self.all_y_analitica
    if self.bidimensional or self.tridimensional or self.quadridimensional:
      if 'euller' in self.methods:
        self.report_euller()
        self.output['X(t) Euller'] = self.all_y_euller[:,0]
        self.output['Y(t) Euller'] = self.all_y_euller[:,1]
        if self.tridimensional or self.quadridimensional:
          self.output['Z(t) Euller'] = self.all_y_euller[:,2]
        if self.quadridimensional:
          self.output['W(t) Euller'] = self.all_y_euller[:,3]
      if 'heuen' in self.methods:
        self.report_heuen()
        self.output['X(t) Heuen'] = self.all_y_heuen[:,0]
        self.output['Y(t) Heuen'] = self.all_y_heuen[:,1]
        if self.tridimensional or self.quadridimensional:
          self.output['Z(t) Heuen'] = self.all_y_heuen[:,2]
        if self.quadridimensional:
          self.output['W(t) Heuen'] = self.all_y_heuen[:,3]
      if 'rk4' in self.methods:
        self.report_rk4()
        self.output['X(t) Runge-Kutta 4ª Ordem'] = self.all_y_rk4[:,0]
        self.output['Y(t) Runge-Kutta 4ª Ordem'] = self.all_y_rk4[:,1]
        if self.tridimensional or self.quadridimensional:
          self.output['Z(t) Runge-Kutta 4ª Ordem'] = self.all_y_rk4[:,2]
        if self.quadridimensional:
          self.output['W(t) Runge-Kutta 4ª Ordem'] = self.all_y_rk4[:,3]
      if self.analytic_solution:
        self.analitica()
        self.output['Y(t) Analítica'] = self.all_y_analitica[:,0]
        if self.bidimensional or self.tridimensional or self.quadridimensional:
          self.output['X(t) Analítica'] = self.all_y_analitica[:,0]
          self.output['Y(t) Analítica'] = self.all_y_analitica[:,1]
        if self.tridimensional or self.quadridimensional:
          self.output['Z(t) Analítica'] = self.all_y_rk4[:,2]
        if self.quadridimensional:
          self.output['W(t) Analítica'] = self.all_y_rk4[:,3]
      
    return self.output


  def plot_output(self, variable=None):
    if self.bidimensional:
      yaxis = "X(t), Y(t)"
    elif self.tridimensional:
      yaxis = "X(t), Y(t), Z(t)"
    elif self.quadridimensional:
      yaxis = "X(t), Y(t), Z(t), W(t)"
    else:
      yaxis = "Y(t)"

    fig = go.Figure()
    if variable == None:
      for column in self.output.columns[1:]:
          fig.add_trace(go.Scatter(x=self.output['Passo'], y=self.output[column],
                      mode='lines',
                      name=column))
    else:
      for column in self.output.columns[1:]:
          if column[:1] == variable:
            fig.add_trace(go.Scatter(x=self.output['Passo'], y=self.output[column],
                        mode='lines',
                        name=column))

      
    fig.update_layout(title_text='Resultados por Método', title_x=0.5,\
                xaxis_title='t', yaxis_title=f'{yaxis}',\
                height = 400, width = 600, font={'size':10})

    fig.show()
 