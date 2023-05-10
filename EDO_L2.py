import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

class Edo:

  def __init__(self, y_inicial, dom, a,b, f_function, methods=[], \
               analytic_function=None, all_error=[], output=pd.DataFrame([]), \
              analytic_solution=False, bidimensional=False):
    self.y = y_inicial
    self.dom = dom
    self.a = a
    self.b = b
    self.methods = methods
    self.f_function = f_function
    self.analytic_solution = analytic_solution
    self.analytic_function = analytic_function
    self.bidimensional = bidimensional

    self.h = (self.b - self.a)/self.dom
    self.size = self.dom-1
    self.all_x = np.arange(self.a, self.b, self.h)
    self.all_y = []
    self.all_y_euller =[self.y]
    self.all_y_analitica = []
    self.yh = [self.y]
    self.y_pm = [self.y]
    self.y_newton = [self.y]

    self.all_y_euller_2d = None
    self.all_y_rk4_2d = None
    self.all_error = all_error
    self.output = output

  def f(self, x, y):
    return self.f_function(x, y)

  def report_euller(self):
    for i in range(len(self.all_x)):
      if len(self.all_y_euller)<=self.size:
        self.all_y_euller.append(self.all_y_euller[i]+self.h*self.f(self.all_x[i],self.all_y_euller[i]))
        
    if self.analytic_solution:
      self.all_y_analitica = list(map(self.analitica, self.all_x))
      self.all_error = abs(np.array(self.all_y_euller)-np.array(self.all_y_analitica))
  

  def report_euller_implicito(self):
    for i in np.arange(0, len(self.all_x)):
      if len(self.yh)<=self.size:
          self.k1 = self.h*self.f(self.all_x[i],self.all_y_euller[i])
          self.k2 = self.h*self.f(self.all_x[i+1],self.all_y_euller[i]+self.k1)
          self.yh.append(self.yh[i] + ((self.k1+self.k2)/2))

  def report_ponto_medio(self):
    for i in np.arange(0, len(self.all_x)):
      if len(self.y_pm)<=self.size:
          self.k1_pm = self.f(self.all_x[i],self.all_y_euller[i])
          self.k2_pm = self.f(self.all_x[i]+(self.h)/2,self.all_y_euller[i]+self.k1_pm*(self.h/2))
          self.y_pm.append(self.y_pm[i] + self.k2_pm*self.h)
          
  def report_newton(self):
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

  def report_euller_2d(self):
    self.all_y_euller_2d = np.zeros((len(self.all_x), len(self.y)))
    k = np.zeros((len(self.all_x), len(self.y)))
    self.all_y_euller_2d[0] = self.y

    for i in range(len(self.all_x)-1):
      k[i] = self.f(self.all_x[i], self.all_y_euller_2d[i])
      self.all_y_euller_2d[i+1] = self.h*k[i]+self.all_y_euller_2d[i]


  def report_rk4_2d(self):
    self.all_y_rk4_2d = np.zeros((len(self.all_x), len(self.y)))
    self.all_y_rk4_2d[0] = self.y
    for i in range(len(self.all_x)-1):
      k1 = self.h*np.array((self.f(self.all_x[i], self.all_y_rk4_2d[i])))
      k2 = self.h*np.array(self.f(self.all_x[i]+self.h/2, self.all_y_rk4_2d[i]+k1/2))
      k3 = self.h*np.array(self.f(self.all_x[i]+self.h/2, self.all_y_rk4_2d[i]+k2/2))
      k4 = self.h*np.array(self.f(self.all_x[i]+self.h, self.all_y_rk4_2d[i]+k3))
      self.all_y_rk4_2d[i+1] = self.all_y_rk4_2d[i] + (k1 + 2*k2 +2*k3 + k4)/6

  def f_analitica(self, x, y):
    return self.analytic_function(x, y)

  def analitica(self):
    self.all_y_analitica = np.zeros((len(self.all_x), len(self.y)))
    k = np.zeros((len(self.all_x), len(self.y)))
    self.all_y_analitica[0] = self.y
    for i in range(len(self.all_x)-1):
      k[i] = self.f_analitica(self.all_x[i], self.all_y_analitica[i])
      self.all_y_analitica[i+1] = k[i]
      

  @property
  def df_output(self):
    self.output['Passo'] = self.all_x
    if self.bidimensional == False:
      if 'euller' in self.methods:
        self.report_euller()
        self.output['Y_Euller'] = self.all_y_euller
      if 'euller_implicito' in self.methods:
        self.report_euller_implicito()
    if self.bidimensional:
      if 'euller' in self.methods:
        self.report_euller_2d()
        self.output['Xt_Euller'] = self.all_y_euller_2d[:,0]
        self.output['Yt_Euller'] = self.all_y_euller_2d[:,1]
      if 'rk4' in self.methods:
        self.report_rk4_2d()
        self.output['Xt_RK4'] = self.all_y_rk4_2d[:,0]
        self.output['Yt_RK4'] = self.all_y_rk4_2d[:,1]

    if self.analytic_solution:
      self.analitica()
      self.output['Xt_Analitica'] = self.all_y_analitica[:,0]
      self.output['Yt_Analitica'] = self.all_y_analitica[:,1]

    return self.output


  def plot_output(self, variable=None):
    fig = go.Figure()
    if self.analytic_solution:
      if variable=="C" or variable==None:
        fig.add_trace(go.Scatter(x=self.output['Passo'], y=self.output['Xt_Analitica'],
                          mode='lines',
                          name='X(t) Analitica'))
      if variable=="N" or variable==None:
        fig.add_trace(go.Scatter(x=self.output['Passo'], y=self.output['Yt_Analitica'],
                    mode='lines',
                    name='Y(t) Analitica'))
      
    if variable=="X" or variable==None:
      fig.add_trace(go.Scatter(x=self.output['Passo'], y=self.output['Xt_Euller'],
                        mode='lines',
                        name='X(t) Euller'))
      fig.add_trace(go.Scatter(x=self.output['Passo'], y=self.output['Xt_RK4'],
                  mode='lines',
                  name='X(t) Runge-Kutta 4ª Ordem'))
    if variable=="Y" or variable==None:
      fig.add_trace(go.Scatter(x=self.output['Passo'], y=self.output['Yt_Euller'],
                        mode='lines',
                        name='Y(t) Euller',
                        line = {'color': '#341f97'}))
      fig.add_trace(go.Scatter(x=self.output['Passo'], y=self.output['Yt_RK4'],
                        mode='lines',
                        name='Y(t) Runge-Kutta 4ª Ordem',
                        line = {'color': '#FF914D'}))
    
    fig.update_layout(title_text='Resultados por Método', title_x=0.5,\
                      xaxis_title='t', yaxis_title=f'{variable}(t)',\
                      height = 400, width = 600, font={'size':10})
    fig.show()
 
  def plot_error_distribution(self):
    if self.analytic_solution:
      fig = px.histogram(self.output, x="Passo", y="Erro", nbins=self.dom,\
                        color_discrete_sequence=['indianred'])
      fig.update_layout(title_text='Distribuição de Erro', title_x=0.5,\
                        xaxis_title='Passo', yaxis_title='Erro',\
                        height = 400, width = 600, font={'size':10})
      fig.show()
