import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

class Edo:

  def __init__(self, y_inicial, dom, a,b, f_function, \
               analytic_function=None, all_error=[], output=pd.DataFrame([]), \
              analytic_solution=False, bidimensional=False):
    self.y = y_inicial
    self.dom = dom
    self.a = a
    self.b = b
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
    if self.analytic_solution:
        self.report_euller_2d()
        self.report_rk4_2d()
        self.analitica()
        self.output = pd.DataFrame({'Passo':self.all_x, 'Ct_Euller': self.all_y_euller_2d[:,0], \
                                  'Nt_Euller': self.all_y_euller_2d[:,1],'Ct_RK4':self.all_y_rk4_2d[:,0], \
                                  'Nt_RK4': self.all_y_rk4_2d[:,1],\
                                  'Ct_Analitica': self.all_y_analitica[:,0], \
                                  'Nt_Analitica':self.all_y_analitica[:,1]})

        return self.output

    if self.bidimensional and self.analytic_solution==False:
      self.report_euller_2d()
      self.report_rk4_2d()
      self.output = pd.DataFrame({'Passo':self.all_x, 'Ct_Euller': self.all_y_euller_2d[:,0], \
                                  'Nt_Euller': self.all_y_euller_2d[:,1],'Ct_RK4':self.all_y_rk4_2d[:,0], \
                                  'Nt_RK4': self.all_y_rk4_2d[:,1]})

      return self.output


  def plot_output(self, variable=None):
    fig = go.Figure()
    if self.analytic_solution:
      if variable=="C" or variable==None:
        fig.add_trace(go.Scatter(x=self.output['Passo'], y=self.output['Ct_Analitica'],
                          mode='lines',
                          name='C(t) Analitica'))
      if variable=="N" or variable==None:
        fig.add_trace(go.Scatter(x=self.output['Passo'], y=self.output['Nt_Analitica'],
                    mode='lines',
                    name='N(t) Analitica'))
      
    if variable=="C" or variable==None:
      fig.add_trace(go.Scatter(x=self.output['Passo'], y=self.output['Ct_Euller'],
                        mode='lines',
                        name='C(t) Euller'))
      fig.add_trace(go.Scatter(x=self.output['Passo'], y=self.output['Ct_RK4'],
                  mode='lines',
                  name='C(t) Runge-Kutta 4ª Ordem'))
    if variable=="N" or variable==None:
      fig.add_trace(go.Scatter(x=self.output['Passo'], y=self.output['Nt_Euller'],
                        mode='lines',
                        name='N(t) Euller',
                        line = {'color': '#341f97'}))
      fig.add_trace(go.Scatter(x=self.output['Passo'], y=self.output['Nt_RK4'],
                        mode='lines',
                        name='N(t) Runge-Kutta 4ª Ordem',
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