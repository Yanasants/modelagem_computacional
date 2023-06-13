import numpy as np
import plotly.graph_objects as go

# Gráfico 01
params = {"gk":36, "gna":120, "gl":0.3, "vk":12, "vna":-115, "vl":-10.613, "cm":1}
# Gráfico 02
params_2 = {"gk":24.27, "gna":120, "gl":0.3, "vk":12, "vna":-115, "vl":-10.613, "cm":1}

# N PARAMS
def alpha_n(V):
    alpha_n = 0.01*(V + 10)/(np.exp((V + 10)/10) - 1)
    return  alpha_n

def beta_n(V):
    beta_n = 0.125*np.exp(V/80)
    return beta_n

# M PARAMS
def alpha_m(V):
    alpha_m = 0.1*(V + 25)/(np.exp((V + 25)/10) - 1)
    return alpha_m

def beta_m(V):
    beta_m = 4*np.exp(V/18)
    return beta_m

# H PARAMS
def alpha_h(V):
    alpha_h = 0.07*np.exp(V/20)
    return  alpha_h

def beta_h(V):
    beta_h = 1.0 / (np.exp((V + 30.0) / 10.0) + 1.0)
    return beta_h

""" # CORRENTES
def i_Na(V, m, h):
    return params["gna"]*m**3 * h * (V - params["vna"])

def i_K(V, n):
    return params["gk"]* n**4 * (V - params["vk"])

def i_L(V):
    return params["gl"]* (V - params["vl"]) """

def dV_dt(V, m, h, n, I_app):
    return (I_app - params["gk"]*n**4*(V - params["vk"]) - params["gna"]*m**3*h* (V - params["vna"]) -  params["gl"]*(V - params["vl"]))/params["cm"]

## Calculo da Equcao 3 com u e (n, m, h). Onde alpha_u = alpha(V) e Beta_u = Beta_u(V)
def dn_dt(V, n):
    return alpha_n(V)*(1.0 - n) - beta_n(V)*n

def dm_dt(V, m):
    return alpha_m(V)*(1.0 - m) - beta_m(V)*m

def dh_dt(V, h):
    return alpha_h(V)*(1.0 - h) - beta_h(V)*h

def calc_current_density(V, m, h, n):
    # Ionica individual
    I_Na = params["gna"]*(m**3)*h*(V - params["vna"])
    I_K = params["gk"]*(n**4)*(V - params["vk"])
    I_L = params["gl"] * (V - params["vl"])
    
    # Calculo da Corrente Ionica Total
    I_total = I_Na + I_K + I_L
    return I_total

def I_function(V, m, h, n, I_app):
    return params["cm"]* dV_dt(V, m, h, n, I_app) + calc_current_density(V, m, h, n)

def plot_output(variable=None, title=None, dfs=None):
    fig = go.Figure()
    for df_key in dfs.keys():
        for column in dfs[df_key].output.columns[1:]:
            if column[:1] == variable:
                fig.add_trace(go.Scatter(x=dfs[df_key].output['Passo'], y=dfs[df_key].output[column],
                            mode='lines',
                            name=f'{column} | {df_key}'))
            
    fig.update_layout(title_text=f'{title}', title_x=0.5,\
                    xaxis_title='t', yaxis_title=f'{variable}(t)',\
                    height = 400, width = 600, font={'size':10})
    fig.show()