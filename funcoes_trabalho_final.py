import numpy as np
import plotly.graph_objects as go

# Gráfico 01
params = {"gk":36, "gna":120, "gl":0.3, "vk":12, "vna":-115, "vl":-10.613, "cm":0.1}
# Gráfico 02
params_2 = {"gk":24.27, "gna":120, "gl":0.3, "vk":12, "vna":-115, "vl":-10.613, "cm":1}
# Gráfico 03
params_3 = {"gk":24.27, "gna":71.755, "gl":0.3, "vk":12, "vna":-115, "vl":-10.613, "cm":1}
# Gráfico 04
params_4 = {"gk":24.27, "gna":71.828, "gl":0.3, "vk":12, "vna":-115, "vl":-10.613, "cm":1}




# N PARAMS
def alpha_n(V):
    alpha_n = 0.01*(V + 10)/(np.exp((V + 10)/10) - 1)
    return  alpha_n

def beta_n(V):
    beta_n = 0.125*np.exp(V/80)
    return beta_n

def dn_dt(V, n):
    return alpha_n(V)*(1.0 - n) - beta_n(V)*n



# M PARAMS
def alpha_m(V):
    alpha_m = 0.1*(V + 25)/(np.exp((V + 25)/10) - 1)
    return alpha_m

def beta_m(V):
    beta_m = 4*np.exp(V/18)
    return beta_m

def dm_dt(V, m):
    return alpha_m(V)*(1.0 - m) - beta_m(V)*m

def dm_params_fixo(V, m):
    return 3.82*(1.0 - m) - 0.15*m




# H PARAMS
def alpha_h(V):
    alpha_h = 0.07*np.exp(V/20)
    return  alpha_h

def beta_h(V):
    beta_h = 1.0 / (np.exp((V + 30.0) / 10.0) + 1.0)
    return beta_h

def dh_dt(V, h):
    return alpha_h(V)*(1.0 - h) - beta_h(V)*h

def dh_params_fixo(V, h):
    return 1.19*h 




def dV_dt(V, n, m, h, I_app, graph=1):
    if graph==1:
        return (I_app - params["gk"]*n**4*(V - params["vk"]) - params["gna"]*m**3*h* (V - params["vna"]) -  params["gl"]*(V - params["vl"]))/params["cm"]
    elif graph==2:
        return (I_app - params_2["gk"]*n**4*(V - params_2["vk"]) - params_2["gna"]*m**3*h* (V - params_2["vna"]) -  params_2["gl"]*(V - params_2["vl"]))/params_2["cm"]
    elif graph==3:
        return (I_app - params_3["gk"]*n**4*(V - params_3["vk"]) - params_3["gna"]*m**3*h* (V - params_3["vna"]) -  params_3["gl"]*(V - params_3["vl"]))/params_3["cm"]
    elif graph==4:
        return (I_app - params_4["gk"]*n**4*(V - params_4["vk"]) - params_4["gna"]*m**3*h* (V - params_4["vna"]) -  params_4["gl"]*(V - params_4["vl"]))/params_4["cm"]





def calc_current_density(V, n, m, h):
    # Ionica individual
    I_Na = params["gna"]*(m**3)*h*(V - params["vna"])
    I_K = params["gk"]*(n**4)*(V - params["vk"])
    I_L = params["gl"] * (V - params["vl"])
    
    # Calculo da Corrente Ionica Total
    I_total = I_Na + I_K + I_L
    return I_total


def gk_function(V, n, m, h):
    Ik = params_2["gk"]*(n**4)*(V - params_2["vk"])
    gk = Ik/(V-params_2["vk"])*(n**4)
    return gk

def gna_function(V, n, m, h):
    I_Na = params_2["gna"]*(m**3)*h*(V - params_2["vna"])
    gna = I_Na/(V-params_2["vna"])*(m**3)*h
    return gna


def I_function_total(V, n, m,h):
    return params["cm"]* V + calc_current_density(V, n, m, h)

def plot_output(variable=None, title=None, dfs=None):
    fig = go.Figure()
    for df_key in dfs.keys():
        for column in dfs[df_key].output.columns[1:]:
            if column[:1] == variable:
                fig.add_trace(go.Scatter(x=dfs[df_key].output['Passo'], y=dfs[df_key].output[column],
                            mode='lines',
                            name=f'{column} | {df_key}'))
            
    fig.update_layout(title_text=f'{title}', title_x=0.5,\
                    xaxis_title='Tempo (ms)', yaxis_title=f'{variable}(t)',\
                    height = 400, width = 600, font={'size':10})
    fig.show()

def plot_condutancias_rk4(variable=None, dfs=None):
    fig = go.Figure()
    for df_key in dfs.keys():
        if variable=="gk":
            var_k = gk_function(V=dfs[df_key].V, m=dfs[df_key].m, n=dfs[df_key].n, h=dfs[df_key].h)
            element = "potássio"
        elif variable=="gna":
            var_k = gna_function(V=dfs[df_key].V, m=dfs[df_key].m, n=dfs[df_key].n, h=dfs[df_key].h)
            element = "sódio"

        fig.add_trace(go.Scatter(x=dfs[df_key].output['Passo'], y=var_k,
                    mode='lines',
                    name=f'{variable} {df_key}'))
        
    fig.update_layout(title_text=f'Condutância do {element}', title_x=0.5,\
                    xaxis_title='Tempo (ms)', yaxis_title=f'{variable} (t)',\
                    height = 400, width = 600, font={'size':10})
    fig.show()

def plot_condutancias(variable=None, dfs=None):
    fig = go.Figure()
    for df_key in dfs.keys():
        if variable=="gk":
            var_k = gk_function(V=dfs[df_key].all_V, m=dfs[df_key].all_m, n=dfs[df_key].all_n, h=dfs[df_key].all_h)
            element = "potássio"
        elif variable=="gna":
            var_k = gna_function(V=dfs[df_key].all_V, m=dfs[df_key].all_m, n=dfs[df_key].all_n, h=dfs[df_key].all_h)
            element = "sódio"

        fig.add_trace(go.Scatter(x=dfs[df_key].output['Passo'], y=var_k,
                    mode='lines',
                    name=f'{variable} {df_key}'))
        
    fig.update_layout(title_text=f'Condutância do {element}', title_x=0.5,\
                    xaxis_title='Tempo (ms)', yaxis_title=f'{variable} (t)',\
                    height = 400, width = 600, font={'size':10})
    fig.show()