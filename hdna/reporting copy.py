import math
import pandas as pd
import numpy as np 
import plotly.graph_objects as go
import plotly.express as px

def upperapprox(num):
    order = math.floor(math.log10(num))
    return math.ceil(num/(10**order))*10**order

def lowerapprox(num):
    order = math.floor(math.log10(num))
    return math.floor(num/(10**order))*10**order

def num_ticks(lower, upper, tick_mag):
    return math.ceil((upper - lower + 1) / tick_mag)

def valplot(data, name, writepath = None, theme = 'dark'):
    
    THM = themetemplates(theme, 'scatter')

    sdata = data.sort_values('experimental')
    # sdata.set_index('seq', inplace=True, drop=False)
    try: sdata.drop('Unnamed: 0', axis=1)
    except KeyError: pass 
    X = list(sdata['experimental'].astype(float))
    Y = list(sdata['computational'].astype(float))
    S = list(sdata['sequences'])
    I = list(sdata['index'])

    lowbound = lowerapprox(min(Y))
    topbound = upperapprox(max(Y))

    XLINE = np.linspace(0, topbound)

    trace1 = go.Scatter(
        x = X, # My list of values for 'x'
        y = Y, # My list of values for 'y'
        marker = dict(color = THM['colordots']),
        mode = 'markers',
        name = '',
        customdata = [f'{i}_{s}' for i, s in zip(I, S)],
        hovertemplate="""emp = %{x:.3e}
                     <br>mod = %{y:.3e}
                     <br>seq:  %{customdata} </b>"""
    )
    trace2 = go.Scatter(
        x = XLINE,
        y = XLINE,
        mode = 'lines',
        marker = dict(color = THM['colorline']),
        line = dict(dash = 'dash'),
        name = 'bisector'
    )
    layout = go.Layout(
        template=THM['template'],
        title = f'Scatterplot for {name}',
        xaxis_title="empirical rates ",
        yaxis_title="computed rates",
        showlegend=False,
        autosize = False,
        width = 600,
        height = 600,
        margin = dict(
            b = 50,
            t = 50,
            l = 50,
            r = 50,
            pad = 0
        ),
        xaxis = dict(
            tickmode = 'array',
            tickvals = np.linspace(lowbound, topbound, num_ticks(lowbound, topbound, 1e6)),
            showgrid = True
            ),
        yaxis = dict(
            tickmode = 'array',
            tickvals = np.linspace(lowbound, topbound, num_ticks(lowbound, topbound, 1e6)),
            showgrid = True
        )
    )
    dados = [trace1, trace2]
    fig = go.Figure(data = dados, layout = layout)
    
    fig.update_xaxes(exponentformat="e", titlefont={'size': 22})
    fig.update_yaxes(exponentformat="e", titlefont={'size': 22})
    
    fig.update_layout(xaxis_range=[0,topbound])
    fig.update_layout(yaxis_range=[0,topbound])

    if writepath == None:
        fig.show()
        return fig 
    else:
        PATH = f'{writepath}/{name}'
        fig.write_html(f"{PATH}.html")

def valplot23(data, name, writepath = None, theme = 'dark'):
    
    THM = themetemplates(theme, 'scatter')

    sdata = data.sort_values('experimental')
    # sdata.set_index('seq', inplace=True, drop=False)
    try: sdata.drop('Unnamed: 0', axis=1)
    except KeyError: pass 
    X = list(sdata['experimental'].astype(float))
    Y3 = list(sdata['threedim'].astype(float))
    Y210 = list(sdata[sdata['length'] == 10]['twodim'].astype(float))
    Y214 = list(sdata[sdata['length'] == 14]['twodim'].astype(float))
    S = list(sdata['sequences'])
    I = list(sdata['index'])

    lowbound = lowerapprox(min(min(Y3), min(Y210), min(Y214)))
    topbound = upperapprox(max(max(Y3), max(Y210), max(Y214)))

    XLINE = np.linspace(0, topbound)

    trace3 = go.Scatter(
        x = X, # My list of values for 'x'
        y = Y3, # My list of values for 'y'
        marker = dict(color = THM['colordots']),
        mode = 'markers',
        name = '3D',
        customdata = [f'{i}_{s}' for i, s in zip(I, S)],
        hovertemplate="""emp = %{x:.3e}
                     <br>mod = %{y:.3e}
                     <br>seq:  %{customdata} </b>"""
    )
    trace210 = go.Scatter(
        x = X, # My list of values for 'x'
        y = Y210, # My list of values for 'y'
        marker = dict(color = '#636EFA'),
        mode = 'markers',
        name = '2D L=10',
        customdata = [f'{i}_{s}' for i, s in zip(I, S)],
        hovertemplate="""emp = %{x:.3e}
                     <br>mod = %{y:.3e}
                     <br>seq:  %{customdata} </b>"""
    )
    trace214 = go.Scatter(
        x = X, # My list of values for 'x'
        y = Y214, # My list of values for 'y'
        marker = dict(color = '#EE6650'),
        mode = 'markers',
        name = '2D L=14',
        customdata = [f'{i}_{s}' for i, s in zip(I, S)],
        hovertemplate="""emp = %{x:.3e}
                     <br>mod = %{y:.3e}
                     <br>seq:  %{customdata} </b>"""
    )
    traceline = go.Scatter(
        x = XLINE,
        y = XLINE,
        mode = 'lines',
        marker = dict(color = THM['colorline']),
        line = dict(dash = 'dash'),
        name = '3D identity'
    )
    layout = go.Layout(
        template=THM['template'],
        title = f'Scatterplot for {name}',
        xaxis_title="empirical rates ",
        yaxis_title="computed rates",
        showlegend=False,
        autosize = False,
        width = 600,
        height = 600,
        margin = dict(
            b = 50,
            t = 50,
            l = 50,
            r = 50,
            pad = 0
        ),
        xaxis = dict(
            tickmode = 'array',
            tickvals = np.linspace(lowbound, topbound, num_ticks(lowbound, topbound, 1e6)),
            showgrid = True
            ),
        yaxis = dict(
            tickmode = 'array',
            tickvals = np.linspace(lowbound, topbound, num_ticks(lowbound, topbound, 1e6)),
            showgrid = True
        )
    )
    dados = [trace3, trace210, trace214, traceline]
    fig = go.Figure(data = dados, layout = layout)
    
    fig.update_xaxes(exponentformat="e", titlefont={'size': 22})
    fig.update_yaxes(exponentformat="e", titlefont={'size': 22})
    
    fig.update_layout(xaxis_range=[0,topbound])
    fig.update_layout(yaxis_range=[0,topbound])

    
    fig.update_layout(legend=dict(
        orientation="h",
        entrywidth=70,
        yanchor="bottom",
        y=1.02,
        xanchor="right",
        x=1),
        showlegend=True,
        scattermode = 'group')

    if writepath == None:
        fig.show()
        return fig 
    else:
        PATH = f'{writepath}/{name}'
        fig.write_html(f"{PATH}.html")

def histotime(data, fit_gamma, runtime, exp=None, mod=None, seq=None, nbins=150, name='timehist', writepath=None, theme='light'):
    
    if theme == 'dark':
        THM = {'template': 'plotly_dark',
                    'colorfit': 'coral',
                    'colorbin': 'white',
                    'fitwidth': 3}
    if theme == 'light':
        THM = {'template': 'ggplot2',
                    'colorfit': 'royalblue',
                    'colorbin': 'tomato',
                    'fitwidth': 2.5}
    
    nbins = nbins
    X = np.linspace(0,runtime,500)
    data1 = go.Scatter(
                x=X, 
                y=fit_gamma.pdf(X),
                line = dict(
                    color = THM['colorfit'],
                    width = THM['fitwidth']))
    data2 = go.Histogram(
                x=data, 
                histnorm='probability density', 
                nbinsx=nbins,
                marker = dict(color = THM['colorbin']))
    
    layout = go.Layout(
                title="""seq %-s<br>mod: %-.2e<br>exp: %-.2e""" % (name.split('_')[0], float(mod), float(exp)),
                xaxis_title="simulation time",
                yaxis_title="# occurrences",
                template=THM['template'],
                xaxis=dict(tickformat=".0e"),
                yaxis=dict(tickformat=".0e"))

    fig = go.Figure(data=[data1, data2], layout=layout)
    fig.update_traces(opacity=0.95)

    if writepath == None:
        fig.show()
    else:
        PATH = f'{writepath}/{name}'
        fig.write_html(f"{PATH}.html")
    


def themetemplates(choice, kind):
    if kind == 'scatter' or 'percent':
        if choice == 'dark':
            return {'template': 'plotly_dark',
                    'colordots': 'tan',
                    'colorline': 'coral',
                    'linewidth': 2}
        if choice == 'light':
            return {'template': 'ggplot2',
                    'colordots': 'tomato',
                    'colorline': 'royalblue',
                    'linewidth': 2}
        

def percomplot(fpts, writepath=None, name=None, theme='dark'):

    THM = themetemplates(theme, 'percent')

    fig = px.ecdf(fpts, ecdfnorm='percent')
    fig.update_layout(
        title="percent completion by time",
        xaxis_title="simulation time",
        yaxis_title="percent",
        yaxis_ticksuffix = "%",
        xaxis=dict(tickformat=".0e"),
        showlegend=False,
        template=THM['template']
        )
    
    fig.data[0].line.color = THM['colorline']
    fig.data[0].line.width = THM['linewidth']

    if writepath == None:
        fig.show()
    else:
        PATH = f'{writepath}/{name}'
        fig.write_html(f"{PATH}.html")



class Tee(object):
    def __init__(self, *files):
        self.files = files
    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush() # If you want the output to be visible immediately
    def flush(self):
        for f in self.files:
            f.flush()
    def close(self):
        for f in self.files:
            f.close()