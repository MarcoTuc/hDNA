import pandas as pd
import numpy as np 
import plotly.graph_objects as go

def valplot(data, name, writepath = None, theme = 'dark'):
    
    if theme == 'dark':
        THM = {'template': 'plotly_dark',
               'colordots': 'tan',
               'colorline': 'coral'}
    if theme == 'light':
        THM = {'template': 'ggplot2',
               'colordots': 'tomato',
               'colorline': 'royalblue'}

    sdata = data.sort_values('expvalue')
    # sdata.set_index('seq', inplace=True, drop=False)
    try: sdata.drop('Unnamed: 0', axis=1)
    except KeyError: pass 
    X = list(sdata['expvalue'].astype(float))
    Y = list(sdata['computed'].astype(float))
    S = list(sdata['seq'])
    I = list(sdata.index)


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
        x = X,
        y = X,
        mode = 'lines',
        marker = dict(color = THM['colorline']),
        line = dict(dash = 'dash'),
        name = 'bisector'
    )
    layout = go.Layout(
        template=THM['template'],
        title = f'Scatterplot for {name}',
        showlegend=False,
        autosize = False,
        width = 1000,
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
            showgrid = True
            ),
        yaxis = dict(
            tickmode = 'array',
            showgrid = True
        )
    )
    dados = [trace1, trace2]
    fig = go.Figure(data = dados, layout = layout)
    
    if writepath == None:
        fig.show()
    else:
        PATH = f'{writepath}/{name}'
        fig.write_image(f"{PATH}.png")
        fig.write_html(f"{PATH}.html")

