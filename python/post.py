from math import sqrt
import pyams.yams as yams
import plotly.graph_objects as go
import uuid
import numpy as np
from math import atan2

def add_span_to_plotly_fig(fig, g: yams.MeridionalGrid, i:int, value_name:str = 'Vm', scale = 1.):
    nj = g.nCols()
    x = []
    y = []
    l_tot = g(i,nj-1).l
    for j in range(nj):
        gp = g(i,j)
        y.append(gp.l / l_tot)
        try:
            x.append(getattr(gp, value_name)*scale)
        except AttributeError as e:
            if value_name == 'alf':
                x.append(atan2(gp.Vu,gp.Vm)*scale)
            else:
                raise e

    fig.add_trace(
        go.Scatter(x=x, y=y, name=f'{value_name} i = {i}') 
    )

def add_mg_to_plotly_fig(fig, g: yams.MeridionalGrid, value_name:str = 'Vm', zmin =None, zmax=None, scale = 1., rg = 287.04):
    
    nj = g.nCols()
    ni = g.nRows()

    a = []
    b = []
    x = []
    y = []
    z = []

    name = str(uuid.uuid1())

    for j in range(nj):
        for i in range(ni):
            a.append(i)
            b.append(j)
            gp = g(i,j)
            x.append(gp.x)
            y.append(gp.y)
            try:
                z.append(getattr(gp, value_name)*scale)
            except AttributeError as e:
                if value_name == 'alf':
                    z.append(atan2(gp.Vu,gp.Vm)*scale)
                if value_name == 'M':
                    z.append(sqrt(gp.Vm*gp.Vm+gp.Vu*gp.Vu)/sqrt(gp.gam*rg*gp.Ts)*scale)
                if value_name == 'Mm':
                    z.append(gp.Vm/sqrt(gp.gam*rg*gp.Ts)*scale)
                if value_name == 'Mu':
                    z.append(gp.Vu/sqrt(gp.gam*rg*gp.Ts)*scale)
                else:
                    raise e

    fig.add_trace(go.Carpet(
        a = a,
        b = b,
        x = x,
        y = y,
        aaxis = dict(
            smoothing = 0,
            minorgridcount = 0,
            type = 'linear',
            showticklabels = "none",
            gridcolor  = 'white',
        ),
        baxis = dict(
            smoothing = 0,
                minorgridcount = 0,
            type = 'linear',
            showticklabels = "none",
            gridcolor  = 'white',
        ),
        opacity = 0.3,
        carpet=name,
    ))

    fig.add_trace(go.Contourcarpet(
        a = a,
        b = b,
        z = z,
        contours = dict(
            showlabels=True,
            showlines=True,
        ),
        line = dict(
            width = 0.5,
            smoothing = 0.3,
        ),
        colorbar = dict(
            len = 0.6,
            y = 0.25,
            title=value_name,
        ),
        carpet = name,
        zmin = zmin,
        zmax = zmax,
    ))

    fig.update_yaxes(
        scaleanchor = "x",
        scaleratio = 1,
    )

    # fig.update_layout(
    #     template='plotly',
    #     coloraxis = {
    #         'colorbar':{'title':value_name},
    #         'colorscale':'portland',
    #     }
    # )

def grid_plotly_fig(g: yams.MeridionalGrid, value_name:str = 'Vm', width=1500, height=1000, scale = 1.):
    fig = go.Figure( )
    add_mg_to_plotly_fig(fig, g, value_name, scale= scale)

    fig.update_layout(
        width=width,
        height=height,
    )
    return fig

def case_plotly_fig(solver_case: yams.SolverCase, value_name:str = 'Vm', width=1500, height=1000, scale = 1.):
    fig = go.Figure( )
    g = solver_case.gi.g
    add_mg_to_plotly_fig(fig, solver_case.gi.g, value_name, scale= scale, rg =solver_case.gi.R)
    
    nj = solver_case.gi.nj

    for bld_info in solver_case.bld_info_lst:
        x_le = []
        y_le = []
        x_te = []
        y_te = []
        i1 = bld_info.i1
        i2 = bld_info.i2
        for j in range(nj):
            gp = g(i1,j)
            x_le.append(gp.x)
            y_le.append(gp.y)
            gp = g(i2,j)
            x_te.append(gp.x)
            y_te.append(gp.y)
        fig.add_trace(
            go.Scatter(x=x_le, y=y_le, line=dict(color='firebrick', width=4, dash='dash'), mode='lines')
        )
        fig.add_trace(
            go.Scatter(x=x_te, y=y_te, line=dict(color='royalblue', width=4, dash='dash'), mode='lines')
        )

    fig.update_layout(
        width=width,
        height=height,
        # template='plotly',
        # coloraxis = {
        #     'colorbar':{'title':value_name},
        #     'colorscale':'portland',
        # }
    )

    # fig.update_layout(
    #     template='plotly',
    #     coloraxis = {
    #         'colorbar':{'title':value_name},
    #         'colorscale':'portland',
    #     }
    # )
    return fig

def residuals_plotly_fig(solver_case: yams.SolverCase, width=800, height=600):
    n =len(solver_case.log.delta_pos_moy)
    x = np.linspace(0, n-1, n)
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=x, y=solver_case.log.delta_pos_moy,
            mode='lines',
            name='Average relative error')
    )
    fig.add_trace(
        go.Scatter(
            x=x, y=solver_case.log.delta_pos_max,
            mode='lines',
            name='Maximal relative error',
            line=dict( dash='dot')
        )
    )
    fig.update_yaxes(type="log")
    fig.update_layout(
        width=width,
        height=height,
    )
    return fig

def couturier_plotly_fig(solver_case: yams.SolverCase, iB:int, value_name:str = 'Vm', width=1500, height=1000, scale = 1. ):
    fig = go.Figure()
    g = solver_case.gi.g
    i1 = solver_case.bld_info_lst[iB].i1
    i2 = solver_case.bld_info_lst[iB].i2
    add_span_to_plotly_fig(fig, g, i1, value_name, scale)
    add_span_to_plotly_fig(fig, g, i2, value_name, scale)
    fig.update_layout(
        width=width,
        height=height,
    )
    return fig