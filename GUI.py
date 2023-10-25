import PySimpleGUI as sg
import numpy as np
import matplotlib.pyplot as plt
from labellines import labelLine, labelLines
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import os
import threading
import time
from matplotlib import use as use_agg

from datetime import datetime, timedelta
import matplotlib.pylab as pylab
from main import *
params = {'legend.fontsize': 'x-small',
          'figure.figsize': (7, 5),
         'axes.labelsize': 'x-small',
         'axes.titlesize':'x-small',
         'xtick.labelsize':'x-small',
         'ytick.labelsize':'x-small'}
pylab.rcParams.update(params)

#### Initialization of useful data structures


energie_moy=[]
x=[]
k=0
config_dict={}
init_param=[]

plotting=False

date=datetime.now()


def Simulation(window,n_step):
    i=0
    global energie_moy
    global date
    global x
    global config_dict
    global k
    global simulation
    while simulation and k<n_step:
        liste = create_liste(N)
        for l in liste:
            old_state = l
            new_state = choose_new_state(config_dict, l, n_cut)[1]
            proba(old_state, new_state, Ex, Ey, Ez, config_dict)
        print(config_dict)
        e = 0
        for cle, valeur in config_dict.items():
            e += mp.sqrt(valeur[0]**2+valeur[1]**2+valeur[2]**2)/N
        x += [k]
        energie_moy += [e]
        #if (date-datetime.now()).total_seconds()>1:
         #   date=datetime.now()
        window.write_event_value('-THREAD-',  (threading.current_thread().name, i))
        time.sleep(0.1)
        k+=1

def pack_figure(graph, figure):
    canvas = FigureCanvasTkAgg(figure, graph.Widget)
    plot_widget = canvas.get_tk_widget()
    plot_widget.pack(side='top', fill='both', expand=1)
    return plot_widget

def plot_figure_first(index):
    fig = plt.figure(index)         # Active an existing figure
    ax = plt.gca()                  # Get the current axes
    ax.cla()                        # Clear the current axes
    ax.set_title(f"Mean energy")
    ax.set_xlabel("Step")
    ax.set_ylabel("Energy")
    ax.grid()
    ax.tick_params(axis='x', labelrotation = 20)
    plt.plot(x,energie_moy)            
    fig.canvas.draw()   

def plot_figure_second(index):
    fig = plt.figure(index)         # Active an existing figure
    ax = plt.gca()                  # Get the current axes
    ax.cla()                        # Clear the current axes
    ax.set_title(f"Momentum space")
    ax.set_xlabel('Nx')
    ax.set_ylabel('Ny')
    ax.set_zlabel('Nz')
    ax.grid()
    ax.tick_params(axis='x', labelrotation = 20)
    xu,yu,zu=[],[],[]
    xd,yd,zd=[],[],[]
    print(config_dict)
    for cle, valeur in config_dict.items():
        if valeur[3]==1:
            xu.append(valeur[0])
            yu.append(valeur[1])
            zu.append(valeur[2])
        if valeur[3]==-1:
            xd.append(valeur[0])
            yd.append(valeur[1])
            zd.append(valeur[2])
    print(zu,zd)
    ax.set_zlim3d(0,2)
    plt.scatter(xu,yu,zu, label='Up', marker='o',color='b')   
    plt.scatter(xd,yd,zd, label='Down', marker='^',color='r') 
    plt.legend()        
    fig.canvas.draw()   
  


use_agg('TkAgg')

dimframe=[ [sg.Text('Length x'), sg.Input(default_text='0.000001',key='-X-'), sg.Text("(m)")],[sg.Text('Length y'), sg.Input(default_text='0.000001',key='-Y-'), sg.Text("(m)")],[sg.Text('Length z'), sg.Input(default_text='0.000001',key='-Z-'), sg.Text("(m)")]]
DimFrame=sg.Frame("Dimensions of the box", dimframe, title_color='red')
layout = [[sg.Text("Welcome to Fermi Dirac simulator")],[sg.Text("Number of particles:"), sg.Input(default_text='10',key='-PART-', enable_events=True)], [sg.Text('Temperature', size =(15, 1)), sg.Input(default_text='1',key='-TEMP-', enable_events=True), sg.Text("(K)")],[sg.Text('Number of steps', size =(15, 1)), sg.Input(default_text='1000',key='-STEPS-', enable_events=True)],
        [DimFrame],[sg.Button("Start simulation", key='-START-')]]

main_window = sg.Window("Fermi Dirac MC Simulation", layout, finalize=True)


while True:
    window, event, values = sg.read_all_windows(timeout=1000)
    print(event)
    # End program if user closes window or
    # presses the OK button
    if event == "Exit" or event == sg.WIN_CLOSED or event == 'OK' or event=='-PARAM-':
        simulation=False
        plotting=False
        if window==main_window:
            break
        else:
            window.close()
    if event == '-START-':
        #try:
        simulation=True
        n_steps, N, T, Lx, Ly, Lz= float(values['-STEPS-']),int(values['-PART-']), float(values['-TEMP-']), float(values['-X-']),float(values['-Y-']),float(values['-Z-'])
        init_param=init(T,Lx,Ly,Lz,N)
        config_dict=init_states(N)
        n_cut=init_param[-1]
        Ex=init_param[0]
        Ey=init_param[1]
        Ez=init_param[2]
        param_column=[[sg.Text(f"Number of particles: {int(N)}")],[sg.Text(f"Temperature: {N} K")],[sg.Text(f"Dimensions: {Lx} x {Ly} x {Lz} m^3")],[sg.Text(f"State cut: {n_cut}")]]
        col=sg.Column(param_column)
        col_graph=sg.Column([[sg.Graph((4, 3), (0, 0), (4, 3), key='Graph1')], [sg.Graph((4, 3), (0, 0), (4, 3), key='Graph2')]])
        param_frame=sg.Frame("Parameters",[[col]], title_color='red')
        graph_frame=sg.Frame("Plots:",[[col_graph]],title_color='red')
        simul_layout=[[param_frame,graph_frame],[sg.Button("Change parameters", key='-PARAM-')]]
        simul_window=sg.Window('Simulation', simul_layout, finalize=True)
        #-----------------------------------------------------------------------
        #  Create graphs
        #-----------------------------------------------------------------------
        plotting=True
        graph1 = simul_window['Graph1']
        plt.ioff()                          
        fig1 = plt.figure(1,figsize=(4,3))
        ax1 = plt.subplot(111) 
        pack_figure(graph1, fig1)   
        plot_figure_first(1)
        graph2 = simul_window['Graph2']                        
        fig2 = plt.figure(2,figsize=(4,3))
        ax2 = fig2.add_subplot(111, projection='3d') 
        pack_figure(graph2, fig2)   
        plot_figure_second(2)
        monitoring_thread = threading.Thread(target=Simulation, args=(window,n_steps), daemon=True)
        monitoring_thread.start()
        
        #except: 
            #popup_layout=[[sg.Text("Please enter suitable parameters!", text_color='red')],[sg.Button("OK")]]
           # popup_wind=sg.Window('Warning', popup_layout, finalize=True, element_justification='c')
    if event == '-THREAD-' or event==sg.TIMEOUT_EVENT:
        if plotting:
            plot_figure_first(1)
            plot_figure_second(2)
