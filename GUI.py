import PySimpleGUI as sg
import numpy as np
import matplotlib.pyplot as plt
import tkinter
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
Fermi_Dirac=[[0,0]]
plotting=False

date=datetime.now()


def Simulation(window,n_step):
    t=0
    global energie_moy
    global date
    global x
    global config_dict
    global k
    global simulation
    global N
    while simulation and k<n_step:
        liste = create_liste(N)
        for l in liste:
            old_state = l
            new_state = choose_new_state(config_dict, l, n_cut)[1]
            proba(old_state, new_state, Ex, Ey, Ez, config_dict)
        #print(config_dict)
        e = 0
        E=[]
        for i in range (0,N):
            E.append(get_energy(i, config_dict,Ex,Ey,Ez))
        x += [k]
        E=np.array(E)
        energie_moy += [np.mean(E)]
        for i in range(0,len(Fermi_Dirac)):
            Fermi_Dirac[i][1]=k/(k+1)*Fermi_Dirac[i][1]+len(E[E==Fermi_Dirac[i][0]])/(k+1)
            E=E[E!=Fermi_Dirac[i][0]]

        treated_energy=[]

        for j in range (0,len(E)):
            if E[j] not in treated_energy:
                Fermi_Dirac.append([E[j],len(E[E==E[j]])/(k+1)])
                treated_energy.append(E[j])
        window.write_event_value('-THREAD-',  (threading.current_thread().name, t))
        time.sleep(0.1)
        k+=1

def pack_figure(graph, figure):
    canvas = FigureCanvasTkAgg(figure, graph.Widget)
    plot_widget = canvas.get_tk_widget()
    plot_widget.pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=True)
    return plot_widget

def plot_figure_first(index):
    global x
    global energie_moy
    fig = plt.figure(index)         # Active an existing figure
    ax = plt.gca()                  # Get the current axes
    ax.cla()                        # Clear the current axes
    ax.set_title(f"Mean energy")
    ax.set_xlabel("Step")
    ax.set_ylabel("Energy")
    ax.grid()
    ax.tick_params(axis='x', labelrotation = 20)
    plt.tight_layout()
    plt.plot(x,energie_moy)            
    fig.canvas.draw()   

def plot_figure_second(index):
    fig = plt.figure(index)         # Active an existing figure
    ax = plt.gca()                  # Get the current axes
    ax.cla()                        # Clear the current axes
    ax.set_title("Momentum space")
    ax.set_xlabel('Nx')
    ax.set_ylabel('Ny')
    ax.set_zlabel('Nz')
    #ax.tick_params(axis='x', labelrotation = 20)
    xu,yu,zu=[],[],[]
    xd,yd,zd=[],[],[]
    #print(config_dict)
    for i in range (0, len(config_dict)):
        valeur=config_dict[f'{i}']
        if valeur[3]==1:
            xu.append(valeur[0])
            yu.append(valeur[1])
            zu.append(valeur[2])
        if valeur[3]==-1:
            xd.append(valeur[0])
            yd.append(valeur[1])
            zd.append(valeur[2])
    print(zu)
    plt.scatter(xu,yu,zu, label='Up', marker='o',color='b')   
    plt.scatter(xd,yd,zd, label='Down', marker='^',color='r') 
    plt.legend()        
    fig.canvas.draw()   
  
def plot_figure_third(index):
    global Fermi_Dirac
    fig = plt.figure(index)         # Active an existing figure
    ax = plt.gca()                  # Get the current axes
    ax.cla()                        # Clear the current axes
    ax.set_title(f"Fermi_Dirac distribution")
    ax.set_xlabel("Energy")
    ax.set_ylabel("Number of particles")
    ax.grid()
    ax.tick_params(axis='x', labelrotation = 20)
    plt.tight_layout()
    Fermi_Dirac_array = np.array(Fermi_Dirac)
    plt.scatter(Fermi_Dirac_array[:,0],Fermi_Dirac_array[:,1])            
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
        print(10*T)
        init_param=init(T,Lx,Ly,Lz,N)
        config_dict=init_states(N)
        n_cut=init_param[-1]
        print(n_cut)
        Ex=init_param[0]
        Ey=init_param[1]
        Ez=init_param[2]
        param_column=[[sg.Text(f"Number of particles: {int(N)}")],[sg.Text(f"Temperature: {T} K")],[sg.Text(f"Dimensions: {Lx} x {Ly} x {Lz} m^3")],[sg.Text(f"State cut: {n_cut}")],[sg.Graph((4, 3), (0, 0), (4, 3), key='Graph3')]]
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
        ax2 = fig2.add_subplot(projection='3d') 
        pack_figure(graph2, fig2)   
        plot_figure_second(2)
        graph3 = simul_window['Graph3']
        fig3 = plt.figure(3,figsize=(4,3))
        ax3 = plt.subplot(111) 
        pack_figure(graph3, fig3)   
        plot_figure_third(3)
        monitoring_thread = threading.Thread(target=Simulation, args=(window,n_steps), daemon=True)
        monitoring_thread.start()
        
        #except: 
            #popup_layout=[[sg.Text("Please enter suitable parameters!", text_color='red')],[sg.Button("OK")]]
           # popup_wind=sg.Window('Warning', popup_layout, finalize=True, element_justification='c')
    if event == '-THREAD-' or event==sg.TIMEOUT_EVENT:
        if plotting:
            plot_figure_first(1)
            plot_figure_second(2)
            plot_figure_third(3)
