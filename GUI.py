import PySimpleGUI as sg
import numpy as np
import matplotlib.pyplot as plt
from labellines import labelLine, labelLines
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import os
import threading
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

def pack_figure(graph, figure):
    canvas = FigureCanvasTkAgg(figure, graph.Widget)
    plot_widget = canvas.get_tk_widget()
    plot_widget.pack(side='top', fill='both', expand=1)
    return plot_widget

def plot_figure_first(index):
    fig = plt.figure(index)         # Active an existing figure
    ax = plt.gca()                  # Get the current axes
    ax.cla()                        # Clear the current axes
    ax.set_title(f"HV channels 1-6")
    ax.set_xlabel("Time")
    ax.set_ylabel("V")
    ax.grid()
    ax.tick_params(axis='x', labelrotation = 20)
    for i in range (0,6):
        plt.plot(times, volts[i], label=f'Ch{i}') 
    try:
        labelLines(plt.gca().get_lines(), zorder=2.5)
    except: pass            
    fig.canvas.draw()   

use_agg('TkAgg')

dimframe=[ [sg.Text('Length x'), sg.Input(default_text='0.000001',key='-X-'), sg.Text("(m)")],[sg.Text('Length y'), sg.Input(default_text='0.000001',key='-Y-'), sg.Text("(m)")],[sg.Text('Length z'), sg.Input(default_text='0.000001',key='-Z-'), sg.Text("(m)")]]
DimFrame=sg.Frame("Dimensions of the box", dimframe, title_color='red')
layout = [[sg.Text("Welcome to Fermi Dirac simulator")],[sg.Text("Number of particles:"), sg.Input(default_text='100',key='-PART-', enable_events=True)], [sg.Text('Temperature', size =(15, 1)), sg.Input(default_text='100',key='-TEMP-', enable_events=True), sg.Text("(K)")],
        [DimFrame],[sg.Button("Start simulation", key='-START-')]]

main_window = sg.Window("Fermi Dirac MC Simulation", layout, finalize=True)


while True:
    window, event, values = sg.read_all_windows()
    print(event)
    # End program if user closes window or
    # presses the OK button
    if event == "Exit" or event == sg.WIN_CLOSED or event == 'OK' or event=='-PARAM-':
        if window==main_window:
            break
        else:
            window.close()
    if event == '-START-':
        try:
            N, T, Lx, Ly, Lz= float(values['-PART-']), float(values['-TEMP-']), float(values['-X-']),float(values['-Y-']),float(values['-Z-'])
            init_param=init(T,Lx,Ly,Lz)
            n_cut=init_param[-1]
            Ex=init_param[0]
            Ey=init_param[1]
            Ez=init_param[2]
            param_column=[[sg.Text(f"Number of particles: {int(N)}")],[sg.Text(f"Temperature: {N} K")],[sg.Text(f"Dimensions: {Lx} x {Ly} x {Lz} m^3")],[sg.Text(f"State cut: {n_cut}")]]
            col=sg.Column(param_column)
            param_frame=sg.Frame("Parameters",[[col]], title_color='red')
            simul_layout=[[param_frame],[sg.Button("Change parameters", key='-PARAM-')]]
            simul_window=sg.Window('Simulation', simul_layout, finalize=True)
            
        except: 
            popup_layout=[[sg.Text("Please enter suitable parameters!", text_color='red')],[sg.Button("OK")]]
            popup_wind=sg.Window('Warning', popup_layout, finalize=True, element_justification='c')

 

