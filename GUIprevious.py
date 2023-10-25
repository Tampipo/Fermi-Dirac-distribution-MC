import PySimpleGUI as sg
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import os
import matplotlib
import matplotlib.dates as mdates
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




layout = [[sg.Text("Welcome to Mo")],[sg.Text("Number of days to scrap"), sg.Input(key='-DAYS-', enable_events=True), sg.Button("Scrap data")], [sg.Text('Keyword', size =(15, 1)), sg.InputText(key='-IN-', enable_events=True)],
        [sg.Text('Number of articles', size =(15, 1)), sg.Input(key='-NUM-', enable_events=True), sg.Button("Search")]]

# Create the window
main_window = sg.Window("Fermi Dirac MC Simulation", layout, finalize=True)
titles=['','','','','','','','']
# Create an event loop
max=10
ndays=10

while True:
    window, event, values = sg.read_all_windows()
    print(event)
    # End program if user closes window or
    # presses the OK button
    if event == "Exit" or event == sg.WIN_CLOSED:
        if window==main_window:
            break
        else:
            window.close()
    if event == "Scrap data":
        layout_data = [
            [sg.Text("Welcome to the data scrapper...")],
            [sg.Text("Number of days to scrap"), sg.Input(key='-DAYS-', enable_events=True), sg.Button("Scrap data")],
            [sg.Image('C:\\Users\\Tanguy\\Documents\\Code\\ArXiVScrapper\\temp.png',expand_x=True, expand_y=True )],
            [sg.Input(key='-FIELD-', enable_events=True), sg.Button("New field")]
        ]
        window_data = sg.Window(
            "Data",
            layout_data,
            finalize=True,
        )
    if event == "New field":
        keys.append(values["-FIELD-"])
        PlotData(ndays)
        if window==window_data:
            window.close()
        layout_data = [
            [sg.Text("Welcome to the data scrapper...")],
            [sg.Text("Number of days to scrap"), sg.Input(key='-DAYS-', enable_events=True), sg.Button("Scrap data")],
            [sg.Image('C:\\Users\\Tanguy\\Documents\\Code\\ArXiVScrapper\\temp.png',expand_x=True, expand_y=True )],
            [sg.Input(key='-FIELD-', enable_events=True), sg.Button("New field")]
        ]
        window_data = sg.Window(
            "Data",
            layout_data,
            finalize=True,
        )
    if event == "Search":
        search = arxiv.Search(
        query = values["-IN-"],
        max_results = int(values["-NUM-"]),
        sort_by = arxiv.SortCriterion.SubmittedDate
        )
        titles=[]
        ids=[]
        urls=[]
        abstracts=[]
        authors=[]
        names=[]
        dates=[]
        max=int(values["-NUM-"])
        for result in search.results():
            titles.append(result.title)
            urls.append(result.pdf_url)
            abstracts.append(result.summary)
            ids.append(result.entry_id.split("/")[-1])
            authors.append(result.authors)
            dates.append(str(result.published).split(' ')[0])
        if len(titles)<int(values["-NUM-"]):
            max=len(titles)
        for i in range (0,max):
            string=''
            for name in authors[i]:
                string+=str(name).strip('arxiv.Result.Author(').strip(')')+' & '
            names.append(string[0:len(string)-2])
        search_layout=[[sg.Text(values['-IN-'], background_color='white', text_color='black')]]
        column=[]
        for i in range (0,max):
            if len(titles[i])>60:
                title=titles[i][0:60]+'...'
            else:
                title=titles[i]
            if len(names[i])>40:
                try:
                    name=names[i].split('&')[0].split(' ')[0][0]+'. '+names[i].split('&')[0].split(' ')[1]+' & Co'
                except:
                    name=names[i][0:40]
            else :
                name=names[i]
            column=column+[[sg.Text(title, text_color='green'),sg.Text(name,text_color='red'),sg.Text(dates[i]), sg.Button("See", key=f"DET{i}")]]
        search_layout=search_layout+[[sg.Column(column, scrollable=True,  vertical_scroll_only=True)]]
        search_window=sg.Window("Search results", search_layout, finalize=True)
    for i in range (0,max):
        if event == f"DET{i}":
            article_layout=[[sg.Text("Publication date: " +dates[i])],[sg.Text("Title:", text_color='green')],[sg.Text(titles[i])],[sg.Text("Authors:",text_color='red')],[sg.Text(names[i])], [sg.Text("Abstract")], [sg.Text(abstracts[i], background_color='white', text_color='black')], [sg.Button("Download", key=f"DL{i}")]]
            article_window=sg.Window("Article abstract", article_layout, finalize=True)
        if event == f"DL{i}":
            paper = next(arxiv.Search(id_list=[ids[i]]).results())
            paper.download_pdf(dirpath="C:\\Users\\Tanguy\\Documents\\Cours\\Cours de phyiques\\ArXiv\\")
            files=os.listdir("C:\\Users\\Tanguy\\Documents\\Cours\\Cours de phyiques\\ArXiv\\")
            for f in files:
                if ids[i] in f:
                    os.startfile("C:\\Users\\Tanguy\\Documents\\Cours\\Cours de phyiques\\ArXiv\\"+f)




