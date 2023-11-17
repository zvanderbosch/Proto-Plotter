import os
import matplotlib, sys
import numpy as np
import matplotlib.pyplot as plt

from glob import glob
from tkinter import *
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
from matplotlib.figure import Figure

matplotlib.use('TkAgg')

"""
A Python GUI built using TKinter/TKK for visualizing
photometric data for young stellar objects (YSOs) and
fitting them with simple one or two component blackbody 
models.

Simply run this script from the command line, but you
must also have the "object_fluxes.csv" file in the 
data directory for this script to work properly.

Author: Zach VAnderbosch
Date: 2018-09-17
"""

# Define some physical constants
h    = 6.626 * 1e-27  ## Planck constant (cgs)
c    = 2.998 * 1e10   ## Speed of light (cgs)
kb   = 1.381 * 1e-16  ## Boltzmann constant (cgs)
jy   = 1e-23          # Conversion from Jy to cgs units

# Define the wavelength range
lstart = -1.0  # In units of micrometers
lfinis = 4.0   # In units of micrometers
larr = np.logspace(lstart,lfinis,200)

# Load in photometric data
wave_vals = np.asarray([0.445,0.551,0.658,0.801,1.25,1.65,
	                    2.2,3.6,4.5,5.8,8.0,24.0,70.0])
freq_vals = (c / wave_vals) * 1e4
path = os.getcwd() + '/'
fname = glob('data/object_fluxes.csv')
flux = np.loadtxt(path+fname[0],dtype=float,
	              usecols=(1,2,3,4,5,6,7,8,9,10,11,12,13),
	              delimiter=",")
obj = np.loadtxt(path+fname[0], dtype=str, 
				 usecols=0, delimiter=",")

# Convert flux values from Jansky to CGS units
flux_cgs = (flux*1e-3) * jy * freq_vals

def _bb_func(lamb,T):
	global c, h, kb
	# Concert wavelength values to centimeters
	norm = 1e-24  # ballpark distance normalization
	lamb_cm = lamb * 1e-4
	b_lam  = 4.0 * lamb * norm * (2.0*h*(c**2)/(lamb_cm**5)) * \
			 ((np.exp((h*c)/(lamb_cm*kb*T)) - 1.0)**(-1.0))
	return b_lam

def _plot_data():
	global wave_vals, flux_vals, ipd, m4

	if ipd == 0:
		if object_chosen.get() == '-- Choose an Object --':
			m4, = a.plot([],[],ls='None',marker='s',ms=6,
						 mew=1,mec='k',mfc='w',zorder=0)
		else:
			flux_vals = flux_cgs[np.where(obj == object_chosen.get())]
			m4, = a.plot(wave_vals,flux_vals,ls='None', marker='s',
				         ms=6,mew=1,mec='k',mfc='w',zorder=0)
			yu = a.get_ylim()[1]
			if 5.0*max(flux_vals[0]) > yu:
				a.set_ylim(1e-15,5.0*max(flux_vals[0]))
		ipd = 1
	else:
		if object_chosen.get() == '-- Choose an Object --':
			m4.set_data([],[])
		else:
			flux_vals = flux_cgs[np.where(obj == object_chosen.get())]
			m4.set_data(wave_vals,flux_vals)
			yu = a.get_ylim()[1]
			if 5.0*max(flux_vals[0]) > yu:
				a.set_ylim(1e-15,5.0*max(flux_vals[0]))


def _plot():
	global a, larr, T, check2, ip, m1, m2, m3

	# Grab values from entry boxes
	T1 = float(temp_entry1.get())
	T2 = float(temp_entry2.get())
	#Grab values from sliders
	N1 = float(var4.get())
	N2 = float(var5.get())

	# Generate arrays of x and y values
	yvals1 = _bb_func(larr,T1) * 10.0**(N1)
	yvals2 = _bb_func(larr,T2) * 10.0**(N2) * 1e2
	ycomb = yvals1 + yvals2

	# Create each plot
	if ip == 0:
		if check2.get() == 1:
			m1, = a.plot(larr,yvals1,ls='-',lw=1.5,c='C0')
			m2, = a.plot(larr,yvals2,ls='-',lw=1.5,c='C3')
			m3, = a.plot(larr,ycomb ,ls='--',lw=1.5,c='k')
			a.set_ylim(1e-15,5.0*max(ycomb))
		else:
			m1, = a.plot(larr,yvals1,ls='-',lw=1.5,c='C0')
			m2, = a.plot([],[],ls='-' ,lw=1.5,c='C3')
			m3, = a.plot([],[],ls='--',lw=1.5,c='k')
			a.set_ylim(1e-15,5.0*max(yvals1))
		_plot_data()
		ip = 1

	# Update data for each plot
	else:
		if check2.get() == 1:
			m1.set_data(larr,yvals1)
			m2.set_data(larr,yvals2)
			m3.set_data(larr,ycomb)
			a.set_ylim(1e-15,5.0*max(ycomb))
		else:
			m1.set_data(larr,yvals1)
			m2.set_data([],[])
			m3.set_data([],[])
			a.set_ylim(1e-15,5.0*max(yvals1))
		_plot_data()
	dataPlot.draw()
	return

def _check_switch():
	global check2

	if check2.get() == 1:
		temp_entry2.configure(state=NORMAL)
		scroll_disc.configure(state=NORMAL, troughcolor='coral1', 
			                  sliderrelief='raised')
		_plot()
	else:
		temp_entry2.configure(state=DISABLED)
		scroll_disc.configure(state=DISABLED, troughcolor='White', 
			                  sliderrelief='flat')
		_plot()
	return

# Command tied 
def _set_object(eventObject):
	#print(object_chosen.get())
	_plot()
	return

# Command tied to "Clear Data" button
def _clear_data():
	object_chosen.set('-- Choose an Object --')
	_plot()
	return

# Command tied to sliders
def _update_value(eventObject):
	_plot()
	return


#####################################################
###############  Configure the GUI  #################

# Create master GUI
window = Tk()
window.title("Proto Plot 3000")
window.configure(bg='#a1dbcd')

# Create Figure Widget
f = Figure(figsize=(9.0,6.0), dpi=100)
a = f.add_subplot(111)
ip = 0
ipd = 0

dataPlot = FigureCanvasTkAgg(f, master=window)
dataPlot.get_tk_widget().grid(row=1, column=0)
f.set_canvas(dataPlot)
toolbar = Frame(window)
toolbar.grid(column=0, row=0)

# Configure the plot to look so good
a.grid(color='k', linestyle='-', linewidth=0.5, alpha=0.3)
a.set_xscale('log')
a.set_yscale('log')
a.set_xlabel(r'Wavelength $(\mu m)$',fontsize=14)
a.set_ylabel(r'Flux Density (cgs)',fontsize=14)
a.set_xlim(0.1,500)
a.minorticks_on()
a.tick_params(which='minor',direction='out',length=3,width=1.0)
a.tick_params(which='major',direction='out',length=5,width=1.1,labelsize=12)
[sp.set_linewidth(1.3) for sp in a.spines.values()]

# Define Styles for Entry Boxes and Sliders
style_temp1 = ttk.Style()
style_temp2 = ttk.Style()
style_block = ttk.Style()
style_temp1.configure("Temp1.TEntry", foreground="black", background="light blue")
style_temp2.configure("Temp2.TEntry", foreground="black", background="coral1")
style_temp2.map("Temp2.TEntry",foreground=[('disabled','silver')],background=[('disabled','white')])

# drop down menu for region and object selection
ttk.Label(toolbar, text="Data Selection").grid(column=0, row=1, sticky="S")
object_chosen = ttk.Combobox(toolbar, width=17)
object_chosen.bind('<<ComboboxSelected>>',_set_object)
object_chosen.grid(column=0, row=2, sticky="W")
object_chosen['values'] = list(np.unique(obj))
object_chosen.set('-- Choose an Object --')

# Add some labels
ttk.Label(toolbar, text="Temperature (Kelvin)").grid(column=1, row=1, sticky="E")
ttk.Label(toolbar, text="Intensity").grid(column=1, row=2, sticky="E")
ttk.Label(toolbar, text="Central Object Blackbody").grid(column=2, row=0, sticky="S")

# Entry Box for Central Object Temperature
temp_entry1 = ttk.Entry(toolbar, style="Temp1.TEntry")
temp_entry1.grid(column=2, row=1)
temp_entry1.insert(END,'3500') # Default Value

# Entry Box for Disc Temperature
temp_entry2 = ttk.Entry(toolbar, style="Temp2.TEntry")
temp_entry2.grid(column=3, row=1)
temp_entry2.insert(END,'300') # Default Value

# slider for Central Object Intensity
var4 = DoubleVar()
scroll_cobj = Scale(toolbar, orient="horizontal", length=180, 
                  showvalue=0, digits=4, from_=-3, to=3, resolution=0.1,
                  variable=var4, bg='#E4E4E4', troughcolor='light blue',
                  command=_update_value)
scroll_cobj.set(0.0)
scroll_cobj.grid(column=2, row=2)

# slider for disc Intensity
var5 = DoubleVar()
scroll_disc = Scale(toolbar, orient="horizontal", length=180, 
                  showvalue=0, digits=4, from_=-10, to=10, resolution=0.1,
                  variable=var5, bg='#E4E4E4', troughcolor='coral1',
                  command=_update_value)
scroll_disc.set(0.0)
scroll_disc.grid(column=3, row=2)


# Checkbox for adding in second blackbody
check2 = IntVar()
check_bb2 = ttk.Checkbutton(toolbar, text="Add Second Blackbody", variable=check2, command=_check_switch)
check_bb2.grid(column=3, row=0, sticky="S")
_check_switch()

# button to clear data
action = ttk.Button(toolbar, text="Clear Data", command=_clear_data)
action.grid(column=4, row=1)

# button to update plot
action = ttk.Button(toolbar, text="Update Plot", command=_plot)
action.grid(column=4, row=2)


#####################################################

# Configure frame and window background colors
color_back='#E4E4E4'
toolbar.configure(background=color_back)
window.configure(background=color_back)

# Configure grid elements for window resizing
window.grid_rowconfigure(1, weight=1)
window.grid_columnconfigure(0, weight=1)

#toolbar.grid(row=0, column=0, sticky="ew")
toolbar.grid(row=0, column=0)
dataPlot.get_tk_widget().grid(row=1, column=0, sticky="nsew")

window.mainloop()