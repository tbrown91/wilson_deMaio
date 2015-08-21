import numpy as np
from bokeh.plotting import figure, show, output_file
from bokeh.models import Range1d
import os.path

x_vec = []
mean_vec = [];
sd_vec = [];
if (os.path.exists("out_0")):
	f_0 = open('out_0','r')
	time_0 = []
	x_vec.append(0)
	for line in f_0:
		time_0.append(float(line))
	f_0.close()
	mean_vec.append(np.mean(time_0));
	sd_vec.append(np.std(time_0));
if (os.path.exists("out_0.001")):
	f_0001 = open('out_0.001','r')
	time_0001 = []
	x_vec.append(0.001)
	for line in f_0001:
		time_0001.append(float(line))
	f_0001.close()
	mean_vec.append(np.mean(time_0001));
	sd_vec.append(np.std(time_0001));
if (os.path.exists("out_0.01")):
	f_001 = open('out_0.01','r')
	time_001 = []
	x_vec.append(0.01)
	for line in f_001:
		time_001.append(float(line))
	f_001.close()
	mean_vec.append(np.mean(time_001));
	sd_vec.append(np.std(time_001));
if (os.path.exists("out_0.1")):
	f_01 = open('out_0.1','r')
	time_01 = []
	x_vec.append(0.1)
	for line in f_01:
		time_01.append(float(line))
	f_01.close()
	mean_vec.append(np.mean(time_01));
	sd_vec.append(np.std(time_01));
if (os.path.exists("out_0.2")):
	f_02 = open('out_0.2','r')
	time_02 = []
	x_vec.append(0.2)
	for line in f_02:
		time_02.append(float(line))
	f_02.close()
	mean_vec.append(np.mean(time_02));
	sd_vec.append(np.std(time_02));
if (os.path.exists("out_0.5")):
	f_05 = open('out_0.5','r')
	time_05 = []
	x_vec.append(0.5)
	for line in f_05:
		time_05.append(float(line))
	f_05.close()
	mean_vec.append(np.mean(time_05));
	sd_vec.append(np.std(time_05));
if (os.path.exists("out_1")):
	f_1 = open('out_1','r')
	time_1 = []
	x_vec.append(1)
	for line in f_1:
		time_1.append(float(line))
	f_1.close()
	mean_vec.append(np.mean(time_1));
	sd_vec.append(np.std(time_1));
if (os.path.exists("out_1.5")):
	f_15 = open('out_1.5','r')
	time_15 = []
	x_vec.append(1.5)
	for line in f_15:
		time_15.append(float(line))
	f_15.close()
	mean_vec.append(np.mean(time_15));
	sd_vec.append(np.std(time_15));
if (os.path.exists("out_2")):
	f_2 = open('out_2','r')
	time_2 = []
	x_vec.append(2)
	for line in f_2:
		time_2.append(float(line))
	f_2.close()
	mean_vec.append(np.mean(time_2));
	sd_vec.append(np.std(time_2));

output_file('time_plot.html')
p = figure(title='Time taken for SimMLST to simulate ARG', width=800, height=800, title_text_font_size='8pt')
p.xaxis.axis_label = 'Recombination rate (%)'
p.yaxis.axis_label = 'Time to simulate ARG (seconds, error bar 1SD)'

p.circle(x_vec, mean_vec, color='blue')
x_err = []
mean_err = []
for x,y,yerr in zip(x_vec,mean_vec,sd_vec):
	x_err.append((x, x))
	mean_err.append((y - yerr, y + yerr))

p.multi_line(x_err, mean_err, color='blue')
show(p)
