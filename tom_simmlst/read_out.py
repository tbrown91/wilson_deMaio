import numpy as np
from bokeh.plotting import figure, show, output_file
from bokeh.models import Range1d

f_0 = open('out_0','r')
material_0 = []
for line in f_0:
	material_0.append(float(line))
f_0.close()
f_0001 = open('out_0.001','r')
material_0001 = []
for line in f_0001:
	material_0001.append(float(line))
f_0001.close()
f_001 = open('out_0.01','r')
material_001 = []
for line in f_001:
	material_001.append(float(line))
f_001.close()
f_01 = open('out_0.1','r')
material_01 = []
for line in f_01:
	material_01.append(float(line))
f_01.close()
f_02 = open('out_0.2','r')
material_02 = []
for line in f_02:
	material_02.append(float(line))
f_02.close()
f_05 = open('out_0.5','r')
material_05 = []
for line in f_05:
	material_05.append(float(line))
f_05.close()
f_1 = open('out_1','r')
material_1 = []
for line in f_1:
	material_1.append(float(line))
f_1.close()
f_15 = open('out_1.5','r')
material_15 = []
for line in f_15:
	material_15.append(float(line))
f_15.close()
f_2 = open('out_2','r')
material_2 = []
for line in f_2:
	material_2.append(float(line))
f_2.close()

mean_vec = [];
mean_vec.append(np.mean(material_0));
mean_vec.append(np.mean(material_0001));
mean_vec.append(np.mean(material_001));
mean_vec.append(np.mean(material_01));
mean_vec.append(np.mean(material_02));
mean_vec.append(np.mean(material_05));
mean_vec.append(np.mean(material_1));
mean_vec.append(np.mean(material_15));
mean_vec.append(np.mean(material_2));
sd_vec = [];
sd_vec.append(np.std(material_0));
sd_vec.append(np.std(material_0001));
sd_vec.append(np.std(material_001));
sd_vec.append(np.std(material_01));
sd_vec.append(np.std(material_02));
sd_vec.append(np.std(material_05));
sd_vec.append(np.std(material_1));
sd_vec.append(np.std(material_15));
sd_vec.append(np.std(material_2));
x_vec = [0,0.001,0.01,0.1,0.2,0.5,1,1.5,2];

output_file('MRCA_materialPlot.html')
p = figure(title='Average amount of ancestral material remaining in Clonal Frame (100 simulations)', width=800, height=800, title_text_font_size='16pt')
p.xaxis.axis_label = 'Recombination rate (%)'
p.yaxis.axis_label = 'Fraction of ancestral material in MRCA of Clonal Frame (error bar 1SD)'

p.circle(x_vec, mean_vec, color='blue')
x_err = []
mean_err = []
for x,y,yerr in zip(x_vec,mean_vec,sd_vec):
	x_err.append((x, x))
	mean_err.append((y - yerr, y + yerr))

p.multi_line(x_err, mean_err, color='blue')
p.set(x_range=Range1d(-0.5, 2.5), y_range=Range1d(-0.1, 1.1))
show(p)
