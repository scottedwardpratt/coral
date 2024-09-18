import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
from pylab import *
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
sformatter=ScalarFormatter(useOffset=True,useMathText=True)
sformatter.set_scientific(True)
sformatter.set_powerlimits((-2,3))

#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}
plt.rc('font', **font)
plt.rc('text', usetex=False)
plt.figure(figsize=(6,5))
fig = plt.figure(1)
ax = fig.add_axes([0.15,0.13,0.8,0.8])
#y00=-100
#y01=40
#y02=-5
#a00=0.75
#a01=1.5
#a02=4.5

colors=['red','green','blue','cyan','violet']
   
mydata = np.loadtxt('../results/Rinv2.txt',skiprows=1,unpack=True)
q2=mydata[0]
cf2=mydata[1]
iplot=0
plt.plot(q2,cf2,linestyle='-',linewidth=2,color='r',marker=None,label='$R_{\\rm inv}=2$ fm')

mydata = np.loadtxt('../results/Rinv3.txt',skiprows=1,unpack=True)
q3=mydata[0]
cf3=mydata[1]
iplot=0
plt.plot(q3,cf3,linestyle='-',linewidth=2,color='g',marker=None,label='$R_{\\rm inv}=3$ fm')

mydata = np.loadtxt('../results/Rinv4.txt',skiprows=1,unpack=True)
q4=mydata[0]
cf4=mydata[1]
iplot=0
plt.plot(q4,cf4,linestyle='-',linewidth=2,color='b',marker=None,label='$R_{\\rm inv}=4$ fm') 

mydata = np.loadtxt('../results/Rinv5.txt',skiprows=1,unpack=True)
q5=mydata[0]
cf5=mydata[1]
iplot=0
plt.plot(q5,cf5,linestyle='-',linewidth=2,color='cyan',marker=None,label='$R_{\\rm inv}=5$ fm')

mydata = np.loadtxt('../results/Rinv5.txt',skiprows=1,unpack=True)
q6=mydata[0]
cf6=mydata[1]
iplot=0
plt.plot(q6,cf6,linestyle='-',linewidth=2,color='orange',marker=None,label='$R_{\\rm inv}=6$ fm')

#data_x=[1.0,2.0,3.0]
#data_y=[-4.5,-20.0,-27.5]
#plt.scatter(data_x,data_y,color='k',marker='*',s=140,zorder=10)

#data_x=[1.0,2.0,3.0]
#data_y=[-37.5,-52.5,-64.0]
#plt.scatter(data_x,data_y,color='k',marker='*',s=140,zorder=10)

#plt.plot(x,z,linestyle=linestyles[1],linewidth=2,color='k',markersize=8, marker=markerstyles[3], markerfacecolor='r', markeredgecolor=colors[3])

#plt.semilogy(x,y)

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,220,40), minor=False)
ax.set_xticklabels(np.arange(0,220,40), minor=False, family='serif')
ax.set_xticks(np.arange(0,220,10), minor=True)
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
plt.xlim(0.0,200)

ax.set_yticks(np.arange(0,2,0.5), minor=False)
ax.set_yticklabels(np.arange(0,2,0.5), minor=False, family='serif')
ax.set_yticks(np.arange(0,2,0.1), minor=True)
plt.ylim(0,1.15)
#ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
#ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
ax.yaxis.set_major_formatter(sformatter)

plt.xlabel('$q$ [MeV/$c$]', fontsize=18, weight='normal')
plt.ylabel('$C(q)$ ',fontsize=18,labelpad=0)

legend(loc="lower right")

#text(80,0.47,"pd femtoscopy",fontsize=24)

#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('cf_pd_vsR.pdf',format='pdf')
os.system('open -a Preview cf_pd_vsR.pdf')
#plt.show()
quit()
