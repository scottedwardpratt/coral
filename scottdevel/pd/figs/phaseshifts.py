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
plt.figure(figsize=(7,10))
fig = plt.figure(1)

MP=938.272
MN=939.565
MD=MP+MN-2.224575
mu=MP*MD/(MP+MD)
colors=['red','green','blue','cyan','violet','orange']

######## 
mydata = np.loadtxt('phaseshifts.txt',skiprows=1,unpack=True)
ke=mydata[1]
q=mydata[0]
delta_s12=mydata[2]
for i in range(0,40):
  if delta_s12[i] < -90.0:
    delta_s12[i]=delta_s12[i]+180.0
delta_s32=mydata[3]
delta_p12=mydata[4]
for i in range(0,40):
  if delta_p12[i] < -90.0:
    delta_p12[i]=delta_p12[i]+180.0
delta_p32=mydata[5]
for i in range(0,40):
  if delta_p32[i] < -90.0:
    delta_p32[i]=delta_p32[i]+180.0
delta_d12=mydata[6]
for i in range(0,40):
  if delta_d12[i] < -90.0:
    delta_d12[i]=delta_d12[i]+180.0
delta_d32=mydata[7]
for i in range(0,40):
  if delta_d32[i] < -90.0:
    delta_d32[i]=delta_d32[i]+180.0

#######

ax = fig.add_axes([0.14,0.06,0.83,0.46])

qq=arange(5.0,195.999,5.0)

ddeltadq_s12=np.zeros(39,dtype='float')
old=delta_s12[0]
for i in range(0,39,1):
  ddeltadq_s12[i]=(delta_s12[i+1]-old)/5.0
  old=delta_s12[i+1]
  
ddeltadq_s32=np.zeros(39,dtype='float')
old=delta_s32[0]
for i in range(0,39,1):
  ddeltadq_s32[i]=(delta_s32[i+1]-old)/5.0
  old=delta_s32[i+1]
  

ddeltadq_p12=np.zeros(39,dtype='float')
old=delta_p12[0]
for i in range(0,39,1):
  ddeltadq_p12[i]=(delta_p12[i+1]-old)/5.0
  old=delta_p12[i+1]
  
ddeltadq_p32=np.zeros(39,dtype='float')
old=delta_p32[0]
for i in range(0,39,1):
  ddeltadq_p32[i]=(delta_p32[i+1]-old)/5.0
  old=delta_p32[i+1]


ddeltadq_d12=np.zeros(39,dtype='float')
old=delta_d12[0]
for i in range(0,39,1):
  ddeltadq_d12[i]=(delta_d12[i+1]-old)/5.0
  old=delta_d12[i+1]
  
ddeltadq_d32=np.zeros(39,dtype='float')
old=delta_d32[0]
for i in range(0,39,1):
  ddeltadq_d32[i]=(delta_d32[i+1]-old)/5.0
  old=delta_d32[i+1]

plt.plot(qq,ddeltadq_s12,linestyle='-',linewidth=2,color=colors[0],marker=None,label='$L=0,S=1/2$')
plt.plot(qq,ddeltadq_s32,linestyle='-',linewidth=2,color=colors[1],marker=None,label='$L=0,S=3/2$')
plt.plot(qq,ddeltadq_p12,linestyle='-',linewidth=2,color=colors[2],marker=None,label='$L=1,S=1/2$')
plt.plot(qq,ddeltadq_p32,linestyle='-',linewidth=2,color=colors[3],marker=None,label='$L=1,S=3/2$')
plt.plot(qq,ddeltadq_d12,linestyle='-',linewidth=2,color=colors[4],marker=None,label='$L=2,S=1/2$')
plt.plot(qq,ddeltadq_d32,linestyle='-',linewidth=2,color=colors[5],marker=None,label='$L=2,S=3/2$')

ax.set_xticks(np.arange(0,300,50), minor=False)
ax.set_xticklabels(np.arange(0,300,50), minor=False, family='serif')
ax.set_xticks(np.arange(0,300,10), minor=True)
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
plt.xlim(0.0,200)
plt.ylim(-2,2)

plt.ylabel('$d\delta/dq$ [deg/MeV]', fontsize=18, weight='normal')
plt.xlabel('$q$ [MeV/$c$]',fontsize=18,labelpad=0)
legend(loc=[0.65,0.01])




#####

ax = fig.add_axes([0.14,0.52,0.83,0.46])


plt.plot(q,delta_s12,linestyle='-',linewidth=2,color=colors[0],marker=None)
plt.plot(q,delta_s32,linestyle='-',linewidth=2,color=colors[1],marker=None)
plt.plot(q,delta_p12,linestyle='-',linewidth=2,color=colors[2],marker=None)
plt.plot(q,delta_p32,linestyle='-',linewidth=2,color=colors[3],marker=None)
plt.plot(q,delta_d12,linestyle='-',linewidth=2,color=colors[4],marker=None)
plt.plot(q,delta_d32,linestyle='-',linewidth=2,color=colors[5],marker=None)

sdata=np.loadtxt('../phaseshifts/data/S/deltabar.txt',skiprows=1,unpack=True)

ke_s=sdata[0]
q_s=mu*sqrt(2.0*ke_s/MP);
del_s12=sdata[1]
del_s32=sdata[2]
plt.scatter(q_s,del_s12,color=colors[0],marker='o',label='$L=0,S=1/2$')
plt.scatter(q_s,del_s32,color=colors[1],marker='o',label='$L=0,S=3/2$')

pdata=np.loadtxt('../phaseshifts/data/P/deltabar.txt',skiprows=1,unpack=True)
ke_p=pdata[0]
q_p=mu*sqrt(2.0*ke_p/MP);
del_p12=pdata[1]
del_p32=pdata[2]
plt.scatter(q_p,del_p12,color=colors[2],marker='o',label='$L=1,S=1/2$')
plt.scatter(q_p,del_p32,color=colors[3],marker='o',label='$L=1,S=3/2$')

ddata=np.loadtxt('../phaseshifts/data/D/deltabar.txt',skiprows=1,unpack=True)
ke_d=ddata[0]
q_d=mu*sqrt(2.0*ke_d/MP);
del_d12=ddata[1]
del_d32=ddata[2]
plt.scatter(q_d,del_d12,color=colors[4],marker='o',label='$L=2,S=1/2$')
plt.scatter(q_d,del_d32,color=colors[5],marker='o',label='$L=2,S=3/2$')


ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,300,50), minor=False)
#ax.set_xticklabels(np.arange(0,300,50), minor=False, family='serif')
ax.set_xticklabels([])
ax.set_xticks(np.arange(0,300,10), minor=True)
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
plt.xlim(0.0,200)

ax.set_yticks(np.arange(-225,225,45), minor=False)
ax.set_yticklabels(np.arange(-225,225,45), minor=False, family='serif')
ax.set_yticks(np.arange(-225,225,15), minor=True)
plt.ylim(-165,50)
#ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
#ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
ax.yaxis.set_major_formatter(sformatter)

#plt.xlabel('$q$ [MeV/$c$]', fontsize=18, weight='normal')
plt.xlabel(None)
plt.ylabel('$\delta$ [degrees]',fontsize=18,labelpad=0)

legend(loc=[0.01,0.01])



#text(120,0.37,"$R_{\\rm inv}=3~{\\rm fm}$",fontsize=24)

#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('delta_pd.pdf',format='pdf')
os.system('open -a Preview delta_pd.pdf')
#plt.show()
quit()
