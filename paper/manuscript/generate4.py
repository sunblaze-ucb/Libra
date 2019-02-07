import numpy as np
import matplotlib.pyplot as plt
from pylab import *

N = 16

dashes = [20,20]

time1 = (0.070766, 0.204950, 1.022308, 6.749271, 53.593483)
x=(16,32,64,128,256)

time2 = (0.4746755,2.519966,20.826396, 192.246441)
x2 = (16,32,64,128)

time3 = (192.246441, 1755.425)
x3 = (128, 256)

time8 = (0.3841,2.9149,24.1450843,196.3208,1545.2795)
x8 = (16,32,64,128, 256)

time10 = (50, 200, 800, 3200, 12000)
x10 = (16,32,64,128, 256)

time12 = (0.986328125,7.890625,63.125)
x12 = (16,32,64)

time13 = (63.125,505,4040)
x13 = (64,128, 256)

#time4 = (1.81724,2.426296,3.304848,5.082816,10.77412,20,32.162384)
#x4 = (2**10, 2**11, 2**12, 2**13, 2**14, 2**15, 50000)

#time5 = (32.162384, 40, 80, 160)
#x5 = (50000, 2**16, 2**17,2**18)

#time6 = (2.349896,2.845056,4.117528,8.995952,14.687584,32.821136)
#x6 = (2**10, 2**11, 2**12, 2**13, 2**14, 2**15)

#time7 = (32.821136, 64, 128, 256)
#x7 = (2**15, 2**16,2**17,2**18)

#time8 = (1.81724,2.426296,3.304848,5.082816,10.77412,20,32.162384)
#x8 = (2**10, 2**11, 2**12, 2**13, 2**14, 2**15, 50000)

#time9 = (32.162384, 40, 80, 160)
#x9 = (50000, 2**16, 2**17,2**18)

#time10 = (2.349896,2.845056,4.117528,8.995952,14.687584,32.821136)
#x10 = (2**10, 2**11, 2**12, 2**13, 2**14, 2**15)

#time11 = (32.821136, 64, 128, 256)
#x11 = (2**15, 2**16,2**17,2**18)

fig, ax = plt.subplots()

#plt.scatter(x1,Dynamic1, s=500, marker='o',facecolor='k')
#plt.scatter(x2,Dynamic2, s=500, marker='o',edgecolor='k',linewidth='3', facecolor='w', hatch='////')
rects1 = ax.plot(x, time1, color='b',linewidth=5,marker='o',markersize=30,fillstyle='full')
rects2 = ax.plot(x2, time2, color='r',linewidth=5,marker='^',markersize=30,fillstyle='full')
rects3 = ax.plot(x3, time3,color='r',linestyle='--',linewidth=5, marker='^',markersize=30)
#rects4 = ax.plot(x4, time4,color='y',linewidth=5,marker='s',markersize=30,fillstyle='full')
#rects5 = ax.plot(x5, time5,color='y',linestyle='-',linewidth=5)
#rects6 = ax.plot(x6, time6,color='k',linewidth=5,marker='d',markersize=30,fillstyle='full')
#rects7 = ax.plot(x7, time7,color='k',linestyle='-',linewidth=5)
rects8 = ax.plot(x8, time8,color='m',linewidth=5,marker='p',markersize=30,fillstyle='full')
#rects9 = ax.plot(x9, time9,color='m',linestyle='-',linewidth=5)
rects10 = ax.plot(x10, time10,color='g',linewidth=5,marker='v',markersize=30,fillstyle='full')
#rects11 = ax.plot(x11, time11,color='g',linestyle='-',linewidth=5)
#rects4 = ax.plot(x, time4,color='k',linewidth=5,marker='s',markersize=30,fillstyle='full')
#rects3 = ax.plot(x, planar,color='r',linewidth=3,marker='s',markersize=20,fillstyle='full')
rects12 = ax.plot(x12, time12, color='c',linewidth=5,marker='>',markersize=30,fillstyle='full')
rects13 = ax.plot(x13, time13, color='c',linestyle = '--', linewidth=5,marker='>',markersize=30,fillstyle='full')
#rects3[0].set_dashes(dashes)
#rects5[0].set_dashes(dashes)
#rects7[0].set_dashes(dashes)

ax.set_yscale('log')
ax.set_xscale('log')

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=26)
ax.xaxis.set_ticks((8,16,32,64,128,256, 512))
ax.set_xticklabels( (' ','$2^4$','$2^5$','$2^6$','$2^7$','$2^8$') )
ax.yaxis.set_ticks((0.01, 0.1, 1, 10, 100, 1000, 10000))
ax.set_yticklabels(('$10^{-2}$','$10^{-1}$','$10^0$','$10^1$','$10^2$','$10^3$', '$10^4$'),ha='left')
#plt.rcParams.update({'legend.labelspacing':0.25})

plt.xlim((14, 300))
plt.ylim((0.01,50000))
#plt.xlim((1000, 10000))
#plt.ylim((10,60000))


#leg = ax.legend( (rects1[0],rects2[0],rects3[0],rects4[0]), ('PP Linear 1', 'PP Linear 2', 'Linear C++', 'Linear Tensorflow') ,loc='upper left', borderpad=0.2,bbox_to_anchor=[0.001, 1.1],fontsize = 45)
#leg = ax.legend( (rects1[0],rects2[0],rects4[0],rects6[0]), ('Ours', 'vnTinyRAM', 'Buffet (ptr chase)', 'Buffet (KMP)'),loc='upper left', borderpad=0.1,bbox_to_anchor=[-0.03, 1.07],fontsize = 45)
ax.yaxis.grid(True)
ax.xaxis.grid(True)



#plt.xlim((5, 1000000))
ax.tick_params(axis='x', pad=15)
ax.tick_params(axis='y', pad=125)

#ax.text(150000, 6, '$2\\times10^5$', fontsize=50)


for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(70) 
for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(70)
ax.set_ylabel('prover time(s)', fontsize=60)
ax.set_xlabel('\#matrix column', fontsize=60)

plt.subplots_adjust(left=0.19, bottom=0.22, top=0.97, right=0.98)
figManager = plt.get_current_fig_manager()
#figManager.window.showMaximized()
fig.set_size_inches(14,10)
plt.savefig('fig4.pdf')

#plt.show()