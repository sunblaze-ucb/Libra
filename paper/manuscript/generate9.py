import numpy as np
import matplotlib.pyplot as plt
from pylab import *

N = 16

dashes = [20,20]

time1 = (0.4777413333,0.5001243333,0.5244003333,0.5470016667,0.5721683333,0.6033633333,0.644886,0.7066146667)
x=(2,4,8,16,32, 64, 128,256)


time2 = (0.901627, 1.010725,1.080326, 1.525268, 2.000217, 3.559008, 4.936292, 9.872230)
x2=(2,4,8,16,32, 64, 128,256)

time6 = (0.031, 0.032, 0.033, 0.034, 0.036, 0.038, 0.043, 0.044)
x6=(2,4,8,16,32, 64, 128,256)

time8 = (0.0036, 0.0036, 0.0036, 0.0036, 0.0036, 0.0036, 0.0036, 0.0036)
x8=(2,4,8,16,32, 64, 128,256)

time10 = (20, 40, 100, 200, 400, 900)
x10=(2,4,8,16,32, 64)

time4 = (0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 2.5, 4)
x4=(2,4,8,16,32, 64, 128,256)

time12 = (0.08926391604,0.2082824708,0.4463195802,0.9223937991,1.874542237)
x12=(2,4,8,16,32)

time13 = (1.874542237,3.778839112,7.58743286,15.20462037)
x13=(32, 64, 128, 256)
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
#rects3 = ax.plot(x3, time3,color='r',linestyle='--',linewidth=5, marker='^',markersize=30)
rects4 = ax.plot(x4, time4,color='y',linewidth=5,marker='s',markersize=30,fillstyle='full')
#rects5 = ax.plot(x5, time5,color='y',linestyle='-',linewidth=5)
rects6 = ax.plot(x6, time6,color='k',linewidth=5,marker='d',markersize=30,fillstyle='full')
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
ax.xaxis.set_ticks((2, 4, 8,16,32,64,128,256, 512))
ax.set_xticklabels( ('$2^1$','$2^2$', '$2^3$','$2^4$','$2^5$','$2^6$','$2^7$','$2^8$') )
ax.yaxis.set_ticks((0.001, 0.01, 0.1, 1, 10, 100, 1000))
ax.set_yticklabels(('$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^0$','$10^1$','$10^2$', '$10^3$'),ha='left')
#plt.rcParams.update({'legend.labelspacing':0.25})

plt.xlim((1.7, 300))
plt.ylim((0.001,5000))
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
ax.set_ylabel('verification time(s)', fontsize=60)
ax.set_xlabel('\#leave', fontsize=60)

plt.subplots_adjust(left=0.19, bottom=0.22, top=0.97, right=0.98)
figManager = plt.get_current_fig_manager()
#figManager.window.showMaximized()
fig.set_size_inches(14,10)
plt.savefig('C:/Users/star/Desktop/papers/fastZKP/paper/manuscript/fig9.pdf')

plt.show()