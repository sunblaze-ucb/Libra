import numpy as np
import matplotlib.pyplot as plt
from pylab import *

N = 16

dashes = [20,20]

time1 = (4.194444, 8.133423, 15.984551, 31.870675, 63.069480, 125.883671, 251.380656, 504.545425)
x=(2,4,8,16,32, 64, 128,256)


time2 = (15.664643,22.975494,35.133959, 66.515923, 124.087195, 249.319634, 503.562050, 1041.573771)
x2=(2,4,8,16,32, 64, 128,256)

time8 = (1.7935, 4.1527, 8.2123, 20.5963, 41.2982, 85.8453, 173.8858, 360.0756)
x8=(2,4,8,16,32, 64, 128,256)

time6 = (300, 800, 1800, 3500, 8000, 15000, 30000)
x6=(2,4,8,16,32, 64, 128)

time10 = (150, 400, 1000, 3000, 6500, 13000)
x10=(2,4,8,16,32, 64)

time4 = (2, 5, 9, 15, 25, 25, 80, 400)
x4=(2,4,8,16,32, 64, 128, 256)

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

#rects3[0].set_dashes(dashes)
#rects5[0].set_dashes(dashes)
#rects7[0].set_dashes(dashes)

ax.set_yscale('log')
ax.set_xscale('log')

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=26)
ax.xaxis.set_ticks((2, 4, 8,16,32,64,128,256, 512))
ax.set_xticklabels( ('$2$','$4$', '$8$','$16$','$32$','$64$','$128$', '$256$') )
ax.yaxis.set_ticks((1, 10, 100, 1000, 10000))
ax.set_yticklabels(('$10^0$','$10^1$','$10^2$','$10^3$', '$10^4$'),ha='left')
#plt.rcParams.update({'legend.labelspacing':0.25})

plt.xlim((1, 512))
plt.ylim((1,50000))
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
ax.set_ylabel('prove time, seconds', fontsize=50)
ax.set_xlabel('Number of leaves in Merkle tree', fontsize=50)

plt.subplots_adjust(left=0.17, bottom=0.22, top=0.97, right=0.95)
figManager = plt.get_current_fig_manager()
#figManager.window.showMaximized()
fig.set_size_inches(16,10)
plt.savefig('fig6.pdf')

#plt.show()