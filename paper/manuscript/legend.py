import pylab
import matplotlib.pyplot as plt

fig = pylab.figure()
figlegend = pylab.figure(figsize=(40,1.4))
ax = fig.add_subplot(111)
rects1 = ax.plot(range(12), color='b',linewidth=5,marker='o',markersize=15,fillstyle='full')
rects2 = ax.plot(range(12), color='r',linewidth=5,marker='^',markersize=15,fillstyle='full')
rects3 = ax.plot(range(12), color='g',linewidth=5,marker='v',markersize=15,fillstyle='full')
rects4 = ax.plot(range(12), color='y',linewidth=5,marker='s',markersize=15,fillstyle='full')
rects5 = ax.plot(range(12), color='m',linewidth=5,marker='p',markersize=15,fillstyle='full')
rects6 = ax.plot(range(12), color='k',linewidth=5,marker='d',markersize=15,fillstyle='full')
rects7 = ax.plot(range(12), color='c',linewidth=5,marker='>',markersize=15,fillstyle='full')

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=30)

#leg = ax.legend( (rects1[0],rects2[0],rects4[0],rects6[0]), ('Ours', 'vnTinyRAM', 'Buffet (ptr chase)', 'Buffet (KMP)'),loc='upper left', borderpad=0.1,bbox_to_anchor=[-0.03, 1.07],fontsize = 45)
#lines = ax.plot(range(10), pylab.randn(10), range(10), pylab.randn(10))
figlegend.legend((rects1[0], rects2[0], rects3[0], rects4[0], rects5[0], rects6[0], rects7[0]), ('Ours', 'Hyrax', 'Bulletproofs', 'Ligero', 'libSNARK', 'libSTARK', 'Aurora'), 'center', ncol = 4, fontsize = 25)
fig.show()
figlegend.show()
figlegend.savefig('legend.pdf') 