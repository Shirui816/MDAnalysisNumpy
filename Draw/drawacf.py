from pylab import *
from sys import argv

norm = Normalize(vmin = 0, vmax = int(argv[1])+1)

cmap = get_cmap('jet', int(argv[1])+1)

cls = cm.ScalarMappable(norm=norm, cmap=cmap)

sm = cm.ScalarMappable(cmap=cmap, norm=norm)

sm._A = []


for i in range(1,int(argv[1])+1):
    data = loadtxt("%s.dat" % (i))
    data = data.T
    t=data[0]
    acf=data[1]
    plot(t, acf, color=cls.to_rgba(i, alpha=0.75), ls='-', lw=2)


for tick in axes().xaxis.get_major_ticks():
    tick.label.set_fontsize(24)
    tick.label.set_fontsize(24)

for tick in axes().yaxis.get_major_ticks():
    tick.label.set_fontsize(24)
    tick.label.set_fontsize(24)
tick_params(pad = 15)

k=int(argv[1])+1
cb = colorbar(sm)
labels = arange(0,k,1)
loc    = labels + .5
cb.set_ticks(loc)
cb.set_ticklabels(labels)
cb.ax.yaxis.label.set_fontsize(24)



#show()
tight_layout()

savefig('acf.png', dpi=(330))

