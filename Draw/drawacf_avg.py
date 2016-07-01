from pylab import *
from sys import argv

data = loadtxt("acf.log")
data = data.T
t=data[0]
acf=data[1]
plot(t, acf, alpha=.75, lw=2, ls='-')


for tick in axes().xaxis.get_major_ticks():
    tick.label.set_fontsize(24)
    tick.label.set_fontsize(24)

for tick in axes().yaxis.get_major_ticks():
    tick.label.set_fontsize(24)
    tick.label.set_fontsize(24)
tick_params(pad = 15)


#show()
tight_layout()

savefig('acf_avg.png', dpi=(330))

