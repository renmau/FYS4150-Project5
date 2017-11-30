import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress

def plot_formatting(fam='serif',fam_font='Computer Modern Roman',font_size=18,tick_size=18):
	""" you get to define what font and size of xlabels and axis ticks you"""
	"""like, if you want bold text or not.								  """
	
	plt.rc('text',usetex=True)
	axis_font={'family': fam,'serif':[fam_font],'size':font_size}
	plt.rc('font',**axis_font)
	plt.rc('font',weight ='bold')
	plt.xticks(fontsize=tick_size)
	plt.yticks(fontsize=tick_size)

plot_formatting()

#m = np.genfromtxt('d_alpha_05_lambda_0_1e5trans_1e3sim.txt')
#m = np.genfromtxt('a_1e5trans_1e3sim.txt')
#m = np.genfromtxt('a_gamma_0_alpha_0_lambda_0_1e5trans_1e3sim_.txt')
m = np.genfromtxt('e_gamma_1_alpha_2_lambda_0_1e5trans_1e3sim.txt')


#Gibbs distribution:
beta=1./np.mean(m)
alpha=2.
a=1.	# find a from calculating slope of simulated result
wm=beta*np.exp(-beta*m)
binsize=0.1
weights = np.ones_like(m)/float(len(m))
n = np.max(m)

counts,edges=np.histogram(m*beta,bins=np.arange(np.amin(m*beta),np.amax(m*beta+binsize),binsize),weights=weights)

centers = (edges[:-1] + edges[1:])/2.

#a=linregress(centers,counts)[0] #something is wrong here
#a=linregress(counts,centers)[0]
print 'slope =', a
wmm=a*(m*beta)**(-1-alpha)
plt.hist(m*beta,bins=np.arange(np.amin(m*beta),np.amax(m*beta+binsize),binsize), weights=weights)

print 'max before norm =',np.amax(wmm)
print np.argmax(wmm)
wmm /=np.amax(wmm)
print wmm
print 'max after norm =', np.amax(wmm)

#plt.plot(m*beta,wmm)
#plt.plot(centers,counts)
plt.loglog(m*beta,wmm)
plt.loglog(centers,counts)
plt.ylabel('P(m)')
plt.xlabel('money $m$')
plt.tight_layout()
plt.show()
