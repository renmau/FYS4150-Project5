import matplotlib.pyplot as plt
import numpy as np

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
m = np.genfromtxt('e_gamma_1_alpha_1_lambda_0_1e5trans_1e3sim.txt')
print m

#Gibbs distribution:
beta=1./np.mean(m)
alpha=0.5
a=1.
wm=beta*np.exp(-beta*m)
wmm=a*m**(-1-alpha)

weights = np.ones_like(m)/float(len(m))

n = np.max(m)
#counts,edges=np.histogram(m,normed=True)
counts,edges=np.histogram(m,bins=np.arange(np.amin(m),np.amax(m)+n/10.,n/10.),normed=True)
centers = (edges[:-1] + edges[1:])/2.0


#plt.hist(m,bins=np.arange(np.amin(m),np.amax(m)+n/20.,n/20.),normed=True)
#plt.plot(m,wm)


#plt.plot(m,wmm)
plt.loglog(centers,counts)
#plt.semilogy(m,wm)
plt.show()
