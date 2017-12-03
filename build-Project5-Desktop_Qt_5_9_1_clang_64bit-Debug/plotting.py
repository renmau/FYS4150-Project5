import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress
from scipy.special import gamma

def plot_formatting(fam='serif',fam_font='Computer Modern Roman',font_size=14,tick_size=14):
	""" you get to define what font and size of xlabels and axis ticks you"""
	"""like, if you want bold text or not.								  """
	
	plt.rc('text',usetex=True)
	axis_font={'family': fam,'serif':[fam_font],'size':font_size}
	plt.rc('font',**axis_font)
	plt.rc('font',weight ='bold')
	plt.xticks(fontsize=tick_size)
	plt.yticks(fontsize=tick_size)


def savings():
	m = np.genfromtxt('a_lambda0_alpha0_gamma0_2e5trans_1e3sim.txt')
	m1 = np.genfromtxt('c_lambda025_alpha0_gamma0_1e6trans_1e3sim.txt')
	m2 = np.genfromtxt('c_lambda05_alpha0_gamma0_2e6trans_1e3sim.txt')
	m3 = np.genfromtxt('c_lambda09_alpha0_gamma0_1e7trans_1e3sim.txt')
	m_list=np.array([m,m1,m2,m3])

	lmbda = np.array([0,0.25,0.5,0.9])
	beta=np.array([1./np.mean(mm) for mm in m_list])
	# ------- reduced money array, m/<m>
	x = np.array([m_list[i]*beta[i] for i in range(len(beta))])
	
	#-------- fitted curves
	n = 1+ 3*lmbda/(1-lmbda)
	an = n**n/gamma(n)
	Pn_list = np.array([an[i]*x[i]**(n[i]-1)*np.exp(-n[i]*x[i]) for i in range(len(m_list))])

	#-------Gibbs distribution
	wmG=np.array([beta[i]*np.exp(-beta[i]*m_list[i]) for i in range(len(lmbda))])
	binsize=0.05
	weights = np.ones_like(x[0])/float(len(x[0])) # to normalize histogram

	#--------------------- a) and b)
	plot_formatting()
	plot_a=0
	if plot_a==1:
		counts,edges=np.histogram(x[0],np.arange(np.amin(x[0]),np.amax(x[0]+binsize),binsize),weights=weights)
		centers = (edges[:-1] + edges[1:])/2.
		plt.plot(centers,counts,label='numerical result')
		#plt.plot(mr,Pn_list[i]*binsize,label='best fit')
		plt.plot(x[0],wmG[0]*binsize/beta[0],label='Gibbs')
		plt.ylabel('$P(m)$')
		plt.xlabel(r'money $m/\langle m\rangle$')
		plt.legend()
		plt.tight_layout()
		plt.show()
		# ------------- log plot
		plt.plot(centers,counts, label='$\lambda=$'+str(lmbda[i]))
		plt.plot(x[0],wmG[0]*binsize/beta[0],label='Gibbs')
		plt.yscale('log')
		plt.ylabel('$P(m)$')
		plt.xlabel(r'money $m/\langle m\rangle$')
		plt.legend()
		plt.tight_layout()
		plt.show()

	#------------------------ c)
	plot_c=0
	if plot_c==1:
		ax = plt.gca()
		for i in range(len(beta)):
			mr = x[i]
			counts,edges=np.histogram(mr,np.arange(np.amin(mr),np.amax(mr+binsize),binsize),weights=weights)
			#plt.hist(mr,bins=np.arange(np.amin(mr),np.amax(mr+binsize),binsize), label='MC hist',weights=weights)
			centers = (edges[:-1] + edges[1:])/2.

			color=next(ax._get_lines.prop_cycler)['color']
			plot_label ='$\lambda=$'+str(lmbda[i])
			#plot_label ='numberical $\lambda=$'+str(lmbda[i])
			plt.plot(centers,counts,'o',markersize=3,color=color, label=plot_label)
			plt.plot(mr,Pn_list[i]*binsize,color=color)
			#plt.plot(mr,wmG[i]*binsize/beta[i],label='Gibbs $\lambda=$ '+str(lmbda[i]))
		
		#plt.yscale('log')
		plt.ylabel('$P(m)$')
		plt.xlabel(r'money $m/\langle m\rangle$')
		plt.legend()
		plt.tight_layout()
		plt.show()

	#------------------------ c) Tail end power law comparison

	#-------Power law distribution for tail ends
	# FIX THIS BIT, FIND PROPER POWER LAWS, PARETO MAYBE?
	wmP = np.array([5*x[i]**(-1-0.7) for i in range(len(lmbda))])
	plot_power=1
	if plot_power==1:
		ax = plt.gca()
		for i in range(2,3):
			mr = x[i]
			counts,edges=np.histogram(mr,np.arange(np.amin(mr),np.amax(mr+binsize),binsize),weights=weights)
			#plt.hist(mr,bins=np.arange(np.amin(mr),np.amax(mr+binsize),binsize), label='MC hist',weights=weights)
			centers = (edges[:-1] + edges[1:])/2.
			color=next(ax._get_lines.prop_cycler)['color']
			plt.plot(centers[10:],counts[10:],'o',markersize=3,color=color, label='$\lambda=$'+str(lmbda[i]))
			plt.plot(mr[260:],Pn_list[i,260:]*binsize,color=color)
			plt.plot(mr,wmP[i]*binsize-binsize)
		plt.yscale('log')
		plt.xscale('log')
		plt.ylabel('$P(m)$')
		plt.xlabel(r'money $m/\langle m\rangle$')
		plt.legend()
		plt.tight_layout()
		plt.show()


#a=linregress(centers,counts)[0] #something is wrong here
#a=linregress(counts,centers)[0]
#print 'slope =', a
def neighbors():
	m = np.genfromtxt('a_lambda0_alpha0_gamma0_2e5trans_1e3sim.txt')
	m1 = np.genfromtxt('d_lambda0_alpha05_gamma0_1e6trans_1e3sim.txt')
	m2 = np.genfromtxt('d_lambda0_alpha1_gamma0_2e6trans_1e3sim.txt')
	m3 = np.genfromtxt('d_lambda0_alpha15_gamma0_2e6trans_1e3sim.txt')
	m4 = np.genfromtxt('d_lambda0_alpha2_gamma0_2e6trans_1e3sim.txt')

	lmbda = np.array([0,0.25,0.5,0.9])
	alpha = np.array([0,0.5,1,1.5,2.])
	m_list=np.array([m,m1,m2,m3,m4])

	beta=np.array([1./np.mean(mm) for mm in m_list])
	x = np.array([m_list[i]*beta[i] for i in range(len(beta))])
	binsize=0.1
	weights = np.ones_like(x[0])/float(len(x[0]))

	a=1.	# find a from calculating slope of simulated result

	wmm=np.array([a*(m[i]*beta[i])**(-1-alpha[i]) for i in range(len(alpha))])


	#print 'max before norm =',np.amax(wmm)
	#print np.argmax(wmm)
	wmm /=np.amax(wmm)
	#print wmm
	#print 'max after norm =', np.amax(wmm)


	#----------------------d)
	plot_formatting()
	ax = plt.gca()
	for i in range(2,5):
		mr = x[i]
		counts,edges=np.histogram(mr,np.arange(np.amin(mr),np.amax(mr+binsize),binsize),weights=weights)
		#plt.hist(mr,bins=np.arange(np.amin(mr),np.amax(mr+binsize),binsize), label='MC hist',weights=weights)
		centers = (edges[:-1] + edges[1:])/2.
		color=next(ax._get_lines.prop_cycler)['color']
		plt.plot(centers,counts,':',color=color, label='$\lambda=$'+str(alpha[i]))
		#plt.plot(mr,Pn_list[i]*binsize,color=color)
	plt.yscale('log')
	plt.xscale('log')
	plt.ylabel('$P(m)$')
	plt.xlabel(r'money $m/\langle m\rangle$')
	plt.legend()
	plt.tight_layout()
	plt.show()

	#plt.hist(mr,bins=np.arange(np.amin(mr),np.amax(mr+binsize),binsize), label='MC hist',weights=weights)
	#plt.plot(mr,wm*5,label='Gibbs')

	#plt.plot(centers,counts,label='MC')
	#plt.loglog(m*beta,wmm)	
	#plt.loglog(centers,counts)
	#plt.ylabel('P(m)')
	#plt.xlabel('money $m$')
	#plt.legend()
	#plt.tight_layout()
	#plt.show()

#savings()
neighbors()