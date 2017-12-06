import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress
from scipy.special import gamma

def plot_formatting(fam='serif',fam_font='Computer Modern Roman',font_size=18,tick_size=18):
	""" you get to define what font and size of xlabels and axis ticks you"""
	"""like, if you want bold text or not.								  """
	
	plt.rc('text',usetex=True)
	axis_font={'family': fam,'serif':[fam_font],'size':font_size}
	plt.rc('font',**axis_font)
	plt.rc('font',weight ='bold')
	plt.xticks(fontsize=tick_size)
	plt.yticks(fontsize=tick_size)


def savings():
	m = np.genfromtxt('a_lambda0_alpha0_gamma0_2e5trans.txt')
	m1 = np.genfromtxt('c_lambda025_alpha0_gamma0_1e6trans.txt')
	m2 = np.genfromtxt('c_lambda05_alpha0_gamma0_2e6trans.txt')
	m3 = np.genfromtxt('c_lambda09_alpha0_gamma0_1e7trans.txt')
	m_list=np.array([m,m1,m2,m3])

	lmbda = np.array([0,0.25,0.5,0.9])
	beta=np.array([1./np.mean(mm) for mm in m_list])
	#beta=np.array([np.mean(mm)/np.mean(mm) for mm in m_list])
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
		plt.plot(centers,counts, label='numerical result')
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
	#wmP = np.array([(x[i]**(-5.6)*10**(0.5)) for i in range(len(lmbda))])
	plot_power=1
	if plot_power==1:
		ax = plt.gca()
		for i in range(3,4):
			mr = x[i]
			counts,edges=np.histogram(mr,np.arange(np.amin(mr),np.amax(mr+binsize),binsize),weights=weights)
			#plt.hist(mr,bins=np.arange(np.amin(mr),np.amax(mr+binsize),binsize), label='MC hist',weights=weights)
			centers = (edges[:-1] + edges[1:])/2.
			color=next(ax._get_lines.prop_cycler)['color']
			plt.plot(centers[0:],counts[0:],'o',markersize=3,color=color, label='$\lambda=$'+str(lmbda[i]))
			#slope, intercept, r_value, p_value, std_err = linregress(np.log10(centers[25:]), np.log10(counts[25:]))
			#print slope
			#plt.plot(mr[260:],Pn_list[i,260:]*binsize,color=color)
			wmP = 10**(1.1)*mr**(-12)
			#plt.plot(mr,wmP[i]*binsize-binsize,label='Power law $\lambda=$'+str(lmbda[i]))
			plt.plot(mr,wmP*binsize,label='Power law fit')
		plt.yscale('log')
		plt.xscale('log')
		plt.ylabel('$\log_{10}P(m)$')
		plt.xlabel(r'money $\log_{10}m/\langle m\rangle$')
		plt.legend()
		plt.tight_layout()
		plt.show()

def neighbors_alpha():
	#m = np.genfromtxt('a_lambda0_alpha0_gamma0_2e5trans_1e3sim.txt')

	# the _ at the end of the file name means avg_m=1.
	m1 = np.genfromtxt('d_lambda0_alpha05_gamma0_2e6trans_.txt')
	m2 = np.genfromtxt('d_lambda0_alpha1_gamma0_2e6trans_.txt')
	m3 = np.genfromtxt('d_lambda0_alpha15_gamma0_2e6trans_.txt')
	m4 = np.genfromtxt('d_lambda0_alpha2_gamma0_2e6trans_.txt')

	c1 = np.genfromtxt('d_lambda0_alpha05_gamma0_2e6trans_1000_.txt')
	c2 = np.genfromtxt('d_lambda0_alpha1_gamma0_2e6trans_1000_.txt')
	c3 = np.genfromtxt('d_lambda0_alpha15_gamma0_2e6trans_1000_.txt')
	c4 = np.genfromtxt('d_lambda0_alpha2_gamma0_2e6trans_1000_.txt')

	# ----------- no of agents
	N1 = 500
	N2 = 1000
	lmbda = np.array([0,0.25,0.5,0.9])
	alpha = np.array([0.5,1,1.5,2.])
	m_list = np.array([m1,m2,m3,m4]) 	# 500 agents
	c_list = np.array([c1,c2,c3,c4])	# 1000 agents

	beta=np.array([1./np.mean(mm) for mm in m_list])

	# to reproduce plots in article, x = m_list, binsize ~1-25
	x1 = np.array([m_list[i]*beta[i] for i in range(len(beta))])
	x2 = np.array([c_list[i]*beta[i] for i in range(len(beta))])
	binsize1 = 20.
	binsize2 = 20.
	weights1 = np.ones_like(x1[0])/float(len(x1[0]))
	weights2 = np.ones_like(x2[0])/float(len(x2[0]))

	# This need fixing
	a=1.	# find a from calculating slope of simulated result

	wmm1=np.array([a*(m_list[i]*beta[i])**(-1-alpha[i]) for i in range(len(alpha))])
	wmm2=np.array([a*(c_list[i]*beta[i])**(-1-alpha[i]) for i in range(len(alpha))])


	#print 'max before norm =',np.amax(wmm)
	#print np.argmax(wmm)
	wmm1 /=np.amax(wmm1)
	wmm2 /=np.amax(wmm2)
	#print wmm
	#print 'max after norm =', np.amax(wmm)


	#----------------------d)
	plot_formatting()
	ax = plt.gca()
	for i in range(len(m_list)):
		mr = m_list[i]#x1[i]
		counts,edges=np.histogram(mr,np.arange(np.amin(mr),np.amax(mr+binsize1),binsize1),weights=weights1)
		#plt.hist(mr,bins=np.arange(np.amin(mr),np.amax(mr+binsize),binsize), label='MC hist',weights=weights1)
		centers = (edges[:-1] + edges[1:])/2.
		color=next(ax._get_lines.prop_cycler)['color']
		plt.plot(centers,counts,color=color, label=r'$\alpha=$'+str(alpha[i]))
		#plt.plot(mr,Pn_list[i]*binsize,color=color)
	plt.annotate('N=$ '+np.str(N1), xy=(0.65, 0.05), xycoords='axes fraction')
	plt.yscale('log')
	plt.xscale('log')
	plt.ylabel('$P(m)$')
	plt.xlabel(r'money $m/\langle m\rangle$')
	plt.legend()
	plt.grid(b=True, which='minor', alpha=0.2)
	plt.tight_layout()
	plt.show()

	for i in range(len(c_list)):
		mr = c_list[i]#x2[i]
		counts,edges=np.histogram(mr,np.arange(np.amin(mr),np.amax(mr+binsize1),binsize1),weights=weights2)
		#plt.hist(mr,bins=np.arange(np.amin(mr),np.amax(mr+binsize),binsize), label='MC hist',weights=weights2)
		centers = (edges[:-1] + edges[1:])/2.
		color=next(ax._get_lines.prop_cycler)['color']
		plt.plot(centers,counts,color=color, label=r'$\alpha=$'+str(alpha[i]))
		#plt.plot(mr,Pn_list[i]*binsize,color=color)
	plt.yscale('log')
	plt.xscale('log')
	plt.ylabel('$P(m)$')
	plt.xlabel(r'money $m/\langle m\rangle$')
	######## get a text box with N stated
	plt.grid(b=True, which='minor', alpha=0.2)
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

# This is not fixed in anyway, just copypasted the general structure.
def neighbors_gamma():
	# the _ at the end of the file name means avg_m=1.
	m = np.genfromtxt('e_lambda0_alpha0_gamma0_2e5trans_.txt')
	m1 = np.genfromtxt('e_lambda0_alpha05_gamma0_2e6trans_.txt')
	m2 = np.genfromtxt('e_lambda0_alpha1_gamma0_2e6trans_.txt')
	m3 = np.genfromtxt('e_lambda0_alpha15_gamma0_2e6trans_.txt')
	m4 = np.genfromtxt('e_lambda0_alpha2_gamma0_2e6trans_.txt')
	m4 = np.genfromtxt('e_lambda0_alpha20_gamma0_2e6trans_.txt') # for alpha>>1

	lmbda = np.array([0,0.25,0.5,0.9])
	alpha = np.array([0.5,1,1.5,2.,20])
	gamma = np.array([0,1.,2.,3.,4.])
	m_list=np.array([m,m1,m2,m3,m4])

	beta=np.array([1./np.mean(mm) for mm in m_list])
	# to reproduce plots in article, x = m_list, binsize ~1-10
	x = np.array([m_list[i]*beta[i] for i in range(len(beta))])
	binsize=0.1
	weights = np.ones_like(x[0])/float(len(x[0]))

	''' THIS NEED FIXING
	a=1.	# find a from calculating slope of simulated result

	wmm=np.array([a*(m[i]*beta[i])**(-1-alpha[i]) for i in range(len(alpha))])


	#print 'max before norm =',np.amax(wmm)
	#print np.argmax(wmm)
	wmm /=np.amax(wmm)
	#print wmm
	#print 'max after norm =', np.amax(wmm)
	'''

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

savings()
#neighbors_alpha()
#neighbors_gamma()
