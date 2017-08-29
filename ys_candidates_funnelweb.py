# TODO: make code more clear

import numpy as np
import scipy.stats #import binned_statistic_2d
import math
import numpy.ma as ma
import pickle
import matplotlib.pyplot as plt
import copy
#~ import matplotlib.colormap as cp

import young_star_candidates as ys

# Convergence point
# This is solar apex; replace it with young stars
#~ from astropy.coordinates import SkyCoord
#~ from astropy import units as u
#~ c = SkyCoord(l=50.2*u.degree, b=24.7*u.degree, frame='galactic')
#~ RA_CP=c.icrs.ra.value
#~ DEC_CP=c.icrs.dec.value

# Convergence point coordinates
# This is Solar apex:
RA_CP=270.0
DEC_CP=30.0

# Read the data
#~ pkl_file = open('gaia_2mass_for_funnelweb.pkl', 'rb')
pkl_file = open('gaia_2mass_for_funnelweb_with_UCAC5.pkl', 'rb')
d = pickle.load(pkl_file)
pkl_file.close()

# Local libraries
#~ import pehta
#~ query='select ra, de, pmra, pmdec, pmRA_u, pmDE_u, w1mpro, w2mpro, phot_g_mean_mag, parallax, ewirt from master where ra is not null and de is not null and w1mpro is not null and de<100'
#query='select ra, de, pmra, pmdec, pmRA_u, pmDE_u, k_m, phot_g_mean_mag, ewirt from master where ra is not null and de is not null and k_m is not null and de<100'
#~ d=pehta.read_mysql(query=query, DB='young_stars')

result=[]
for x in d:
	# pmra = mu_alpha* (meaning that cos(delta) has already been taken into account)
	
	# Determine angle theta between the two coordinate systems
	theta=ys.determine_theta(alpha_cp=RA_CP, delta_cp=DEC_CP, alpha=x['ra'], delta=x['de'])
	
	# Rotate vector for angle theta (rotation is done for -theta)
	#~ rotated_vector=ys.rotate_vector(vec=np.array([x['pmra'], x['pmdec']]), theta=theta)

	try:
		if np.isfinite(x['pmRA_u']):
			rotated_vector=ys.rotate_vector(vec=np.array([x['pmRA_u'], x['pmDE_u']]), theta=theta)
	except:
		try:
			if np.isfinite(x['pmra']):
				rotated_vector=ys.rotate_vector(vec=np.array([x['pmra'], x['pmdec']]), theta=theta)
		except:
			continue

	pm_parallel=rotated_vector[0]
	pm_perp=rotated_vector[1]
	
	# Is this ok?
	pm_parallel /= np.sin(np.deg2rad(90.0-theta))
	
	
	# Reduced proper motions
	# 1e-3: convert from mas/yr to arcsec/yr
	#~ gmag=x['phot_g_mean_mag']
	#~ gmag=x['g_psf']
	
	gmag=x['phot_g_mean_mag']
	m1=gmag
	m2=x['w1mpro']
	#~ m2=x['k_m']
	
	
	#~ H_parallel = gmag+5.0+5.0*np.log10(pm_parallel*1e-3)
	#~ H_perp = gmag+5.0+5.0*np.log10(pm_perp*1e-3)
	
	# Take absolute values of proper motion
	H_parallel = gmag+5.0+5.0*np.log10(np.abs(pm_parallel)*1e-3)
	H_perp = gmag+5.0+5.0*np.log10(np.abs(pm_perp)*1e-3)
	H = gmag+5.0+5.0*np.log10(np.sqrt(pm_parallel**2+pm_perp**2)*1e-3)
	
	try:
		if x['ewirt']>-10:
			ewirt=x['ewirt']
	except:
		ewirt=-10

	try:
		G=gmag+5.0-5.0*np.log10(1.0/x['parallax'])
	except:
		G=999
	
	#~ result.append([pm_parallel, pm_perp, x['ra'], x['de'], H_parallel, H_perp, gmag-x['k_m']])
	result.append([pm_parallel, pm_perp, x['ra'], x['de'], H_parallel, H_perp, m1-m2, H, G, ewirt])

result=np.array(result)
print('Number of stars:', len(result))
C={'pm_parallel': 0, 'pm_perp': 1, 'ra': 2, 'de': 3, 'H_parallel': 4, 'H_perp': 5, 'mag': 6, 'H': 7, 'G': 8, 'EWirt': 9}

def known_moving_groups():
	# From Malo et al. 20..
	TODO=True

def reduced_proper_motion_diagrams(result, C):

	def shift_edges(e=False): # to plot histograms
		d=e[1]-e[0]
		e2=[x+d/2.0 for x in e]
		return e2[:-1]

	xbins=100
	ybins=100
	#~ xlim=[-2, 8]
	xlim=[-2, 6]
	ylim=[-20, 30]
	x_label='gmag-W1'
	
	colorbar_index='pm_parallel'
	
	cbi={'pm_parallel': 0, 'pm_perp': 1}
	amp={'pm_parallel': 400, 'pm_perp': 150}

	vmin=-amp[colorbar_index]
	vmax=amp[colorbar_index]
	pointsize=1

	result=np.array(sorted(result, key=lambda x: x[cbi[colorbar_index]]))
	
	# Subplots
	nx=4
	ny=2

	fig=plt.figure()


	# Diagram of reduced proper motion (proper motion size including both components)
	ax0=fig.add_subplot(ny, nx, 1)

	hist, xedges, yedges, binnumber = scipy.stats.binned_statistic_2d(result[:,7], result[:,6], result[:,6], statistic='count', bins=[ybins, xbins], range=[ylim, xlim])
	hist=ma.array(hist, mask=hist<0.000001)
	hist=np.log10(hist)
	#~ shift edges
	xedges=shift_edges(xedges)
	yedges=shift_edges(yedges)
	cs = ax0.contour(yedges, xedges, hist, colors='k')
	ax0.set_xlabel(x_label)
	ax0.set_ylabel('H')
	ax0.set_xlim(xlim[0], xlim[1])
	ax0.set_ylim(ylim[1], ylim[0])

	#RAVE
	result2=copy.deepcopy(result)
	result2=np.array(sorted(result2, key=lambda x: x[-1]))
	
	cb=ax0.scatter(result2[:,6], result2[:,7], c=result2[:,-1], vmin=0, vmax=1.5, s=10, edgecolor='none', cmap=plt.cm.bone_r)
	c=plt.colorbar(cb)
	c.set_label('EWirt [A]')
	

	# Contours H_parallel
	hist, xedges, yedges, binnumber = scipy.stats.binned_statistic_2d(result[:,4], result[:,6], result[:,6], statistic='count', bins=[ybins, xbins], range=[ylim, xlim])
	hist=ma.array(hist, mask=hist<0.000001)
	histPA=np.log10(hist)
	#~ X, Y = np.meshgrid(yedges, xedges)
	#~ shift edges
	xedgesPA=shift_edges(xedges)
	yedgesPA=shift_edges(yedges)
	
	
	# Contours H_perp
	hist, xedges, yedges, binnumber = scipy.stats.binned_statistic_2d(result[:,5], result[:,6], result[:,6], statistic='count', bins=[ybins, xbins], range=[ylim, xlim])
	hist=ma.array(hist, mask=hist<0.000001)
	histPE=np.log10(hist)
	#~ X, Y = np.meshgrid(yedges, xedges)
	#~ shift edges
	xedgesPE=shift_edges(xedges)
	yedgesPE=shift_edges(yedges)

	# H_parallel vs G-K
	ax=fig.add_subplot(ny, nx, 2)
	#~ cb=ax.scatter(result[:,6], result[:,4], c=result[:,2], vmin=0, vmax=200, s=5, edgecolor='none', cmap=plt.cm.RdYlGn)
	cb=ax.scatter(result[:,6], result[:,4], c=result[:,cbi[colorbar_index]], vmin=vmin, vmax=vmax, s=pointsize, edgecolor='none', cmap=plt.cm.RdYlGn)
	c=plt.colorbar(cb)
	c.set_label(colorbar_index)
	ax.set_xlabel(x_label)
	ax.set_ylabel('H parallel')
	ax.set_xlim(xlim[0], xlim[1])
	ax.set_ylim(ylim[1], ylim[0])
	cs = ax.contour(yedgesPA, xedgesPA, histPA, colors='k')

	# H_perp vs G-K
	ax=fig.add_subplot(ny, nx, 6)
	cb=ax.scatter(result[:,6], result[:,5], c=result[:,cbi[colorbar_index]], vmin=vmin, vmax=vmax, s=pointsize, edgecolor='none', cmap=plt.cm.RdYlGn)
	c=plt.colorbar(cb)
	c.set_label(colorbar_index)
	ax.set_xlabel(x_label)
	ax.set_ylabel('H perpendicular')
	ax.set_xlim(xlim[0], xlim[1])
	ax.set_ylim(ylim[1], ylim[0])
	cs = ax.contour(yedgesPE, xedgesPE, histPE, colors='k')


	# Reverse the point order
	result1=copy.deepcopy(result)
	result1=np.array(sorted(result1, key=lambda x: x[cbi[colorbar_index]], reverse=True))
	#~ result=np.array(sorted(result, key=lambda x: x[0], reverse=True))
	#~ result=result1

	# H_parallel vs G-K
	ax=fig.add_subplot(ny, nx, 3)
	#~ cb=ax.scatter(result[:,6], result[:,4], c=result[:,2], vmin=0, vmax=200, s=5, edgecolor='none', cmap=plt.cm.RdYlGn)
	cb=ax.scatter(result1[:,6], result1[:,4], c=result1[:,cbi[colorbar_index]], vmin=vmin, vmax=vmax, s=pointsize, edgecolor='none', cmap=plt.cm.RdYlGn)
	c=plt.colorbar(cb)
	c.set_label(colorbar_index)
	ax.set_xlabel(x_label)
	ax.set_ylabel('H parallel')
	ax.set_xlim(xlim[0], xlim[1])
	ax.set_ylim(ylim[1], ylim[0])
	cs = ax.contour(yedgesPA, xedgesPA, histPA, colors='k')

	# H_perp vs G-K
	ax=fig.add_subplot(ny, nx, 7)
	cb=ax.scatter(result1[:,6], result1[:,5], c=result1[:,cbi[colorbar_index]], vmin=vmin, vmax=vmax, s=pointsize, edgecolor='none', cmap=plt.cm.RdYlGn)
	c=plt.colorbar(cb)
	c.set_label(colorbar_index)
	ax.set_xlabel(x_label)
	ax.set_ylabel('H perpendicular')
	ax.set_xlim(xlim[0], xlim[1])
	ax.set_ylim(ylim[1], ylim[0])
	cs = ax.contour(yedgesPE, xedgesPE, histPE, colors='k')


	# RAVE
	# H_parallel vs G-K
	ax=fig.add_subplot(ny, nx, 4)
	#~ cb=ax.scatter(result[:,6], result[:,4], c=result[:,2], vmin=0, vmax=200, s=5, edgecolor='none', cmap=plt.cm.RdYlGn)
	cb=ax.scatter(result2[:,6], result2[:,4], c=result2[:,-1], vmin=0, vmax=1.5, s=10, edgecolor='none', cmap=plt.cm.bone_r)
	c=plt.colorbar(cb)
	c.set_label('EWirt [A]')
	ax.set_xlabel(x_label)
	ax.set_ylabel('H parallel')
	ax.set_xlim(xlim[0], xlim[1])
	ax.set_ylim(ylim[1], ylim[0])
	cs = ax.contour(yedgesPA, xedgesPA, histPA, colors='k')

	# H_perp vs G-K
	ax=fig.add_subplot(ny, nx, 8)
	cb=ax.scatter(result2[:,6], result2[:,5], c=result2[:,-1], vmin=0, vmax=1.5, s=10, edgecolor='none', cmap=plt.cm.bone_r)
	c=plt.colorbar(cb)
	c.set_label('EWirt [A]')
	ax.set_xlabel(x_label)
	ax.set_ylabel('H perpendicular')
	ax.set_xlim(xlim[0], xlim[1])
	ax.set_ylim(ylim[1], ylim[0])
	cs = ax.contour(yedgesPE, xedgesPE, histPE, colors='k')






# Color-magnitude diagram
	fig=plt.figure()
	ax=fig.add_subplot(221)

	hist, xedges, yedges, binnumber = scipy.stats.binned_statistic_2d(result[:,C['G']], result[:,C['mag']], result[:,6], statistic='count', bins=[ybins, xbins], range=[[2, 25], xlim])
	hist=ma.array(hist, mask=hist<0.000001)
	hist=np.log10(hist)
	#~ shift edges
	xedges=shift_edges(xedges)
	yedges=shift_edges(yedges)
	cs = ax.contour(yedges, xedges, hist, colors='k')

	keyword='pm_parallel'

	result4=copy.deepcopy(result)
	result4=np.array(sorted(result4, key=lambda x: x[C[keyword]]))
	cb=ax.scatter(result4[:,C['mag']], result4[:,C['G']], c=result4[:,C[keyword]], vmin=-amp[keyword], vmax=amp[keyword], s=1, edgecolor='none', cmap=plt.cm.RdYlGn)
	c=plt.colorbar(cb)
	c.set_label(keyword)
	ax.set_xlabel(x_label)
	ax.set_ylabel('G')
	ax.set_xlim(xlim[0], xlim[1])
	ax.set_ylim(25, 5)


	ax=fig.add_subplot(222)
	cs = ax.contour(yedges, xedges, hist, colors='k')

	result5=copy.deepcopy(result)
	result5=np.array(sorted(result5, key=lambda x: x[C[keyword]], reverse=True))
	cb=ax.scatter(result5[:,C['mag']], result5[:,C['G']], c=result5[:,C[keyword]], vmin=-amp[keyword], vmax=amp[keyword], s=1, edgecolor='none', cmap=plt.cm.RdYlGn)
	c=plt.colorbar(cb)
	c.set_label(keyword)
	ax.set_xlabel(x_label)
	ax.set_ylabel('G')
	ax.set_xlim(xlim[0], xlim[1])
	ax.set_ylim(25, 5)


	ax=fig.add_subplot(223)
	cs = ax.contour(yedges, xedges, hist, colors='k')
	
	keyword='pm_perp'

	result6=copy.deepcopy(result)
	result6=np.array(sorted(result6, key=lambda x: x[C[keyword]]))
	cb=ax.scatter(result6[:,C['mag']], result6[:,C['G']], c=result6[:,C[keyword]], vmin=-amp[keyword], vmax=amp[keyword], s=1, edgecolor='none', cmap=plt.cm.RdYlGn)
	c=plt.colorbar(cb)
	c.set_label(keyword)
	ax.set_xlabel(x_label)
	ax.set_ylabel('G')
	ax.set_xlim(xlim[0], xlim[1])
	ax.set_ylim(25, 5)


	ax=fig.add_subplot(224)
	cs = ax.contour(yedges, xedges, hist, colors='k')

	result7=copy.deepcopy(result)
	result7=np.array(sorted(result7, key=lambda x: x[C[keyword]], reverse=True))
	cb=ax.scatter(result7[:,C['mag']], result7[:,C['G']], c=result7[:,C[keyword]], vmin=-amp[keyword], vmax=amp[keyword], s=1, edgecolor='none', cmap=plt.cm.RdYlGn)
	c=plt.colorbar(cb)
	c.set_label(keyword)
	ax.set_xlabel(x_label)
	ax.set_ylabel('G')
	ax.set_xlim(xlim[0], xlim[1])
	ax.set_ylim(25, 5)

	plt.show()


def other_plots(): # TODO



	# H_perp vs H_parallel
	fig=plt.figure()
	ax=fig.add_subplot(111)
	result=np.array(sorted(result, key=lambda x: x[6]))
	cb=ax.scatter(result[:,4], result[:,5], c=result[:,6], s=5, edgecolor='none', cmap=plt.cm.RdYlGn, vmin=0.5, vmax=3)
	c=plt.colorbar(cb)
	c.set_label('G-K')
	ax.set_xlabel('H parallel')
	ax.set_ylabel('H perpendicular')
	#~ ax.set_xlim(-1, 5)
	#~ ax.set_ylim(20, -15)




	amp=100




	# Parallel component: distribution in the sky
	fig=plt.figure()
	ax=fig.add_subplot(211)
	cb=ax.scatter(result[:,2], result[:,3], c=result[:,0], vmin=-amp, vmax=amp, s=5, edgecolor='none')
	c=plt.colorbar(cb)
	c.set_label('pm parallel [mas/yr]')
	ax.set_xlabel('RA')
	ax.set_ylabel('DEC')

	#~ plt.show()
	#~ exit(0)

	result=np.array(sorted(result, key=lambda x: x[0], reverse=True))
	ax=fig.add_subplot(212)
	cb=ax.scatter(result[:,2], result[:,3], c=result[:,0], vmin=-amp, vmax=amp, s=5, edgecolor='none')
	c=plt.colorbar(cb)
	c.set_label('pm parallel [mas/yr]')
	ax.set_xlabel('RA')
	ax.set_ylabel('DEC')

	# Perp component
	result=np.array(sorted(result, key=lambda x: x[1]))
	fig=plt.figure()
	ax=fig.add_subplot(211)
	cb=ax.scatter(result[:,2], result[:,3], c=result[:,1], vmin=-amp, vmax=amp, s=5, edgecolor='none')
	c=plt.colorbar(cb)
	c.set_label('pm perp [mas/yr]')
	ax.set_xlabel('RA')
	ax.set_ylabel('DEC')

	result=np.array(sorted(result, key=lambda x: x[1], reverse=True))
	ax=fig.add_subplot(212)
	cb=ax.scatter(result[:,2], result[:,3], c=result[:,1], vmin=-amp, vmax=amp, s=5, edgecolor='none')
	c=plt.colorbar(cb)
	c.set_label('pm perp [mas/yr]')
	ax.set_xlabel('RA')
	ax.set_ylabel('DEC')

	# PM parallel vs PM perpendicular
	fig=plt.figure()
	ax=fig.add_subplot(111)
	cb=ax.scatter(result[:,1], result[:,0], s=5)
	ax.set_xlabel('pm perp [mas/yr]')
	ax.set_ylabel('pm parallel / sin(theta) [mas/yr]')
	ax.axhline(y=0, c='k')
	ax.axvline(x=0, c='k')
	ax.set_xlim(-100, 100)
	ax.set_ylim(-1000, 1000)

	# PM parallel vs PM perpendicular vs color





	# G-K histogram
	#~ fig=plt.figure()
	#~ ax=fig.add_subplot(111)
	#~ ax.hist(result[:,6], bins=100)

	plt.show()

def histograms(result):
	fig=plt.figure()
	ax=fig.add_subplot(111)
	ax.hist(result[:,0], bins=np.linspace(-200, 200, 1000))
	ax.axvline(x=0, c='k')
	ax2=ax.twinx()
	ax2.hist(result[:,0], bins=np.linspace(-200, 200, 1000), histtype='step', cumulative=True, normed=1, color='k')
	ax2.axhline(y=0.5, c='k')
	ax2.axhline(y=0.9, c='k')
	ax2.axhline(y=0.1, c='k')
	#~ [pm_parallel, pm_perp, x['ra'], x['de'], H_parallel, H_perp, m1-m2, H, ewirt]
	
	
	fig=plt.figure()
	ax=fig.add_subplot(111)
	ax.hist(result[:,1], bins=np.linspace(-100, 100, 1000))
	ax.axvline(x=0, c='k')
	ax2=ax.twinx()
	ax2.hist(result[:,1], bins=np.linspace(-100, 100, 1000), histtype='step', cumulative=True, normed=1, color='k')
	ax2.axhline(y=0.5, c='k')
	ax2.axhline(y=0.9, c='k')
	ax2.axhline(y=0.1, c='k')
	#~ [pm_parallel, pm_perp, x['ra'], x['de'], H_parallel, H_perp, m1-m2, H, ewirt]
	
	
	
	plt.show()

if __name__ == "__main__":
	reduced_proper_motion_diagrams(result, C)
	#~ histograms(result)
