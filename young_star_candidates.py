import numpy as np
import math
import pickle
import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord
from astropy import units as u

# Convergence point
# This is solar apex; replace it with young stars
#~ c = SkyCoord(l=50.2*u.degree, b=24.7*u.degree, frame='galactic')
#~ RA_CP=c.icrs.ra.value
#~ DEC_CP=c.icrs.dec.value

RA_CP=270.0
DEC_CP=30.0
DEC_CP=0.0

# Read the data
pkl_file = open('gaia_2mass_for_funnelweb.pkl', 'rb')
d = pickle.load(pkl_file)
pkl_file.close()


def determine_theta(alpha_cp=None, delta_cp=None, alpha=None, delta=None):
	""" Determine angle between equatorial system and system alligned with direction of motion of nearby young stars
	
	Args:
	alpha_cp: [deg] right ascension of convergence point
	delta_cp: [deg] declination of convergence point
	alpha: [deg] right ascension of a star
	delta: [deg] declination of a star
	
	Return: theta [deg]: angle between the two systems
	"""
	
	# Convert to radians
	d=np.deg2rad(delta)
	a=np.deg2rad(alpha)
	
	dcp=np.deg2rad(delta_cp)
	acp=np.deg2rad(alpha_cp)
	
	# Angular distance A between the star and convergence point CP
	cos_A = np.sin(dcp)*np.sin(d) + np.cos(dcp)*np.cos(d)*np.cos(acp-a)
	sin_A = np.sqrt(1.0-cos_A**2)
	
	# Spherical triangle with the following vertices: north celestial pole, a star (alpha, delta) and convergence point (alpha_cp, delta_cp)
	cos_theta = (np.sin(dcp) - cos_A*np.sin(d)) / (sin_A*np.cos(d))
	theta=np.arccos(cos_theta)
	theta=np.rad2deg(theta)

	
	sin_theta=np.sin(acp-a)/sin_A*np.cos(dcp) # TODO: np.sin(acp-a) sign???
	theta2=np.arcsin(sin_theta)
	theta2=np.rad2deg(theta2)
	
	theta = np.rad2deg(math.atan2(sin_theta, cos_theta))
	
	#~ print (cos_theta, sin_theta, theta, theta2, 180.0-theta2, np.rad2deg(math.atan2(sin_theta, cos_theta)))

	return theta

def components_in_a_new_system(theta=None, pmRA=None, pmDE=None):
	""" Determine components of proper motion vector in a new coordinate system
	
	"""
	beta = np.arctan(pmDE/pmRA) - theta
	pm=np.sqrt(pmRA**2+pmDE**2)
	pm_parallel=pm*np.cos(beta)
	pm_perp=pm*np.sin(beta)
	return pm_parallel, pm_perp

def rotate_vector(vec=None, theta=None):
	""" Rotate 2D vector for angle theta
	
	Args:
	vec: 2D vector (np.array)
	theta: rotation angle [deg]
	
	return: rotated vector vec2
	
	"""

	theta=np.deg2rad(theta)
	R=np.array([[np.sin(theta), np.cos(theta)], [-np.cos(theta), np.sin(theta)]]) # rotation matrix
	vec2=R.dot(vec)
	return vec2

result=[]
for x in d:
	# pmra = mu_alpha* (meaning that cos(delta) has already been taken into account)
	
	# Determine angle theta between the two coordinate systems
	theta=determine_theta(alpha_cp=RA_CP, delta_cp=DEC_CP, alpha=x['ra'], delta=x['de'])
	
	# Method 1: Components in a new coordinate system
	#~ pm_parallel, pm_perp = components_in_a_new_system(theta=theta, pmRA=x['pmra'], pmDE=x['pmdec'])
	
	# Method 2: Rotate vector for angle theta (rotation is done for -theta)
	rotated_vector=rotate_vector(vec=np.array([x['pmra'], x['pmdec']]), theta=theta)
	pm_parallel=rotated_vector[0]
	pm_perp=rotated_vector[1]
	
	result.append([pm_parallel, pm_perp, x['ra'], x['de']])
result=np.array(sorted(result, key=lambda x: x[0]))

amp=100

# Parallel component: distribution in the sky
fig=plt.figure()
ax=fig.add_subplot(211)
cb=ax.scatter(result[:,2], result[:,3], c=result[:,0], vmin=-amp, vmax=amp, s=5)
c=plt.colorbar(cb)
c.set_label('pm parallel [mas/yr]')
ax.set_xlabel('RA')
ax.set_ylabel('DEC')

result=np.array(sorted(result, key=lambda x: x[0], reverse=True))
ax=fig.add_subplot(212)
cb=ax.scatter(result[:,2], result[:,3], c=result[:,0], vmin=-amp, vmax=amp, s=5)
c=plt.colorbar(cb)
c.set_label('pm parallel [mas/yr]')
ax.set_xlabel('RA')
ax.set_ylabel('DEC')

# Perp component
result=np.array(sorted(result, key=lambda x: x[1]))
fig=plt.figure()
ax=fig.add_subplot(211)
cb=ax.scatter(result[:,2], result[:,3], c=result[:,1], vmin=-amp, vmax=amp, s=5)
c=plt.colorbar(cb)
c.set_label('pm perp [mas/yr]')
ax.set_xlabel('RA')
ax.set_ylabel('DEC')

result=np.array(sorted(result, key=lambda x: x[1], reverse=True))
ax=fig.add_subplot(212)
cb=ax.scatter(result[:,2], result[:,3], c=result[:,1], vmin=-amp, vmax=amp, s=5)
c=plt.colorbar(cb)
c.set_label('pm perp [mas/yr]')
ax.set_xlabel('RA')
ax.set_ylabel('DEC')

# PM parallel vs PM perpendicular
fig=plt.figure()
ax=fig.add_subplot(111)
cb=ax.scatter(result[:,1], result[:,0], s=5)
ax.set_xlabel('pm perp [mas/yr]')
ax.set_ylabel('pm parallel [mas/yr]')
ax.axhline(y=0, c='k')
ax.axvline(x=0, c='k')

plt.show()
