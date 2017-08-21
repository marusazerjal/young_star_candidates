import numpy as np
import math
import random
import matplotlib.pyplot as plt

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

def distance_between_two_points_in_the_sky(alpha1=None, delta1=None, alpha2=None, delta2=None):
	a1=np.deg2rad(alpha1)
	d1=np.deg2rad(delta1)
	a2=np.deg2rad(alpha2)
	d2=np.deg2rad(delta2)
	
	cos_A = np.sin(d1)*np.sin(d2) + np.cos(d1)*np.cos(d2)*np.cos(a1-a2)
	A=np.arccos(cos_A) # [0, 180] in radians
	#~ A=np.rad2deg(A)
	return A

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
	
	# Angular distance A between the star and convergence point CP (spherical cosine law)
	cos_A = np.sin(dcp)*np.sin(d) + np.cos(dcp)*np.cos(d)*np.cos(acp-a)
	sin_A = np.sqrt(1.0-cos_A**2)
	
	# Spherical triangle with the following vertices: north celestial pole, a star (alpha, delta) and convergence point (alpha_cp, delta_cp)
	cos_theta = (np.sin(dcp) - cos_A*np.sin(d)) / (sin_A*np.cos(d))
	
	# The same spherical triangle, sine law
	#~ sin_theta=np.sin(acp-a)/sin_A*np.cos(dcp) # TODO: np.sin(acp-a) sign???
	sin_theta=np.sin(acp-a)/sin_A*np.cos(dcp) # TODO: np.sin(acp-a) sign???

	theta = np.rad2deg(math.atan2(sin_theta, cos_theta))
	
	#~ print (cos_theta, sin_theta, theta, theta2, 180.0-theta2, np.rad2deg(math.atan2(sin_theta, cos_theta)))

	return theta
def determine_theta2(alpha_cp=None, delta_cp=None, alpha=None, delta=None):
	""" Determine angle between equatorial system and system alligned with direction of motion of nearby young stars
	
	Args:
	alpha_cp: [deg] right ascension of convergence point
	delta_cp: [deg] declination of convergence point
	alpha: [deg] right ascension of a star
	delta: [deg] declination of a star
	
	Return: theta [deg]: angle between the two systems
	"""
	
	# Distances between the points
	# Distance between CP and a star
	d1=distance_between_two_points_in_the_sky(alpha1=RA_CP, delta1=DEC_CP, alpha2=alpha, delta2=delta)
	# Distance between CP and North pole
	d2=distance_between_two_points_in_the_sky(alpha1=RA_CP, delta1=DEC_CP, alpha2=RA_CP, delta2=90)
	# Distance between a star and North pole
	d3=distance_between_two_points_in_the_sky(alpha1=alpha, delta1=delta, alpha2=alpha, delta2=90)
	
	cos_theta=(np.cos(d2)-np.cos(d1)*np.cos(d3)) / (np.sin(d1)*np.sin(d3))
	
	# Convert to radians
	d=np.deg2rad(delta)
	a=np.deg2rad(alpha)
	
	dcp=np.deg2rad(delta_cp)
	acp=np.deg2rad(alpha_cp)
	
	
	# The same spherical triangle, sine law
	sin_theta=np.sin(acp-a)/np.sin(d1)*np.sin(d2) # TODO: np.sin(acp-a) sign???

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


### TESTS ##################
random.seed(a=9)

def test_theta_with_random_points_all_over_the_sky():
	RA_CP=270.0
	DEC_CP=30.0

	N=50000

	result=[]
	for x in range(N):
		# Random point in the sky
		alpha=random.uniform(0, 360)
		delta=random.uniform(-90, 90)
		
		# Determine angle theta between the two coordinate systems
		theta=determine_theta(alpha_cp=RA_CP, delta_cp=DEC_CP, alpha=alpha, delta=delta)
		
		result.append([alpha, delta, theta])

	result=np.array(result)

	# Theta: distribution in the sky
	fig=plt.figure()
	ax=fig.add_subplot(211)
	cb=ax.scatter(result[:,0], result[:,1], c=result[:,2], s=5)
	c=plt.colorbar(cb)
	c.set_label('theta')
	ax.set_xlabel('RA')
	ax.set_ylabel('DEC')

	ax=fig.add_subplot(212)
	cb=ax.scatter(result[:,0], result[:,1], c=(np.cos([np.deg2rad(x) for x in result[:,2]])), s=5)
	c=plt.colorbar(cb)
	c.set_label('cos(theta)')
	ax.set_xlabel('RA')
	ax.set_ylabel('DEC')
	
	# Points at RA approx RA_CP to see if they have values 0 or +/- 180 or something else
	test=result[np.abs(result[:,0]-RA_CP)<4,:]
	#~ test=result[np.abs(result[:,0]-300)<0.5,:]
	fig=plt.figure()
	ax=fig.add_subplot(111)
	ax.scatter(test[:,1], test[:,2])
	ax.set_xlabel('Dec')
	ax.set_ylabel('Theta [deg]')
	
	# Points at constant DEC
	test=result[np.abs(result[:,1]-0)<0.5,:]
	fig=plt.figure()
	ax=fig.add_subplot(111)
	ax.scatter(test[:,0], test[:,2])
	ax.set_xlabel('RA')
	ax.set_ylabel('Theta [deg]')
	
	plt.show()


if __name__ == "__main__":
	test_theta_with_random_points_all_over_the_sky()
