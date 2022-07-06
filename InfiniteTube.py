import numpy as np
import scipy.special as ssf
import matplotlib.pyplot as plt

# define parameters of the shell
# density
rho_s = 7850
# elasticity modulus
E = 200e9
# Poission's ration
sigma = 0.3
# radius of the neutral surface
a = 0.5
# thickness
h= 0.005
# longitudinal wave velocity
cp = np.sqrt(E/(rho_s*(1-sigma**2)))
# thickness factor
beta = np.sqrt(h**2/(12*a**2))

# define parameters of water
# density
rho = 1000
# speed of sound
c = 1500

# define parameters of EM excitations
# amplitude
F0 = 100000
# order
n = 0
# frequency
f = 50


# derivative of the Hankel function of the first kind
dHm = lambda n, z: n*ssf.hankel1(n, z)/z-ssf.hankel1(n+1,z)

# plot mode
mode = 1 # plot impedance and vibration curve
# mode = 2 # plot sound pressure contour


if (mode == 1) :
  f = np.arange(1, 3000, 100)
  omega = 2. * np.pi * f
  k = omega/c
  omega = 2. * np.pi * f

  Omega = omega*a/cp
  Omega1 = np.sqrt( 0.5*(1+n**2+beta**2*n**4+np.sqrt((1+n**2+beta**2*n**4)**2-4*beta**2*n**6)) )
  Omega2 = np.sqrt( 0.5*(1+n**2+beta**2*n**4-np.sqrt((1+n**2+beta**2*n**4)**2-4*beta**2*n**6)) )

  Zp = rho*c*omega*ssf.hankel1(n,k*a)/dHm(n,k*a)
  Zm = -E*h/(a**2*(1-sigma**2)) * (Omega**2-Omega1**2)*(Omega**2-Omega2**2)/(Omega**2-n**2)

  w1 = F0/(Zp+Zm)
  w2 = F0/Zm

  fig, (ax1, ax2) = plt.subplots(2)
  ax1.plot(f, np.absolute(Zp), label='Zp')
  ax1.plot(f, Zm, label='Zm')
  ax2.plot(f, np.absolute(w1), label='w with acoustic impedance')
  ax2.plot(f, np.absolute(w2), label='w without acoustic impedance')
  ax1.legend()
  ax2.legend()
  ax2.set_yscale('log')
  plt.show()
elif (mode == 2) :
  n = 2
  F0 = 100
  f = 1000

  omega = 2. * np.pi * f
  k = omega/c
  omega = 2. * np.pi * f

  Omega = omega*a/cp
  Omega1 = np.sqrt( 0.5*(1+n**2+beta**2*n**4+np.sqrt((1+n**2+beta**2*n**4)**2-4*beta**2*n**6)) )
  Omega2 = np.sqrt( 0.5*(1+n**2+beta**2*n**4-np.sqrt((1+n**2+beta**2*n**4)**2-4*beta**2*n**6)) )

  Zp = rho*c*omega*ssf.hankel1(n,k*a)/dHm(n,k*a)
  Zm = -E*h/(a**2*(1-sigma**2)) * (Omega**2-Omega1**2)*(Omega**2-Omega2**2)/(Omega**2-n**2)

  w = F0/(Zp+Zm)
  #print(np.absolute(w))
  u0 = -1j*omega*w

  azimuths = np.radians(np.linspace(0, 360, 72))
  zeniths = np.linspace(0.5, 10, 100)
  r, theta = np.meshgrid(zeniths, azimuths)

  P = 1j*rho*c*u0/dHm(n, k*a)*np.cos(n*theta)*ssf.hankel1(n, k*r)
  print('Max sound pressure ',np.amax(np.absolute(P)))

  fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
  Pre = plt.contourf(theta, r, np.absolute(P), levels=100, cmap='jet')
  ax.set_rmax(10)
  ax.set_rmin(0)
  plt.colorbar()
  plt.show()

















