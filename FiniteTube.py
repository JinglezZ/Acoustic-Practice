import numpy as np
import scipy.special as ssf
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

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
# length
L = 0.2
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
F0 = 1000
# order
n = 2
# frequency
f = 3000


# derivative of the Hankel function of the first kind
dHm = lambda n, z: n*ssf.hankel1(n, z)/z-ssf.hankel1(n+1,z)
j0 = lambda x : np.sin(x)/x

omega = 2. * np.pi * f
k = omega/c
omega = 2. * np.pi * f

Omega = omega*a/cp
Omega1 = np.sqrt( 0.5*(1+n**2+beta**2*n**4+np.sqrt((1+n**2+beta**2*n**4)**2-4*beta**2*n**6)) )
Omega2 = np.sqrt( 0.5*(1+n**2+beta**2*n**4-np.sqrt((1+n**2+beta**2*n**4)**2-4*beta**2*n**6)) )

Zp = rho*c*omega*ssf.hankel1(n,k*a)/dHm(n,k*a)
Zm = -E*h/(a**2*(1-sigma**2)) * (Omega**2-Omega1**2)*(Omega**2-Omega2**2)/(Omega**2-n**2)

w = F0/(Zp+Zm)
print(np.absolute(w))
u0 = -1j*omega*w

(n, m) = (50, 50)

# Meshing a unit sphere according to n, m 
theta = np.linspace(np.pi/6, 0.99 * np.pi, num=n, endpoint=False)
phi = np.linspace(np.pi * (-0.5 + 1./(m+1)), np.pi*0.5, num=m, endpoint=False)
theta, phi = np.meshgrid(theta, phi)
theta, phi = theta.ravel(), phi.ravel()
#theta = np.append(theta, [0.]) # Adding the north pole...
#phi = np.append(phi, [np.pi*0.5])
mesh_x, mesh_y = ((np.pi*0.5 - phi)*np.cos(theta), (np.pi*0.5 - phi)*np.sin(theta))
triangles = mtri.Triangulation(mesh_x, mesh_y).triangles
x, y, z = np.cos(phi)*np.cos(theta), np.cos(phi)*np.sin(theta), np.sin(phi)

# Defining a custom color scalar field
R = 10
P = -1j*2*L*omega*u0*rho*np.exp(1j*k*R)/(np.pi*k*R*np.sin(theta))#*(j0(k*np.cos(theta))*(-1j)**(n+1)*np.cos(n*phi))/(dHm(n, k*a*np.sin(theta)))
vals = np.absolute(P)
print(vals)
colors = np.mean(vals[triangles], axis=1)

# Plotting
fig = plt.figure()
ax = fig.gca(projection='3d')
cmap = plt.get_cmap('jet')
triang = mtri.Triangulation(x, y, triangles)
collec = ax.plot_trisurf(triang, z, cmap=cmap, shade=False, linewidth=0.)
collec.set_array(colors)
collec.autoscale()
plt.show()