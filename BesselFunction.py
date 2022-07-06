from matplotlib.pyplot import legend
import scipy.special as ssf
import numpy as np
import matplotlib.pylab as plt

# Choose the kind of function to be plotted
# kind = 1  # Bessel function of the first kind
# kind = 2  # Bessel function of the second kind
# kind = 3  # Hankel function of the first kind
kind = 4  # Hankel function of the second kind


z = np.linspace(2, 20, 1000)

fig, ax = plt.subplots()
for i in np.arange(0,4,1):
  if kind == 1 :
    ax.plot(z, ssf.jv(i, z), label='m='+str(i))
    ax.set_title('Bessel function of the first kind')
  elif kind == 2 :
    ax.plot(z, ssf.yv(i, z), label='m='+str(i))
    ax.set_title('Bessel function of the second kind')
  elif kind == 3 :
    ax.plot(z, ssf.hankel1(i, z).imag, label='m='+str(i))
    ax.set_title('Hankel function of the first kind')
  elif kind == 4 :
    ax.plot(z, ssf.hankel2(i, z).imag, label='m='+str(i))
    ax.set_title('Hankel function of the second kind')
ax.legend()

plt.show()