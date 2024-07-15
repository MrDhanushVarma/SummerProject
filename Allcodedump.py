# idk
# --------------------__--------------___--_--_--------Fourier SPACE ______-------------__---__-__-__-_--_----_--_----------------------------------


# 1d lattice 

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import eig

N = 1000
q = np.arange(-np.pi/2, np.pi/2, 1*np.pi/N)
h = 0.00025
L = 5

def M(q,h):
  m = []
  for i in range(L):
    c = []
    for j in range(L):
      if i == j:
        c.append((2*(np.cos((q+(j-2)*np.pi)*h) - 1))/(h**2))
      elif i - j == 1 or j - i == 1:
        c.append(1/2)
      else:
        c.append(0)
    m.append(c)

  eigen_values,v = eig(m)
  eigen_values.tolist()
  eigen_values.sort()

  return eigen_values


for i in range(L):
  band = []
  for j in q:
    band.append(np.abs((M(j,h)[i])))
  plt.plot(q,band)

plt.show()

# --------------------__--------------___--_--_--------# --------------------__--------------___--_--_--------

# 2d sq potential

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import eig
from matplotlib import cbook, cm
from matplotlib.colors import LightSource

G_11 = 1 * 2* np.pi
G_12 = 0* 2* np.pi
G_21 = 0* 2* np.pi
G_22 = 2* np.pi


def V(x,y):
  return np.cos(G_11*x + G_12*y )*np.cos(G_22*y + G_21*x)




# Define the grid of points
x = np.linspace(-1, 1, 100)
y = np.linspace(-1, 1, 100)
X, Y = np.meshgrid(x, y)

# Calculate Z using the function
Z = V(X, Y)

# Create the contour plot
contour = plt.contour(X, Y, Z, levels=200, cmap='viridis')
plt.colorbar(contour)
plt.title('V')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')



plt.tight_layout()
plt.show()

# --------------------__--------------___--_--_--------# --------------------__--------------___--_--_--------
# 2d sq lattice
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import eig
from matplotlib import cbook, cm
from matplotlib.colors import LightSource


N = 100
q_x = np.arange(-np.pi, np.pi, 2*np.pi/N)
q_y = np.arange(-np.pi, np.pi, 2*np.pi/N)


h = 0.00025


def M(q_x,q_y,h):
  m = np.zeros((9, 9), dtype=complex)
  for i in range(9):
    for j in range(9):
      if i == j and i < 3:
        m[i][j] = (2*np.cos((q_x-1*np.pi)*h) + 2*np.cos((q_y+(1*(i-1)*np.pi))*h) - 4)/h**2
      elif i == j and i<6 and i>2:
        m[i][j] = (2*np.cos((q_x)*h) + 2*np.cos((q_y+(1*(i-4)*np.pi))*h) - 4)/h**2
      elif i == j and i<9 and i>5:
        m[i][j] = (2*np.cos((q_x+1*np.pi)*h) + 2*np.cos((q_y+(1*(i-7)*np.pi))*h) - 4)/h**2


  m[0][4] = -1/4
  m[1][3] = -1/4
  m[1][5] = -1/4
  m[2][4] = -1/4
  m[3][1] = -1/4
  m[4][0] = -1/4
  m[4][2] = -1/4
  m[4][6] = -1/4
  m[4][8] = -1/4
  m[5][1] = -1/4
  m[5][7] = -1/4
  m[6][4] = -1/4
  m[7][3] = -1/4
  m[7][5] = -1/4
  m[8][4] = -1/4
  m[3][7] = -1/4


  eigen_values, _ = eig(-1*m)
  eigen_values = np.real((eigen_values))
  eigen_values.tolist()
  eigen_values.sort()

  return eigen_values


mesh = np.zeros((len(q_x), len(q_y), 9))

for i in range(len(q_x)):
    for j in range(len(q_y)):
        mesh[i, j] = M(q_x[i], q_y[j], h)

fig, axs = plt.subplots(3, 3, figsize=(15, 15))

for i, ax in enumerate(axs.flat):
  if i == 1:
    CS = ax.contourf(q_x, q_y, mesh[:, :, i+1], cmap='hot')
    ax.set_title(f'E{i+1}')
    ax.set_xlabel('$q_x$')
    ax.set_ylabel('$q_y$')
    cbar = fig.colorbar(CS)
    cbar.ax.set_ylabel('Energy')
  elif i == 2:
    CS = ax.contourf(q_x, q_y, mesh[:, :, i-1], cmap='hot')
    ax.set_title(f'E{i+1}')
    ax.set_xlabel('$q_x$')
    ax.set_ylabel('$q_y$')
    cbar = fig.colorbar(CS)
    cbar.ax.set_ylabel('Energy')
  elif i != 1 and i != 2:
    CS = ax.contourf(q_x, q_y, mesh[:, :, i], cmap='hot')
    ax.set_title(f'E{i+1}')
    ax.set_xlabel('$q_x$')
    ax.set_ylabel('$q_y$')
    cbar = fig.colorbar(CS)
    cbar.ax.set_ylabel('Energy')




plt.tight_layout()
plt.show()


import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import eig

# Parameters
N = 100
q_x = np.linspace(-np.pi, np.pi, N)
q_y = np.linspace(-np.pi, np.pi, N)
h = 0.00025

# Define the function M
def M(q_x,q_y,h):
  m = np.zeros((9, 9), dtype=complex)
  for i in range(9):
    for j in range(9):
      if i == j and i < 3:
        m[i][j] = (2*np.cos((q_x-1*np.pi)*h) + 2*np.cos((q_y+(1*(i-1)*np.pi))*h) - 4)/h**2
      elif i == j and i<6 and i>2:
        m[i][j] = (2*np.cos((q_x)*h) + 2*np.cos((q_y+(1*(i-4)*np.pi))*h) - 4)/h**2
      elif i == j and i<9 and i>5:
        m[i][j] = (2*np.cos((q_x+1*np.pi)*h) + 2*np.cos((q_y+(1*(i-7)*np.pi))*h) - 4)/h**2

  m[0][4] = -1/4
  m[1][3] = -1/4
  m[1][5] = -1/4
  m[2][4] = -1/4
  m[3][1] = -1/4
  m[4][0] = -1/4
  m[4][2] = -1/4
  m[4][6] = -1/4
  m[4][8] = -1/4
  m[5][1] = -1/4
  m[5][7] = -1/4
  m[6][4] = -1/4
  m[7][3] = -1/4
  m[7][5] = -1/4
  m[8][4] = -1/4
  m[3][7] = -1/4

  eigen_values, _ = eig(-m)
  eigen_values = np.real((eigen_values))
  eigen_values.tolist()
  eigen_values.sort()

  return eigen_values

# Create the meshgrid
q_x_grid, q_y_grid = np.meshgrid(q_x, q_y)

# Compute eigenvalues for the meshgrid
mesh = np.zeros((len(q_x), len(q_y), 9))

for i in range(len(q_x)):
    for j in range(len(q_y)):
        mesh[i, j] = M(q_x[i], q_y[j], h)

# Plotting the surface plots for each eigenvalue
fig = plt.figure(figsize=(15, 15))

for i in range(9):
    ax = fig.add_subplot(3, 3, i+1, projection='3d')
    ax.plot_surface(q_x_grid, q_y_grid, mesh[:, :, i], cmap='viridis')
    ax.set_title(f'E{i+1}')
    ax.set_xlabel('$q_x$')
    ax.set_ylabel('$q_y$')
    ax.set_zlabel('$E$')

plt.tight_layout()
plt.show()

# --------------------__--------------___--_--_--------# --------------------__--------------___--_--_--------
# 2d hexagonal lattice

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import eig
from matplotlib import cbook, cm
from matplotlib.colors import LightSource

G_11 = 1 * 2* np.pi
G_12 = (-1/3**0.5)* 2* np.pi
G_21 = 0* 2* np.pi
G_22 = (2/ 3**0.5)* 2* np.pi


def V(x,y):
  return np.cos(G_11*x + G_12*y )*np.cos(G_22*y + G_21*x)




# Define the grid of points
x = np.linspace(-1, 1, 100)
y = np.linspace(-1, 1, 100)
X, Y = np.meshgrid(x, y)

# Calculate Z using the function
Z = V(X, Y)

# Create the contour plot
contour = plt.contour(X, Y, Z, levels=200, cmap='viridis')
plt.colorbar(contour)
plt.title('V')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')



plt.tight_layout()
plt.show()


import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import eig
from matplotlib import cbook, cm
from matplotlib.colors import LightSource


N = 200
q_x = np.arange(-2*np.pi, 2*np.pi, 4*np.pi/N)
q_y = np.arange(-2*np.pi, 2*np.pi, 4*np.pi/N)


h = 0.00025


def M(q_x,q_y,h):
  m = np.zeros((9, 9), dtype=complex)
  for i in range(9):
    for j in range(9):
      if i == j and i < 3:
        m[i][j] = (2*np.cos((q_x-2*np.pi)*h) + 2*np.cos((q_y+(4*(i-1)*np.pi/(3**0.5))+(2*np.pi/(3**0.5)))*h) - 4)/h**2
      elif i == j and i<6 and i>2:
        m[i][j] = (2*np.cos((q_x)*h) + 2*np.cos((q_y+(4*(i-4)*np.pi/(3**0.5)))*h) - 4)/h**2
      elif i == j and i<9 and i>5:
        m[i][j] = (2*np.cos((q_x+2*np.pi)*h) + 2*np.cos((q_y+(4*(i-7)*np.pi/(3**0.5))-(2*np.pi/(3**0.5)))*h) - 4)/h**2
  m[0][4] = -1/4
  m[1][3] = -1/4
  m[1][5] = -1/4
  m[2][4] = -1/4
  m[3][1] = -1/4
  m[4][0] = -1/4
  m[4][2] = -1/4
  m[4][6] = -1/4
  m[4][8] = -1/4
  m[5][1] = -1/4
  m[5][7] = -1/4
  m[6][4] = -1/4
  m[7][3] = -1/4
  m[7][5] = -1/4
  m[8][4] = -1/4
  m[3][7] = -1/4
  eigen_values,v = eig(-m)
  eigen_values =  np.real((eigen_values))
  eigen_values.tolist()
  eigen_values.sort()
  # print(m)
  return eigen_values


mesh = np.zeros((len(q_x), len(q_y), 9))

for i in range(len(q_x)):
    for j in range(len(q_y)):
        mesh[i, j] = M(q_x[i], q_y[j], h)

fig, axs = plt.subplots(3, 3, figsize=(15, 15))

for i, ax in enumerate(axs.flat):
  if i == 1:
    CS = ax.contourf(q_x, q_y, mesh[:, :, i+1], cmap='hot')
    ax.set_title(f'E{i+1}')
    ax.set_xlabel('$q_x$')
    ax.set_ylabel('$q_y$')
    cbar = fig.colorbar(CS)
    cbar.ax.set_ylabel('Energy')
  elif i == 2:
    CS = ax.contourf(q_x, q_y, mesh[:, :, i-1], cmap='hot')
    ax.set_title(f'E{i+1}')
    ax.set_xlabel('$q_x$')
    ax.set_ylabel('$q_y$')
    cbar = fig.colorbar(CS)
    cbar.ax.set_ylabel('Energy')
  elif i != 1 and i != 2:
    CS = ax.contourf(q_x, q_y, mesh[:, :, i], cmap='hot')
    ax.set_title(f'E{i+1}')
    ax.set_xlabel('$q_x$')
    ax.set_ylabel('$q_y$')
    cbar = fig.colorbar(CS)
    cbar.ax.set_ylabel('Energy')


plt.tight_layout()
plt.show()


import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import eig

# Parameters
N = 100
q_x = np.linspace(-np.pi, np.pi, N)
q_y = np.linspace(-np.pi, np.pi, N)
h = 0.00025

# Define the function M
def M(q_x, q_y, h):
  m = np.zeros((9, 9), dtype=complex)
  for i in range(9):
    for j in range(9):
      if i == j and i < 3:
        m[i][j] = (2*np.cos((q_x-2*np.pi)*h) + 2*np.cos((q_y+(4*(i-1)*np.pi/(3**0.5))+(2*np.pi/(3**0.5)))*h) - 4)/h**2
      elif i == j and i<6 and i>2:
        m[i][j] = (2*np.cos((q_x)*h) + 2*np.cos((q_y+(4*(i-4)*np.pi/(3**0.5)))*h) - 4)/h**2
      elif i == j and i<9 and i>5:
        m[i][j] = (2*np.cos((q_x+2*np.pi)*h) + 2*np.cos((q_y+(4*(i-7)*np.pi/(3**0.5))-(2*np.pi/(3**0.5)))*h) - 4)/h**2

  m[0][4] = -1/4
  m[1][3] = -1/4
  m[1][5] = -1/4
  m[2][4] = -1/4
  m[3][1] = -1/4
  m[4][0] = -1/4
  m[4][2] = -1/4
  m[4][6] = -1/4
  m[4][8] = -1/4
  m[5][1] = -1/4
  m[5][7] = -1/4
  m[6][4] = -1/4
  m[7][3] = -1/4
  m[7][5] = -1/4
  m[8][4] = -1/4
  m[3][7] = -1/4

  eigen_values, _ = eig(-m)
  eigen_values = np.real((eigen_values))
  eigen_values.tolist()
  eigen_values.sort()

  return eigen_values

# Create the meshgrid
q_x_grid, q_y_grid = np.meshgrid(q_x, q_y)

# Compute eigenvalues for the meshgrid
mesh = np.zeros((len(q_x), len(q_y), 9))

for i in range(len(q_x)):
    for j in range(len(q_y)):
        mesh[i, j] = M(q_x[i], q_y[j], h)

# Plotting the surface plots for each eigenvalue
fig = plt.figure(figsize=(15, 15))

for i in range(9):
    ax = fig.add_subplot(3, 3, i+1, projection='3d')
    ax.plot_surface(q_x_grid, q_y_grid, mesh[:, :, i], cmap='viridis')
    ax.set_title(f'E{i+1}')
    ax.set_xlabel('$q_x$')
    ax.set_ylabel('$q_y$')
    ax.set_zlabel('$E$')

plt.tight_layout()
plt.show()



# --------------------__--------------___--_--_--------REAL SPACE ______-------------__---__-__-__-_--_----_--_----------------------------------
# 1d lattice(working)

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import eig

N = 50
q = np.linspace(-1,1,num=N)

x_values = np.linspace(-np.pi,np.pi,num=N)


def H(q_x):



  h=x_values[1]-x_values[0]
  Dx=-np.eye(N)+np.diagflat(np.ones(N-1),1)
  Dx[-1,0] = 1

  Dx = Dx / h




  # D=np.diagflat(np.ones(n_grid-1),1)
  Dx[-1,-1]=Dx[0,0]

  D2x=Dx.dot(-Dx.T)
  D2x[-1,-1]=D2x[0,0]




#100x100 matrix
  V_values = 0.5*np.cos(2*np.pi*x_values*h)
  V = np.diag(V_values)

  T3 = np.identity(N)

  # H = T_1(q) + T_2(q) + V(q) + Q(q)
  H =   -(D2x) + q_x*1j*Dx  + V  + (q_x**2)*T3*h
  eigen_values, _ = eig(H)
  eigen_values = (eigen_values)
  eigen_values.tolist()
  eigen_values.sort()
  return eigen_values


for i in range(N):
  if i<6:
    band = []
    for j in q:
      band.append(H(j)[i])
    plt.plot(q,band)

plt.show()

# --------------------__--------------___--_--_--------# --------------------__--------------___--_--_--------
# 2d sq lattice(not working)


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from numpy.linalg import eig


%matplotlib inline
sns.set_style("white")


n_grid=5
N = n_grid**2

k_x = np.arange(-2, 2, 0.1)
k_y = np.arange(-2, 2, 0.1)
x_values = np.linspace(-np.pi,np.pi,num=N)
y_values = np.linspace(-np.pi,np.pi,num=N)

def H(q_x,q_y):



  h=x_values[1]-x_values[0]
  Dx=-np.eye(N)+np.diagflat(np.ones(N-1),1)
  Dx[-1,0] = 1
  Dx = Dx / h




  # D=np.diagflat(np.ones(n_grid-1),1)
  # Dx[-1,-1]=Dx[0,0]

  D2x=Dx.dot(-Dx.T)
  # D2x[-1,-1]=D2x[0,0]



  hy=(y_values[1]-y_values[0])*n_grid
  Dy=-np.eye(N)
  for i in range(N):
    for j in range(N):
      if i - j == -n_grid:
        Dy[i][j] = 1

  Dy[-1,n_grid] = 1

  Dy = Dy / hy
  # print(D)

  # Dy[-1,-1]=Dy[0,0]



  D2y=Dy.dot(-Dy.T)
  # D2y[-1,-1]=D2y[0,0]

#100x100 matrix
  V_values = np.cos(2*np.pi*x_values*h)*np.cos(2*np.pi*y_values*hy)
  V = np.diag(V_values)

  T3 = np.identity(N)

  # H = T_1(q) + T_2(q) + V(q) + Q(q)
  H =   -(D2y+D2x) - q_x*1j*Dx - q_y*1j*Dy  + V  + (q_y**2+q_x**2)*T3
  eigen_values, _ = eig(H)
  eigen_values = (eigen_values)
  eigen_values.tolist()
  eigen_values.sort()
  return eigen_values

q_x_grid, q_y_grid = np.meshgrid(k_x, k_y)


mesh = np.zeros((len(k_x), len(k_y), N))


for i in range(len(k_x)):
    for j in range(len(k_y)):
      mesh[i, j] = H(k_x[i],k_y[j])

fig, axs = plt.subplots(3, 3, figsize=(15, 15))

for i, ax in enumerate(axs.flat):
  CS = ax.contourf(k_x, k_y, mesh[:, :, i], cmap='hot')
  ax.set_title(f'E{i+1}')
  ax.set_xlabel('$q_x$')
  ax.set_ylabel('$q_y$')
  cbar = fig.colorbar(CS)
  cbar.ax.set_ylabel('Energy')

plt.tight_layout()
plt.show()

# --------------------__--------------___--_--_--------# --------------------__--------------___--_--_--------
# Hofstader butterfly(working)

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import eig
from fractions import Fraction



def M(E,a,v,m1):
  m = np.identity(2)
  m[0][0] =E - 2*np.cos((2*np.pi*m1*a)-v)
  m[1][0] = 1
  m[0][1] = -1
  m[1][1] = 0
  return m

def Q(E,a,v,m1):
  x = []
  for i in range(50):
    n = Fraction(f'{a}').limit_denominator(i+1)
    x.append(n)

  y = [x[i] for i in range(len(x)) if abs((x[i].numerator/x[i].denominator)-a) < 0.001]
  q = y[0].denominator
  p = y[0].numerator

  # v = 1/(2*q)

  r = M(E,a,v,m1)
  for i in range(q-1):
    r = r@M(E,a,v,m1-(i+1))

  return r

a = np.arange(0,1,0.01)
eps = np.arange(-4,4,0.01)
k_y = np.arange(0,2*np.pi,1)

# print(a)

alpha = []
energies = []
for v in k_y :
  for j in eps:
    for i in a:
      try:
        if abs(np.trace(Q(j,i,v,9))) <= 4:
          alpha.append(i)
          energies.append(j)
      except IndexError:
          # print("Hey QT, its an irrational number.")
          pass


plt.scatter(energies,alpha)
plt.xlim(-4,4)

plt.show()

# --------------------__--------------___--_--_--------# --------------------__--------------___--_--_--------
# chern no. 2d sq lattice(working)

import numpy as np

# Parameters
M = 1.0
N = 100  # Grid size
kx, ky = np.meshgrid(np.linspace(-np.pi, np.pi, N), np.linspace(-np.pi, np.pi, N))

# Hamiltonian components
d_x = np.sin(kx)
d_y = np.sin(ky)
d_z = M + np.cos(kx) + np.cos(ky)
d_norm = np.sqrt(d_x**2 + d_y**2 + d_z**2)

# Berry curvature calculation
def berry_curvature(kx, ky, d_x, d_y, d_z, d_norm):
    # Partial derivatives
    d_dz_dkx = -np.sin(kx)
    d_dz_dky = -np.sin(ky)

    d_dz_dkx = -np.sin(kx)
    d_dz_dky = -np.sin(ky)

    # Berry curvature
    F_xy = (
        (d_dz_dkx * d_dz_dky)
        / d_norm**3
    )
    return F_xy

# Compute Berry curvature
F_xy = berry_curvature(kx, ky, d_x, d_y, d_z, d_norm)

# Integrate over the Brillouin zone
C = np.sum(F_xy) * (2 * np.pi / N)**2 / (2 * np.pi)
print(f"Chern number: {round(C)}")


