from matplotlib import pyplot as plt
import numpy as np
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import animation

fig = plt.figure(figsize=(10, 3))
Plot = p3.Axes3D(fig)
Plot.view_init(25,75)

def gen(n):
    phi = 0
    while phi < 2*np.pi:
        yield np.array([np.cos(phi), np.sin(phi), phi])
        phi += 2*np.pi/n


def setup(data,PlotE):
    PlotE = Plot.scatter(data[0, 0:1], data[1, 0:1], data[2, 0:1], animated=True)
    #PlotE.set_3d_properties(data[2, :num])

def update(num, data,PlotE):
    PlotE.set_offsets(data[:2, :num])
    #PlotE.set_3d_properties(data[2, :num])


N = 100
data = np.array(list(gen(N))).T

PlotE = Plot.scatter(data[0, 0:1], data[1, 0:1], data[2, 0:1], animated=True)
# Setting the axes properties
Plot.set_xlim3d([-1.0, 1.0])
Plot.set_xlabel('X')

Plot.set_ylim3d([-1.0, 1.0])
Plot.set_ylabel('Y')

Plot.set_zlim3d([0.0, 10.0])
Plot.set_zlabel('Z')

ani = animation.FuncAnimation(fig, update, N, fargs=(data, PlotE), interval=10000/N, blit=False)
ani.save('matplot003.gif', writer='imagemagick')
plt.show()