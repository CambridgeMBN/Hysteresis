#import matplotlib.pyplot as plt
import matplotlib, sys
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

arpath  = 'Arrow_Data/'

field = '-1'
temp = '6'
if len(sys.argv) > 1:
    file = sys.argv[1]
else:
    file = 'Arrows' + '_' + field + '_' + temp
path = arpath + file
dat = np.loadtxt(path)
#dat = np.loadtxt(arpath + 'Arrows49_10K')

ardat = np.delete(dat, [0,len(dat)-1])

fig = plt.figure()
ax = fig.gca(projection='3d')

class Arrow3D(FancyArrowPatch):
    
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

cols = (['k'] * 17) + (['r'] * 26) + (['k'] * 18)

for i in range(len(ardat)):
    a = Arrow3D([0,0.5*np.sin(ardat[i])],[0,0.5*np.cos(ardat[i])],[i,i],
                mutation_scale=20, arrowstyle="-", color=cols[i])
    ax.add_artist(a)

b = Arrow3D([0,0], [0,-0.75], [-5,-5], arrowstyle='-|>', mutation_scale=20, color='k')
ax.add_artist(b)

ax.set_xlim(-1,1)
ax.set_ylim(-1,1)
ax.set_zlim(0,70)

ax.view_init(26, 0) #elev, azim

plt.axis('off')

#plt.savefig('Arrow_Plots/Arrow_3D.pdf', bbox_inches='tight', transparent=True)
plt.savefig('Arrow_Plots/Arrow_%s_%s.pdf' %(field, temp), bbox_inches='tight', transparent=True)

#plt.show()
