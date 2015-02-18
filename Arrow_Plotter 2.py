import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

class Arrow3D(FancyArrowPatch):
    
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

arpath  = 'Arrow_Data/'

#fields = ['-49', '-1', '-6', '49', '1', '6']
#temps = ['6', '30', '37', '40', '50']

fields = ['32.5']
temps = ['6.40']

cols = (['k'] * 18) + (['r'] * 26) + (['k'] * 18)

for temp in temps:

    print '-----------------\nTemp = %s\n________________' %temp
    
    for field in fields:

        print 'Field = %s' %field

        dat = np.loadtxt(arpath + 'Arrows' + '_' + field + '_' + temp)
        ardat = np.delete(dat, [0,len(dat)-1])

        fig = plt.figure()
        ax = fig.gca(projection='3d')

        z = Arrow3D([0,0], [0,0], [-5,100], mutation_scale = 20, arrowstyle = '->')
        x = Arrow3D([-0.6, 0.6], [0,0], [0,0], mutation_scale = 20, arrowstyle = '->')
        y = Arrow3D([0,0], [0.6,-0.6], [0,0], mutation_scale = 20, arrowstyle = '->')
        #ax.add_artist(z)
        ax.add_artist(x)
        ax.add_artist(y)

        for i in range(len(ardat)):
            a = Arrow3D([0,0.5*np.sin(ardat[i])],[0,0.5*np.cos(ardat[i])],[1.5*i,1.5*i],
                        mutation_scale=20, arrowstyle="-", color=cols[i])
            ax.add_artist(a)

        #  Direction trace
        #z_plots = range(len(ardat))
        #plt.plot(0.5*np.sin(ardat), 0.5*np.cos(ardat), z_plots)

        # Field direction arrow
        #b = Arrow3D([0,0], [0,-0.5], [-5,-5], arrowstyle='-|>', mutation_scale=20, color='k')
        #ax.add_artist(b)

        ax.set_xlim(-1,1)
        ax.set_ylim(-1,1)
        ax.set_zlim(0,60)

        #ax.view_init(56, 71) #elev, azim

        plt.axis('off')

        #plt.savefig('Arrow_Plots/Arrow_3D.pdf', bbox_inches='tight', transparent=True)
        plt.savefig('Arrow_Plots/Arrow_%s_%s.pdf' %(field, temp), bbox_inches='tight', transparent=True)
        plt.show()
        plt.clf()
        plt.close()
        print 'Done'
        #plt.show()
