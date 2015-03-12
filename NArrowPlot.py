from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import math

plt.close("all")

def imp(file):
    
    with open(''.join(['Arrow_Data/', file])) as f:
        alist = [float(line.rstrip()) for line in f]
    x = []
    y = []
    z = []
    j = 0
    print 'alist: ', alist
    for i in range(0, len(alist) -1):
        print 'phi: ', i
        phi = alist[i]
        x.append(math.cos(phi))
        y.append(math.sin(phi))
        z.append(j)
        j -= 1
    return x,y,z

x,y,z = imp('Arrows_-50_6.40')
x2, y2, z2 = imp('Arrows_-5_6.40')
# print alist



fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(x, y, z, c='r', marker='o')
# ax.scatter(x2, y2, z2, c='b', marker='x')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()
plt.close()