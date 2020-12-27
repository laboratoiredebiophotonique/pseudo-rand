#!/usr/bin/python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as pl
import numpy as np
import scipy.spatial as sp


    def voronoi(self):
        import sys

        eps = sys.float_info.epsilon

    # Select pts inside the box
        i = in_box()

    # Mirror points
        center = self.pts[i, :]
        left   = np.copy(center)
        right  = np.copy(center)
        down   = np.copy(center)
        up     = np.copy(center)

        left[:, 0]  = self.box[0] - (left[:, 0] - self.box[0])
        right[:, 0] = self.box[1] + (self.box[1] - right[:, 0])
        down[:, 1]  = self.box[2] - (down[:, 1] - self.box[2])
        up[:, 1]    = self.box[3] + (self.box[3] - up[:, 1])

        points = np.append(center, np.append(np.append(left, right, axis=0), np.append(down, up, axis=0), axis=0), axis=0)
    # Compute Voronoi
        vor = sp.Voronoi(points)
    # Filter regions
        regions = []
        for region in vor.regions:
            flag = True
            for index in region:
                if index == -1:
                    flag = False
                    break
                else:
                    x = vor.vertices[index, 0]
                    y = vor.vertices[index, 1]
                    if not(self.box[0] - eps <= x and x <= self.box[1] + eps and
                           self.box[2] - eps <= y and y <= self.box[3] + eps):
                        flag = False
                        break
            if region != [] and flag:
                regions.append(region)
        self.filtered_points   = center
        self.filtered_regions  = regions

        return vor



    def centroid_region(self, vertices):
        A = 0
        C_x = 0
        C_y = 0
        for i in range(0, len(vertices) - 1):
            s = (vertices[i, 0] * vertices[i + 1, 1] - vertices[i + 1, 0] * vertices[i, 1])
            A = A + s
            C_x = C_x + (vertices[i, 0] + vertices[i + 1, 0]) * s
            C_y = C_y + (vertices[i, 1] + vertices[i + 1, 1]) * s
        A = 0.5 * A
        C_x = (1.0 / (6.0 * A)) * C_x
        C_y = (1.0 / (6.0 * A)) * C_y
        return np.array([[C_x, C_y]])









# a ajouter myCVT = CVT(pts, bounding_box)



vor = myCVT.voronoi()

fig = pl.figure()
ax = fig.gca()

# Plot points
ax.plot(myCVT.filtered_points[:, 0], myCVT.filtered_points[:, 1], 'b.')


# Plot ridges
for region in myCVT.filtered_regions:
    vertices = vor.vertices[region + [region[0]], :]
    ax.plot(vertices[:, 0], vertices[:, 1], color='lightgray')

# Compute and plot centroids
centroids = []
for region in myCVT.filtered_regions:
    vertices = vor.vertices[region + [region[0]], :]
    centroid = myCVT.centroid_region(vertices)
    centroids.append(list(centroid[0, :]))
    ax.plot(centroid[:, 0], centroid[:, 1], 'r.')


ax.set_xlim([-0.1, 1.1])
ax.set_ylim([-0.1, 1.1])
pl.savefig("bounded_voronoi.png")

#pl.savefig("voronoi.png")

pl.show()
