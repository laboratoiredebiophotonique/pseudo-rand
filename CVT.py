# -*- coding: utf-8 -*-

import numpy as np
from scipy.spatial import cKDTree


def in_box(pts, box):
    x_in = np.logical_and(box[0] <= pts[:, 0], pts[:, 0] <= box[1])
    y_in = np.logical_and(box[2] <= pts[:, 1], pts[:, 1] <= box[3])
    index = np.logical_and(x_in, y_in)
    return pts[index, :].copy()


class PseudoRand:

    def __init__(self, nbpt=100, datafile=None, initial='Random', boundary='Free', delta=0.1):
        self.box = np.array([0.0, 1.0, 0.0, 1.0])  # INTEREST Box
        self.initial = initial
        self.boundary = boundary
        self.delta = delta
        self.d = 1.0/np.sqrt(nbpt)
        self.nbpt = nbpt
        self.datafile = datafile
        self.dist_mean = np.array([])
        self.dist_sigma = np.array([])
        self.pts = np.zeros((nbpt, 2))
        self.site = np.zeros((nbpt, 2))
        self.site_new_x = np.zeros((nbpt, 2))
        self.site_new_y = np.zeros((nbpt, 2))

        self.s = np.zeros((nbpt, 2))
        self.set_initial()

    def set_nbpt(self, nbpt=100):
        self.nbpt = nbpt
        self.set_grid()

    def set_datafile(self, datafile=None):
        self.datafile = datafile

    def set_initial(self, initial=None):
        if initial:
            self.initial = initial

        if self.initial == 'Random':
            self.set_initial_random()
        elif self.initial == 'Square':
            self.set_initial_square()
        elif self.initial == 'Hexagonal_Compact':
            self.set_initial_hexagonal()
        elif self.initial == 'From_File':
            self.set_initial_from_file()
        else:
            self.set_initial_random()
        self.set_site()
        self.set_grid()

    def set_initial_random(self):
        self.pts = np.random.rand(self.nbpt, 2)

    def set_initial_from_file(self):
        self.pts = np.genfromtxt(self.datafile, delimiter=',')
        self.set_nbpt(len(self.pts))

    def set_initial_square(self):
        self.pts = np.empty([0, 2])
        x0 = self.d / 2
        y0 = self.d / 2
        nb = np.sqrt(self.nbpt).astype(int)

        for id_c in range(nb):
            for id_l in range(nb):
                self.pts = np.append(self.pts, [[id_c * self.d + x0, id_l * self.d + y0]], axis=0)

    def set_initial_hexagonal(self):
        self.pts = np.empty([0, 2])
        dx = self.d
        dy = np.sqrt(3) * self.d / 2
        nb_x = np.ceil(1.0 / dx).astype('int') - 1
        nb_y = np.ceil(1.0 / dy).astype('int') - 2
        x0 = dx / 2
        y0 = dy / 2

        for id_c in range(nb_x):
            for id_l in range(nb_y):
                x = (2 * id_c + id_l % 2) * x0 + x0/2
                y = id_l * dy + y0
                self.pts = np.append(self.pts, [[x, y]], axis=0)

    def set_boundary(self, boundary):
        self.boundary = boundary
        self.set_site()
        self.set_grid()

    def set_delta(self, delta=0.1):
        self.delta = delta
        self.set_site()

    def set_grid(self):
        delta = self.delta
        step_s = 20 * int(np.sqrt(self.nbpt))
        if self.boundary == 'Free':
            sy, sx = np.mgrid[0.0: 1.0: step_s * 1j, 0.0: 1.0: step_s * 1j]
        else:
            sy, sx = np.mgrid[-delta: 1.0 + delta: step_s * 1j, -delta: 1.0 + delta: step_s * 1j]
        self.s = np.c_[sx.ravel(), sy.ravel()]

    def set_site(self):
        delta = self.delta
        if self.boundary == 'Free':
            self.site = self.pts.copy()
        elif self.boundary == 'Mirror':
            top_left = [-1.0, 1.0] * in_box(self.pts, np.array([0.0, delta, 1.0 - delta, 1.0]))
            top_left[:, 1] = 2.0 - top_left[:, 1]
            top = in_box(self.pts, np.array([0.0, 1.0, 1.0 - delta, 1.0]))
            top[:, 1] = 2.0 - top[:, 1]
            top_right = [2.0, 2.0] + [-1.0, -1.0] * in_box(self.pts, np.array([1.0 - delta, 1.0, 1.0 - delta, 1.0]))
            left = [-1.0, 1.0] * in_box(self.pts, np.array([0.0, delta, 0.0, 1.0]))
            right = in_box(self.pts, np.array([1.0 - delta, 1.0, 0.0, 1.0]))
            right[:, 0] = 2.0 - right[:, 0]
            bottom_left = [-1.0, -1.0] * in_box(self.pts, np.array([0.0, delta, 0.0, delta]))
            bottom = [1.0, -1.0] * in_box(self.pts, np.array([0.0, 1.0, 0.0, delta]))
            bottom_right = [1.0, -1.0] * in_box(self.pts, np.array([1.0 - delta, 1.0, 0.0, delta]))
            bottom_right[:, 0] = 2.0 - bottom_right[:, 0]
            self.site = np.vstack((top_left, top, top_right,
                                   left, self.pts, right,
                                   bottom_left, bottom, bottom_right))
        elif self.boundary == 'Periodic':
            top_left = [-1.0, 1.0] + in_box(self.pts, np.array([1 - delta, 1.0, 0.0, delta]))
            top = [0.0, 1.0] + in_box(self.pts, np.array([0.0, 1.0, 0.0, delta]))
            top_right = [1.0, 1.0] + in_box(self.pts, np.array([0.0, delta, 0.0, delta]))
            left = [1.0, 0.0] + in_box(self.pts, np.array([0.0, delta, 0.0, 1.0]))
            right = [-1.0, 0.0] + in_box(self.pts, np.array([1.0 - delta, 1.0, 0.0, 1.0]))
            bottom_left = [-1.0, -1.0] + in_box(self.pts, np.array([1 - delta, 1.0, 1 - delta, 1.0]))
            bottom = [0.0, -1.0] + in_box(self.pts, np.array([0.0, 1.0, 1.0 - delta, 1.0]))
            bottom_right = [1.0, -1.0] + in_box(self.pts, np.array([0.0, delta, 1 - delta, 1.0]))
            self.site = np.vstack((top_left, top, top_right,
                                   right, self.pts, left,
                                   bottom_left, bottom, bottom_right))

    def get_centroid(self):
        """
            cKDTree() create a kd-tree of Voronoi sites
            cKDTree.query() get closest Voronoi site neighbour of each grid point
            bincount() calculate the histogram distribution
        """
        nbpt = len(self.site)
        tree = cKDTree(self.site)
        d, k = tree.query(self.s)
        m = np.bincount(k, minlength=nbpt)
        self.site_new_x = np.bincount(k, weights=self.s[:, 0], minlength=nbpt)
        self.site_new_y = np.bincount(k, weights=self.s[:, 1], minlength=nbpt)

        for i in range(0, nbpt):
            if m[i] > 0:
                self.site_new_x[i] = self.site_new_x[i] / float(m[i])
                self.site_new_y[i] = self.site_new_y[i] / float(m[i])

        tree = cKDTree(self.pts)
        d, k = tree.query(self.pts, 2)

        self.dist_mean = np.append(self.dist_mean, np.mean(d[:, 1]))
        self.dist_sigma = np.append(self.dist_sigma, np.std(d[:, 1]))

    def update(self):
        self.site[:, 0] = self.site_new_x
        self.site[:, 1] = self.site_new_y
        self.pts = in_box(self.site, self.box)
        self.set_site()

    def iteration(self):
        self.get_centroid()

    def save_pts(self, datafile):
        np.savetxt(datafile, self.pts, delimiter=',')
