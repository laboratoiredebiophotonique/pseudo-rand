# -*- coding: utf-8 -*-
import numpy as np


class PairCorrelation:
    def __init__(self, datafile, L=1):
        self.L = L
        self.ptS = np.genfromtxt(datafile, delimiter=',')
        self.ptS -= [0.5, 0.5]

        r_in = 0.2
        r_out = 0.5
        r_min = 1e-6
        r_max = r_out - r_in
        dr = 1 / (20 * np.sqrt(len(self.ptS)))  # same as s_step in CVT.py
        self.radius = np.arange(r_min, r_max + dr, dr)

        d2 = np.sqrt(self.ptS[:, 0] ** 2 + self.ptS[:, 1] ** 2)

        self.o_pt = self.ptS[np.where(d2 < r_out)]  # pts to evaluate
        self.i_pt = self.ptS[np.where(d2 < r_in)]  # pts for measurement

        distance = np.asarray(
            [np.sqrt((pt - self.o_pt)[:, 0] ** 2 + (pt - self.o_pt)[:, 1] ** 2) for pt in self.i_pt]).flatten()

        distance = np.sort(distance)
        distance = distance[np.where(distance < self.radius[-1] + dr)]

        self.y = np.float64(
            np.asanyarray([len(np.where(distance < r + dr)[0]) - len(np.where(distance < r)[0]) for r in self.radius]))

        self.ptS *= self.L
        self.o_pt *= self.L
        self.i_pt *= self.L
        self.radius *= self.L
        dr *= self.L

        self.y /= 2 * np.pi * self.radius * dr * len(distance)
