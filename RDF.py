# -*- coding: utf-8 -*-
import numpy as np
from scipy.spatial import cKDTree


def rdf_2d(pts, nb_bin=100):
    r_in = 0.2
    r_out = 0.5
    r_min = 1e-6
    r_max = r_out - r_in
    dr = 0.01  # 1.0 / (20*np.sqrt(nb_bin))
    print(20 * np.sqrt(nb_bin))
    radius = np.arange(r_min, r_max + dr, dr)
    distance = np.sqrt(pts[:, 0] ** 2 + pts[:, 1] ** 2)
    pts_in = pts[np.where(distance <= r_in)]  # pts for measurement
    pts_out = pts[np.where(distance <= r_out)]  # pts to evaluate
    nb_pts_in = pts_in.shape[0]
    density = pts.shape[0]

    # Compute pairwise correlation for each inner particle pts_in
    pairwise = np.asarray([np.sqrt((pt - pts_out)[:, 0] ** 2 + (pt - pts_out)[:, 1] ** 2) for pt in pts_in]).flatten()
    pairwise = pairwise[np.where(pairwise <= radius[-1])]
    g, bins = np.histogram(pairwise, bins=radius)

    normalization = 2 * np.pi * radius[:-1] * dr * nb_pts_in * density

    rdf = g / normalization
    rdf[0] = 1.0
    return rdf, radius[:-1]


def rdf2d(particles, dr, rho=None, eps=1e-15):
    particles = particles - np.min(particles, axis=0)
    min_x, min_y = np.min(particles, axis=0)
    max_x, max_y = np.max(particles, axis=0)

    # dimensions of box
    w, h = (max_x - min_x), (max_y - min_y)

    r_max = (np.min([w, h]) / 2) * 0.8
    radii = np.arange(dr, r_max, dr)
    g_r = np.zeros(shape=(len(radii)))

    nb_pts = len(particles)
    if not rho:
        rho = nb_pts / (w * h)  # number density

    # create a KDTree for fast nearest-neighbor lookup of particles
    tree = cKDTree(particles)

    for r_idx, r in enumerate(radii):
        # find all particles that are at least r + dr away from the edges of the box
        valid_id = (particles[:, 0] - (r + dr) >= min_x) & (particles[:, 0] + (r + dr) <= max_x) \
                     & (particles[:, 1] - (r + dr) >= min_y) & (particles[:, 1] + (r + dr) <= max_y)
        valid_particles = particles[valid_id]

        # compute n_i(r) for valid particles.
        for particle in valid_particles:
            n = tree.query_ball_point(particle, r + dr - eps, return_length=True) \
                - tree.query_ball_point(particle, r, return_length=True)
            g_r[r_idx] += n

        # normalize
        n_valid = len(valid_particles)
        shell_vol = np.pi * ((r + dr) ** 2 - r ** 2)
        g_r[r_idx] /= n_valid * shell_vol * rho

    return radii, g_r
