
import numpy as np
from scipy.optimize import linear_sum_assignment as hungarian
import math
import random
import dionysus as dion


def fpd_cluster(data, c, verbose=False, max_iter=5, p=1, max_range=10, T=100):
    # Compute topological fuzzy clusters of a collection of point clouds
    #
    # inputs: list of datasets data, number of cluster centres c,
    #         verbose True/False, maximum iterations max_iter,
    #         dimension p-PH  p, Rips parameter max_range, hyper-parameter T
    #         (if unsure of value for max_range or T, set as the furthest distance between two points)
    # outputs: membership values r, cluster centres M

    diagrams = []
    for i in range(len(data)):
        filtration = dion.fill_rips(data[i], p+1, max_range)
        homology = dion.homology_persistence(filtration)
        diagram = dion.init_diagrams(homology, filtration)
        diagrams.append(diagram[p])

    diagrams = reformat_diagrams(diagrams, T)
    r, M = pd_fuzzy(diagrams, c, verbose, max_iter)

    return r, M


def calc_W(D, M):
    # calculate pairwise Wasserstein distances
    # inputs: array of diagrams D, array of centres M
    # returns: array W with w[j][k] = W_2(D_j, M_k)

    n = len(D)
    c = len(M)

    # W[j,k] = W_2(D_j, M_k)
    W = np.zeros((n, c))

    for j in range(n):
        for k in range(c):
            wass = calc_wasserstein(D[j], M[k])
            if wass != 0:
                W[j][k] = wass
            else:
                W[j][k] = 0.001
    return W


def calc_r(W):
    # calculate membership values
    # inputs: array of Wasserstein distances W
    # returns: array r of membership values r[j][k]

    n = np.shape(W)[0]
    c = np.shape(W)[1]
    r = np.zeros((n, c))

    for j in range(n):
        for k in range(c):
            sum = 0
            for l in range(c):
                sum += W[j][k] / W[j][l]
            r[j][k] = 1 / sum

    return r


def add_diagonals(D):
    # adds diagonals to diagrams so distance is well defined
    # diagonal denoted -1
    # inputs: diagrams D
    # returns: diagrams D

    n = len(D)
    # c = len(M)

    # find m = max number of off-diagonal points
    m = 0
    for j in range(n):
        if len(D[j]) > m:
            m = len(D[j])
    """for k in range(c):
        if len(M[k]) > m:
            m = len(M[k])"""

    # add diagonals so every diagram has m features
    for j in range(n):
        for i in range(m - len(D[j])):
            D[j].append([-1, -1])
    """for k in range(c):
        for i in range(m - len(M[k])):
            M[k].append([-1, -1])"""

    return D


def calc_cost_matrix(Dj, Mk):
    # calculates the cost matrix for optimal transport problem
    # inputs: diagram Dj, cluster centre Mk
    # returns: cost matrix c

    m = len(Dj)
    if m != len(Mk):
        exit("Incompatible diagram size in calc_cost_matrix: " + str(len(Dj)) + " and " + str(len(Mk)))

    c = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            # both off-diagonal
            if Dj[i][0] != -1 and Mk[j][0] != -1:
                c[i][j] = (Dj[i][0]-Mk[j][0])**2 + (Dj[i][1]-Mk[j][1])**2
            # only Dj[i] off-diagonal
            elif Dj[i][0] != -1 and Mk[j][0] == -1:
                c[i][j] = ((Dj[i][1] - Dj[i][0]) * 1/math.sqrt(2))**2
            # only Mk[j] off-diagonal
            elif Dj[i][0] == -1 and Mk[j][0] != -1:
                c[i][j] = ((Mk[j][1] - Mk[j][0]) * 1/math.sqrt(2))**2

    return c


def calc_wasserstein(Dj, Mk):
    # calculates the 2-Wasserstein L2 distance between two diagrams
    # inputs: diagram Dj, centre Mk
    # returns: W_2(Dj, Mk)

    m = len(Dj)
    c = calc_cost_matrix(Dj, Mk)
    X = hungarian(c)
    total = 0
    for i in range(m):
        total += c[X[0][i]][X[1][i]]
    return math.sqrt(total)


def calc_frechet_mean(D, r, k, verbose):
    # computes the weighted frechet mean of D with weights r[.][k]
    # inputs: diagrams D, membership values r, centre index k, verbose
    # returns: weighted frechet mean y, optimal pairings x

    n = len(D)
    m = len(D[0])
    # initialise to random diagram in D
    random.seed(0)
    M_update = D[random.randint(0, n-1)]

    # first run to find matching
    matching = []
    for j in range(n):
        c = calc_cost_matrix(M_update, D[j])
        x_indices = hungarian(c)
        matching.append(x_indices)

    # loop until stopping condition is found
    counter2 = 0

    while True:
        counter2 += 1

        # update matched points
        x = np.zeros((n, m, 2))
        for j in range(n):
            for i in range(m):
                index = matching[j][1][i]
                x[j][i] = D[j][index]

        # generate y to return
        y = np.zeros((m, 2))

        # loop over each point
        for i in range(m):
            # calculate w and w_\Delta
            r2_od = 0
            r2x_od = [0, 0]
            for j in range(n):
                if x[j][i][0] != -1:
                    r2_od += r[j][k]**2
                    r2x_od[0] += r[j][k]**2 * x[j][i][0]
                    r2x_od[1] += r[j][k]**2 * x[j][i][1]

            # if all points are diagonals
            if r2_od == 0:
                # then y[i] is a diagonal
                y[i] = [-1, -1]

            # else proceed
            else:
                w = [r2x_od[0]/r2_od, r2x_od[1]/r2_od]
                w_delta = [(w[0]+w[1])/2, (w[0]+w[1])/2]

                r2_d = 0
                r2_w_delta = [0, 0]
                for j in range(n):
                    if x[j][i][0] == -1:
                        r2_d += r[j][k] ** 2
                        r2_w_delta[0] += r[j][k]**2 * w_delta[0]
                        r2_w_delta[1] += r[j][k]**2 * w_delta[1]

                # calculate weighted mean
                y[i][0] = (r2x_od[0] + r2_w_delta[0]) / (r2_od + r2_d)
                y[i][1] = (r2x_od[1] + r2_w_delta[1]) / (r2_od + r2_d)

        old_matching = matching.copy()
        matching = []
        for j in range(n):
            c = calc_cost_matrix(y, D[j])
            x_indices = hungarian(c)
            matching.append(x_indices)

        comparison = (np.array(matching) == np.array(old_matching))
        if comparison.all():
            if verbose:
                print("      Frechet iterations for M_" + str(k) + ": " + str(counter2))
            return y, x


def init_clusters(D, c):
    # initialise cluster centres to Frechet mean of two diagrams
    # inputs: diagrams D, number of clusters c
    # outputs: initialised cluster cenres M

    M = []
    ones = np.ones((len(D)+1, c+1))
    for i in range(c):
        diagram, _ = calc_frechet_mean([D[i], D[i+1]], ones, i, verbose=False)
        M.append(diagram)

    return M


def J(r, D, M):
    # computes the cost function J
    # inputs: membership values r, diagrams D, cluster centres M
    # returns: clustering cost

    W = calc_W(D, M)

    n = np.shape(W)[0]
    c = np.shape(W)[1]

    sum = 0
    for j in range(n):
        for k in range(c):
            sum += r[j][k]**2 * W[j][k]**2

    return sum


def reformat_diagrams(D, T=100):
    # reformat from dionysus to custom structure
    # replace d=inf with hyper-parameter T
    # inputs: dionysus diagrams D, hyper-parameter T
    # returns: reformatted diagrams D

    D_new = []
    for i in range(len(D)):
        D_temp = []
        for p in D[i]:
            if np.isinf(p.death):
                D_temp.append([p.birth, T])
            else:
                D_temp.append([p.birth, p.death])
        D_new.append(D_temp)

    return D_new


def pd_fuzzy(D, c, verbose=False, max_iter=5):
    # Performs Fuzzy c-Means Clustering on Persistence Diagrams
    #
    # inputs: persistence diagrams D, number of cluster centres c, verbose, maximum iterations max_iter
    # outputs: membership values r, cluster centres M

    D = add_diagonals(D)
    M = init_clusters(D, c)

    n = len(D)
    m = len(D[0])

    counter = 0
    while counter < max_iter:
        if verbose:
            print("Fuzzy iteration: " + str(counter))
        counter += 1

        # update membership values
        W = calc_W(D, M)
        r = calc_r(W)

        if verbose:
            J_temp = J(r, D, M)
            print(" -- Update r -- ")
            print("   J(r, M) = " + str(J_temp))

            print(" -- Update M -- ")

        x = np.zeros((c, n, m, 2))
        # update cluster centres
        for k in range(c):
            M[k], x[k] = calc_frechet_mean(D, r, k, verbose)

        # compute J
        # J_prev = J_new
        # J_new = J(r, D, M)
        if verbose:
            J_new = J(r, D, M)
            print("   J(r, M) = " + str(J_new) + "\n")

    return r, M
