
import numpy as np
import matplotlib.pyplot as plt
import timeit
from sklearn import datasets
import random
import clustering
import dionysus as dion


def gen_data(n_samples, seed=0, noise=0.05):

    # print("\nGenerating data...\n")

    np.random.seed(seed)
    random.seed(seed)

    data = []
    data.append(np.random.normal(size=(n_samples, 2), scale=0.5))
    data.append(0.9 * np.random.normal(size=(n_samples, 2), scale=0.5))
    data.append(datasets.make_circles(n_samples=n_samples, factor=0.99, noise=noise, random_state=seed)[0])
    data.append(datasets.make_circles(n_samples=n_samples, factor=0.99, noise=noise, random_state=seed + 1)[0])

    return data


def save_data(data, i):
    file = "../data/timing/data" + str(i) + ".d2"
    f = open(file, 'a')

    for k in range(len(data)):
        f.write("2\n")
        f.write(str(len(data[k])) + "\n")
        weight = 1/len(data[k])

        line = (str(weight) + " ") * len(data[k])
        f.write(line + "\n")
        for j in range(len(data[k])):
            f.write(str(data[k][j][0]) + " " + str(data[k][j][1]) + "\n")


def compute_diagrams(data):

    diagrams = []
    for i in range(len(data)):
        # print("Processing data: " + str(i))
        filtration = dion.fill_rips(data[i], 2, 3.0)
        homology = dion.homology_persistence(filtration)
        diagram = dion.init_diagrams(homology, filtration)
        diagrams.append(diagram[1])
    # print()

    return diagrams


def time_fcm(num_points, c=2, max_iter=5, num_repeats=5):

    data = gen_data(num_points)
    # save_data(data, num_points)

    def packaged():
        diagrams = compute_diagrams(data)
        diagrams = clustering.reformat_diagrams(diagrams)
        clustering.pd_fuzzy(diagrams, 2, verbose=False, max_iter=max_iter)

    return round(timeit.timeit(packaged, number=num_repeats) * (1000 / (max_iter * num_repeats)), 4)


def check_timing(number_points, max_iter=5, num_repeats=5, write_file=False):

    for i in range(len(number_points)):
        time = time_fcm(number_points[i], max_iter, num_repeats)
        print('Total points: ' + str(int(4*number_points[i])) + ", time (ms): " + str(time))

        if write_file:
            with open('../data/timing_data/timings/pdcluster.csv', 'a') as the_file:
                the_file.write(str(number_points[i]) + ", " + str(time) + '\n')


if __name__ == '__main__':
    # need to generate datasets
    # 2 noise, 2 rings
    # increase number of points in each one
    # time how long it takes to cluster

    # number of points per dist from 25, 50, 100, 125, ..., 250.. giving total points 100 -> 1000

    number_points = np.arange(25, 275, 25)
    check_timing(number_points, write_file=False)
