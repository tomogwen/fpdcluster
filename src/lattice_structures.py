
import clustering
import numpy as np
import matplotlib.pyplot as plt
import dionysus as dion
import file_utils


def plot_diagram(dgm, ax=False, show=False, labels=False, line_style=None, pt_style=None, lims=False):

    # taken from Dionysus2 package

    line_kwargs = {}
    pt_kwargs = {}
    if pt_style is not None:
        pt_kwargs.update(pt_style)
    if line_style is not None:
        line_kwargs.update(line_style)

    inf = float('inf')

    if ax==False:
        ax = plt.axes()
        ax.set_aspect('equal', 'datalim')

    if lims==False:
        min_birth = min(p.birth for p in dgm if p.birth != inf)
        max_birth = max(p.birth for p in dgm if p.birth != inf)
        min_death = min(p.death for p in dgm if p.death != inf)
        max_death = max(p.death for p in dgm if p.death != inf)

    else:
        min_birth = lims[0]
        max_birth = lims[1]
        min_death = lims[2]
        max_death = lims[3]

    ax.set_aspect('equal', 'datalim')

    min_diag = min(min_birth, min_death)
    max_diag = max(max_birth, max_death)
    ax.scatter([p.birth for p in dgm], [p.death for p in dgm], **pt_kwargs)
    ax.plot([min_diag, max_diag], [min_diag, max_diag], **line_kwargs)

    ax.set_xlabel('birth')
    ax.set_ylabel('death')


def plot_all(data, diagrams):

    fig = plt.figure(figsize=(20, 10))

    for i in range(len(data)):
        num = 241 + i

        ax = plt.subplot(num)
        plt.scatter(data[i][:, 0], data[i][:, 1])

        ax = plt.subplot(num + 4)
        plot_diagram(diagrams[i], ax, lims=[0, 1.5, 0, 1.75])

    fig.suptitle("Datasets with corresponding persistence diagrams")
    plt.show()


def compute_diagrams(data, k_filt=3, to_return=2):

    diagrams = []
    for i in range(len(data)):
        # print("Processing data: " + str(i) + "\n")
        filtration = dion.fill_rips(data[i], k_filt, 3.0)
        homology = dion.homology_persistence(filtration)
        diagram = dion.init_diagrams(homology, filtration)

        diagrams.append(diagram[to_return])

    return diagrams


def plot_clusters(M):
    plt.scatter(M[0].T[0], M[0].T[1], c='r', label='Rings')
    plt.scatter(M[1].T[0], M[1].T[1], c='b', label='Noise')
    plt.xlim([0, 1.5])
    plt.ylim([0, 1.75])
    plt.plot([0.1, 1.2], [0.1, 1.2])
    plt.legend()
    plt.title("Persistence Diagram Cluster Centres")

    plt.show()


def add_noise(data, noise):
    for i in range(len(data)):
        for j in range(len(data[i])):
            data[i][j][0] += noise * np.random.normal()
            data[i][j][1] += noise * np.random.normal()
            data[i][j][2] += noise * np.random.normal()

    return data


def reflect(data, axis=0):
    datax = data.T[0]
    datay = data.T[1]
    dataz = data.T[2]

    if axis == 0:
        data = np.array([datay, datax, dataz]).T
    elif axis == 1:
        data = np.array([datax, dataz, datay]).T

    return data


def translate(data, amount):
    return data + amount


def rotate(data, degrees, axis=0):
    theta = np.radians(degrees)
    c, s = np.cos(theta), np.sin(theta)
    if axis == 0:
        R = np.array(((1, 0, 0), (0, c, -s), (0, s, c)))
    elif axis == 1:
        R = np.array(((c, 0, s), (0, 1, 0), (-s, 0, c)))
    else:
        R = np.array(((c, -s, 0), (s, c, 0), (0, 0, 1)))

    return np.dot(data, R.T)


def apply_transformations(data1, data2, transforms):
    data = np.array([data2, data1, data1, data1, data2, data2])

    for i in range(len(transforms)):
        if transforms[i] == 0:
            data[1] = reflect(data[1], axis=0)
            data[2] = reflect(data[2], axis=1)
            data[4] = reflect(data[4], axis=0)
            data[5] = reflect(data[5], axis=1)

        if transforms[i] == 1:
            data[1] = rotate(data[1], degrees=180, axis=0)
            data[2] = rotate(data[2], degrees=180, axis=1)
            data[4] = rotate(data[4], degrees=180, axis=0)
            data[5] = rotate(data[5], degrees=180, axis=1)

        if transforms[i] == 2:
            data[1] = translate(data[1], 5)
            data[2] = translate(data[2], -5)
            data[4] = translate(data[4], 5)
            data[5] = translate(data[5], -5)

    return data


if __name__ == '__main__':

    for i in range(2):

        seed = 0
        data = []
        np.random.seed(0)

        if i == 0:
            data1 = np.genfromtxt('data/cubic_structures/bcc/bcc.csv', delimiter=',')
            data2 = np.genfromtxt('data/cubic_structures/fcc/fcc.csv', delimiter=',')
            print("\n--- CUBIC STRUCTURES DATA ---\n")
        else:
            data1 = np.genfromtxt('data/carbon_allotropes/dia/dia.csv', delimiter=',')
            data2 = np.genfromtxt('data/carbon_allotropes/ths3/ths3.csv', delimiter=',')
            print("\n--- CARBON ALLOTROPES DATA ---\n")

        # transformations key
        # 0 - reflection
        # 1 - rotation
        # 2 - translation
        transformations = [[], [0], [1], [2]]

        for i in range(len(transformations)):
            print("Case: " + str(i) + ", Transformations: " + str(transformations[i]))
            data = apply_transformations(data1, data2, transformations[i])

            # file_utils.save_as_d2(data, i)

            diagrams = compute_diagrams(data, k_filt=3, to_return=2)
            diagrams_cluster = clustering.reformat_diagrams(diagrams)

            r, M = clustering.pd_fuzzy(diagrams_cluster, 2, verbose=False, max_iter=5)
            # rearrange r to match data order for interpretation
            r[[0, 3], :] = r[[3, 0], :]

            if i == 0:
                t = 'no transformation'
            elif i == 1:
                t = 'reflection'
            elif i == 2:
                t = 'rotation'
            elif i == 3:
                t = 'translation'
            else:
                t = ''
            print("Membership values for " + t + ":")
            print(r)


            c1 = r[0][0] > 0.5 and r[1][0] > 0.5 and r[2][0] > 0.5 and r[3][1] > 0.5 and r[4][1] > 0.5 and r[5][1] > 0.5
            c2 = r[0][1] > 0.5 and r[1][1] > 0.5 and r[2][1] > 0.5 and r[3][0] > 0.5 and r[4][0] > 0.5 and r[5][0] > 0.5

            if c1 or c2:
                print("Clustered Succesfully\n")
            else:
                print("Clustered Unsuccesfully\n")
