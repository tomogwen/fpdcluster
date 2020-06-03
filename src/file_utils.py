
# file utils used for processing data

def xyz_to_csv():
    # used to convert xyz files exported from Vesta to csv files

    file_read = '../data/iron_structures/Fe_mp-150_conventional_standard2.xyz'
    file_write = '../data/iron_structures/fcc.csv'
    atom_to_replace = 'Fe'

    f_r = open(file_read, 'r')
    lines = f_r.readlines()
    lines = lines[2:]

    f_w = open(file_write, 'a')

    for i in range(len(lines)):
        print(lines[i])
        line = lines[i]
        line = line.replace('\n', '')
        line = line.replace(' ' + str(atom_to_replace) + '    ', '')
        line = line.replace('    ', ',')

        f_w.write(line)
        f_w.write('\n')


def save_as_d2(data, i):
    # saves data as d2 format used by https://github.com/bobye/WBC_Matlab

    file = "../data/iron_structures/transformed_data/data" + str(i) + ".d2"
    f = open(file, 'a')

    for k in range(len(data)):
        f.write("3\n")
        f.write(str(len(data[k])) + "\n")
        weight = 1/len(data[k])

        line = (str(weight) + " ") * len(data[k])
        f.write(line + "\n")
        for j in range(len(data[k])):
            f.write(str(data[k][j][0]) + " " + str(data[k][j][1]) + " " + str(data[k][j][2]) + "\n")
