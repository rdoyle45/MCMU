import heapq
from operator import itemgetter
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
from os import makedirs
from os.path import isdir, dirname

helium_names = ['Elastic',
                'He \\rightarrow \mathrm{2^3S_1}',
                'He \\rightarrow \mathrm{2^1S_0}',
                'He \\rightarrow \mathrm{2^3P_{0,1}}',
                'He \\rightarrow \mathrm{2^1P_1}',
                'He \\rightarrow \mathrm{3^3S_1}',
                'He \\rightarrow \mathrm{3^1S_0}',
                'He \\rightarrow \mathrm{3^3P_{0,1}}',
                'He \\rightarrow \mathrm{3^3D_{1,3}}',
                'He \\rightarrow \mathrm{3^1D_2}',
                'He \\rightarrow \mathrm{3^1P_1}',
                'He \\rightarrow \mathrm{4^3S_1}',
                'He \\rightarrow \mathrm{4^1S_0}',
                'He \\rightarrow \mathrm{4^3P_{0,1}}',
                'He \\rightarrow \mathrm{4^1D_2}',
                'He \\rightarrow \mathrm{4^3D_{1,3}}',
                'He \\rightarrow \mathrm{4^1F_3}',
                'He \\rightarrow \mathrm{4^3F_{2,3,4}}',
                'He \\rightarrow \mathrm{4^1P_1}',
                'He \\rightarrow \mathrm{5^3S_1}',
                'He \\rightarrow \mathrm{5^1S_0}',
                'He \\rightarrow \mathrm{5^3P_{0,1}}',
                'He \\rightarrow \mathrm{51^D_2}',
                'He \\rightarrow \mathrm{5^3D_{1,3}}',
                'He \\rightarrow \mathrm{5^1S_3}',
                'He \\rightarrow \mathrm{5^3F_{2,3,4}}',
                'He \\rightarrow \mathrm{5^1P_1}',
                'He \\rightarrow \mathrm{6^3S_1}',
                'He \\rightarrow \mathrm{6^1S_0}',
                'He \\rightarrow \mathrm{6^3P_{0,1}}',
                'He \\rightarrow \mathrm{6^3D_{1,3}}',
                'He \\rightarrow \mathrm{6^1D_2}',
                'He \\rightarrow \mathrm{6^1P_1}',
                'He \\rightarrow \mathrm{7^3S_1}',
                'He \\rightarrow \mathrm{7^1S_0}',
                'He \\rightarrow \mathrm{7^3P_{1,0}}',
                'He \\rightarrow \mathrm{7^3D_{1,3}}',
                'He \\rightarrow \mathrm{7^1D_2}',
                'He \\rightarrow \mathrm{7^1P_1}',
                'He \\rightarrow \mathrm{N3S Sum}',
                'He \\rightarrow \mathrm{N1S Sum}',
                'He \\rightarrow \mathrm{N3P Sum}',
                'He \\rightarrow \mathrm{N1D Sum}',
                'He \\rightarrow \mathrm{N3D Sum}',
                'He \\rightarrow \mathrm{8^1P_1}',
                'He \\rightarrow \mathrm{9^1P_1}',
                'He \\rightarrow \mathrm{10^1P_1}',
                'He \\rightarrow \mathrm{11^1P_1}',
                'He \\rightarrow \mathrm{12^1P_1}',
                'He \\rightarrow \mathrm{N1P Sum}',
                'He \\rightarrow \mathrm{He+}'
                ]

nitrogen_names = ['Effective',
                  'Rotational',
                  '\mathrm{N_2} \\rightarrow Vib. 0.29 eV',
                  '\mathrm{N_2} \\rightarrow Vib. 0.291 eV',
                  '\mathrm{N_2} \\rightarrow Vib. 0.59 eV',
                  '\mathrm{N_2} \\rightarrow Vib. 0.88 eV',
                  '\mathrm{N_2} \\rightarrow Vib. 1.17 eV',
                  '\mathrm{N_2} \\rightarrow Vib. 1.47 eV',
                  '\mathrm{N_2} \\rightarrow Vib. 1.76 eV',
                  '\mathrm{N_2} \\rightarrow Vib. 2.06 eV',
                  '\mathrm{N_2} \\rightarrow Vib. 2.35 eV',
                  '\mathrm{N_2} \\rightarrow \mathrm{A ^3\\Sigma} 6.17 eV',
                  '\mathrm{N_2} \\rightarrow \mathrm{A ^3\\Sigma} 7 eV',
                  '\mathrm{N_2} \\rightarrow \mathrm{B ^3\\Pi}',
                  '\mathrm{N_2} \\rightarrow \mathrm{w ^3\\Delta}',
                  '\mathrm{N_2} \\rightarrow \mathrm{A ^3\\Sigma} 7.8 eV',
                  '\mathrm{N_2} \\rightarrow \mathrm{B\' ^3\\Sigma}',
                  '\mathrm{N_2} \\rightarrow \mathrm{a\' ^1\\Sigma}',
                  '\mathrm{N_2} \\rightarrow \mathrm{a ^1\\Pi}',
                  '\mathrm{N_2} \\rightarrow \mathrm{w ^1\\Delta}',
                  '\mathrm{N_2} \\rightarrow \mathrm{C ^3\\Delta}',
                  '\mathrm{N_2} \\rightarrow \mathrm{E ^3\\Sigma}',
                  '\mathrm{N_2} \\rightarrow \mathrm{a\" ^1\\Sigma}',
                  '\mathrm{N_2} \\rightarrow \mathrm{Ex. Sum}',
                  '\mathrm{N_2} \\rightarrow \mathrm{N_2+} 15.6 eV',
                  '\mathrm{N_2} \\rightarrow \mathrm{N_2+ B^2\\Sigma}'
                ]

tex_names = {'He':helium_names, 'N2':nitrogen_names}


def plot_barh(data, title, indexes, filename, cross_sections):

    pos = np.arange(len(data)) + 0.5

    colour = []
    bar_label = []

    for item in data:
        if item > 0:
            colour.append('g')
        else:
            colour.append('r')

    for num in indexes:
        bar_label.append(r'$' + cross_sections[num] + '$')

    fig = plt.figure()
    ax = fig.add_subplot(111)
    #ax.xaxis.set_major_locator(MaxNLocator(symmetric=True))
    if data:
        data_max = max(data)
        ax.set_xlim(xmax=data_max*1.01)
        ax.set_xticks(np.linspace(0, data_max*1.01, 7))
    ax.barh(pos, data, color=colour)
    plt.title(title)
    ax.set_yticks(pos + 0.0)
    ax.set_yticklabels(bar_label)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(12)
    fig.tight_layout()
    plt.savefig(filename)
    plt.close()
    return


def extract_data(data_line, num_bars):

    data = data_line.split()

    data = [float(item) for item in data]
    E_over_N = data.pop(0)

    if E_over_N < 1:
        E_over_N = str(E_over_N).replace(".", "_")
    else:
        E_over_N = str(int(E_over_N))

    index_list = list(map(itemgetter(0), heapq.nlargest(num_bars, enumerate(data), key=itemgetter(1))))

    top_data = []
    usuable_indexes = []
    for i in index_list:
        if data[i] != 0:
            top_data.append(data[i])
            usuable_indexes.append(i)

    return E_over_N, usuable_indexes, top_data


def create_barchart(stats_file, species, coef_index, filename, num_bars):

    data_dir = dirname(stats_file)
    graph_dir = data_dir + "/Bar_charts"
    if not isdir(graph_dir):
        makedirs(graph_dir)
    filename = graph_dir + "/" + filename

    data_file = open(stats_file, 'r')
    lines = data_file.readlines()

    if coef_index != 0:
        coef_sort = coef_index
    else:
        coef_sort = coef_index

    data_file.close()

    coef_data = lines[coef_sort*21:(coef_sort+1)*21]  # Extract coef data

    cross_sections = tex_names[species]

    if coef_index in [0, 1, 2]:
        title = coef_data[0].strip()    # Extract coefs name
    elif coef_index == 3:
        title = r'$' + cross_sections[0] + '$' + " Collision"
    else:
        title = "Rate Coef. for " + r'$' + cross_sections[coef_index-3] + '$'

    for data_line in coef_data[4:11]:

        E_over_N, index_list, top_data = extract_data(data_line, num_bars)
        print(E_over_N)
        new_filename = filename + "_" + E_over_N

        plot_barh(top_data, title, index_list, new_filename, cross_sections)

    return








