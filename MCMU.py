from os import getcwd, chdir, makedirs, listdir, remove
from os.path import join, exists
from numpy import random as rd, empty, format_float_scientific as f_f, zeros, nonzero, divide, sqrt, seterr
from math import log, ceil
from tempfile import NamedTemporaryFile as NTF
from subprocess import Popen, DEVNULL
from matplotlib import pyplot as plt
from shutil import rmtree
import time
import logging
from random import randint, shuffle
from scipy import stats
import glob
from concurrent import futures
from multiprocessing import cpu_count
from sys import argv
from progress.bar import Bar

__author__ = "Rhys Doyle"
__copyright__ = "Copyright 2018, LTP Uncertainty Calculator"
__credits__ = ["Prof. Miles Turner", "Dublin City University"]
__maintainer__ = "Rhys Doyle"
__email__ = "rhys.doyle45@mail.dcu.ie"
__status__ = "Testing"


def run_program():

    with open(str(argv[1]), 'r') as f:

        lines = f.readlines()

    program = lines[9].split()[2]
    subject = lines[10].split()[2]
    data = lines[11].split()[2]
    uncert_set = lines[12].split()[2::]
    runs = int(lines[13].split()[2])
    species = lines[14].split()[2]
    num_EN_values = int(lines[15].split()[2])
    bolsig = lines[16].split()[2]

    try:
        num_cpus = int(lines[17].split()[2])
    except ValueError:
        num_cpus = cpu_count()

    p_values = int(lines[20].split()[2])
    if program == '1':
        Monte_carlo(subject, data, uncert_set, runs, species, num_EN_values, bolsig, num_cpus)
    elif program == '2':
        Morris(subject, data, uncert_set, runs, species, num_EN_values, bolsig, num_cpus, p_values)
    else:
        print("A valid program code was not written in your input file. Please check and try again.")

    return

"""

Monte-Carlo Simulation Code:

*** This code was created to determine the possible uncertainties associated with plasma transport coefficients
    using the experimental uncertainties associated with the cross-sections used ***

"""


# Run function for Monte-Carlo simulations
def Monte_carlo(subject, data, uncert_set, runs, species, num_EN_values, bolsig, num_cpus):

    start = time.time()

    og_cwd = getcwd()

    cwd = changedir(og_cwd, subject)

    cwd2 = join(cwd, 'Monte-Carlo')

    # Checking if directory exists, if yes remove and replace with clean directory
    if exists(cwd2):

        rmtree(cwd2)
        # Making and navigating to a new directory to store the new files
        makedirs(cwd2)

    elif not exists(cwd2):

        makedirs(cwd2)

    logging.basicConfig(filename=join(cwd2, "program.log"), filemode='w', level=logging.DEBUG,
                        format='%(asctime)s %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p')

    linelist = data_set(data)

    startsec, dashsec  = search_set(linelist)

    com_format(startsec, dashsec, linelist, data)

    linelist = data_set(data)

    startsec, dashsec  = search_set(linelist)

    commsec = comment_sec(startsec, dashsec, linelist)

    lognorm_set = lognorm(uncert_set, runs)

    cs_data = find_num(startsec, dashsec, linelist)

    inputdir = data_edit(cs_data, lognorm_set, data, cwd, dashsec, runs, commsec, cwd2)

    logging.info('Data files edited')

    outputdir = bolsig_file(inputdir, cwd2, og_cwd, species, num_EN_values)

    logging.info('BOLSIG run files created')
    logging.info('BOLSIG- running')

    bolsig_minus(cwd2, bolsig, num_cpus)

    logging.info('BOLSIG- completed')

    red_e_field, gen_mat, e_loss_dict, rate_dict, name_rate, name_energy, titlelist, num_e = output_file_read(outputdir, startsec, runs, num_EN_values)

    stats = gen_stats(gen_mat, runs)

    rate_stats = rate_cof(rate_dict, name_rate, runs, num_e)

    e_loss_stats = energy_loss(e_loss_dict, name_energy, runs, num_e)

    logging.info('Stats Generated')

    stats_file(red_e_field, stats, rate_stats, e_loss_stats, name_rate, name_energy, titlelist)

    logging.info('Stats file created')

#    graphing_data(red_e_field, stats, rate_stats, e_loss_stats, name_rate, name_energy)

    end = time.time()

    print(end - start)


# Changing directory for data sets
def changedir(path, subject):

    # Changing main directory to above location
    chdir(join(path, subject))
    # Displaying new main directory
    cwd = getcwd()

    return cwd


# Importing data file to program
def data_set(data):
    # Opening data file in a Read/Write capacity
    with open(data, "r+", encoding='utf8') as datafile:

        # split the file content into lines and save a list
        linelist = datafile.read().splitlines()

    return linelist


# Searching Data set for locations of COMMENT lines of different Cross Sections
def search_set(linelist):

    # Section words
    items = ['EXCITATION', 'ELASTIC', 'IONIZATION', 'EFFECTIVE', 'ROTATIONAL', 'ATTACHMENT']

    # Find how many lines in the data set
    x = len(linelist)
    # Initial value
    y = 0

    # Lists to save the locations of the start of each data set and their corresponding COMMENT section
    startsec = []
    dashsec = []

    while y <= x:

        # Calling each line individually
        line = str(linelist[y:y + 1])

        # Looking for the beginning of a data set
        if any(word in line for word in items):
            startsec.append(y)
            y += 1
        # Looking for the Column section
        elif line.find('---------------') != -1:
            dashsec.append(y)
            y += 1
        else:
            y += 1

    return startsec, dashsec


# FUnction to search for Comment sections and add those that do not have them
def com_format(startsec, dashsec, linelist, data):
    # Variable to store the location of lines containing param.:
    paramline = []
    p_location = 0

    # Loop to check is the data set has comment sections
    for i in range(0, (len(startsec))):
        x = startsec[i]
        y = dashsec[i * 2]

        while x <= y:
            line = str(linelist[x:x + 1])

            # checking line for the word comment
            if line.find('COMMENT:') != -1:
                x += 1
                break
            # If comment isn't present checking for param.:
            elif line.find('PARAM.:') != -1:
                p_location = x
                x += 1
            # Continuation of loop
            else:
                x += 1

            # If no COMMENT: was found this will store the location of PARAM.:
            if x == y:
                paramline.append(p_location)
            else:
                continue

    # Writing changes to the existing file
    out_file = []
    for i in range(0, (len(linelist))):

        # This will add a comment line after param.: to any data set without one
        if i in paramline:
            out_file.append(linelist[i] + '\nCOMMENT:\n')
        else:
            out_file.append(linelist[i] + '\n')

    with open(data, 'w', encoding='utf-8') as f:
        f.writelines(out_file)

    return

# Function to collect all comment sections
def comment_sec(startsec, dashsec, linelist):
    # Array to store comment line locations
    commsec = []

    # Loop to find and store the comment lines
    for i in range(0, (len(startsec))):

        x = startsec[i]
        y = dashsec[i * 2]

        while x <= y:
            line = str(linelist[x:x + 1])

            # Looking for comment in the selected line
            if line.find('COMMENT:') != -1:
                # If the word comment was found storing the line number
                commsec.append(x)
                x += 1
                break
            else:
                x += 1

    return commsec


# Generating a set amount of numbers from a lognormal distribution of specified parameters
def lognorm(uncert_set, num_values):
    # Distribution Parameters ( MEan and Variance)
    mean = 1
    lognorm_set = []

    for i in range(len(uncert_set)):
        var = (float(uncert_set[i]))**2

        # Calcualting std from a specified formula
        sigma = sqrt(log(1 + (var / (mean ** 2))))
        # Calcualting mean from a specified formula
        mu = log(mean ** 2 / sqrt(var + (mean ** 2)))

        # Gemerating the random numbers
        uncert = rd.lognormal(mu, sigma, num_values)
        lognorm_set.append(uncert)

    return lognorm_set


# Finding and extracting each individual data set
def find_num(startsec, dashsec, linelist):
    # Dictionary to collect the data for each Cross Section
    cs_data = {}

    progressbar = Bar('Extracting Data ', max=len(startsec), suffix='%(percent)d%%')

    # Iteratively adding to the above dictionary
    for i in range(len(startsec)):
        cs_name = linelist[startsec[i] + 1]
        cs_start = dashsec[i * 2] + 1
        cs_end = dashsec[(i * 2) + 1]

        cs_set = linelist[cs_start:cs_end]

        cs_data[(str(cs_name) + ", No." + str(i))] = cs_set

        progressbar.next()
    progressbar.finish()

    return cs_data


# Using the generated random numbers to create an equal amount of perturbed data files
def data_edit(cs_data, lognorm_set, oldfile, cwd, dashsec, num_uncert, comsec, cwd2):
    # Importing the list of cross-section data from the dictionary created earlier
    cross_values = list(cs_data.values())
    # Creating a second variable of the same list to reset the first when altered
    reset = cross_values.copy()

    chdir(cwd2)

    cwd3 = join(cwd2, 'Input_Files')
    makedirs(cwd3)

    progressbar = Bar('Generating Datasets ', max=num_uncert, suffix='%(percent)d%%')

    # For loop for each random number
    for j in range(num_uncert):

        # Reseting the variables
        cross_values = reset.copy()

        # Opening the currently used data file in order to read write equivalent files
        with open(join(cwd, oldfile), 'r', encoding='utf8') as old_data:
            # read lines into a list
            oldlinelist = old_data.read().splitlines()

        # Iterating over each set of Cross section data
        for i in range(len(cross_values)):
            # Selecting 1 set of data at a time
            dataset = cross_values[i].copy()

            # Iterating over each cross section individually
            for k in range(len(dataset)):
                # Splitting the Energy and Cross-Section numbers
                E, cs = dataset[k].split('\t')

                # Perturbing the Cross section number by the jth generated number
                cs = (float(cs) * lognorm_set[i][j])
                # Creating a new data line with the new cross section number
                newdata = (E + '\t' + str(cs))

                # Updating the list
                dataset[k] = newdata

            # Updating the entire data set
            cross_values[i] = dataset.copy()

            # Adding the newly perturbed data lines to the main list
            for x, y in zip(range((dashsec[i * 2] + 1), (dashsec[(i * 2) + 1])), range(len(dataset))):

                oldlinelist[x] = dataset[y]

            oldlinelist[comsec[i]] = ('COMMENT: ' + str(lognorm_set[i][j]))

        # Creating a new file for each uncertainty perturbation
        new_file = NTF(mode='w', dir=cwd3, suffix='.txt', delete=False, encoding="utf8")

        # and write everything back into their own new file
        for item in oldlinelist:
            new_file.write("%s\n" % item)

        new_file.close()
        progressbar.next()
    progressbar.finish()

    return cwd3


# Creating the input file for Bolsig minus
def bolsig_file(inputdir, cwd, og_cwd, species_name, num_EN_values):
    # Making a directory for the run files for Bolsig
    makedirs('BOLSIG_Run_Files')

    # Opening the template runfile
    file = open(join(og_cwd, 'bolsig_run_script.dat'), "r+", encoding='utf8')

    # Splitting the lines of the file into a readable list
    linelist = file.read().splitlines()

    # Initial Conditions
    inputlocation = 0
    outputlocation = 0
    species = 0
    num_EN_loc = 0

    outputdir = str(cwd) + "/Output_Files"
    makedirs(outputdir)

    # Loop to find were in the file the input and output filenames go
    for i in range(len(linelist)):

        # Finding the Input filename
        if linelist[i].find('/ File') != -1:
            inputlocation = i
        elif linelist[i].find('/ Species') != -1:
            species = i
        elif linelist[i].find('/ Output File') != -1:
            outputlocation = i
        elif linelist[i].find('/ Number') != -1:
            num_EN_loc = i

    dir_list = listdir(inputdir)
    progressbar = Bar('Creating Run Files ', max=len(dir_list), suffix='%(percent)d%%')

    # Creating a run file for each data file
    for inputfile in dir_list:

        runfilename = ('run_' + inputfile)
        runfile = ('BOLSIG_Run_Files/' + runfilename)

        # Opening the runfile in a writable capacity
        with open(runfile, 'w', encoding='utf8') as f:

            # Updating the file with the respective input and output filenames
            # Along with every other line in the template file
            for j in range(len(linelist)):

                if j == inputlocation:
                    f.write(('\"' + inputdir + "/" + inputfile + '\"' + '    / File\n'))
                elif j == outputlocation:
                    f.write(('\"' + outputdir + '/output_' + inputfile + '\"' + '      / Output File\n'))
                elif j == species:
                    f.write((species_name + '    / Species\n'))
                elif j == num_EN_loc:
                    f.write((str(num_EN_values) + '    / Number\n'))
                else:
                    f.write(linelist[j] + '\n')

        progressbar.next()
    progressbar.finish()

    return outputdir


# Function to Run BOLSIG-
def bolsig_minus(cwd, bolsig, num_cpus):

    # Gets list of runfiles
    infile_list = glob.glob(str(cwd)+"/BOLSIG_Run_Files/run*")
    print(infile_list)

    progressbar = Bar('Running BOLSIG- ', max=len(infile_list), suffix='%(percent)d%%')

    # Function that will execute BOLSIG-
    def solver(infiles):
        proc = Popen([bolsig, infiles], stdout=DEVNULL)
        proc.wait()
        progressbar.next()

    with futures.ThreadPoolExecutor(max_workers=num_cpus) as execute_solver:
        execute_solver.map(solver, infile_list)

    progressbar.finish()
    return

# Read all output files and extract data
def output_file_read(outputdir, startsec, num_values, num_e):

    # List of data set title to be searched for
    titlelist = ['Electric field / N (Td)', 'Mobility *N (1/m/V/s)',
                 'Diffusion coefficient *N (1/m/s)',
                 'Energy mobility *N (1/m/V/s)', 'Energy diffusion coef. D*N (1/m/s)',
                 'Total collision freq. /N (m3/s)', 'Momentum frequency /N (m3/s)',
                 'Total ionization freq. /N (m3/s)', 'Townsend ioniz. coef. alpha/N (m2)',
                 'Inelastic energy loss coefficient (eV m3/s)',
                 'Rate coefficient (m3/s)', 'Energy loss coefficient (eV m3/s)']

    # Array to hold E/N values
    red_e_field = []

    # Setting up holding mattrices for each set of coefficients
    mobil = empty((num_e, num_values))
    dif_cof = empty((num_e, num_values))
    e_mob = empty((num_e, num_values))
    e_dif_cof = empty((num_e, num_values))
    tot_col_freq = empty((num_e, num_values))
    p_freq = empty((num_e, num_values))
    tot_ion_freq = empty((num_e, num_values))
    town_ion_cof = empty((num_e, num_values))
    inelastic_loss_cof = empty((num_e, num_values))
    rate_matrix = empty((num_e, len(startsec)))
    e_loss_matrix = empty((num_e, len(startsec)))
    name_rate = []
    name_energy = []
    rate_dict = {}
    e_loss_dict = {}

    outfile_list = glob.glob(outputdir + "/output*")

    # Loop to extract the outputted data from each file
    for outputfile, i in zip(outfile_list, range(len(outfile_list))):

        # Initial conditions for later
        r = 0
        e = 0
        name_rate.clear()
        name_energy.clear()

        # Opening each data file to be checked
        with open(outputfile, 'r', encoding='utf-8') as f:
            # extracting each line from the file and listing it to be searched
            linelist = f.read().splitlines()

        # Iterating over the list of lines
        for j in range(len(linelist)):
            # CHecking for the phrase specified above
            if linelist[j].find(titlelist[0]) != -1:
                # collecting all the subsequent data for that set
                for k in range(1, (num_e + 1)):
                    e_n, unused = linelist[(j + k)].split('\t')
                    # Saving E/N values
                    red_e_field.append(e_n)
            elif linelist[j].find(titlelist[1]) != -1:
                for k in range(1, (num_e + 1)):
                    e_n, mob = linelist[(j + k)].split('\t')
                    mobil[(k - 1), i] = mob
            # As above
            elif linelist[j].find(titlelist[2]) != -1:
                for k in range(1, (num_e + 1)):
                    e_n, dif = linelist[(j + k)].split('\t')
                    dif_cof[(k - 1), i] = dif
            # As above
            elif linelist[j].find(titlelist[3]) != -1:
                for k in range(1, (num_e + 1)):
                    e_n, e_mobility = linelist[(j + k)].split('\t')
                    e_mob[(k - 1), i] = e_mobility
            # As above
            elif linelist[j].find(titlelist[4]) != -1:
                for k in range(1, (num_e + 1)):
                    e_n, e_diffusion = linelist[(j + k)].split('\t')
                    e_dif_cof[(k - 1), i] = e_diffusion
            # As above
            elif linelist[j].find(titlelist[5]) != -1:
                for k in range(1, (num_e + 1)):
                    e_n, total_collision = linelist[(j + k)].split('\t')

                    tot_col_freq[(k - 1), i] = total_collision
            # As above
            elif linelist[j].find(titlelist[6]) != -1:
                for k in range(1, (num_e + 1)):
                    e_n, mom_freq = linelist[(j + k)].split('\t')
                    p_freq[(k - 1), i] = mom_freq
            # As above
            elif linelist[j].find(titlelist[7]) != -1:
                for k in range(1, (num_e + 1)):
                    e_n, total_ion = linelist[(j + k)].split('\t')
                    tot_ion_freq[(k - 1), i] = total_ion
            # As above
            elif linelist[j].find(titlelist[8]) != -1:
                for k in range(1, (num_e + 1)):
                    e_n, townsend_ion = linelist[(j + k)].split('\t')
                    town_ion_cof[(k - 1), i] = townsend_ion
            # As above
            elif linelist[j].find(titlelist[9]) != -1:
                for k in range(1, (num_e + 1)):
                    e_n, inelastic = linelist[(j + k)].split('\t')
                    inelastic_loss_cof[(k - 1), i] = inelastic
            # As above
            elif linelist[j].find(titlelist[10]) != -1:
                # Taking the name of the data set
                name_rate.append(linelist[j - 1].rstrip())
                for k in range(1, (num_e + 1)):
                    e_n, rate = linelist[(j + k)].split('\t')
                    # Saving the data to a collective matrix
                    rate_matrix[(k - 1), r] = rate
                r += 1
            # As above
            elif linelist[j].find(titlelist[11]) != -1:
                name_energy.append(linelist[j - 1].rstrip())
                for k in range(1, (num_e + 1)):
                    e_n, e_loss = linelist[(j + k)].split('\t')
                    e_loss_matrix[(k - 1), e] = e_loss
                e += 1

        # Saving each set of rate and energy loss coefficients to collective dictionary
        rate_dict[i] = rate_matrix.copy()
        e_loss_dict[i] = e_loss_matrix.copy()

    gen_mat = [mobil, dif_cof, e_mob, e_dif_cof, tot_col_freq, p_freq, tot_ion_freq, town_ion_cof,
               inelastic_loss_cof]

    # Reaturning all releveant mattrices and dictionaries
    return red_e_field, gen_mat, e_loss_dict, rate_dict, name_rate, name_energy, titlelist, num_e


# Calculate stats for each mobility set
def gen_stats(gen_mat, num_values):

        # Resultant matrix for the calculated values
    result = []

    for x in range(len(gen_mat)):

        y = gen_mat[x]
        result_single = []

        # Iterate over each E/N value
        for i in range(len(y[:, 0])):
            # Initial values
            m_0 = y[i, 0]
            s_0 = 0
            std_mean = 0
            std_pop = 0

            # Iterate for each values in the E/N set
            for j in range(1, num_values):
                m = m_0 + (y[i, j] - m_0) / (j + 1)
                s = s_0 + ((y[i, j] - m_0) * (y[i, j] - m))

                m_0 = m
                s_0 = s

                std_mean = sqrt(s_0 / ((num_values ** 2) - num_values))
                std_pop = sqrt(s_0 / (num_values - 1))

            result_single.append([m_0, std_mean, std_pop])

        result.append(result_single)

    return result


# Calculation of rate coefficient stats
def rate_cof(rate, names, num_values, num_e):
    # Retrieving just the values from the dictionary
    rate_matrix = list(rate.values())

    # Setting up a subsequent dictionary
    rates_mean_std = {}

    # Iterating over the number of types of rate coefficients
    for i in range(len(names)):

        # GEneral saving matrix for each data type
        result = []

        for j in range(0, num_e):

            # Extracting the initial value for each E/N value from each data type
            m_0 = rate_matrix[0][j, i]
            s_0 = 0
            std_mean = 0
            std_pop = 0

            # Iterate over number of data values
            for k in range(1, num_values):
                m = m_0 + ((rate_matrix[k][j, i]) - m_0) / (k + 1)

                s = s_0 + ((rate_matrix[k][j, i]) - m_0) * (rate_matrix[k][j, i] - m)

                m_0 = m
                s_0 = s

                std_mean = sqrt(s_0 / ((num_values ** 2) - num_values))
                std_pop = sqrt(s_0 / (num_values - 1))

            # Final calculation and saving of the mean and std values for each E/N value
            result.append([m_0, std_mean, std_pop])

        rates_mean_std[names[i]] = result.copy()

    return rates_mean_std


# As described in the above function but for energy loss
def energy_loss(e_loss, names, num_values, num_e):
    e_matrix = list(e_loss.values())
    energy_loss_mean_std = {}

    for i in range(len(names)):

        result = []

        for j in range(0, num_e):

            m_0 = e_matrix[0][j, i]
            s_0 = 0
            std_mean = 0
            std_pop = 0

            for k in range(1, num_values):
                m = m_0 + ((e_matrix[k][j, i]) - m_0) / (k + 1)

                s = s_0 + ((e_matrix[k][j, i]) - m_0) * (e_matrix[k][j, i] - m)

                m_0 = m
                s_0 = s

                std_mean = sqrt(s_0 / ((num_values ** 2) - num_values))
                std_pop = sqrt(s_0 / (num_values - 1))

            # Final calculation and saving of the mean and std values for each E/N value
            result.append([m_0, std_mean, std_pop])

        energy_loss_mean_std[names[i]] = result.copy()

    return energy_loss_mean_std


# Saving all calculated statistics to a single text file
def stats_file(red_e_field, all_stats, rate_stats, e_loss_stats, names_rate, names_energy, titlelist):

    # Opening said text file
    with open('data_stats.txt', 'w') as f:

        f.write('Below is the mean and standard deviation of the mean \nfor a range of data types listed, '
                'such data was \nobtained by using BOLSIG- and inputed data file.\n')

        for i in range(len(all_stats)):

            stats = all_stats[i]
            f.write('\n' + str(titlelist[i+1]) + '\n' + str(titlelist[0]) + '        Mean           STD of the Mean   '
                                                                            '   STD of the Pop.  Percentage Error\n')

            for j in range(len(stats[:])):

                if stats[j][0] != 0:

                    perc_error = (float(stats[j][2]) / float(stats[j][0])) * 100

                else:

                    perc_error = 0

                f.write('         ' + str(red_e_field[j]) + '             ' +
                        f_f(float(stats[j][0]), precision=6, unique=False) + '         ' +
                        f_f(float(stats[j][1]), precision=6, unique=False) + '         ' +
                        f_f(float(stats[j][2]), precision=6, unique=False) + '        ' +
                        str("%.3f" % perc_error) + '%\n')

        for i in range(len(names_rate)):

            f.write('\n' + names_rate[i] + ' - Rate Coefficient\n' + str(titlelist[0]) +
                    '        Mean           STD of the Mean      STD of the Pop.  Percentage Error\n')

            rate_data = rate_stats[names_rate[i]]

            for j in range(len(rate_data[:])):

                if rate_data[j][0] != 0:

                    perc_error = (float(rate_data[j][2])/float(rate_data[j][0]))*100

                else:

                    perc_error = 0

                f.write('         ' + str(red_e_field[j]) + '             ' +
                        f_f(float(rate_data[j][0]), precision=6, unique=False) + '         ' +
                        f_f(float(rate_data[j][1]), precision=6, unique=False) + '         ' +
                        f_f(float(rate_data[j][2]), precision=6, unique=False) + '        ' +
                        str("%.3f" % perc_error) + '%\n')

        for i in range(len(names_energy)):

            f.write('\n' + names_energy[i] + ' - Energy Loss Coefficient\n' + str(titlelist[0]) +
                    '        Mean           STD of the Mean      STD of the Pop.  Percentage Error\n')

            e_loss_data = e_loss_stats[names_energy[i]]

            for j in range(len(e_loss_data[:])):

                if e_loss_data[j][0] != 0:

                    perc_error = (float(e_loss_data[j][2]) / float(e_loss_data[j][0])) * 100

                else:

                    perc_error = 0

                f.write('         ' + str(red_e_field[j]) + '             ' +
                        f_f(float(e_loss_data[j][0]), precision=6, unique=False) + '         ' +
                        f_f(float(e_loss_data[j][1]), precision=6, unique=False) + '         ' +
                        f_f(float(e_loss_data[j][2]), precision=6, unique=False) + '        ' +
                        str("%.3f" % perc_error) + '%\n')


# Graph each set of statistics and save
def graphing_data(red_e_field, stats, rate_stats, e_loss_stats, names_rate, names_energy):

    # Data type names
    image_name = ['Mobility', 'Diffusion coefficient',
                  'Energy mobility', 'Energy diffusion coef.', 'Total collision freq.', 'Momentum frequency',
                  'Total ionization freq.', 'Townsend ioniz. coef.', 'Inelastic energy loss coefficient']

    # Correspoding units
    units = [' *N (1/m/V/s)', ' *N (1/m/s)', ' *N (1/m/V/s)', ' D*N (1/m/s)', ' /N (m3/s)', ' /N (m3/s)',
             ' /N (m3/s)', '  alpha/N (m2)', ' (eV m3/s)']

    # X-axis values are constant for all data types
    x = [float(i) for i in red_e_field[0:50]]

    # Iterating through all data set types to produce a graph of each
    for i in range(len(stats)):
        # Taking a the y-values and corresponding error for each x-value
        y = [item[0] for item in stats[i]]
        yerr = [item[2] for item in stats[i]]

        # Opening figure for the graph
        plt.figure()
        # X and Y axis labels
        plt.xlabel('Reduced Electric Field, E/N (Td)')
        plt.ylabel(image_name[i] + units[i])
        # Changing to log scales
        plt.yscale('log')
        plt.xscale('log')
        # Plotting errorbar graphs
        plt.errorbar(x, y, yerr=yerr, ecolor='k')
        plt.tight_layout()
        plt.savefig(image_name[i] + '.png')
        plt.close()

    # Name and units for Rate and Energy loss coefficients
    r_and_e = ['Rate coefficient', 'Energy loss coefficient']
    units_re = [' (m3/s)', ' (eV m3/s)']

    # Rate Coefficient
    plt.figure()
    plt.xlabel('Reduced Electric Field, E/N (Td)')
    plt.ylabel(r_and_e[0] + units_re[0])
    plt.yscale('log')
    plt.xscale('log')

    for j in range(len(names_rate)):
        y = [item[0] for item in rate_stats[names_rate[j]]]
        yerr = [item[2] for item in rate_stats[names_rate[j]]]

        plt.errorbar(x, y, yerr=yerr)
        plt.tight_layout()

    # plt.legend(names_rate)
    plt.savefig(r_and_e[0] + '.png')
    plt.close()

    # Energy loss Coeffiecient
    plt.figure()
    plt.xlabel('Reduced Electric Field, E/N (Td)')
    plt.ylabel(r_and_e[1] + units_re[1])
    plt.yscale('log')
    plt.xscale('log')

    for j in range(len(names_energy)):
        y = [item[0] for item in e_loss_stats[names_energy[j]]]
        yerr = [item[2] for item in e_loss_stats[names_energy[j]]]

        plt.errorbar(x, y, yerr=yerr)
        plt.tight_layout()

    # plt.legend(names_rate)
    plt.savefig(r_and_e[1] + '.png')
    plt.close()

    # for k in range(len(names_energy)):
    #  plt.figure()
    # plt.title(names_energy[k])
    # plt.xlabel('Reduced Electric Field, E/N (Td)')
    # plt.ylabel(r_and_e[1] + units_re[1])
    # plt.yscale('log')
    # plt.xscale('log')
    # y = [item[0] for item in e_loss_stats[names_energy[k]]]
    # yerr = [item[2] for item in e_loss_stats[names_energy[k]]]
    # plt.errorbar(x, y, yerr=yerr, ecolor='k')
    # plt.savefig(r_and_e[1] + str(k) + '.png')
    # plt.close()

    mob = [item[0] for item in stats[0]]
    dif = [item[0] for item in stats[1]]
    mob_err = [item[2] for item in stats[0]]
    dif_err = [item[2] for item in stats[1]]
    dif_mob = []
    dif_mob_err = []

    for i in range(len(stats[0])):

        dif_mob.append(dif[i]/mob[i])

        dif_mob_err.append(dif_mob[i]*sqrt((mob_err[i]/mob[i])**2 + (dif_err[i]/dif[i])**2))

    plt.figure()
    plt.xlabel('Reduced Electric Field, E/N (Td)')
    plt.ylabel(r'$D_T/\mu, (V)$')
    plt.yscale('log')
    plt.xscale('log')
    plt.errorbar(x, dif_mob, yerr=dif_mob_err, ecolor='k')
    plt.tight_layout()
    plt.savefig('Dif_mobility_ratio.png')
    plt.close()


"""

Morris Method Code:

*** This code was created to assess the sensitivity plasma transport coefficients have to variations in their
    corresponding cross-sections ***

"""


def Morris(subject, data, uncert_set, runs, species, num_EN_values, bolsig, num_cpus, p_values):

    start = time.time()

    og_cwd = getcwd()

    cwd = changedir(og_cwd, subject)

    cwd2 = join(cwd, 'Morris')

    # Checking if directory exists, if yes remove and replace with clean directory
    if exists(cwd2):

        rmtree(cwd2)
        # Making and navigating to a new directory to store the new files
        makedirs(cwd2)

    elif not exists(cwd2):

        makedirs(cwd2)

    logging.basicConfig(filename=join(cwd2, "program.log"), filemode='w', level=logging.DEBUG,
                        format='%(asctime)s %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p')

    linelist = data_set(data)

    startsec, dashsec = search_set(linelist)

    com_format(startsec, dashsec, linelist, data)

    linelist = data_set(data)

    startsec, dashsec  = search_set(linelist)

    commsec = comment_sec(startsec, dashsec, linelist)

    dataset_name = read_energy(startsec, linelist)

    sample_set, delta, num_params = run_morris(uncert_set, p_values, runs)

    sample_set = cdf(uncert_set, sample_set)

    logging.info('Morris trajectories created - working directories created')

    cs_data = find_num(startsec, dashsec, linelist)

    inputfiles = data_edit_morris(cs_data, sample_set, data, cwd, dashsec, commsec, cwd2)

    logging.info('Input files created')

    inputdir = bolsig_file_morris(inputfiles, cwd2, og_cwd, species, num_EN_values)

    logging.info('BOLSIG run files created')

    bolsig_minus_morris(inputdir, bolsig, num_cpus)

    logging.info('All BOLSIG simulations complete')

    mean, std, energy = morris_stats(inputfiles, inputdir, sample_set, delta, num_EN_values, num_params)

    logging.info('Statistics completed')

    print_stats(mean, std, energy, cwd2, dataset_name)

    logging.info('Stats file created')
    logging.info('Normalising Stats')

    normal_stats(mean, std, energy, dataset_name)

    logging.info('Stats Normalised')

    end = time.time()

    print(end-start)


# Reading the energy of each collision for information later
def read_energy(startsec, linelist):

    data_energy = []
    for x in startsec:
        data_energy.append(linelist[x+2])

    return data_energy


# Compiling all matrices and calculation Morris trajectory matrix
def run_morris(uncert_set, num_levels, runs):

    num_params = len(uncert_set)
    delta = num_levels/2

    # Empty array to contain trajectories
    initial_traj = zeros(((num_params + 1), num_params))
    sample_set = {}

    progressbar = Bar('Generating Trajectories  ', max=runs, suffix='%(percent)d%%')

    # Calculating trajectories
    for iterations in range(runs):

        # Initial and input parameters for randomisation
        traj = initial_traj.copy()
        route = list(range(num_params))
        shuffle(route)

        # Setting up first location
        for i in range(num_params):
            traj[0][i] = randint(0, num_levels-1)

        # Changing one parameter at a time to create a random walk
        for j in range(num_params):
            # Copying previous location
            for m in range(num_params):
                traj[j+1][m] = traj[j][m]

            # Defining which parameter needs to be changed
            rnd_loc = route[j]
            # Changing the parameter
            if traj[j+1][rnd_loc] < delta:
                traj[j+1][rnd_loc] += delta
            else:
                traj[j+1][rnd_loc] -= delta

        # Normalising values to exclusively between 0 and 1
        traj += 0.5
        traj = divide(traj, num_levels)
        sample_set[iterations] = traj

        progressbar.next()
    progressbar.finish()

    return sample_set, delta, num_params


def cdf(uncert_set, sample_set):

    progressbar = Bar('Converting Trajectories ', max=len(sample_set), suffix='%(percent)d%%')

    # Distribution Parameters ( MEan and Variance)
    mean = 1
    for key in sample_set:
        sample = sample_set[key]

        for i in range(len(sample)):
            for j in range(len(uncert_set)):
                x = sample[i, j].copy()
                var = float(uncert_set[j])**2

                # Calcualting std from a specified formula
                sigma = sqrt(log(1 + (var / (mean ** 2))))
                # Calcualting mean from a specified formula
                mu = log(mean ** 2 / sqrt(var + (mean ** 2)))

                # Percetange Point Function to convert Morris sample values to perturbation values
                sample[i, j] = stats.lognorm.ppf(x, sigma, mu)
        progressbar.next()
    progressbar.finish()

    return sample_set


# Using the generated random numbers to create an equal amount of perturbed data files
def data_edit_morris(cs_data, sample_set, oldfile, cwd, dashsec, commsec, cwd2):
    # Importing the list of cross-section data from the dictionary created earlier
    cross_values = list(cs_data.values())

    chdir(cwd2)

    cwd3 = join(cwd2, 'Input_Files')
    makedirs(cwd3)
    new_file_database = {}

    progressbar = Bar('Generating Datasets ', max=len(sample_set), suffix='%(percent)d%%')

    for key in sample_set:
        new_dir_name = ('File_set_' + str(key))
        cwd4 = join(cwd3, new_dir_name)
        makedirs(cwd4)

        # Opening the currently used data file in order to read write equivalent files
        with open(join(cwd, oldfile), 'r', encoding='utf8') as old_data:
            # read lines into a list
            oldlinelist = old_data.read().splitlines()

        sample = sample_set[key]
        new_file_list = []

        for i in range(len(sample)):
            cvals = cross_values.copy()

            for j in range(len(sample[0])):
                # Selecting 1 set of data at a time
                dataset = cvals[j].copy()
                # Iterating over each cross section individually
                for k in range(len(dataset)):
                    # Splitting the Energy and Cross-Section numbers
                    E, cs = dataset[k].split('\t')

                    # Perturbing the Cross section number by the jth generated number
                    new_cs = (float(cs) * sample[i, j])
                    # Creating a new data line with the new cross section number
                    newdata = (E + '\t' + str(new_cs))
                    # Updating the list
                    dataset[k] = newdata

                # Updating the entire data set
                cvals[j] = dataset.copy()
                # Adding the newly perturbed data lines to the main list
                for x, y in zip(range((dashsec[j * 2] + 1), (dashsec[(j * 2) + 1])), range(len(dataset))):
                    oldlinelist[x] = dataset[y]

                oldlinelist[commsec[j]] = ('COMMENT: ' + str(sample[i, j]))

            # Creating a new file for each uncertainty perturbation
            new_file = NTF(mode='w', dir=cwd4, suffix='.txt', delete=False, encoding="utf8")
            new_file_list.append(new_file.name.replace((cwd4+str("/")), ""))

            # and write everything back into their own new file
            for item in oldlinelist:
                new_file.write("%s\n" % item)

            new_file.close()

        new_file_database[new_dir_name] = new_file_list
        progressbar.next()
    progressbar.finish()

    return new_file_database


# Creating the input file for Bolsig minus
def bolsig_file_morris(filename, cwd, og_cwd, species_name, num_EN_values):

    # Opening the template runfile
    file = open(join(og_cwd, 'bolsig_run_script.dat'), "r+", encoding='utf8')
    # Splitting the lines of the file into a readable list
    linelist = file.read().splitlines()

    # Initial Conditions
    inputlocation = 0
    outputlocation = 0
    species = 0
    num_EN_loc = 0

    # Loop to find were in the file the input and output filenames go
    for i in range(len(linelist)):

        # Finding the Input filename
        if linelist[i].find('/ File') != -1:
            inputlocation = i
        # Finding the output file name
        elif linelist[i].find('/ Species') != -1:
            species = i
        # Finding the output file name
        elif linelist[i].find('/ Output File') != -1:
            outputlocation = i
        elif linelist[i].find('/ Number') != -1:
            num_EN_loc = i

    input_dirs = list(filename.keys())
    inputdir = str(cwd) + '/Input_Files/'

    progressbar = Bar('Creating Run Files ', max=len(input_dirs), suffix='%(percent)d%%')

    # Creating a run file for each data file
    for x in range(len(input_dirs)):

        chdir(join(inputdir, input_dirs[x]))
        cwd = getcwd()
        file_list = listdir(cwd)

        for inputfile in file_list:
            runfile = ('run_' + inputfile)

            # Opening the runfile in a writable capacity
            with open(runfile, 'w', encoding='utf8') as f:
                # Updating the file with the respective input and output filenames
                # Along with every other line in the template file
                for j in range(len(linelist)):
                    if j == inputlocation:
                        f.write(('\"' + cwd + "/" + inputfile + '\"' + '    / File\n'))
                    elif j == outputlocation:
                        f.write(('\"' + cwd + '/output_' + inputfile + '\"' + '      / Output File\n'))
                    elif j == species:
                        f.write((species_name + '    / Species\n'))
                    elif j == num_EN_loc:
                        f.write((str(num_EN_values) + '    / Number\n'))
                    else:
                        f.write(linelist[j] + '\n')
        progressbar.next()
    progressbar.finish()

    return inputdir


# Function to Run BOLSIG-
def bolsig_minus_morris(inputdir, bolsig, num_cpus):

    dir_list = listdir(inputdir)
    logging.info('BOLSIG- calculations have begun')

    progressbar = Bar('Running BOLSIG- ', max=len(dir_list), suffix='%(percent)d%%')

    # Function that will execute BOLSIG-
    def solver(infiles):
        proc = Popen([bolsig, infiles], stdout=DEVNULL)
        proc.wait()

    for item in dir_list:

        logging.debug('BOLSIG now running in ' + item)
        runfiles = glob.glob(inputdir + item + "/run*")

        # This will map the solver function to multiple processors
        with futures.ThreadPoolExecutor(max_workers=num_cpus) as execute_solver:
            execute_solver.map(solver, runfiles)
        progressbar.next()

    progressbar.finish()
    return


def morris_stats(file_list, inputdir, sample_set, delta, num_EN_val, num_params):

    logging.debug('Mobility statistical calculations have begun')

    # Initial Matrices and values
    mob_start = 0
    diff_start = 0
    ion_start = 0
    excite_rate = []
    ion_rate = []

    mean = {}
    std = {}
    energy = zeros(num_EN_val)

    # Interating through each created directory in Input Files
    for keys, num in zip(file_list, sample_set):

        # Entering each directory and setting up list of input files and corresponding morris trajectories
        chdir(join(inputdir, keys))
        input_files = file_list[keys]
        sample = sample_set[num]

        # Iterate over the one less than the number of input file, so as to take the the changes between each
        # trajectory point
        for i in range(len(input_files)-1):

            # Setting up two files being tested
            filename_1 = 'output_' + input_files[i]
            filename_2 = 'output_' + input_files[i+1]

            logging.debug('Computing statistics between ' + filename_1 + ' and ' + filename_2 + ' within directory '
                          + keys)

            # Calculating from trajectory matrix which parameter changed and in which direction
            sample[i, :] -= sample[i + 1, :]
            row_sum = sample.sum(axis=1)
            value_loc = nonzero(sample[i])

            # Opening files in question for data to be extracted
            file_1 = open(filename_1, 'r', encoding='utf8').read().splitlines()
            file_2 = open(filename_2, 'r', encoding='utf8').read().splitlines()

            # Searching the first file for the location of mobility data
            # This is only done once as each file is of the exact same format
            if mob_start == 0:
                for j in range(len(file_1)):
                    if 'Mobility' in file_1[j]:
                        mob_start = j
                    elif 'Diffusion coefficient' in file_1[j]:
                        diff_start = j
                    elif 'Total ionization freq.' in file_1[j]:
                        ion_start = j
                    elif 'Excitation' in file_1[j] and 'Rate coefficient' in file_1[j+1]:
                        excite_rate.append(j+1)
                    elif 'Ionization' in file_1[j] and 'Rate coefficient' in file_1[j+1]:
                        ion_rate.append(j+1)
                    elif 'Energy loss coefficient' in file_1[j+1]:
                        break
                    else:
                        continue

            for stat, stat_name in zip([mob_start, diff_start, ion_start, excite_rate, ion_rate], ['mob', 'dif', 'ion', 'excite_rate', 'ion_rate']):

                if isinstance(stat, list):
                    stat_list = stat
                    stat_list_name = stat_name
                    for i in range(len(stat_list)):
                        stat_name = stat_list_name + "_{}".format(i)
                        stat = stat_list[i]

                        # Iterating through each line of mobility data
                        for k in range(num_EN_val):

                            # extracting energy and mobility values
                            E, stat_1 = file_1[stat + k + 1].split('\t')
                            E, stat_2 = file_2[stat + k + 1].split('\t')

                            # The set of files must act as the data that sets the initial values
                            if num == 0:

                                mean[stat_name] = zeros((num_EN_val, num_params))
                                std[stat_name] = zeros((num_EN_val, num_params))

                                stat_mean = mean[stat_name]
                                stat_std = std[stat_name]

                                # The direction of trajectory change determines the formula used
                                if row_sum[i] < 0:

                                    # Elementary Effect calculation
                                    EE = (float(stat_2)-float(stat_1))/delta
                                    EE = abs(EE)

                                    # Setting the initial elements of the matrices
                                    stat_mean[k, value_loc] = EE
                                    stat_std[k, value_loc] = 0
                                    energy[k] = E

                                else:

                                    EE = (float(stat_1) - float(stat_2)) / delta
                                    EE = abs(EE)

                                    stat_mean[k, value_loc] = EE
                                    stat_std[k, value_loc] = 0
                                    energy[k] = E

                            # For all data after the first data
                            else:

                                stat_mean = mean[stat_name]
                                stat_std = std[stat_name]

                                if row_sum[i] < 0:

                                    EE = (float(stat_2) - float(stat_1)) / delta
                                    EE = abs(EE)

                                    # Calculating moving averages and standard deviations as data is read
                                    m = stat_mean[k, value_loc] + (EE - stat_mean[k, value_loc]) / (num + 1)
                                    s = stat_std[k, value_loc] + ((EE - stat_mean[k, value_loc]) * (EE - m))

                                    # Redefined mean and std values
                                    stat_mean[k, value_loc] = m
                                    stat_std[k, value_loc] = s

                                else:

                                    EE = (float(stat_1) - float(stat_2)) / delta
                                    EE = abs(EE)

                                    m = stat_mean[k, value_loc] + (EE - stat_mean[k, value_loc]) / (num + 1)
                                    s = stat_std[k, value_loc] + ((EE - stat_mean[k, value_loc]) * (EE - m))

                                    stat_mean[k, value_loc] = m
                                    stat_std[k, value_loc] = s
                else:
                    # Iterating through each line of mobility data
                    for k in range(num_EN_val):

                        # extracting energy and mobility values
                        E, stat_1 = file_1[stat + k + 1].split('\t')
                        E, stat_2 = file_2[stat + k + 1].split('\t')

                        # The set of files must act as the data that sets the initial values
                        if num == 0:

                            mean[stat_name] = zeros((num_EN_val, num_params))
                            std[stat_name] = zeros((num_EN_val, num_params))

                            stat_mean = mean[stat_name]
                            stat_std = std[stat_name]

                            # The direction of trajectory change determines the formula used
                            if row_sum[i] < 0:

                                # Elementary Effect calculation
                                EE = (float(stat_2) - float(stat_1)) / delta
                                EE = abs(EE)

                                # Setting the initial elements of the matrices
                                stat_mean[k, value_loc] = EE
                                stat_std[k, value_loc] = 0
                                energy[k] = E

                            else:

                                EE = (float(stat_1) - float(stat_2)) / delta
                                EE = abs(EE)

                                stat_mean[k, value_loc] = EE
                                stat_std[k, value_loc] = 0
                                energy[k] = E

                        # For all data after the first data
                        else:

                            stat_mean = mean[stat_name]
                            stat_std = std[stat_name]

                            if row_sum[i] < 0:

                                EE = (float(stat_2) - float(stat_1)) / delta
                                EE = abs(EE)

                                # Calculating moving averages and standard deviations as data is read
                                m = stat_mean[k, value_loc] + (EE - stat_mean[k, value_loc]) / (num + 1)
                                s = stat_std[k, value_loc] + ((EE - stat_mean[k, value_loc]) * (EE - m))

                                # Redefined mean and std values
                                stat_mean[k, value_loc] = m
                                stat_std[k, value_loc] = s

                            else:

                                EE = (float(stat_1) - float(stat_2)) / delta
                                EE = abs(EE)

                                m = stat_mean[k, value_loc] + (EE - stat_mean[k, value_loc]) / (num + 1)
                                s = stat_std[k, value_loc] + ((EE - stat_mean[k, value_loc]) * (EE - m))

                                stat_mean[k, value_loc] = m
                                stat_std[k, value_loc] = s

    # Final calculation for STD
    for name in std.keys():
        std[name] = sqrt(divide(std[name], (len(sample_set) - 1)))

    return mean, std, energy


def print_stats(mean_stats, std_stats, energy, cwd, dataset):

    chdir(cwd)

    # Creating a file contain the calculated data
    with open('stats.txt', 'w', encoding='utf8') as f:

        coef_names = ['Mobility', '\nDiffusion Coef.', '\nIonisation coef.']
        excite_index = 1
        ion_index = 1

        for keys in mean_stats:
            if 'excite' in keys:
                coef_names.append("\nExcition Rate Coefficient {}".format(excite_index))
                excite_index += 1
            if 'ion_rate' in keys:
                coef_names.append("\nIonization Rate Coefficient {}".format(ion_index))
                ion_index += 1

        # Iterate over all transport coefficients
        for keys, name in zip(mean_stats, coef_names):

            mean = mean_stats[keys]
            std = std_stats[keys]

            f.write(name + '\n\n')

            f.write('Mean\nE/N    ')

            for i in range(len(dataset)):
                f.write(dataset[i] + ' eV    ')

            f.write('\n')

            for i in range(len(mean)):
                f.write(f_f(float(energy[i])) + '    ')

                for j in range(len(dataset)):
                    f.write(f_f(float(mean[i, j]), precision=3, unique=False) + '          ')
                f.write('\n')

            f.write('STD\nE/N    ')

            for i in range(len(dataset)):
                f.write(dataset[i] + ' eV    ')
            f.write('\n')

            for i in range(len(mean)):
                f.write(f_f(float(energy[i])) + '    ')
                for j in range(len(dataset)):

                    f.write(f_f(float(std[i, j]), precision=3, unique=False) + '          ')
                f.write('\n')


def normal_stats(mean_stats, std_stats, energy, dataset):

    # Creating a file contain the calculated data
    with open('norm_stats.txt', 'w', encoding='utf8') as f:

        coef_names = ['Mobility', '\nDiffusion Coef.', '\nIonisation coef.']
        excite_index = 1
        ion_index = 1

        for keys in mean_stats:
            if 'excite' in keys:
                coef_names.append("\nExcition Rate Coefficient {}".format(excite_index))
                excite_index += 1
            if 'ion_rate' in keys:
                coef_names.append("\nIonization Rate Coefficient {}".format(ion_index))
                ion_index += 1

        # Iterate over all transport coefficients
        for keys, name in zip(mean_stats, coef_names):

            mean = mean_stats[keys]
            std = std_stats[keys]

            seterr(divide='ignore', invalid='ignore')
            mean = mean / mean.sum(axis=1)[:, None]
            std = std/std.sum(axis=1)[:, None]

            f.write(name + '\n\n')

            f.write('Mean\nE/N    ')

            for i in range(len(dataset)):
                f.write(dataset[i] + ' eV    ')

            f.write('\n')

            for i in range(len(mean)):

                f.write(f_f(float(energy[i])) + '    ')

                for j in range(len(dataset)):
                    f.write(f_f(float(mean[i, j]), precision=3, unique=False) + '          ')
                f.write('\n')

            f.write('STD\nE/N    ')

            for i in range(len(dataset)):
                f.write(dataset[i] + ' eV    ')

            f.write('\n')

            for i in range(len(mean)):

                f.write(f_f(float(energy[i])) + '    ')

                for j in range(len(dataset)):
                    f.write(f_f(float(std[i, j]), precision=3, unique=False) + '          ')
                f.write('\n')

run_program()
