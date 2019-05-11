from os import getcwd, chdir, makedirs, listdir, remove
from os.path import join, exists
from numpy import random as rd, empty, format_float_scientific as f_f, zeros, nonzero, divide, sqrt, seterr
from math import log
from tempfile import NamedTemporaryFile as NTF
from subprocess import Popen
from matplotlib import pyplot as plt
from shutil import rmtree
import time
import logging
from random import randint, shuffle
from scipy import stats


__author__ = "Rhys Doyle"
__copyright__ = "Copyright 2019, LTP Uncertainty Calculator, MCMU"
__credits__ = ["Prof. Miles Turner", "Dublin City University"]
__maintainer__ = "Rhys Doyle"
__email__ = "rhys.doyle45@mail.dcu.ie"
__status__ = "Build Complete"


def run_program():

    program = input('What program would you like to run? (1:Monte-Carlo, 2:Morris) ')
    i = 0

    while i == 0:

        if program == '1':

            i = 1
            Monte_carlo()

        elif program == '2':

            i = 1
            Morris()

        else:

            print('You did not enter a valid input. Please try again.')
            program = input('What program would you like to run? (1:Monte-Carlo, 2:Morris) ')


"""  

Monte-Carlo Simulation Code:

*** This code was created to determine the possible uncertainties associated with plasma transport coefficients
    using the experimental uncertainties associated with the cross-sections used ***

"""


# Run function for Monte-Carlo simulations
def Monte_carlo():

    start = time.time()

    og_cwd = getcwd()

    cwd, data = changedir(getcwd())

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

    file = data_set(data)

    items = wordlist()

    startsec, dashsec, linelist = search_set(file, items)

    com_format(startsec, dashsec, linelist, data)

    file = data_set(data)

    startsec, dashsec, linelist = search_set(file, items)

    commsec = comment_sec(startsec, dashsec, linelist)

    uncert_set = com_uncert(startsec, linelist)

    lognorm_set, num_values = lognorm(uncert_set)

    cs_data = find_num(startsec, dashsec, linelist)

    inputfiles = data_edit(cs_data, lognorm_set, data, cwd, dashsec, num_values, commsec, cwd2)

    logging.info('Data files edited')

    runfilelist, outputfilelist, num_en = bolsig_file(inputfiles, cwd2, og_cwd)

    logging.info('BOLSIG run files created')
    logging.info('BOLSIG- running')

    outputdir = bolsig_minus(runfilelist, cwd2, og_cwd)

    logging.info('BOLSIG- completed')

    red_e_field, gen_mat, e_loss_dict, rate_dict, name_rate, name_energy, titlelist, num_e = \
        output_file_read(outputfilelist, outputdir, startsec, num_values, num_en)

    stats = gen_stats(gen_mat, num_values)

    rate_stats = rate_cof(rate_dict, name_rate, num_values, num_e)

    e_loss_stats = energy_loss(e_loss_dict, name_energy, num_values, num_e)

    logging.info('Stats Generated')

    stats_file(red_e_field, stats, rate_stats, e_loss_stats, name_rate, name_energy, titlelist)

    logging.info('Stats file created')

    graphing_data(red_e_field, stats, rate_stats, e_loss_stats, name_rate, name_energy)

    end = time.time()

    print(end - start)


# Changing directory for data sets
def changedir(path):
    # requesting data file location within main directory
    subject = str(input("Please enter the path to your data location. "))

    # Changing main directory to above location
    chdir(join(path, subject))
    # Displaying new main directory
    cwd = getcwd()
    print(cwd)

    # Prompt for data filename
    data = str(input("Please enter the full file name. "))

    return cwd, data


# Importing data file to program
def data_set(data):
    # Opening data file in a Read/Write capacity
    file = open(data, "r+", encoding='utf8')

    return file


# Creating list of words to search for
def wordlist():
    # global items

    items = ['EXCITATION', 'ELASTIC', 'IONIZATION', 'EFFECTIVE', 'ROTATIONAL', 'ATTACHMENT']
    print(items)

    return items


# Searching Data set for locations of COMMENT lines of different Cross Sections
def search_set(a, j):
    # split the file content into lines and save a list
    linelist = a.read().splitlines()
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
        if any(word in line for word in j):

            startsec.append(y)

            y += 1

        # Looking for the Column section
        elif line.find('---------------') != -1:

            dashsec.append(y)

            y += 1
        else:

            y += 1

    return startsec, dashsec, linelist


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

    print(out_file)
    with open(data, 'w', encoding='utf-8') as f:
        f.writelines(out_file)


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


# Editing the data file to add the uncertainty on each cross section to the COMMENT section the the CS
def com_uncert(startsec, linelist):

    # Array to store the uncertainties of each data set as inputted below
    uncert_set = []

    # For loop to edit each cross section
    for i in range(0, (len(startsec))):
        # Line data lists
        x = startsec[i]
        uncert = float(input("What is the uncertainty for " + str(linelist[x]) + ", " + str(linelist[x + 1]) + ":"))
        # Adding uncertainty input to data set
        uncert_set.append(uncert)

    return uncert_set


# Generating a set amount of numbers from a lognormal distribution of specified parameters
def lognorm(uncert_set):
    # Distribution Parameters ( MEan and Variance)
    mean = 1
    lognorm_set = []
    num_values = int(input('How many perturbation values would you like? '))

    for i in range(len(uncert_set)):
        var = (uncert_set[i])**2

        # Calcualting std from a specified formula
        sigma = sqrt(log(1 + (var / (mean ** 2))))
        # Calcualting mean from a specified formula
        mu = log(mean ** 2 / sqrt(var + (mean ** 2)))

        # Gemerating the random numbers
        uncert = rd.lognormal(mu, sigma, num_values)
        # uncert = ones(num_values)

        lognorm_set.append(uncert)

    return lognorm_set, num_values


# Finding and extracting each individual data set
def find_num(startsec, dashsec, linelist):
    # Dictionary to collect the data for each Cross Section
    cs_data = {}

    # Iteratively adding to the above dictionary
    for i in range(0, len(startsec)):
        cs_name = linelist[startsec[i] + 1]
        cs_start = dashsec[i * 2] + 1
        cs_end = dashsec[(i * 2) + 1]

        cs_set = linelist[cs_start:cs_end]

        cs_data[(str(cs_name) + ", No." + str(i))] = cs_set

    return cs_data


# Using the generated random numbers to create an equal amount of perturbed data files
def data_edit(cs_data, lognorm_set, oldfile, cwd, dashsec, num_uncert, comsec, cwd2):
    # Importing the list of cross-section data from the dictionary created earlier
    cross_values = list(cs_data.values())
    # Creating a second variable of the same list to reset the first when altered
    reset = cross_values.copy()

    chdir(cwd2)

    cwd3 = join(cwd2, 'Input Files')
    makedirs(cwd3)

    # For loop for each random number
    for j in range(0, num_uncert):

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

    newfilelist = listdir(cwd3)
    print(newfilelist)

    return newfilelist


# Creating the input file for Bolsig minus
def bolsig_file(filename, cwd, og_cwd):
    # Making a directory for the run files for Bolsig
    makedirs('BOLSIG Run Files')

    # Opening the template runfile
    file = open(join(og_cwd, 'bolsig_run_script.dat'), "r+", encoding='utf8')

    # Splitting the lines of the file into a readable list
    linelist = file.read().splitlines()

    # Initial Conditions
    inputlocation = 0
    outputlocation = 0
    species = 0
    num_EN_loc = 0

    inputdir = str(cwd) + '\\Input Files\\'

    # Loop to find were in the file the input and output filenames go
    for i in range(len(linelist)):

        # Finding the Input filename
        if linelist[i].find('/ File ') != -1:

            inputlocation = i

        elif linelist[i].find('/ Species') != -1:

            species = i

        elif linelist[i].find('/ Output File ') != -1:

            outputlocation = i

        elif linelist[i].find('/ Number') != -1:

            num_EN_loc = i

    # Creating an array for the runfile names
    runfilelist = []
    output_filelist = []

    species_name = str(input('What species is being studied? '))
    num_EN_values = int(input('How many E/N values would you like? '))

    # Creating a run file for each data file
    for x in range(len(filename)):

        runfilename = ('run_' + str(x) + '.dat')

        runfile = ('BOLSIG Run Files\\' + runfilename)

        # Opening the runfile in a writable capacity
        with open(runfile, 'w', encoding='utf8') as f:

            # Updating the file with the respective input and output filenames
            # Along with every other line in the template file
            for j in range(len(linelist)):

                if j == inputlocation:

                    f.write(('\"' + inputdir + filename[x] + '\"' + '    / File\n'))

                elif j == outputlocation:

                    output_filename = ('output_' + filename[x])

                    # Updating the list of output file names
                    output_filelist.append(output_filename)

                    f.write((output_filename + '    / Output File\n'))

                elif j == species:

                    f.write((species_name + '    / Species\n'))

                elif j == num_EN_loc:

                    f.write((str(num_EN_values) + '    / Number\n'))

                else:

                    f.write(linelist[j] + '\n')

        # Updating the list with the name of the run file
        runfilelist.append(runfilename)

    return runfilelist, output_filelist, num_EN_values


# Function to Run BOLSIG-
def bolsig_minus(runfilelist, cwd, og_cwd):
    # Create a new output directory
    outdir = 'Output Files'

    makedirs(outdir)
    # Specify the location of BOLSIG-
    bolsig = ("\"" + str(og_cwd) + "\\bolsigminus.exe\"")

    print(bolsig)

    # Specify the output directory
    outputdir = str(cwd) + '\\Output Files'

    # Iteratively create each put file through BOLSIG-
    for i in range(len(runfilelist)):
        # Specify the input file location
        infile = ("\"" + str(cwd) + "\\BOLSIG Run Files\\\"" + str(runfilelist[i]))

        # Temporary file allowing for Popen to run
        f = NTF('w', suffix='.txt', delete=True, encoding='utf-8')

        # Run BOLSIG-
        r = Popen("%s %s" % (bolsig, infile), stdout=f, cwd=outputdir)
        # Wait for BOLSIG- to complete
        Popen.wait(r)

    return outputdir


# Read all output files and extract data
def output_file_read(outputfilelist, outputdir, startsec, num_values, num_e):
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

    # Loop to extract the outputted data from each file
    for i in range(len(outputfilelist)):

        # Initial conditions for later
        r = 0
        e = 0

        name_rate.clear()
        name_energy.clear()

        # Opening each data file to be checked
        with open((outputdir + '\\' + outputfilelist[i]), 'r', encoding='utf-8') as f:

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


def Morris():

    start = time.time()

    og_cwd = getcwd()

    cwd, data = changedir_morris(getcwd())

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

    file = data_set_morris(data)

    items = wordlist_morris()

    startsec, dashsec, linelist = search_set_morris(file, items)

    com_format_morris(startsec, dashsec, linelist, data)

    file = data_set_morris(data)

    startsec, dashsec, linelist = search_set_morris(file, items)

    commsec = comment_sec_morris(startsec, dashsec, linelist)

    uncert_set, dataset_name = com_uncert_morris(startsec, linelist)

    sample_set, delta, num_params = run_morris(uncert_set)

    sample_set = cdf(uncert_set, sample_set)

    logging.info('Morris trajectories created - working directories created')

    cs_data = find_num_morris(startsec, dashsec, linelist)

    inputfiles = data_edit_morris(cs_data, sample_set, data, cwd, dashsec, commsec, cwd2)

    logging.info('Input files created')

    num_EN_val = bolsig_file_morris(inputfiles, cwd2,og_cwd)

    logging.info('BOLSIG run files created')

    bolsig_minus_morris(cwd2, og_cwd)

    logging.info('All BOLSIG simulations complete')

    mean, std, energy, coef_names = morris_stats(inputfiles, sample_set, delta, num_EN_val, num_params)

    logging.info('Statistics completed')

    print_stats(mean, std, energy, cwd2, dataset_name, coef_names)

    logging.info('Stats file created')
    logging.info('Normalising Stats')

    normal_stats(mean, std, energy, dataset_name, coef_names)

    logging.info('Stats Normalised')

    end = time.time()

    print(end-start)


# Changing directory for data sets
def changedir_morris(path):
    # requesting data file location within main directory
    subject = str(input("Please enter the path to your data location. "))

    # Changing main directory to above location
    chdir(join(path, subject))
    # Displaying new main directory
    cwd = getcwd()
    print(cwd)

    # Prompt for data filename
    data = str(input("Please enter the full file name. "))

    return cwd, data


# Importing data file to program
def data_set_morris(data):
    # Opening data file in a Read/Write capacity
    file = open(data, "r+", encoding='utf8')

    return file


# Creating list of words to search for
def wordlist_morris():
    # global items

    items = ['EXCITATION', 'ELASTIC', 'IONIZATION', 'EFFECTIVE', 'ROTATIONAL', 'ATTACHMENT']
    print(items)

    return items


# Searching Data set for locations of COMMENT lines of different Cross Sections
def search_set_morris(a, j):
    # split the file content into lines and save a list
    linelist = a.read().splitlines()
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
        if any(word in line for word in j):

            startsec.append(y)

            y += 1

        # Looking for the Column section
        elif line.find('---------------') != -1:

            dashsec.append(y)

            y += 1
        else:

            y += 1

    return startsec, dashsec, linelist


# Function to search for Comment sections and add those that do not have them
def com_format_morris(startsec, dashsec, linelist, data):
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

    print(out_file)
    with open(data, 'w', encoding='utf-8') as f:
        f.writelines(out_file)

    return


# Function to collect all comment sections
def comment_sec_morris(startsec, dashsec, linelist):
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


# Editing the data file to add the uncertainty on each cross section to the COMMENT section the the CS
def com_uncert_morris(startsec, linelist):

    # Array to store the uncertainties of each data set as inputted below
    uncert_set = []
    energy = []

    # For loop to edit each cross section
    for i in range(0, (len(startsec))):
        # Line data lists
        x = startsec[i]
        # Experimental uncertainty in that cross section
        uncert = float(input("What is the uncertainty for " + str(linelist[x]) + ", " + str(linelist[x + 1]) + ":"))
        # Adding uncertainty input to data set
        uncert_set.append(uncert)

        energy.append(linelist[x+2])

    return uncert_set, energy


# General analysis specific inputs of Morris analysis
def morris_inputs(uncert_set):

    # Number of input parameters
    # num_params = int(input("Please enter the number of parameters describing the system. "))
    num_params = len(uncert_set)
    # Number of discrete values each input can take
    num_levels = int(input("Please enter the number of p-levels. "))
    # Number of iterations of the Morris analysis
    runs = int(input("How many runs? "))

    # Size of discrete step
    delta = num_levels/2

    return num_params, delta, runs, num_levels


# Compiling all matrices and calculation Morris trajectory matrix
def run_morris(uncert_set):

    num_params, delta, runs, num_levels = morris_inputs(uncert_set)

    # Empty array to contain trajectories
    initial_traj = zeros(((num_params + 1), num_params))

    sample_set = {}

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

    return sample_set, delta, num_params


def cdf(uncert_set, sample_set):

    # Distribution Parameters ( MEan and Variance)
    mean = 1

    for key in sample_set:

        sample = sample_set[key]

        for i in range(len(sample)):

            for j in range(len(uncert_set)):
                x = sample[i, j].copy()

                var = uncert_set[j]**2

                # Calcualting std from a specified formula
                sigma = sqrt(log(1 + (var / (mean ** 2))))
                # Calcualting mean from a specified formula
                mu = log(mean ** 2 / sqrt(var + (mean ** 2)))

                # Percetange Point Function to convert Morris sample values to perturbation values
                sample[i, j] = stats.lognorm.ppf(x, sigma, mu)

    print(sample_set)

    return sample_set


# Finding and extracting each individual data set
def find_num_morris(startsec, dashsec, linelist):
    # Dictionary to collect the data for each Cross Section
    cs_data = {}

    # Iteratively adding to the above dictionary
    for i in range(0, len(startsec)):
        cs_name = linelist[startsec[i] + 1]
        cs_start = dashsec[i * 2] + 1
        cs_end = dashsec[(i * 2) + 1]

        cs_set = linelist[cs_start:cs_end]

        cs_data[(str(cs_name) + ", No." + str(i))] = cs_set

    return cs_data


# Using the generated random numbers to create an equal amount of perturbed data files
def data_edit_morris(cs_data, sample_set, oldfile, cwd, dashsec, commsec, cwd2):
    # Importing the list of cross-section data from the dictionary created earlier
    cross_values = list(cs_data.values())

    chdir(cwd2)

    cwd3 = join(cwd2, 'Input Files')
    makedirs(cwd3)
    new_file_database = {}

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

            new_file_list.append(new_file.name.replace((cwd4+str("\\")), ""))

            # and write everything back into their own new file
            for item in oldlinelist:
                new_file.write("%s\n" % item)

            new_file.close()

        new_file_database[new_dir_name] = new_file_list

    print(new_file_database)
    print(getcwd())

    return new_file_database


# Creating the input file for Bolsig minus
def bolsig_file_morris(filename, cwd, og_cwd):

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
        if linelist[i].find('/ File ') != -1:

            inputlocation = i

        # Finding the output file name
        elif linelist[i].find('/ Species') != -1:

            species = i

        # Finding the output file name
        elif linelist[i].find('/ Output File ') != -1:

            outputlocation = i

        elif linelist[i].find('/ Number') != -1:

            num_EN_loc = i

    species_name = str(input('What species is being studied? '))
    num_EN_values = int(input('How many E/N values would you like? '))
    input_dirs = list(filename.keys())
    inputdir = str(cwd) + '\\Input Files\\'

    # Creating a run file for each data file
    for x in range(len(input_dirs)):

        chdir(join(inputdir, input_dirs[x]))
        cwd = getcwd()

        file_list = listdir(cwd)

        for i in range(len(file_list)):

            runfile = ('run_' + str(i) + '.dat')

            # Opening the runfile in a writable capacity
            with open(runfile, 'w', encoding='utf8') as f:

                # Updating the file with the respective input and output filenames
                # Along with every other line in the template file

                for j in range(len(linelist)):

                    if j == inputlocation:

                        f.write(('\"' + cwd + '\\' + file_list[i] + '\"' + '    / File\n'))

                    elif j == outputlocation:

                        output_filename = ('output_' + file_list[i])

                        f.write((output_filename + '    / Output File\n'))

                    elif j == species:

                        f.write((species_name + '    / Species\n'))

                    elif j == num_EN_loc:

                        f.write((str(num_EN_values) + '    / Number\n'))

                    else:

                        f.write(linelist[j] + '\n')

    return num_EN_values


# Function to Run BOLSIG-
def bolsig_minus_morris(cwd, og_cwd):

    # Specify the location of BOLSIG-
    bolsig = ("\"" + str(og_cwd) + "\\bolsigminus.exe\"")

    chdir(join(cwd, 'Input files'))
    cwd = getcwd()
    dir_list = listdir(cwd)

    logging.info('BOLSIG- calculations have begun')

    for item in dir_list:

        # Specify the output directory
        outputdir = str(cwd) + '\\' + item
        print('Current working directory is ' + item)
        logging.debug('BOLSIG now running in ' + item)

        for file in listdir(outputdir):

            if file.startswith('run_'):

                # Specify the input file location
                infile = ("\"" + outputdir + '\\' + file + "\"")

                # Temporary file allowing for Popen to run
                f = NTF('w', suffix='.txt', delete=True, encoding='utf-8')

                # Run BOLSIG-
                r = Popen("%s %s" % (bolsig, infile), stdout=f, cwd=outputdir)
                # Wait for BOLSIG- to complete
                Popen.wait(r)

    return


def morris_stats(file_list, sample_set, delta, num_EN_val, num_params):

    logging.debug('Mobility statistical calculations have begun')

    # Initial Matrices and values
    cwd = getcwd()
    mob_start = 0
    diff_start = 0
    ion_start = 0
    excite_1 = 0
    excite_2 = 0
    excite_3 = 0
    ex_pass = 100000000
    coef_names = ['Mobility', '\nDiffusion Coef.', '\nIonisation coef.', '\n1st Excitation', '\n2nd Excitation',
                  '\n3rd Excitation']

    mean = {}
    std = {}
    energy = zeros(num_EN_val)

    # Interating through each created directory in Input Files
    for keys, num in zip(file_list, sample_set):

        # Entering each directory and setting up list of input files and corresponding morris trajectories
        chdir(join(cwd, keys))
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
                    if file_1[j].find('Mobility') != -1:
                        mob_start = j

                    elif file_1[j].find('Diffusion coefficient') != -1:
                        diff_start = j

                    elif file_1[j].find('Total ionization freq.') != -1:
                        ion_start = j
                        ex_pass = j

                    elif file_1[j].find('C2') != -1 and j > ex_pass:
                        excite_1 = j
                        print(file_1[excite_1])

                    elif file_1[j].find('C3') != -1 and j > ex_pass:
                        excite_2 = j

                    elif file_1[j].find('C4') != -1 and j > ex_pass:
                        excite_3 = j
                        break

                    else:
                        continue

            for stat, stat_name in zip([mob_start, diff_start, ion_start, (excite_1+1), (excite_2+1), (excite_3+1)],
                                       coef_names):

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

    # Final calculation for STD
    for name in coef_names:
        std[name] = sqrt(divide(std[name], (len(sample_set) - 1)))

    return mean, std, energy, coef_names


def print_stats(mean_stats, std_stats, energy, cwd, dataset, coef_names):

    chdir(cwd)

    # Creating a file contain the calculated data
    with open('stats.txt', 'w', encoding='utf8') as f:

        # Iterate over all transport coefficients
        for name in coef_names:

            mean = mean_stats[name]
            std = std_stats[name]

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


def normal_stats(mean_stats, std_stats, energy, dataset, coef_names):

    # Creating a file contain the calculated data
    with open('norm_stats.txt', 'w', encoding='utf8') as f:

        # Iterate over all transport coefficients
        for name in coef_names:

            mean = mean_stats[name]
            std = std_stats[name]

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
