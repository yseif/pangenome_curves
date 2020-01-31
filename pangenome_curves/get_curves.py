
from matplotlib import pyplot as plt
%matplotlib inline
import scipy
import pandas as pd
from glob import glob
from scipy.optimize import curve_fit

def get_pangenome_counts(input_table, bootstrap_no):

    core_rows = defaultdict(dict)
    pan_rows = defaultdict(dict)

    for i in range(bootstrap_no):
        gid_list = np.array(input_table.columns)
        np.random.shuffle(gid_list)

        core = []
        pan = []
        count = 0

        for gid in gid_list:
            if core == []:
                core = set(input_table.loc[input_table[gid] != 0].index)
            else:
                core = set(input_table.loc[input_table[gid] != 0].index) & set(core)
            if pan == []:
                pan = np.array(input_table.loc[input_table[gid] != 0].index)
            else:
                pan = np.union1d(np.array(input_table.loc[input_table[gid] != 0].index), pan)

            core_rows['try_%d'%i][count] = len(core)
            pan_rows['try_%d'%i][count] = len(pan)
            count += 1
            
    core_matrix = pd.DataFrame(core_rows)
    pan_matrix = pd.DataFrame(pan_rows)
    
    pangenome_res = pd.DataFrame({'Core genome (mean)':core_matrix.apply(np.mean, axis = 1), 'Core genome (std)':core_matrix.apply(np.std, axis = 1), 
 'Pan genome (mean)':pan_matrix.apply(np.mean, axis = 1), 'Pan genome (std)':pan_matrix.apply(np.std, axis = 1)})
    return pangenome_res


def plot_pangenome_curves(pangenome_res, ax, plot_std = False):

    xrange = pangenome_res.index

    ax.plot(xrange, pangenome_res['Core genome (mean)'], 'maroon', label = 'Core genome')

    ax.plot(xrange, pangenome_res['Pan genome (mean)'], 'midnightblue', label = 'Pan genome')

    ax.legend()

    if plot_std == True:

        ax.plot(xrange, pangenome_res['Core genome (mean)'] + pangenome_res['Core genome (std)'], 'maroon', linestyle = 'dashed', linewidth = 0.8)
        ax.plot(xrange, pangenome_res['Core genome (mean)'] - pangenome_res['Core genome (std)'], 'maroon', linestyle = 'dashed', linewidth = 0.8)

        ax.plot(xrange, pangenome_res['Pan genome (mean)'] + pangenome_res['Pan genome (std)'], 'midnightblue', linestyle = 'dashed', linewidth = 0.8)
        ax.plot(xrange, pangenome_res['Pan genome (mean)'] - pangenome_res['Pan genome (std)'], 'midnightblue', linestyle = 'dashed', linewidth = 0.8)    


    ax.set_ylabel('No. genes')
    ax.set_xlabel('No. strains')  
    
    return ax
    

def heaps_func(X, k, gamma):
    return k*X[0]**gamma+ X[1]

def return_heaps_estimates(matrix):
    popt, pcov = curve_fit(heaps_func, [range(len(pangenome_res)),pangenome_res['Pan genome (mean)'].iloc[0]], pangenome_res['Pan genome (mean)'])
    stds = np.sqrt(np.diag(pcov))
    heaps_estimates = pd.DataFrame(list({'K estimate':round(popt[0],3), 'Gamma estimate': round(popt[1],3), 'K standard deviation':round(stds[0], 3), 'Gamma standard deviation': round(stds[1], 3)}.items())).set_index(0).T
    return heaps_estimates