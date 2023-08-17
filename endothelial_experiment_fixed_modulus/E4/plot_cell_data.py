import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
# sns.set_theme()
sns.set(rc={'figure.figsize':(11.7,8.27)})


def cell_count(file):
    fig, axs = plt.subplots(ncols=1, nrows=2)
    plt.subplots_adjust(left = 0.07, bottom = 0.075, right=0.800, top=0.96)


    data = pd.read_csv(file, sep=",",usecols = ["mcsteps", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"])
    # print(data)
    # Load an example dataset with long-form data
    # fmri = sns.load_dataset("fmri")
    cell_count_labels = ["Endothelial", "Resting Neutrophils", "Monocytes", "Fibroblasts", 
    "Activated Neutrophils", "NDN Neutrophils","Resting Monocytes",
     "Macrophages I", "Macrophages II", "Myofibroblasts"]
    cols = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
    colors = ["blue", "brown", "cyan", "violet", "red", "pink", "yellow", "orange", "darkblue", "green"]

    # scoped plot need 1, 3, 4, 5, 6, 7, 8, 10
    scoped_cols = ["1", "3", "4", "5", "6", "7", "8", "10"]
    scoped_colors = ["blue", "cyan", "violet", "red", "pink", "yellow", "orange", "green"]
    scoped_labels = ["Endothelial",  "Monocytes", "Fibroblasts", 
    "Activated Neutrophils", "NDN Neutrophils","Resting Monocytes",
     "Macrophages I", "Myofibroblasts"]

    for i in range(0, len(cols)):
    # Plot the responses for different events and regions
        sns.lineplot(x="mcsteps", y= data[cols[i]],
                     data=data, color = colors[i], label = cell_count_labels[i], ax = axs[0])
    for i in range(0, len(scoped_cols)):
        sns.lineplot(x="mcsteps", y= data[scoped_cols[i]],
                         data=data, color = scoped_colors[i], label = scoped_labels[i], ax = axs[1])

    hour_ticks = np.arange(-20, 120, 20)
    # xlabels[1] = 'Testing'
    xlabels = [r"{}".format(i) for i in hour_ticks]
    print(xlabels)

    axs[0].set_xticklabels(xlabels)
    axs[1].set_xticklabels(xlabels)

    axs[1].set_ylim([0, 400])
    axs[0].set_ylabel("Cell count")
    axs[1].set_ylabel("Cell count")
    axs[0].set_xlabel("Hours of simulation")
    axs[1].set_xlabel("Hours of simulation")
    axs[0].legend(loc='center left', bbox_to_anchor=(1, 0.5))
    axs[1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
    fig.savefig("cell_count_data.png", dpi=300)
    plt.show()
    return 


def cyt_concentrations(file):
    fig, axs = plt.subplots(ncols=1, nrows=2)
    plt.subplots_adjust(left = 0.07, bottom = 0.075, right=0.900, top=0.96)


    data = pd.read_csv(file, sep=",",usecols = ["meanconcen","il8mean","il1mean",
        "il6mean","il10mean","tnfmean",
        "tgfmean","il8std","il1std","il6std",
        "il10std","tnfstd","tgfstd"])
    # print(data["meanconcen"][:,0])
    cytokine_labels = ["IL-8", "IL-1","IL-6", "IL-10", "TNF", 
    "TGF"]
    cols = ["il8mean","il1mean",
        "il6mean","il10mean","tnfmean",
        "tgfmean"]
    fill_cols = ["il8std","il1std","il6std",
        "il10std","tnfstd","tgfstd"]
    colors = ["blue", "brown", "violet", "red", "yellow", "orange", "green"]

    # scoped plot need "IL-8", "IL-1","IL-6", "IL-10", "TNF", 
    scoped_cols = ["il8mean","il1mean",
        "il6mean","il10mean","tnfmean"]
    scoped_fill_cols = ["il8std","il1std","il6std",
        "il10std","tnfstd"]
    scoped_colors = ["blue", "cyan", "violet", "red", "pink", "yellow","green"]
    scoped_labels = ["IL-8", "IL-1","IL-6", "IL-10", "TNF"]

    for i in range(0, len(cols)):
    # Plot the responses for different events and regions
        sns.lineplot(x="meanconcen", y= data[cols[i]],
                     data=data, color = colors[i], label = cytokine_labels[i], ax = axs[0])
        

    for i in range(0, len(scoped_cols)):
        sns.lineplot(x="meanconcen", y= data[scoped_cols[i]],
                         data=data, color = scoped_colors[i], label = scoped_labels[i], ax = axs[1])
        # axs[1].fill_between(data["meanconcen"].values, data[scoped_cols[i]].values-scoped_fill_cols[cols[i]].values, 
        #     data[scoped_cols[i]].values+scoped_fill_cols[scoped_cols[i]].values,color = scoped_colors[i], alpha=0.2)

    hour_ticks = np.arange(-20, 120, 20)
    xlabels = [r"{}".format(i) for i in hour_ticks]

    axs[0].set_xticklabels(xlabels)
    axs[1].set_xticklabels(xlabels)

    axs[1].set_ylim([0, 0.5*10**(-9)])
    axs[0].set_ylabel("Concentration")
    axs[1].set_ylabel("Concentration")
    axs[0].set_xlabel("Hours of simulation")
    axs[1].set_xlabel("Hours of simulation")
    axs[0].legend(loc='center left', bbox_to_anchor=(1, 0.5))
    axs[1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
    fig.savefig("cytokine_concentration.png", dpi=300)
    plt.show()
    return


cell_count('combi_code/datafiles/cellcount.txt')
# cyt_concentrations("mean_concentration.txt")
