#!/usr/bin/python
# -*- coding: utf-8 -*-
import pandas
import numpy
import sys
sys.path.append('../')
from core_functions import taxonomy_colors
from sklearn.manifold import MDS
from scipy.spatial import distance
from scipy.stats import pearsonr
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerLine2D
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.lines
import matplotlib
matplotlib.rcParams['font.sans-serif'] = ['FreeSans', ]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['axes.titlepad'] = 4
matplotlib.rcParams['xtick.major.pad']='0.3'
matplotlib.rcParams['ytick.major.pad']='0.3'


class SymHandler(HandlerLine2D):
    def create_artists(self, legend, orig_handle,xdescent, ydescent, width, 
                       height, fontsize, trans):
        xx= self.y_align*height
        return super(SymHandler, self).create_artists(legend, orig_handle,xdescent, 
                     xx, width, height, fontsize, trans)


def mds_worker(filename):
    mds_obj = MDS(dissimilarity='precomputed', n_components=2, n_jobs=4, 
                  random_state=100, eps=1e-6, n_init=10) # -6 10
    df = pandas.read_csv(filename, index_col=0, sep='\t')
    colors = [taxonomy_colors[species_dict[i]] for i in df.index.tolist()]
    coords = mds_obj.fit_transform(df.values)
    return coords, colors


def worker_for_scatter_plot(ax, coords, colors, title, lim):
    ax.grid(True, linestyle='--', alpha=0.5, linewidth=0.3)
    ax.scatter(coords[:,0], coords[:,1], c=colors, s=5, linewidths=0., 
               edgecolor='black', alpha=0.85)
    ax.set_xlim([-lim,lim])
    ax.set_ylim([-lim,lim])
    plt.xticks(scatter_plots_ticks, size=6)
    plt.xlabel('MDS Component 1', size=8, labelpad=1)
    plt.yticks(scatter_plots_ticks, size=6)
    plt.ylabel('MDS Component 2', size=8, labelpad=1)
    plt.title(title, size=10)
    plt.tick_params(labelsize=6, pad=0.5, width=0.3, length=1.5, direction='out')




if __name__ == '__main__':

    '''
    DESCRIPTION:

    This script creates the 2nd figure of the manuscript. It consists of two
    separate subplots. The first one depicts the species of each kingdom
    for each criterion on a 2-dimensional space, constructed using the MDS method.
    This dimensionality reduction methods uses the distances of species to map them
    on an orthogonal space, with respect to the initial distance values. The 
    second subplot concerns the correlation of distances between rRNA and all
    the other criteria. An equal number of species from each kingdom is used (60)
    for the calculation of pearson's correlation, in order to avoid to introduce 
    bias because of the bacterial set, as it contains significantly more species 
    than the other two.
    '''

    filenames = {
        'rRNA':'./rRNAs_analysis/final_distance_matrix.tsv',
        'hsp40':'./HSPs_analysis/hsp40/final_distance_matrix.tsv',
        'hsp70':'./HSPs_analysis/hsp70/final_distance_matrix.tsv',
        'proteostasis_network':'./PN_analysis/final_distance_matrix.tsv'
    }    
    titles = {
        'rRNA': 'Ribosomal Sequences',
        'hsp40': 'HSP40 Sequences',
        'hsp70': 'HSP70 Sequences',
        'proteostasis_network': 'Proteostasis Network'
    }

    '''
    STEP 1: Construct the scatter plots with MDS
    '''
    species_df = pandas.read_csv('./files/species.tsv', sep='\t', index_col=0)
    species_dict = dict(zip(species_df.species_name, species_df.taxonomy))   
    taxonomy_groups = species_df.groupby(by='taxonomy')
    numpy.random.seed(4234)
    selected_species = []
    for taxonomy, tmp_df in taxonomy_groups:
        random_species = numpy.random.choice(tmp_df.species_name, 60)
        selected_species.extend(random_species)
    selected_species = sorted(selected_species)

    # Create the MDS matrices
    mds_coord_matrices = {}
    for method, filename in list(filenames.items()):
        coords, colors = mds_worker(filename)
        mds_coord_matrices.update({method:{'coords':coords, 'colors':colors}})

    # Create the framework of the plot with GridSpec
    fig = plt.figure(1, figsize=(7,4.66), frameon=False)
    grids = GridSpec(22,32, hspace=0.0, wspace=0.0)
    scatter_plots_lim = 0.85
    scatter_plots_ticks = numpy.arange(-0.8, 0.8+0.1, 0.4)    

    # Plot for rRNAs distances
    ax = fig.add_subplot(grids[0:9,0:9])
    worker_for_scatter_plot(ax, mds_coord_matrices['rRNA']['coords'], 
                            mds_coord_matrices['rRNA']['colors'],
                            titles['rRNA'], scatter_plots_lim)
    ax.text(x=-1.23, y=0.9, s='A', fontsize=16, fontweight=500)
    #plt.setp(ax.get_xticklabels(), visible=False)
    #plt.setp(ax.get_yticklabels(), visible=False)

    # Plot for PN distances    
    ax = fig.add_subplot(grids[0:9,12:21])
    worker_for_scatter_plot(ax, mds_coord_matrices['proteostasis_network']['coords'],
                            mds_coord_matrices['proteostasis_network']['colors'],
                            titles['proteostasis_network'], scatter_plots_lim)
    
    # Plot for hsp40 distances
    ax = fig.add_subplot(grids[13:22,0:9])
    worker_for_scatter_plot(ax, mds_coord_matrices['hsp40']['coords'], 
                            mds_coord_matrices['hsp40']['colors'], 
                            titles['hsp40'], scatter_plots_lim)
    
    # Plot for hsp70 distances      
    ax = fig.add_subplot(grids[13:22,12:21])
    worker_for_scatter_plot(ax, mds_coord_matrices['hsp70']['coords'],
                            mds_coord_matrices['hsp70']['colors'], 
                            titles['hsp70'], scatter_plots_lim)
    
    patches= [[], []]
    for taxonomy, color in taxonomy_colors.items():
        patches[0].append(taxonomy.capitalize())
        p = Line2D([0], [0], marker='o', color=color, markerfacecolor=color, 
                   markersize=12, label=taxonomy.capitalize())
        patches[1].append(p)
    SH = SymHandler()
    SH.y_align = 0.6
    leg = plt.legend(handles=patches[1], labels=patches[0], bbox_to_anchor=(0.85,-0.15),
               fontsize=10, ncol=3, frameon=False, fancybox=False,  handlelength=0.75, 
               handletextpad=0.5, handleheight=3, title='', labelspacing=0, 
               borderpad=0.5, handler_map={matplotlib.lines.Line2D: SH})
    leg._legend_box.align = "left"
    plt.setp(leg.get_title(), fontsize=6, multialignment='center')

    '''
    STEP 2: Construct the correlation density plots
    '''
    df = pandas.read_csv(filenames['rRNA'], sep='\t', index_col=0)
    df = df.loc[selected_species, selected_species]
    reference_values = distance.squareform(df.values, force='to_vector')
    
    corr_plots_ticks = list(i/10. for i in range(0,11,2))
    rows = range(23)

    for i, method in enumerate(['proteostasis_network', 'hsp40', 'hsp70']):
    
        filename = filenames[method]
        tmp_df = pandas.read_csv(filename, sep='\t', index_col=0)
        tmp_df = tmp_df.loc[selected_species, selected_species]
        tmp_values = distance.squareform(tmp_df.values, force='to_vector')
        
        corr = pearsonr(reference_values, tmp_values)[0]
        tmp_matrix = numpy.array([reference_values, tmp_values]).transpose()

        pearson_correlation = "% 0.2f" %round(corr, 2)
        cmap1 = LinearSegmentedColormap.from_list("my_colormap", ('#ffffff','#0b3c5d','#072a41','#020c12'), 
                                                  N=25, gamma=1.0)
        contour_axes = fig.add_subplot(grids[rows[1]:rows[6],25:30])
        if i == 0 :
            contour_axes.text(x=-0.85, y=1.4, s='B', fontsize=16, fontweight=500)
            
        sns.kdeplot(tmp_matrix[:,1], tmp_matrix[:,0], cmap=cmap1, n_levels=25, 
                    shade=True, bw_method=0.25, ax=contour_axes)
    
        contour_axes.tick_params(labelsize=0, pad=0.5, width=0, length=0.1, direction='out')
        contour_axes.grid(linestyle='--', alpha=0.5, linewidth=0.3)
        contour_axes.set_frame_on(True)
        contour_axes.text(0.0, 0.95, r"Pearson's r: "+pearson_correlation, fontsize=6)
        contour_axes.set_xlim([-0.1,1.1])
        contour_axes.set_xticks(corr_plots_ticks)
        contour_axes.set_xticklabels(contour_axes.get_xticks(), size=5)
        contour_axes.set_ylim([-0.1,1.1])
        contour_axes.set_yticks(corr_plots_ticks)
        contour_axes.set_yticklabels(contour_axes.get_yticks(), size=5)
    
        t1 = 'Distances of '+titles[method]
        t1 = ' '.join(t1.split(' ')[0:-1]) + '\n' + t1.split(' ')[-1]
        contour_axes.set_xlabel(t1, fontsize=7, labelpad=2)
        t2 = 'Distances of '+ titles['rRNA']
        t2 = ' '.join(t2.split(' ')[0:-1]) + '\n' + t2.split(' ')[-1]
        contour_axes.set_ylabel(t2, fontsize=7, labelpad=2)
    
        y_density_axes = fig.add_subplot(grids[rows[1]:rows[6], 30])
        sns.kdeplot(tmp_matrix[:,0], color="#08304a", shade=True, bw=0.1, 
                    vertical=True, linewidth=1, ax=y_density_axes)
        y_density_axes.axes.axvline(x=0, linewidth=1, color='black')
        y_density_axes.axes.axis('off')
        y_density_axes.set_ylim([-0.1,1.1])
        y_density_axes.set_frame_on(False)
    
        x_density_axes = fig.add_subplot(grids[rows[0], 25:30])
        sns.kdeplot(tmp_matrix[:,1], color="#08304a", shade=True, bw=0.1, 
                    vertical=False, linewidth=1, ax=x_density_axes)
        x_density_axes.axes.axhline(y=0, linewidth=1, color='black')
        x_density_axes.axes.axis('off')
        x_density_axes.set_xlim([-0.1,1.1])
        x_density_axes.set_frame_on(False)
        
        rows = rows[8:]

    plt.savefig('./comparison_of_clustergrams/figure_2.png', 
                dpi=800, format='png', bbox_inches='tight', pad_inches=0.05)

    #plt.savefig('./comparison_of_clustergrams/figure_2.pdf', 
    #            dpi=600, format='pdf', bbox_inches='tight', pad_inches=0.05)


    


