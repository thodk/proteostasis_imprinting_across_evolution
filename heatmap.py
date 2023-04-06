# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas
import numpy
from scipy.spatial import distance
from scipy.cluster import hierarchy
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerLine2D
import matplotlib.lines
from matplotlib import pyplot as plt
import seaborn as sns



class SymHandler(HandlerLine2D):
    def create_artists(self, legend, orig_handle,xdescent, ydescent, width, 
                       height, fontsize, trans):
        xx= self.y_align*height
        return super(SymHandler, self).create_artists(legend, orig_handle,xdescent, 
                     xx, width, height, fontsize, trans)
        

class Clustergram(object):

    def __init__(self, heatmap_data, row_clustering_data=None, row_colors_data=None,
		col_clustering_data=None, col_colors_data=None, x_axis_data=None,
		y_axis_data=None, cmap=None, norm=None, vmax = None, vmin=None, colors_legend_data=None,
		legend_data=None, linecolor='#333333', linewidth=1, figsize=(4,4)):

        self.heatmap_data = heatmap_data
        self.row_clustering_data = row_clustering_data
        self.row_colors_data = row_colors_data
        self.col_clustering_data = col_clustering_data
        self.col_colors_data = col_colors_data
        self.colors_legend_data = colors_legend_data
        self.x_axis_data = x_axis_data
        self.y_axis_data = y_axis_data
        self.cmap = cmap
        self.norm = norm
        self.vmax = vmax        
        self.vmin = vmin
        self.figsize = figsize
        self.legend_data = legend_data
        self.linecolor = linecolor
        self.linewidth = linewidth


    def construction(self):
        
	#if self.heatmap_data['type'] == 'distances':
	    #condensed_matrix = distance.squareform(self.heatmap_data['data'], force='to_vector')
	#else: # features
	    #condensed_matrix = distance.pdist(self.heatmap_data['data'])


        # define the row_dendrogram
        if self.row_clustering_data == None:
            self.row_clustering_data = None
            row_cluster=False
            row_linkage=None
        elif type(self.row_clustering_data['data']) == list: #Z
            row_cluster=True
            row_linkage = self.row_clustering_data['data']
        else:
            row_cluster=True
            if self.heatmap_data['type'] == 'distances':
                tmp_condensed_matrix = distance.squareform(self.row_clustering_data['data'], force='to_vector')
            else:
                tmp_condensed_matrix = distance.pdist(self.row_clustering_data['data'])
            row_linkage = hierarchy.linkage(y=tmp_condensed_matrix, method=self.row_clustering_data['method'],
                                            metric=self.row_clustering_data['metric'], optimal_ordering=True)
        # define row_colors
        if self.row_colors_data == None:
            self.row_colors_data = None
            row_colors = None
        else:
            row_colors = self.row_colors_data['colors']
        
        # define the col_dendrogram
        if self.col_clustering_data == None:
            self.col_clustering_data = None
            col_cluster=False
            col_linkage=None
        elif type(self.col_clustering_data['data']) == list: #Z
            col_cluster=True
            col_linkage = self.col_clustering_data['data']
        else:
            col_cluster=True
            if self.heatmap_data['type'] == 'distances':
                tmp_condensed_matrix = distance.squareform(self.col_clustering_data['data'], force='to_vector')
            else:
                tmp_condensed_matrix = distance.pdist(self.col_clustering_data['data'])
            col_linkage = hierarchy.linkage(y=tmp_condensed_matrix, method=self.col_clustering_data['method'],
                                            metric=self.col_clustering_data['metric'], optimal_ordering=True)
        
        # define col_colors
        if self.col_colors_data == None:
            self.col_colors_data = None
            col_colors = None
        else:
            col_colors = self.col_colors_data['colors']        

        # construct self.clustergrid
        try:
            cbar_kws=self.legend_data['cbar_kws']
        except (TypeError, KeyError):
            cbar_kws = None
                    

        self.clustergrid = sns.clustermap(self.heatmap_data['data'], 
                                          row_cluster=row_cluster, 
                                          row_linkage=row_linkage, 
                                          row_colors=row_colors,
                                          col_cluster=col_cluster, 
                                          col_linkage=col_linkage, 
                                          col_colors=col_colors, 
                                          xticklabels=True if self.x_axis_data !=None else False, 
                                          yticklabels=True if self.y_axis_data !=None else False,
                                          cmap=self.cmap, norm=self.norm,
                                          vmax=self.vmax, vmin=self.vmin,
                                          figsize=self.figsize, linecolor=self.linecolor,
                                          linewidth=self.linewidth, 
                                          cbar_kws=cbar_kws)


    def sizes_normalization(self, sizes):
        sizes = numpy.array(sizes)
        return sizes/sum(sizes)*.8

    
    def set_coords(self):
        coords = ['x','y','w','h']
        step=0.005

        if self.row_clustering_data == None:
            row_den_coords = dict(zip(coords, numpy.zeros(4)))
        else:
            row_den_coords = dict(zip(coords, self.clustergrid.ax_row_dendrogram.get_position().bounds))
            row_den_coords['w'] = self.row_clustering_data['ratio']
        if self.row_colors_data == None:
            row_col_coords = dict(zip(coords, numpy.zeros(4)))
        else:
            row_col_coords = dict(zip(coords, self.clustergrid.ax_row_colors.get_position().bounds))
            row_col_coords['w'] = self.row_colors_data['ratio']
        if self.col_clustering_data == None:
            col_den_coords = dict(zip(coords, numpy.zeros(4)))
        else:
            col_den_coords = dict(zip(coords, self.clustergrid.ax_col_dendrogram.get_position().bounds))
            col_den_coords['h'] = self.col_clustering_data['ratio']
        if self.col_colors_data == None:
            col_col_coords = dict(zip(coords, numpy.zeros(4)))
        else:
            col_col_coords = dict(zip(coords, self.clustergrid.ax_col_colors.get_position().bounds))
            col_col_coords['h'] = self.col_colors_data['ratio']
        heatmap_coords = dict(zip(coords, self.clustergrid.ax_heatmap.get_position().bounds))
        heatmap_coords['w'] = self.heatmap_data['x_ratio']
        heatmap_coords['h'] = self.heatmap_data['y_ratio']
        
        
        x_sizes = self.sizes_normalization([row_den_coords['w'], row_col_coords['w'], heatmap_coords['w']])
        y_sizes = self.sizes_normalization([col_den_coords['h'], col_col_coords['h'], heatmap_coords['h']])

        x_margin = 0.2
        y_margin = 0.2   
        
        # place row dendrogram 
        row_den_coords['x'] = x_margin
        if x_sizes[0] == 0:
            row_den_coords['w'] = 0
        else:
            row_den_coords['w'] = x_sizes[0] - step
        row_den_coords['y'] = 1 - (y_margin + y_sizes[0] + y_sizes[1] + y_sizes[2]) - step
        row_den_coords['h'] = y_sizes[2]
        #if self.row_clustering_data != None:
        self.clustergrid.ax_row_dendrogram.set_position(map(lambda i: row_den_coords[i], coords))

        # place row colors
        row_col_coords['x'] = x_margin + row_den_coords['w'] + step
        row_col_coords['w'] = x_sizes[1]
        row_col_coords['y'] = 1 - (y_margin + y_sizes[0] + y_sizes[1] + y_sizes[2]) - step
        row_col_coords['h'] = y_sizes[2]
        if self.row_colors_data != None:
            self.clustergrid.ax_row_colors.set_position(map(lambda i: row_col_coords[i], coords))
        
        # place col dendrogram
        col_den_coords['x'] = x_margin + row_den_coords['w'] + row_col_coords['w'] + step
        col_den_coords['w'] = x_sizes[2]
        col_den_coords['y'] = 1 - (y_margin + y_sizes[0])
        if y_sizes[0] == 0:
            col_den_coords['h'] = 0
        else:
            col_den_coords['h'] = y_sizes[0] - step
        #if self.col_clustering_data != None:
        self.clustergrid.ax_col_dendrogram.set_position(map(lambda i: col_den_coords[i], coords))
        
        # place col colors
        col_col_coords['x'] = x_margin + row_den_coords['w'] + row_col_coords['w'] + step
        col_col_coords['w'] = x_sizes[2]
        col_col_coords['y'] = 1 - (y_margin + y_sizes[0] + y_sizes[1]) - step
        col_col_coords['h'] = y_sizes[1]
        if self.col_colors_data != None:
            self.clustergrid.ax_col_colors.set_position(map(lambda i: col_col_coords[i], coords))        
        
        heatmap_coords['x'] = x_margin + row_den_coords['w'] + row_col_coords['w'] + step
        heatmap_coords['w'] = x_sizes[2]
        heatmap_coords['y'] = 1 - (y_margin + y_sizes[0] + y_sizes[1] + y_sizes[2]) - step
        heatmap_coords['h'] = y_sizes[2]
        
        self.clustergrid.ax_heatmap.set_position(map(lambda i: heatmap_coords[i], coords))        

        try:
            x = self.legend_data['x']
            y = self.legend_data['y']
            w = self.legend_data['w']
            h = self.legend_data['h']
            self.clustergrid.cax.set_position([x, y, w, h])
            self.clustergrid.cax.set_yticklabels(labels=self.legend_data['labels'], 
                                                 fontdict={'size':self.legend_data['labels_size']})        
            self.clustergrid.cax.set_title(label=self.legend_data['title'], 
                                           fontdict={'size':self.legend_data['title_size']})
        except TypeError:
            pass
        
        
        if self.colors_legend_data != None:
            SH = SymHandler()
            try:
                SH.y_align = self.colors_legend_data['y_align']
            except KeyError:
                SH.y_align = 0.6
            if type(self.colors_legend_data['patches']) == dict:

                leg = self.clustergrid.fig.legend(handles=self.colors_legend_data['patches'][1], 
                                                  labels=self.colors_legend_data['patches'][0], 
                                                  bbox_to_anchor=self.colors_legend_data['bbox_anchor'],
                                                  fontsize=self.colors_legend_data['fontsize'], 
                                                  ncol=1, fancybox=False, framealpha=0, 
                                                  handlelength=self.colors_legend_data['handlelength'], 
                                                  handletextpad=self.colors_legend_data['handletextpad'],
                                                  handleheight=self.colors_legend_data['handleheight'],
                                                  title=self.colors_legend_data['title'],
                                                  handler_map={matplotlib.lines.Line2D: SH})
            else:
                leg = self.clustergrid.fig.legend(handles=self.colors_legend_data['patches'][1], 
                                                  labels=self.colors_legend_data['patches'][0], 
                                                  bbox_to_anchor=self.colors_legend_data['bbox_anchor'],
                                                  fontsize=self.colors_legend_data['fontsize'], 
                                                  ncol=1, fancybox=False, framealpha=0, 
                                                  handlelength=self.colors_legend_data['handlelength'], 
                                                  handletextpad=self.colors_legend_data['handletextpad'],
                                                  handleheight=self.colors_legend_data['handleheight'],
                                                  title=self.colors_legend_data['title'],
                                                  handler_map={matplotlib.lines.Line2D: SH})                
            leg._legend_box.align = "left"
            plt.setp(leg.get_title(), fontsize=self.colors_legend_data['title_size'],
                     multialignment='center')
        
        
    def set_labels(self):
        
        if self.x_axis_data != None:
            labels = self.x_axis_data['labels']
            try:
                fontdict = self.x_axis_data['fontdict']
            except KeyError:
                fontdict = None
            try:
                rotation = self.x_axis_data['rotation']
            except KeyError:
                rotation = None

            sorted_labels = []
            for i in self.clustergrid.ax_heatmap.get_xticklabels():
                index = int(i.get_text())
                label = labels[index]
                sorted_labels.append(label)
            self.clustergrid.ax_heatmap.set_xticklabels(labels=sorted_labels, rotation=rotation, fontdict=fontdict)
            
            try:
                self.x_axis_data['specific_labels_format']
                for index, label_obj in enumerate(self.clustergrid.ax_heatmap.get_xticklabels()):
                    label = label_obj.get_text()
                    try:
                        for key, value in self.x_axis_data['specific_labels_format'][label].items():
                            if key == 'color':
                                self.clustergrid.ax_heatmap.get_xticklabels()[index].set_color(value)
                            elif key == 'weight':
                                self.clustergrid.ax_heatmap.get_xticklabels()[index].set_weight(value)
                            else:
                                pass
                    except KeyError:
                        continue
            except KeyError:
                pass



        if self.y_axis_data != None:
            labels = self.y_axis_data['labels']
            try:
                fontdict = self.y_axis_data['fontdict']
            except KeyError:
                fontdict = None
            try:
                rotation = self.y_axis_data['rotation']
            except KeyError:
                rotation = None

            sorted_labels = []
            for i in self.clustergrid.ax_heatmap.get_yticklabels():
                index = int(i.get_text())
                label = labels[index]
                sorted_labels.append(label)
            self.clustergrid.ax_heatmap.set_yticklabels(labels=sorted_labels, rotation=rotation, fontdict=fontdict)
            try:
                self.y_axis_data['specific_labels_format']
                for index, label_obj in enumerate(self.clustergrid.ax_heatmap.get_yticklabels()):
                    label = label_obj.get_text()
                    try:
                        for key, value in self.y_axis_data['specific_labels_format'][label].items():
                            if key == 'color':
                                self.clustergrid.ax_heatmap.get_yticklabels()[index].set_color(value)
                            elif key == 'weight':
                                self.clustergrid.ax_heatmap.get_yticklabels()[index].set_weight(value)
                            else:
                                pass
                    except KeyError:
                        continue
            except KeyError:
                pass
