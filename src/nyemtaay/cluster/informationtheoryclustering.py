# -*- coding: utf-8 -*-
##############################################################################
## Copyright (c) 2022 Adrian Ortiz-Velez.
## All rights reserved.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer.
##     * Redistributions in binary form must reproduce the above copyright
##       notice, this list of conditions and the following disclaimer in the
##       documentation and/or other materials provided with the distribution.
##     * The names of its contributors may not be used to endorse or promote
##       products derived from this software without specific prior written
##       permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
## ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
## WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL ADRIAN ORTIZ-VELEZ BE LIABLE FOR ANY DIRECT,
## INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
## BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
## LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
## OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
## ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
##############################################################################
#
# Library for clustering based on information theory.
#
###########################################################################
## Reader Implementation

# Libraries

from community import community_louvain
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import networkx as nx

def louvian_clustering(G,metadata): #pass metadata
#    from mpl_toolkits.basemap import Basemap
    
    print('clustering')

    #################### compute the best partition
    partition = community_louvain.best_partition(G)

    # draw the graph
    pos = nx.circular_layout(G)
    
    # draw
    nx.draw_networkx(G,pos,node_size=200)

        
    # color the nodes according to their partition
    cmap = cm.get_cmap('Set3', max(partition.values()) + 1)
    nx.draw_networkx_nodes(G, pos, partition.keys(), node_size=500,
                           cmap=cmap, node_color=list(partition.values()))
    nx.draw_networkx_edges(G, pos, alpha=0.5)
    plt.show()
    
#    # load the karate club graph
#    G = nx.karate_club_graph()

#    # compute the best partition
#    partition = community_louvain.best_partition(G)

#    # draw the graph
#    pos = nx.draw_planar(G)
#    # color the nodes according to their partition
#    cmap = cm.get_cmap('viridis', max(partition.values()) + 1)
#    nx.draw_networkx_nodes(G, pos, partition.keys(), node_size=40,
#                           cmap=cmap, node_color=list(partition.values()))
#    nx.draw_networkx_edges(G, pos, alpha=0.5)
#    plt.show()
    return 

def louvian_clustering_overmap(G,metadata): #pass metadata
    from mpl_toolkits.basemap import Basemap
    
    print('clustering')
    min_lat = metadata['LAT'].min()
    min_long = metadata['LONG'].min()
    max_lat = metadata['LAT'].max()
    max_long = metadata['LONG'].max()
    
    
    m = Basemap(projection='mill',
                llcrnrlat = min_lat,
                llcrnrlon = min_long,
                urcrnrlat = max_lat,
                urcrnrlon = max_long,
                resolution='l')


    plt.title('Analysis plotted over map')

    #################### compute the best partition
    partition = community_louvain.best_partition(G)

    # draw the graph
    pos = nx.spring_layout(G)
    
    # draw
    nx.draw_networkx(G,pos,node_size=200)

        
    # color the nodes according to their partition
    cmap = cm.get_cmap('Set3', max(partition.values()) + 1)
    nx.draw_networkx_nodes(G, pos, partition.keys(), node_size=50,
                           cmap=cmap, node_color=list(partition.values()))
    nx.draw_networkx_labels(G, pos, ax=ax)    
    nx.draw_networkx_edges(G, pos, alpha=0.5)


    m.drawcoastlines()
    m.drawcountries(linewidth=2)
    ##m.drawstates(color='b')
    ##m.drawcounties(color='darkred')
    #m.fillcontinents()
    #m.etopo()
    m.bluemarble()
    plt.show()
    return 
    
    
def dbscan_imp(X):
    import numpy as np

    from sklearn import metrics
    from sklearn.cluster import DBSCAN
    
    db = DBSCAN(eps=0.2, min_samples=1,metric='precomputed').fit(X) 
    labels = db.labels_

    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0) 
    n_noise_ = list(labels).count(-1)

    print("Estimated number of clusters: %d" % n_clusters_)
    print(labels)
    print("Estimated number of noise points: %d" % n_noise_)
    
    return
    
    
    
