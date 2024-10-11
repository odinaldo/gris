# GRIS v1.0
#
# This simple Python script implements the
# Gabbay-Rodrigues Iteration Schema, as described in
# Gabbay, D.M. and Rodrigues, O., Equilibrium States in
# Numerical Argumentation Networks, Logica Universalis,
# pp. 1-63, DOI: 10.1007/s11787-015-0119-7, 2015. 

# If you do not want too see the argumentation frameworks,
# comment ouu the networkx libraries below and the function
# plot_graph

# This program is for illustrative purposes only, has not been
# throroughly tested, and comes with no guarantees
#
# Run it with python3 gris.py [tgffile.tgf]

import os.path
import sys
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout

plt.rcParams["figure.figsize"] = [12, 4]
plt.rcParams["figure.autolayout"] = True

change_threshold = 0.0001
crisp_threshold = 0.01
pos = []

def get_colours(G,pos):
    colour_map=[]
    for node in G:
        strength=G.nodes[node]['values'][pos]
        if ((1-strength) < crisp_threshold):
            node_colour = "green"
        else:
            if (strength <= crisp_threshold):
                node_colour = "red"
            else:
                node_colour = "grey"
        colour_map.append(node_colour)
    return colour_map

def read_tgf(filename,G):
    node_list=[]
    att_list=[]
    if os.path.isfile(filename):
        file_handler = open(filename)
        started = False
        readArgs = True
        error = False
        for line in file_handler:
            if ('#' in line):
                readArgs=False
            else:
                if (readArgs):
                    argLine=line.split()
                    argument=argLine[0]
                    # print(argLine,argument)
                    # break
                    node_list.append(argument)
                    if (len(argLine)>1):
                        try:
                            strength=float(argLine[1])
                        except ValueError:
                            print("\nError in strength of argument \"",argument,"\": \"",argLine[1],"\"",sep="")
                            error=True
                            break
                    else:
                        strength=0.5
                    G.add_node(argument,weight=strength,values=[strength])
                else:
                    attack=line.strip().split(" ")
                    # print("#######",attack)
                    # break
                    if (len(attack)==2):
                        source=attack[0]
                        target=attack[1]
                        if (source in G.nodes()) and (target in G.nodes()):
                            att_list.append((source,target))
                        else:
                            print("Invalid source of target argument in attack line:",attack)
                            error = True
                            break
                    else:
                        error = True
                        break
        file_handler.close()
        # print(node_list)
        G.add_edges_from(att_list)
        return True
    else:
        print(filename, 'does not exist')
        return False
def GR(this_val,att_vals):
    if (len(att_vals) > 0):
        max_att = max(att_vals)
    else:
        max_att = 0
    next_val = (1 - this_val) * min([0.5, 1 - max_att]) + this_val * max([0.5, 1 - max_att])
    return next_val

def plot_graph(G):
    figure, (ax1, ax2, ax3) = plt.subplots(1, 3, width_ratios=[1.25, 2.5, 1.25])
    pos = nx.nx_pydot.graphviz_layout(G, prog="dot")
    ival_pos={}
    ivals={}
    my_labels={}
    for node in pos:
        this_label=node+"  "+"{:.2f}".format(G.nodes[node]['values'][0])
        (x1,y1)=pos[node]
        ax1.text(x1,y1,this_label, fontsize=12, 
                 verticalalignment='center_baseline', horizontalalignment='left') 
    colour_map = get_colours(G,0)
    ax1.set_title("Initial Values")
    nx.draw(G, pos, node_color=colour_map, node_size=500, font_size=12, with_labels=False, ax=ax1)

    ax2.set_title("Value Evolution")
    ax2.set_xticks(range(0, min(max_iter, iter) + 1, 5))
    ax2.grid(visible=True, which='major', axis='both')
    ax2.set_ylabel('value', fontsize=12)
    ax2.set_xlabel('iteration', fontsize=12)
    for key in G.nodes():
        values = G.nodes[key]['values']
        ax2.plot(range(min(max_iter, iter) + 1), values, marker='.', label=key)
    ax2.legend(loc="best")
    
    for node in pos:
        this_label=node+"  "+"{:.2f}".format(G.nodes[node]['values'][-1])
        (x1,y1)=pos[node]
        # print("======>",(x1,y1),(x2,y2),"{:.2f}".format(G.nodes[node]['values'][-1]))
        ax3.text(x1,y1,this_label, fontsize=12, 
                 verticalalignment='center_baseline', horizontalalignment='left') 
    colour_map = get_colours(G,-1)
    ax3.set_title("Final Values")
    nx.draw(G, pos, node_color=colour_map, node_size=500, font_size=12, with_labels=False, ax=ax3)
    plt.show(block=False)

if __name__ == '__main__':
    print("Initial argument values are allowed in TGF file.")
    print("Use empty filename to stop.\n")
    if (len(sys.argv)>1):
        filename=sys.argv[1]
    else:
        filename=input("Enter the name of the TGF file: ").strip()
    while (True):
        G = nx.DiGraph()
        if filename and read_tgf(filename,G):
            max_iter=100
            iter=0
            val_change=True
            stable=False
            stabIter=0
            while (iter < max_iter) and val_change:
                val_change = False
                prevStable = stable
                stable = True
                for node in G.nodes():
                    # Retrieve the value of the node in previous iteration
                    node_val=G.nodes[node]['values'][iter]
                    # Retrieve the value of all of its attackers
                    att_vals=[]
                    for att in G.predecessors(node):
                         att_vals.append(G.nodes[att]['values'][iter])
                    # The node's new value is given by a function of its value and that
                    # of its attackers. Below we use the Gabbay-Rodrigues Equation function:
                    next_val=GR(node_val,att_vals)
                    val_change=val_change or ((next_val-node_val)>change_threshold)
                    stable=stable and (not ((node_val == 0) and (next_val > 0)) or ((node_val == 1) and (next_val < 1)))
                    G.nodes[node]['values'].append(next_val)
                iter+=1
                if not(prevStable) and stable:
                    stabIter=iter-1
            print("\nStable at iteration ",stabIter,". Values after ",iter," iteration(s):",sep="")
            for node in G.nodes():
                values = G.nodes[node]['values']
                print(node, ["%.4f" % member for member in values])
            print("")
            print("Nodes accepted at stable iteration (modulo threshold "+"{:.5f}".format(crisp_threshold)+"):")
            print("[ ",sep="", end="") 
            for node in G.nodes():
                if ((1-G.nodes[node]['values'][stabIter]) < crisp_threshold):
                    print(node," ",sep="",end="")
            print("]\n",sep="") 
            print("Nodes accepted at final iteration (modulo threshold "+"{:.5f}".format(crisp_threshold)+"):")
            print("[ ",sep="", end="")
            for node in G.nodes():
                if ((1-G.nodes[node]['values'][-1]) < crisp_threshold):
                    print(node," ",sep="",end="")
            print("]\n",sep="") 
            plot_graph(G)
        else:
            if not filename:
                break
        filename=input("Enter the name of the TGF file: ").strip()

