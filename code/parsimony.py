import ete3 as ete 
import pandas as pd
import argparse
import random
import sys
import os
import numpy as np
random.seed(0)

def parse():
    parser = argparse.ArgumentParser("reconstruct gain/losses using max parsimony")
    parser.add_argument("tree", help="input phylogenetic tree")
    parser.add_argument("table", help="output")
    parser.add_argument("outdir", help="isolate group specific nodes out")
    parser.add_argument("--tree_is_named", action = "store_true", help = "set if tree has internal node names")
    parser.add_argument("--is_transposed", action = "store_true", help = "set to indicate input is genes (rows) X samples (columns) like in Panaroo presence/absence")
    args = parser.parse_args()
    parsimony(**vars(args))

def parsimony(tree, table, outdir, tree_is_named, is_transposed):
    
    t = ete.Tree(tree, format = 1)
    i = 0
    if not tree_is_named:
        for n in t.traverse("preorder"):
            if not n.is_leaf():
                i += 1
                n.name = "N%s" % i
        t.write(outfile = outdir + "/named_tree.nwk", format = 1)

    table = pd.read_csv(table, index_col = 0, sep = "\t")
    if is_transposed:
        table = table.T
    table.index  = [i.replace(".velvet", "") for i in table.index]
    table.index  = [i.replace(".spades", "") for i in table.index]
    #set all indels to WT for PAO1
    #table.loc["PAO1",] = 0
    #or remove PAO1 node for GPA
    #pao1 = t.search_nodes(name="AE004091")[0]
    #pao1.detach()
    table.fillna(float('inf'), inplace = True)
    node_names =  [n.name for n in t.traverse("preorder")]
    out_table = pd.DataFrame(np.zeros((len(node_names), table.shape[1])), index = node_names, columns = table.columns)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    events = open("%s/events.txt" % outdir, 'w')
    events.write("variant_id\tnode\tparent_node\tnode_state\tparent_node_state\n")
    for c in table.columns:
        if all(pd.isnull(table.loc[:, c])) | all(table.loc[:, c] == np.inf):
            out_table.loc[:, c] = np.nan 
        else:
            node2state = down_pass(t, table.loc[:, c])
            reconstruction = up_pass(node2state, t, c, events) 
            out_table.loc[reconstruction.keys(), c] = list(reconstruction.values())
    events.close()
    out_table.to_csv("%s/ancestral_states.txt" % outdir, sep = "\t")

def up_pass(node2state,tree, c, events):
    reconstruction = {}
    for n in tree.traverse("preorder"):
        #todo consider nas
        if n == tree.get_tree_root():
            states = list(node2state[n.name])
            if len(states) > 1:
                #reconstruction[n.name] = states[random.randrange(len(states))] 
                reconstruction[n.name] = states[0] 
            else:
                reconstruction[n.name] = states.pop()
        else:
            up_state = reconstruction[n.up.name]
            n_states = list(node2state[n.name])
            if up_state in node2state[n.name]:
                reconstruction[n.name] = up_state
            else:
                #pick a random state and record event
                if float('inf') in node2state[n.name]:
                    reconstruction[n.name] = float('nan')
                elif len(n_states) == 0:
                    reconstruction[n.name] = float('nan')
                else:
                    #print(n_states, n.name)
                    #reconstruction[n.name] = n_states[random.randrange(len(n_states))] 
                    reconstruction[n.name] = n_states[0] 
                    events.write("%s\t%s\t%s\t%s\t%s\n" % (c, n.name, n.up.name, reconstruction[n.up.name], reconstruction[n.name]))
    return reconstruction


def down_pass(tree, annotation):
    node2state = {}
    for n in tree.traverse("postorder"):
        #todo consider nas
        if n.is_leaf():
            node2state[n.name] = set([annotation.loc[n.name]])
        else:
            node2state[n.name] = set.intersection(*[node2state[i.name] for i in n.children])
            if len(node2state[n.name]) == 0:
                node2state[n.name] = set.union(*[node2state[i.name] for i in n.children])
            try:
                node2state[n.name].remove(float('inf'))
            except:
                pass
    return node2state

        

    

if __name__ == "__main__":
    parse()
