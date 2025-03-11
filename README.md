# panaroo_team


##### Parsimony ancestral character state reconstruction and joining up genes co-gained/lost into events
You need to install ete3, biopython, pandas and networkx e.g. using conda.
For the reconstruction  run

    python parsimony.py <tree> <gene_presence_absence.Rtab> <outdir> --is_transposed --tree_is_named 

which will create a directory for the output if it doesn't already exist and write a table with the events ```events.txt``` and the ancestral states ```ancestral_states.txt```. If the ```--tree_is_named``` option is not set it will also (re)name the internal nodes of the tree and write that tree to disk ```named_tree.nwk```.

For joining up genes co-gained/lost into events run

    reconstructGeneEvents.py <graph> <recon> <tree> <root_name> <outdir/prefix>

You have to create ```outdir``` first. This will create a number of output files. 
