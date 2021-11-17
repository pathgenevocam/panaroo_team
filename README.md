# panaroo_team

##### Comments

Work in progress--email or slack me if you have any trouble adding to our new repo.
I added separate folders for code and data for basic organisation (feel free to change/expand). 
(?) Aaron: would it be ok for me to host a copy of your data here, or should I leave it in the MRC cloud?

##### Useful commands

Set-up:

    git clone git@github.com:pathgenevocam/panaroo_team.git
    cd panaroo_team
    
    # git checkout caitlin  # replace caitlin w your branch name
    
Push local changes to repo: 

    git add -A
    git commit -m "did something useful"
    git status
    
    # git push origin caitlin  # replace caitlin w your branch name
    git push origin main
##### Parsimony ancestral character state reconstruction and joining up genes co-gained/lost into events
For the reconstruction run

    python parsimony.py <tree> <gene_presence_absence.Rtab> <outdir> --is_transposed --tree_is_named 

which will create a directory for the output if it doesn't already exist and write a table with the events ```events.txt``` and the ancestral states ```ancestral_states.txt```. If the ```--tree_is_named``` option is not set it will also (re)name the internal nodes of the tree and write that tree to disk ```named_tree.nwk```.

For joining up genes co-gained/lost into events run

    reconstructGeneEvents.py <graph> <recon> <tree> <root_name> <outdir/prefix>

You have to create ```outdir``` first. This will create a number of output files. 
