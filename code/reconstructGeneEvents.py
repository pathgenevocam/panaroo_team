#Reconstructs gene changes identified by PastML onto a tree and joins gene gains and losses into an event if they occur along the same branch and are
#adjacent in the pangenome graph
#To run: python3 reconstructGeneEvents.py -g final_graph.gml -r combined_ancestral_states.tab -t tree -o outfilePrefix

import argparse
from Bio import Phylo
import networkx as nx

def labelTreeClades(tree): #Takes a tree and returns the tree with node names as in the gene dictionary
    cladeIterator = 0 #Will be appended with each internal non-root node
    for clade in tree.find_clades(): #Iterate through the clade
        if len(tree.get_path(clade)) == 0: #Check if the clade is the root
            clade.nodeName = "root"
        else:
            if not clade.name: #Check if the clade is an internal node
                clade.nodeName = "n" + str(cladeIterator)
                cladeIterator += 1
            else:
                clade.nodeName = clade.name 
    return(tree)

def getParentNode(tree,clade): #Takes a tree and a clade and returns the identifier of the upstream clade
    if len(tree.get_path(clade)) >= 2: #Check if the clade is one downstream of the root
        upstreamNode = tree.get_path(clade)[-2]
    else:
        upstreamNode = "root"
    return(upstreamNode)

def joinGenes(genes): #Takes a set of genes and combines the adjacent genes
    canJoin = True

    while canJoin == True:
        joined = 0
        for geneA in genes: #Iterate through the genes
            for geneB in genes: #Iterate through the genes
                if geneA != geneB: #Do not compare the same gene
                    for adjacentGene in geneB[1]: #Iterate through the neighbouring genes
                        if adjacentGene in geneA[0]: #Check if the neighbouring gene is in another event
                            geneA[0] = geneA[0] + geneB[0]
                            geneA[1] = geneA[1] + geneB[1]
                            genes.remove(geneB)
                            for geneName in geneA[1]:
                                if geneName in geneA[0]:
                                    geneA[1].remove(geneName)
                            joined += 1
                            break
        if joined == 0:
            canJoin = False
    return(genes)

def convertNodeNames(geneGraph,nodeName): #Converts a given node name to the node number in the graph
    print(geneGraph)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("graph", help = "final_graph.gml from Panaroo")
    parser.add_argument("recon", help = "tabular ancestral character state reconstruction e.g. from parsimony.py")
    parser.add_argument("tree", help = "Rooted tree containing all and only samples used in the pangenome reconstruction, with exactly the same sample names")
    parser.add_argument("root_name", help = "name of the root in the ancestral character state reconstruction table")
    parser.add_argument("out", help = "Output file prefix")
    args = parser.parse_args()

    geneNodeOutFile = open(args.out + ".branch.gene.gains.losses.txt","w")
    geneNodeOutFile.write("Branch\tGenes.gained\tGenes.lost\n")
    labelledTreeOutFile = open(args.out + ".tree.with.node.labels.nex","w")
    eventNodeOutFile = open(args.out + ".branch.event.gains.losses.txt","w")
    eventNodeOutFile.write("Branch\tEvents.gained\tEvents.lost\n")
    eventNodeSummaryOutFile = open(args.out + ".branch.event.number.gains.losses.txt","w")
    eventNodeSummaryOutFile.write("Branch\tEvents.gained\tEvents.lost\n")
    eventSummaryOutFile = open(args.out + ".event.summary.txt","w")
    eventSummaryOutFile.write("Event\tNumber.of.genes\tGenes\tBranch\tGain.or.loss\n")

    geneGraph = nx.read_gml(args.graph) #Import the gene graph

    reconstruction = open(args.recon).readlines() #Import the gene reconstructions
    geneReconstruction = [] #Will be filled with the nodes and tips and their gene presence absence
    processedNodes = []
    geneDictionary = {} #Will be filled with each node and the gene presence absence as a string
    
    geneNames = reconstruction[0].strip().split("\t") #The name of each gene

    for node in reconstruction[1:]: #Iterate through the nodes and remove the second line for the internal nodes
        if node.rstrip().split("\t")[0] not in processedNodes:
            geneReconstruction.append(node)
            processedNodes.append(node.rstrip().split("\t")[0])
            if node.strip().split("\t")[0] != "node": #Append the gene presence absence information for the node
                geneDictionary[node.rstrip().split("\t")[0]] = "".join([str(int(float((i)))) if i != "" else " " for i in node.rstrip().split("\t")[1:]])
    

    nodeGenes = [] #Will be filled with each node and the genes that are gained and lost along the corresponding branch
    
    tree = Phylo.read(args.tree,"newick") #Import the tree
    labelledTree = labelTreeClades(tree)
    Phylo.write(labelledTree,labelledTreeOutFile,"nexus")
    for clade in labelledTree.find_clades(): #Iterate through the clades
        if len(tree.get_path(clade)) != 0: #Do not analyse the root
            upstreamNode = getParentNode(tree,clade)
            if upstreamNode == "root": #Check if the upstream node is the root
                differentGenes = [i + j for i,j in zip(geneDictionary[args.root_name],geneDictionary[clade.nodeName])] #1 or 0 for each gene if the gene is present or absent
                geneGains = [] #Will be filled with the gene names that are gained along the branch
                geneLosses = [] #Will be filled with the gene names that are lost along the branch
                for i,gene in enumerate(differentGenes):
                    if gene == "01": #Check if the gene was gained along the branch
                        geneGains.append(geneNames[i])
                    elif gene == "10": #Check if the gene was lost along the branch
                        geneLosses.append(geneNames[i])
                nodeGenes.append(["root:" + clade.nodeName,geneGains,geneLosses])
                geneNodeOutFile.write("root:" + clade.nodeName + "\t" + ",".join(geneGains) + "\t" + ",".join(geneLosses) + "\n")
            else:
                differentGenes = [i + j for i,j in zip(geneDictionary[upstreamNode.nodeName],geneDictionary[clade.nodeName])] #1 or 0 for each gene if the gene is present or absent
                geneGains = [] #Will be filled with the gene names that are gained along the branch
                geneLosses = [] #Will be filled with the gene names that are lost along the branch
                for i, gene in enumerate(differentGenes):
                    if gene == "01": #Check if the gene was gained along the branch
                        geneGains.append(geneNames[i])
                    elif gene == "10": #Check if the gene was lost along the branch
                        geneLosses.append(geneNames[i])
                nodeGenes.append([upstreamNode.nodeName + ":" + clade.nodeName,geneGains,geneLosses])
                geneNodeOutFile.write(upstreamNode.nodeName + ":" + clade.nodeName + "\t" + ",".join(geneGains) + "\t" + ",".join(geneLosses) + "\n")
    
    eventNumber = 1 #Will be incremented with each gain or loss event

    for branch in nodeGenes: #Iterate through the branches
        branchGainEvents = [] #Will be filled with the gain events along the branch
        branchLossEvents = [] #Will be filled with the loss events along the branch
        gainNeighbours = [] #Will be filled with each gained gene and its neighbours
        lossNeighbours = [] #Will be filled with each lost gene and its neighbours
        gainEvents = [] #Will be filled with the gene gain events
        lossEvents = [] #Will be filled with the gene loss events
        for node in geneGraph: #Iterate through the nodes in the gene graph
            if geneGraph.nodes[node]["name"] in branch[1]: #Check if the node is gained along the branch
                genes = [[geneGraph.nodes[node]["name"]],[]] #Append the gene name and an empty list, the empty list will be filled with the neighbours
                for geneNeighbour in geneGraph.neighbors(node):
                    genes[1].append(geneGraph.nodes[geneNeighbour]["name"])
                gainNeighbours.append(genes)
                gainEvents = joinGenes(gainNeighbours) #Assign to the gene gain events
            elif geneGraph.nodes[node]["name"] in branch[2]: #Check if the node is lost along the branch
                genes = [[geneGraph.nodes[node]["name"]],[]]
                for geneNeighbour in geneGraph.neighbors(node):
                    genes[1].append(geneGraph.nodes[geneNeighbour]["name"])
                lossNeighbours.append(genes)
                lossEvents = joinGenes(lossNeighbours) #Assign to the gene loss events
        totalGainEvents = [] #Will be filled with all of the gain events along the branch
        totalLossEvents = [] #Will be filled with all of the loss events along the branch
        totalGainEventNumbers = [] #Will be filled with each gain event number
        totalLossEventNumbers = [] #Will be filled with each loss event number
        for geneGain in gainEvents: #Iterate through the gene gain events
            totalGainEvents.append(",".join(geneGain[0]))
            totalGainEventNumbers.append("event" + str(eventNumber))
            eventSummaryOutFile.write("event" + str(eventNumber) + "\t" + str(len(geneGain[0])) + "\t" + ",".join(geneGain[0]) + "\t" + branch[0] + "\t" + "gain" + "\n")
            eventNumber += 1
        for geneLoss in lossEvents: #Iterate through the gene loss events
            totalLossEvents.append(",".join(geneLoss[0]))
            totalLossEventNumbers.append("event" + str(eventNumber))
            eventSummaryOutFile.write("event" + str(eventNumber) + "\t" + str(len(geneLoss[0])) + "\t" + ",".join(geneLoss[0]) + "\t" + branch[0] + "\t" + "loss" + "\n")
            eventNumber += 1
        eventNodeOutFile.write(branch[0] + "\t" + "|".join(totalGainEvents) + "\t" + "|".join(totalLossEvents) + "\n")
        eventNodeSummaryOutFile.write(branch[0] + "\t" + ",".join(totalGainEventNumbers) + "\t" + ",".join(totalLossEventNumbers) + "\n")
    
    geneDF = [[["A"],["B","D"]],[["D"],["A"]],[["B"],["A","C","E"]],[["E"],["B","F"]],[["F"],["E","H","I","J"]],[["U"],["T","V"]]]
    
    geneNodeOutFile.close()
    labelledTreeOutFile.close()
