__author__ = 'Quinn'

import msprime, pyslim
import numpy as np
import argparse
from time import perf_counter

parser = argparse.ArgumentParser(description="For F4 simulations convert SLiM tree output into a VCF")
parser.add_argument('--out', help="Output prefix, required", required=True)
parser.add_argument('--input', help="Input tree, required", required=True)
parser.add_argument('--numInds', help="Number of individuals to subset from each pop")
parser.set_defaults(numInds=10)
args = parser.parse_args()

outputPrefix = args.out
inputFile = args.input
subNum = int(args.numInds)

totalStart = perf_counter()
start = perf_counter()
ts = pyslim.load(inputFile).simplify()
end = perf_counter()
#print("Load = " + str(end - start) + " sec")
start = perf_counter()
mut_ts = msprime.mutate(ts, rate=1.2e-8, random_seed=1, keep=True, model=msprime.InfiniteSites(alphabet=1))
#mut_ts = pyslim.SlimTreeSequence(msprime.mutate(ts, rate=1e-5, keep=True))
end = perf_counter()
#print("Mutate = " + str(end - start) + " sec")

numTrees = ts.num_trees
numInds = mut_ts.num_individuals
numPops = mut_ts.num_populations
popList = ["p1","p2","p3","p4","p5"]
popIndDict = {}
for pop in popList:
    popIndDict[pop] = []

for i in range(0,numInds):
    indID = ts.individual(i).id
    indPop = ts.individual(i).population
    popName = popList[indPop]
    if ts.individual(i).time == 0.0:
        popIndDict[popName].append(indID)

subPopDict = {}
for j in popIndDict:
   #print(f"We have {len(popIndDict[j])} individuals in the {j} population.")
   subPopDict[j] = np.random.choice(popIndDict[j], size=subNum, replace=False)

indivlist = []
indivnames = []
with open(outputPrefix +"_sim_individuals.txt", "w") as indfile:
  #indfile.writelines("\t".join(["vcf_label", "tskit_id", "slim_id", "popNum", "popName"]) + "\n")
  for pop in popList:
     for i in subPopDict[pop]:
        indivlist.append(i)
        ind = mut_ts.individual(i)
        vcf_label = pop + "_" + str(ind.id)
        indivnames.append(vcf_label)
        data = [vcf_label, pop, str(ts.individual(i).metadata.pedigree_id)]
        indfile.writelines("\t".join(data) + "\n")

with open(outputPrefix + "_sim_genotypes.vcf", "w") as vcffile:
  mut_ts.write_vcf(vcffile, individuals=indivlist, individual_names=indivnames)

#start = perf_counter()
#p5time = 1
#p3time = 1
#for x in ts.nodes():
#    if x.population == 4:
#        if int(x.time) > p5time:
#            p5time = x.time + 1
#    elif x.population == 3:
#        if int(x.time) > p3time:
#            p3time = x.time + 1

#was_p3 = [x.id for x in ts.nodes() if (x.population == 2 and x.time == p5time)]
#was_p3 = [x.id for x in ts.nodes() if (x.population == 2)]
#end = perf_counter()
#print("Was p3 = " + str(end - start) + " sec")

samp_inds = subPopDict['p5']
start = perf_counter()
indNodes = {}
wNodes = list()
for ind in samp_inds:
    nodes = list(ts.individual(ind).nodes)
    indNodes[ind] = list(ts.individual(ind).nodes)
    for node in nodes:
        #print(ts.nodes(node).population)
        wNodes.append(node)

#samp_nodes = np.concatenate([ind.nodes for ind in ts.individuals() if ind.id in samp_inds]).tolist()
end = perf_counter()
#print("Sample nodes = " + str(end - start) + " sec")

outputName = outputPrefix + "_p5_indAnc.txt"
outFile = open(outputName, 'w')
header = "chr\tstart\tend"
for indName  in subPopDict['p5']:
    header = header + "\tp5_" + str(indName)

outFile.write(header + "\n")
#outFile.close()

print("Number of Trees: " + str(numTrees))
count = 0
outStr = ""
for tree in ts.trees():
    loopStart = perf_counter()
    #start = perf_counter()
    #outFile = open(outputName, 'a+')
    count += 1
    interval = tree.interval
    #print(interval)
    #outFile.write("1")
    outStr += "1"
    for i in interval:
        #outFile.write("\t" + str(i))
        outStr += "\t" + str(i)
    for ind in samp_inds:
        sampGeno = 0
        for indNode in indNodes[ind]:
            pop = tree.population(indNode)
            node = indNode
            countBack = 0
            while pop == 4:
                parentNode = tree.parent(node)
                parentPop = tree.population(parentNode)
                if parentPop == pop:
                    pop = parentPop
                    node = parentNode
                    countBack += 1
                else:
                    pop = parentPop
                    #print(ts.node(parentNode))
            if pop == 2:
                sampGeno += 1
        #outFile.write("\t" + str(sampGeno))
        outStr += "\t" + str(sampGeno)
    #outFile.write("\n")
    outStr += "\n"
    #outFile.close()
    loopEnd = perf_counter()
    #print("Tree parse = " + str(loopEnd - loopStart) + " sec")
    #print("At tree: " + str(count) + " of " + str(numTrees))
    #print("Total run time: " + str(loopEnd - totalStart) + " sec")

#outFile = open(outputName, 'a+')
outFile.write(outStr)
outFile.close()
totalEnd = perf_counter()
print(outputPrefix + " Total run time: " + str(totalEnd - totalStart) + " sec")
