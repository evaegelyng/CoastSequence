# Remember to activate gwf_new!

from gwf import Workflow
import os, sys
import math
from glob import glob

project_name = "CoastSeq"

gwf = Workflow(defaults={"account": "edna"}) 

#Taxonomic assignment with DADA2

###Split fasta file (the nochim one with chimeras removed) into K parts
def splitter(inputFile, K=99):
    inputs = [inputFile]
    outputs = ["tmp/split/split.log.txt"]
    options = {
        'cores': 1,
        'memory': '2g',
        'walltime': '1:00:00'
    }
    spec = '''
    eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
    conda activate metabar_2021
    seqkit split -O tmp/split/ {inputFile} -p {K} -2
    echo "hello" > tmp/split/split.log.txt
    '''.format(inputFile=inputFile, K=K)
    return inputs, outputs, options, spec

#####classify a single k-th file
def assignment(k, outFolder):
    inputFasta = 'tmp/split/DADA2_nochim.part_'+'{:0>3d}'.format(k)+'.fasta'
    inputs = [inputFasta]
    outTax = outFolder + '/assign.' + str(k) + '.silva.named'
    outputs = [
      outTax
    ]
    options = {
        'cores': 2,
        'memory': '32g',
        'walltime': '6:00:00'
    }
    spec = '''
    mkdir -p {out}
    echo "RUNNING THREAD {k} TAX"
    eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
    conda activate metabar_2021
    Rscript scripts/dada2_tax.r {inputFasta} {outTax} 
    echo "DONE THREAD {k}"
    '''.format(out=outFolder, k=k, inputFasta=inputFasta, outTax=outTax)
    return inputs, outputs, options, spec

inputName = "../base/results/DADA2_nochim.otus"

gwf.target_from_template( 'split', splitter(inputFile=inputName) )

parts=glob('tmp/split/DADA2_nochim.part_*.fasta')
K=len(parts)
                                                                
for k in range(1,K+1):
    gwf.target_from_template( 'assign_{}'.format(k), assignment(k=k, outFolder='tmp/assignment') )

### Combine all the small taxonomical classfication files into one large file

input_files = ['tmp/assignment/assign.' + str(k) + '.silva.named' for k in range(1,K+1)]
 
output_files = ['results/assignment.txt']
    
gwf.target(
   name="combine_taxonomy_{}".format(project_name),
   inputs=input_files,
   outputs=output_files,
   cores=1,
   memory="1g",
   walltime="00:10:00",
 ) << """
    head -n1 tmp/assignment/assign.1.silva.named > results/assignment.txt
    for fname in tmp/assignment/assign*.silva.named
    do
        tail -n +2 $fname >> results/assignment.txt
    done
   """                
