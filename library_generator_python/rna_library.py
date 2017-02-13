import sys, os, getopt,popen2, time,shutil
from glob import glob
import ntpath
import numpy as np
import IPython
import subprocess
from multiprocessing import Pool
import logging
import library_utilities as lutil

NTs = "AGCU"
num_vienna_models=6
autofill=3000

tmp_dir = os.path.dirname(os.path.abspath(__file__)) + "/tmp_rna"
if not os.path.exists(tmp_dir): os.mkdir(tmp_dir)
outs_dir = os.path.dirname(os.path.abspath(__file__)) + "/outs"
if not os.path.exists(outs_dir): os.mkdir(outs_dir)
current_dir = os.path.dirname(os.path.abspath(__file__))

def delete_temporary():
    shutil.rmtree(tmp_dir)
    for file in glob(current_dir+"/*.ss"):
        os.remove(file)

def get_property_table():
    no_properties=4
    no_letters=4
    property_table = np.zeros((no_properties, no_letters))
    property_table[0]=[x/float(600) for x in [599, 403, 336,372]]
    property_table[1]=[x/float(100) for x in [91, 194, 96, 63]]
    property_table[2]=[x/float(2.8) for x in[2.8, 1.73,1.99,2.73]]
    property_table[3]=[x/float(0.82) for x in[0.43,0.82,0.58,0.5]]

    return property_table


def rna_subopt_parallel(num_vienna_models,frag_name,property_table):
    #logging.debug('rna_subopt_parallel')
    fragments_file=os.path.join(tmp_dir,frag_name)
    cmd="./bin/RNAsubopt "
    args="-e 1 -s < "+fragments_file +"  | head -10 | sort -n -k 2 | awk '{print $1}'"
    output = popen2.Popen3(cmd + args)

    while output.poll() < 0:
        try:
            output.wait()
            time.sleep(0.001)
        except:
            break

    seq=output.fromchild.readline().strip()
    structure=output.fromchild.readlines()

    fragment_library=""
    if len(structure) >= num_vienna_models:
        fragment_library=analyse_fragment(num_vienna_models,frag_name,seq,structure,property_table)


    return fragment_library


def rna_plot(frag_name,frag_file):
    #logging.debug('rna_plot.')
    cmd="./bin/RNAplot "
    args=" -o xrna < "+frag_file

    output = popen2.Popen3(cmd + args)

    while output.poll() < 0:
        try:
            output.wait()
            time.sleep(0.001)
        except:
            break


    z = np.array([[complex(float(line.split()[2]), float(line.split()[3])) for line in open(frag_name+'_ss.ss','r').readlines() if line.rstrip() and not line.startswith('#')]])
    distances=abs(z.T-z)
    neighbors=np.array([(sum(i<50 for i in distances[j])-1)/float(50) for j in range(0,len(distances))])


    myarray = np.zeros(autofill)
    myarray[:neighbors.shape[0]] = neighbors

    return neighbors

library=[]
def log_result_rnasubopt(result):
    # This is called whenever foo_pool(i) returns a result.
    # result_list is modified only by the main process, not the pool workers.
    if result:
        library.append(result)
    #logging.debug('Result done.')


def analyse_fragment(num_vienna_models,frag_name,seq,structures,property_table):
    #logging.debug('analyse_fragment.')
    models=[]
    #logging.debug('analyse rnaplot.')
    for j in range(0,num_vienna_models):

        frag_file=os.path.join(tmp_dir,frag_name)
        frag_file_handle=open( frag_file ,"w")
        frag_file_handle.writelines(">"+frag_name+"\n"+str(seq)+"\n")

        frag_file_handle.writelines(structures[j])
        frag_file_handle.close()

        models.append(rna_plot(frag_name,frag_file))

    #logging.debug('Score and smooth.')
    scored=lutil.score_sequence(seq,NTs,property_table)
    smoothed=lutil.smooth(scored, "rna", autofill)

    models=models+[np.array(s) for s in smoothed]

    #logging.debug('extend and ft')
    models_extended=list(map(lambda m: lutil.adaptator(m, autofill), models))
    models_compiled=list(map(lutil.ft_coefficients, models_extended))


    #print "\n".join(frag_name+"\t"+"\t".join(map(str, x)) for x in models_compiled)
    return "\n".join(frag_name+"\t"+"\t".join(map(str, x)) for x in models_compiled)


def run_library(input_file):
    input_file_handle=open(input_file,'r')
    fragments_names=[]
    property_table=get_property_table()
    logging.debug('Reading fragments.')


    lines = input_file_handle.readlines()
        # stuff that is empty after strip is not really important for us,
        # most likely an empty line
    lines = filter(lambda x: x.strip(), lines)

    for line in lines:

        frag_name=line.split()[0]
        sequence=line.split()[1].upper().replace("T","U")

        fragments_file=os.path.join(tmp_dir,frag_name)
        fragments_file_handle=open(fragments_file,"w")

        fragments_names.append(frag_name)
        fragments_file_handle.writelines(sequence+'\n')
        fragments_file_handle.close()

    cores=min(10, len(fragments_names))
    pool = Pool(processes=cores)

    logging.debug('Parrallelizing RNAsubopt.')
    for frag_name in fragments_names:
        #rna_subopt_parallel(num_vienna_models,frag_name,property_table)
        pool.apply_async(rna_subopt_parallel, args=(num_vienna_models,frag_name,property_table),callback = log_result_rnasubopt)
        #pool.apply_async(rna_subopt_parallel, args=([fragments_file]),callback = log_result_rnasubopt)


    pool.close()
    pool.join()

    logging.debug('RNAsubopt parallel done.')
    logging.debug('Writting library')

    new_lib=open(os.path.join(outs_dir,ntpath.basename(input_file)+".rna.lib") ,"w")
    new_lib.write("\n".join(library))
    new_lib.close()
