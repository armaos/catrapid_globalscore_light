import sys, os, getopt,popen2, time,shutil
from glob import glob
import ntpath
import numpy as np
import IPython
import subprocess
from multiprocessing import Pool
import logging
import library_utilities as lutil

AAs = "ACDEFGHIKLMNPQRSTVWYX"
autofill=750

tmp_dir = os.path.dirname(os.path.abspath(__file__)) + "/tmp_protein"
if not os.path.exists(tmp_dir): os.mkdir(tmp_dir)
current_dir = os.path.dirname(os.path.abspath(__file__))
outs_dir = os.path.dirname(os.path.abspath(__file__)) + "/outs"
if not os.path.exists(outs_dir): os.mkdir(outs_dir)


def delete_temporary():
    shutil.rmtree(tmp_dir)

def get_property_table():
    return lutil.parse_scales("protein.scales",AAs)[1]

library=[]
def log_result(result):
    # This is called whenever foo_pool(i) returns a result.
    # result_list is modified only by the main process, not the pool workers.
    if result:
        library.append(result)
    #logging.debug('Result done.')

def analyse_fragment(fragment,property_table):

    seq=fragment[1]
    frag_name=fragment[0]
    scored=lutil.score_sequence(seq,AAs,property_table)
    smoothed=lutil.smooth(scored, "protein", autofill)
    models=[np.array(s) for s in smoothed]

    #logging.debug('extend and ft')
    models_extended=list(map(lambda m: lutil.adaptator(m, autofill),models))

    models_compiled=list(map(lutil.ft_coefficients, models_extended))

    #IPython.embed()
    #print "\n".join(frag_name+"\t"+"\t".join(map(str, x)) for x in models_compiled)
    return "\n".join(frag_name+"\t"+"\t".join(map(str, x)) for x in models_compiled)

def run_library(input_file):
    input_file_handle=open(input_file,'r')
    fragments=[]
    property_table=get_property_table()
    logging.debug('Reading fragments.')


    lines = input_file_handle.readlines()
        # stuff that is empty after strip is not really important for us,
        # most likely an empty line
    lines = filter(lambda x: x.strip(), lines)

    for line in lines:
        fragments.append((line.split()[0],line.split()[1]))


    cores=min(10, len(fragments))
    pool = Pool(processes=cores)

    logging.debug('Parrallelizing Fragment Analysis.')
    for fragment in fragments:
        #analyse_fragment(fragment,property_table)
        pool.apply_async(analyse_fragment, args=(fragment,property_table),callback = log_result)
        #pool.apply_async(rna_subopt_parallel, args=([fragments_file]),callback = log_result_rnasubopt)


    pool.close()
    pool.join()

    logging.debug('Writting library')

    new_lib=open(os.path.join(outs_dir,ntpath.basename(input_file)+".protein.lib") ,"w")
    new_lib.write("\n".join(library))
    new_lib.close()
