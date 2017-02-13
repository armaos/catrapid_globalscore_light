import sys, os, getopt,popen2
import numpy as np
import IPython
import logging


NTs = "AGCU"
AAs = "ACDEFGHIKLMNPQRSTVWYX"


def parse_scales(filename, AAs):
    """

    """
    with open(filename, "r") as f:
        lines = f.readlines()
        # stuff that is empty after strip is not really important for us,
        # most likely an empty line
        lines = filter(lambda x: x.strip(), lines)

    # variable definition
    titles = []
    no_properties = len(lines)/2
    # TODO - why are we using the extra line? ---removed now, let's see
    # property_table = np.zeros((no_properties+1, len(AAs)))
    property_table = np.zeros((no_properties, len(AAs)))

    # now go through the file and assign values
    for i, line in enumerate(lines):
        if i % 2 == 0:
            # scale title line
            titles.append(line.strip())
        else:

            # line of numbers
            numbers = map(float, line.split())
            property_table[i/2] = numbers
            # print numbers
            # print len(numbers)

    return (titles, property_table)

def ft_coefficients(vector):
    shape=52
    length=len(vector)
    coefficients=np.zeros(shape)

    vector[0]=0

    for k in range(0,shape):
        c=0
        for i in range(0,length):

            c+=np.sqrt(2/(float(length)))*vector[i]* np.cos((3.14/float(length))*(k+0.5)*(i+0.5))
        coefficients[k]=np.around(c,2)

    return coefficients

library=[]



def adaptator(vector,autofill):
    myarray = np.zeros(autofill)
    myarray[:vector.shape[0]] = vector


    return myarray

def smooth(scale_values, smooth_type, autofill, window_size=7):
    # TODO - make the multipliers change according to the window size
    # multipliers = [0.25,0.5, 0.75, 1, 0.75,0.5, 0.25]

    one_half = window_size/2 + 1
    step = 1.0/one_half
    one_side = np.arange(0, 1, step).tolist()

    multipliers = one_side[1:] + [1] + one_side[::-1][:-1]
    multipliers = [1,1,1,1,1,1,1]
    # print multipliers

    # TODO - define the window sizes better
    scales, length = scale_values.shape
    length = length -window_size + 1

    try:
        # res = np.empty((scales, length))
        def _proc_scale(scale, smooth_type):
            """
                Inner function smoothing individual scales
            """
            l = len(scale)
            arr = np.zeros((window_size, l))
            if smooth_type=="rna":
                first=np.array([scale[0]/float(1.5),(scale[0]+scale[1]+scale[2])/float(3.5),(scale[0]+scale[1]+scale[2]+scale[3]+scale[4])/float(5.5)])
                last=np.array([(scale[-1]+scale[-2]+scale[-3]+scale[-4]+scale[-5])/float(5.5),(scale[-1]+scale[-2]+scale[-3])/float(3.5),scale[-1]/float(1.5)])
            elif smooth_type=="protein":
                first=np.array([scale[0]/float(1),(scale[0]+scale[1]+scale[2])/float(3.0),(scale[0]+scale[1]+scale[2]+scale[3]+scale[4])/float(5.0)])
                last=np.array([(scale[-1]+scale[-2]+scale[-3]+scale[-4]+scale[-5])/float(5.0),(scale[-1]+scale[-2]+scale[-3])/float(3.0),scale[-1]/float(1.0)])
            for i in range(0, window_size):
                # arr[i] = np.roll(scale, i) #old and slower way
                arr[i, i:] = scale[0:l-i]

            # print scale_values.shape, length, window_size
            my_array=np.append(first,np.average(arr, 0, multipliers)[window_size-1:])
            my_array=np.append(my_array,last)

            myarray = np.zeros(autofill)
            myarray[:my_array.shape[0]] = my_array

            return np.around(my_array,2)
            #return np.average(arr, 0, multipliers)[window_size-1:]
        try:
            return np.apply_along_axis(_proc_scale, 1, scale_values, smooth_type)
        except FloatingPointError, fpe:
            print scale_values
            # import code; code.interact(local=locals())
            raise

    except ValueError:
        # TODO - make explicit length check
        # probably the protein is too small
        print "ValueError"
        return scale_values
        # import code; code.interact(local=locals())

def score_sequence(seq, letters, property_table):
    """
        Applies scales to proteins
        Also converts the profile to uppercase for matching
    """
    indices =  map(lambda s: letters.find(s), seq.upper())
    return property_table[:, indices]
