#!/usr/bin/env python

# Simon H. Martin 2020
# simon.martin@ed.ac.uk

# This script accompanies the paper:
# "Signatures of introgression across the allelel frequency spectrum"
# by Simon H. Martin and William Amos

# It is a wrapper for the moments package (Jouganous et al. 2017, https://doi.org/10.1534/genetics.117.200493).

# Simulation models are defined in periods.
# Each period can have a defined length and size for each population

# "--onePopPeriod" refers to periods before any population split.
# If none of these are specified, it is assumed to have constant population size

# "--twoPopPeriod" refers to periods after the first population split, and so on
# for "--threePopPeriod" and "--fourPopPeriod".

#Any number of periods can be specified.

# Output is a frequency spectrum in tabular format.
# The first N columns corrspond to base counts in the N populations
# Final column represents the proportion of sites with that combination of base counts

import numpy as np

import moments

import itertools, argparse, sys

from collections import defaultdict

def moments_model(params, ns):
    
    onePopPeriods,twoPopPeriods,threePopPeriods,fourPopPeriods = params
    
    # steady state for the equilibrium ancestral population
    sts = moments.LinearSystem_1D.steady_state_1D(sum(ns))
    
    fs = moments.Spectrum(sts)
    
    # any changes in the ancestral population size are defined in the onePopPeriods
    for period in onePopPeriods:
        fs.integrate([period["N"]], period["T"])
    
    if len(twoPopPeriods) >= 1:
        # Split the fs
        fs = moments.Manips.split_1D_to_2D(fs, ns[0], sum(ns[1:]))
        
        # any changes in the two resulting populations are defined in the twoPopPeriods
        # Second pop is ancestor of future pops 2, 3 and 4
        for period in twoPopPeriods:
            fs.integrate([period["N1"],period["N2"]], period["T"],
                        m=np.array([[0, period["m12"]],
                                    [period["m21"],0]]))
    
    if len(threePopPeriods) >= 1:
        # Second Split (split ancestor of 2 and 3)
        fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], sum(ns[2:]))
        
        # periods after second split have three pops and many migration rates
        for period in threePopPeriods:
            fs.integrate([period["N1"], period["N2"], period["N3"]], period["T"],
                        m = np.array([[0,            period["m12"], period["m13"]],
                                        [period["m21"],0,             period["m23"]],
                                        [period["m31"],period["m32"],0             ]]))
    
    if len(fourPopPeriods) >= 1:
        # Second Split (split ancestor of 2 and 3)
        fs = moments.Manips.split_3D_to_4D_3(fs, ns[2], ns[3])
        
        # periods after second split have three pops and many migration rates
        for period in fourPopPeriods:
            fs.integrate([period["N1"], period["N2"], period["N3"], period["N4"]], period["T"],
                        m = np.array([[      0,      period["m12"],period["m13"],period["m14"]],
                                        [period["m21"],      0,      period["m23"],period["m24"]],
                                        [period["m31"],period["m32"],      0      ,period["m34"]],
                                        [period["m41"],period["m42"],period["m43"],      0      ]]))
    
    return fs



def sfsArrayToTable(sfs_array, excludeZero=True):
    indices = itertools.product(*[range(i) for i in sfs_array.shape])
    return [list(idx) + [sfs_array[idx]] for idx in indices if not (sfs_array.mask[idx] or (excludeZero and sfs_array[idx]==0))]



#a function to add the necessary parser arguments. This is so that you can import this function in other scripts and it'll automatically add the required arguments
def addSimArgsToParser(parser):
    parser.add_argument("--Nsam", help="Number of samples in each population", type=int, action = "store", nargs="+")
    parser.add_argument("--onePopPeriod", help="Params for one pop period in format 'T=1 N=0.5'", nargs="+", action = "append")
    parser.add_argument("--twoPopPeriod", help="Params for two pop period in format 'T=1 N2=0.5 m12=0.01 etc'",
                        nargs="+", action = "append")
    parser.add_argument("--threePopPeriod", help="Params for three pop period in format 'T=1 N2=0.5 m12=0.01 etc'",
                        nargs="+", action = "append")
    parser.add_argument("--fourPopPeriod", help="Params for four pop period in format 'T=1 N2=0.5 m12=0.01 etc'",
                        nargs="+", action = "append")


def getParamDict(args):
    paramDict={'Nsam': args.Nsam}
    
    paramDict["onePopPeriods"] = []
    paramDict["twoPopPeriods"] = []
    paramDict["threePopPeriods"] = []
    paramDict["fourPopPeriods"] = []
    
    if args.onePopPeriod is not None:
        for period in args.onePopPeriod:
            paramDict["onePopPeriods"].append(defaultdict(int))
            paramDict["onePopPeriods"][-1].update(dict([(x,float(y)) for x,y in [item.split("=") for item in period]]))
            #check N values
            if paramDict["onePopPeriods"][-1]["N"] == 0: paramDict["onePopPeriods"][-1]["N"] = 1
    
    if args.twoPopPeriod is not None:
        for period in args.twoPopPeriod:
            paramDict["twoPopPeriods"].append(defaultdict(int))
            paramDict["twoPopPeriods"][-1].update(dict([(x,float(y)) for x,y in [item.split("=") for item in period]]))
            #check N values
            for N in "N1","N2":
                if paramDict["twoPopPeriods"][-1][N] == 0: paramDict["twoPopPeriods"][-1][N] = 1
    
    if args.threePopPeriod is not None:    
        for period in args.threePopPeriod:
            paramDict["threePopPeriods"].append(defaultdict(int))
            paramDict["threePopPeriods"][-1].update(dict([(x,float(y)) for x,y in [item.split("=") for item in period]]))
            #check N values
            for N in "N1","N2","N3":
                if paramDict["threePopPeriods"][-1][N] == 0: paramDict["threePopPeriods"][-1][N] = 1
    
    if args.fourPopPeriod is not None:
        for period in args.fourPopPeriod:
            paramDict["fourPopPeriods"].append(defaultdict(int))
            paramDict["fourPopPeriods"][-1].update(dict([(x,float(y)) for x,y in [item.split("=") for item in period]]))
            #check N values
            for N in "N1","N2","N3","N4":
                if paramDict["fourPopPeriods"][-1][N] == 0: paramDict["fourPopPeriods"][-1][N] = 1
    
    return paramDict



if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--commandLinesFile", help="File of multiple command lines", action = "store")
    
    addSimArgsToParser(parser)
    
    #args = parser.parse_args("--Nsam 4 4 20 20 --twoPopPeriod T=1 --threePopPeriod T=1 --fourPopPeriod T=1".split())
    args = parser.parse_args()
    
    if args.commandLinesFile:
        with open(args.commandLinesFile, "rt") as commandLinesFile:
            argsList = [parser.parse_args(l.split()) for l in commandLinesFile]
    
    else: argsList = [args]
    
    for _args_ in argsList:
        paramDict = getParamDict(_args_)
        
        #print(paramDict, file=sys.stderr)
        
        sys.stderr.write("\n{} periods in first phase\n".format(len(paramDict["onePopPeriods"])))
        for period in paramDict["onePopPeriods"]:
            sys.stderr.write("    Length: {} x2N generations\n".format(period["T"]))
            #print(period,file=sys.stderr)

        sys.stderr.write("\n{} periods in second phase\n".format(len(paramDict["twoPopPeriods"])))
        for period in paramDict["twoPopPeriods"]:
            sys.stderr.write("    Length: {} x2N generations\n".format(period["T"]))
            #print(period,file=sys.stderr)

        sys.stderr.write("\n{} periods in third phase\n".format(len(paramDict["threePopPeriods"])))
        for period in paramDict["threePopPeriods"]:
            sys.stderr.write("    Length: {} x2N generations\n".format(period["T"]))
            #print(period,file=sys.stderr)
        
        sys.stderr.write("\n{} periods in fourth phase\n".format(len(paramDict["fourPopPeriods"])))
        for period in paramDict["fourPopPeriods"]:
            sys.stderr.write("    Length: {} x2N generations\n".format(period["T"]))
            #print(period,file=sys.stderr)
        
        ###
        sfs = moments_model((paramDict["onePopPeriods"],
                            paramDict["twoPopPeriods"],
                            paramDict["threePopPeriods"],
                            paramDict["fourPopPeriods"]), paramDict["Nsam"])
        
        sfsTable = sfsArrayToTable(sfs.round(6))

        sys.stdout.write("\n".join(["\t".join([str(val) for val in lst]) for lst in sfsTable]) + "\n")
