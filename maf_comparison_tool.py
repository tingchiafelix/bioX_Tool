#/!/usr/local/bin/env python3.7
# -*- coding: utf-8 -*-
"""
MAF files comparison:
    - Based on the chromosome# and position of a variant
    - Generate the common and unique variants from both files
"""
import os
import csv
import argparse
import sys

version = '1.1.02112019'

def get_args():
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument(
        '-f','--frontfile',
        metavar='<frontfile>',
        help='The first MAF file'
    )
    parser.add_argument(
        '-b','--backfile',
        metavar='<backfile>',
        help='The second MAF file'
    )
    parser.add_argument(
        '-l', '--listfile',
        metavar='<list_file>',
        help='A list of MAF files: The MAF files compared to each other\
        should be listed in the same row and separared by tab. e.g.: \
        => same row: \
        022348P.clean.annotated.maf 022348P.exome.clean.annotated.maf'
    )
    ###To Do: add save function in the argument###
    parser.add_argument(
        '-v', '--version',
        action = 'version',
        version = '%(prog)s - v' + version
    )
    args = parser.parse_args()
    if args.listfile:
        return args
    elif (args.frontfile is None) or (args.backfile is None):
        sys.stderr.write('ERROR: You must input two MAF files')
        sys.exit(1)

def dict_fun(ref):
    refile_dict = {}
    with open(ref) as rf:
        for va in rf:
            if not va.startswith('#version'):
                va = va.rstrip().split('\t')
                if not 'Hugo_Symbol' in va[0]:
                    va_key = ','.join(va[4:6])
                    refile_dict[va_key] = 1            
    return refile_dict

def main(frontfile, backfile):

    ### Step1: make front file as ref and compare to back file ###
    comfile_total_n = 0
    common_va_n = 0
    comfile_uni_n = 0
    comfile_total_n_R =0
    common_va_n_R = 0
    comfile_uni_n_R =0
    frontfile_path = os.path.abspath(frontfile)
    frontfile_name = (frontfile_path.split('/')[-1]).split('.')[0]
    backfile_path = os.path.abspath(backfile)
    backfile_name = (backfile_path.split('/')[-1]).split('.')[0]

    reference_file = dict_fun(frontfile)
    with open(backfile) as bf:
        for row in bf:
            if not row.startswith('#version'):
                comfile_total_n +=1
                row = row.rstrip().split('\t')
                new_row = [row[0]] + row[4:7] + [row[8]] + \
                [row[10]] + [row[12]] +[row[15]] + [row[36]] +\
                row[39:42] + row[132:140]
                comfile_key = ','.join(new_row[1:3])
                if 'Hugo_Symbol' in new_row[0]:
                    pass
                    #print ('this is a header:', new_row)
                else:
                    if comfile_key in reference_file:
                        common_va_n += 1
                        #print ('common variants:', new_row)
                    else:
                        comfile_uni_n += 1
                        #print ('Unique variants:', new_row)
    print ('-----------------------------------------------------------------')
    print ('Reference file:', frontfile)
    print ('Compared file:', backfile)
    print ('Total # of variant in compared file:', comfile_total_n-1)
    print ('Common variants in both files:', common_va_n)
    print ('Unique variants in compared file:',comfile_uni_n, '\n')

    ### Step2: make back file as ref and compare to front file ###
    back_ref = dict_fun(backfile)
    with open(frontfile) as ff:
        for line in ff:
            if not line.startswith('#version'):
                comfile_total_n_R += 1
                line = line.rstrip().split('\t')
                new_line = [line[0]] + line[4:7] + [line[8]] + \
                        [line[10]] + [line[12]] +[line[15]] + [line[36]] +\
                        line[39:42] + line[132:140]
                comfile_key_1 = ','.join(new_line[1:3])
                if 'Hugo_Symbol' in new_line[0]:
                    pass
                else:
                    if comfile_key_1 in back_ref:
                        common_va_n_R += 1
                    else:
                        comfile_uni_n_R += 1
    print ('Reference file:', backfile)
    print ('Compared file:', frontfile)
    print ('Total # of variant in compared file:',comfile_total_n_R-1)
    print ('Common variants in both files:', common_va_n_R)
    print ('Unique variants in compared file:',comfile_uni_n_R)
    print ('-----------------------------------------------------------------')

if __name__ == '__main__':
    args = get_args()
    if args.frontfile and args.backfile:
        frontfile = args.frontfile
        backfile = args.backfile
        main(frontfile, backfile)
    else:
        if args.listfile:
            with open(args.listfile) as lf:
                for lf_row in lf:
                    lf_row = lf_row.rstrip().split('\t')
                    frontfile = lf_row[0]
                    backfile = lf_row[1]
                    main(frontfile,backfile)
        else:
            sys.exit(1)
