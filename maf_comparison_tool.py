#/!/usr/bin python3
# -*- coding: utf-8 -*-
"""
MAF files comparision:
    - Based on the chromosome# and position of a variant
    - Generate the common and unique variant from both files
"""
import os
import csv
import argparse
import sys

version = '1.0.02062019'

def get_args():
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument(
        'frontfile',
        metavar='<frontfile>',
        help='The first MAF file'
    )
    parser.add_argument(
        'backfile',
        metavar='<backfile>',
        help='The second MAF file'
    )
### To Do: add one arg for a file with large samples###
    parser.add_argument(
        '-v', '--version',
        action = 'version',
        version = '%(prog)s - v' + version
    )
    args = parser.parse_args()

    return args

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
    print ('Reference MAF file:', frontfile, 'compared to MAF file:', backfile)
    print ('Total # of variant in', backfile,':',comfile_total_n-1)
    print ('Common variants in both MAF files:', common_va_n)
    print ('Unique variants in',backfile,':',comfile_uni_n, '\n')
### To do: save common and unique variants in a file ###

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
    print ('Reference MAF file:', backfile, 'compared to MAF file:', frontfile)
    print ('Total # of variant in', frontfile,':',comfile_total_n_R-1)
    print ('Common variants in both MAF files:', common_va_n_R)
    print ('Unique variants in',frontfile,':',comfile_uni_n_R)


if __name__ == '__main__':
    args = get_args()
    frontfile = args.frontfile
    backfile = args.backfile
    main(frontfile, backfile)

