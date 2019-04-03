#/!/usr/local/bin/env python3.7
# -*- coding: utf-8 -*-
"""
MAF files comparison: Based on the chromosome# and position of a variant and 
generate the statistics table from both files
"""
import os
import csv
import argparse
import sys
from pprint import pprint as pp
import json

version = '1.1.04022019'

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
    parser.add_argument(
        '-o','--savefile',
        metavar='<savefile>',
        help='Save the output result as csv file'
    )
    parser.add_argument(
        '-v', '--version',
        action = 'version',
        version = '%(prog)s - v' + version
    )
    args = parser.parse_args()
    if args.listfile:
        return args
    elif args.frontfile and args.backfile:
        return args
    elif (args.frontfile is None) or (args.backfile is None):
        sys.stderr.write('ERROR: You must input two MAF files or a list of MAF file in tsv'
        )
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

def main(input_file):
    ###Separate the files ###
    input_file = input_file.rstrip().split(',')
    frontfile = input_file[0]
    backfile = input_file[1]

    comfile_total_n = 0
    common_va_n = 0
    comfile_uni_n = 0
    comfile_total_n_R =0
    common_va_n_R = 0
    comfile_uni_n_R =0
    frontfile_path = os.path.abspath(frontfile)
    frontfile_name = str((frontfile_path.split('/')[-1]).split('.')[0])
    backfile_path = os.path.abspath(backfile)
    backfile_name = str((backfile_path.split('/')[-1]).split('.')[0])

    reference_file = dict_fun(frontfile)
    ### Step1-1: make front file as ref and compare to back file ###
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

    ### Step1-2: make back file as ref and compare to front file ###
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

    ### Step 2: putting the data together###
    #header = ['Samples', '#_Common_Var','#_Uni_file1','#_file1','#_Uni_file2','#_file2']
    #print (header)
    final_dict = {}
    if frontfile_name not in final_dict:
        final_dict[frontfile_name] = {}
    final_dict[frontfile_name]['#_Common_Var'] = str(common_va_n)
    final_dict[frontfile_name]['#_Uni_frontfile'] = str(comfile_uni_n)
    final_dict[frontfile_name]['#_frontfile'] = str(comfile_total_n -1)

    if backfile_name in final_dict:
        if common_va_n == common_va_n_R:
            final_dict[frontfile_name]['#_Uni_backfile'] = str(comfile_uni_n_R)
            final_dict[frontfile_name]['#_backfile'] = str(comfile_total_n_R -1) 
        else:
            print ('Something wrong in the script regarding the matched method!')
    else:
        print (backfile_name, 'is not matched up with ', frontfile_name)
   
    ### Step 3: save a csv file ###

    header = ['Samples', '#_Common_Var','#_Uni_frontfile', \
            '#_frontfile','#_Uni_backfile','#_backfile']

    try:
        savefile = input_file[2]
        file_exists = os.path.isfile(savefile)
        with open(savefile,'a') as csvfile:
            if not file_exists:
                csvfile.write(','.join(header) + '\n')
                print (','.join(header))
            headers = {}
            for k, v in final_dict.items():
                tem = []
                for kk, vv in v.items():
                    tem.append(vv)
                    headers[kk] = 1
                new_line = [k] + tem
                csvfile.write(','.join(new_line) + '\n')
                print (','.join(new_line))

    except IndexError:
        header_out = ','.join(header)

        headers = {}
        for k, v in final_dict.items():
            tem = []
            for kk, vv in v.items():
                tem.append(vv)
                headers[kk] = 1
            new_line = [k] + tem
            print (','.join(new_line))


if __name__ == '__main__':
    args = get_args()
    frontfile = args.frontfile
    backfile = args.backfile
    savefile = args.savefile
    listfile = args.listfile

    if frontfile and backfile:
        if savefile:
            front_back_save = ','.join([frontfile]+[backfile]+[savefile])
            main(front_back_save)
        else:
            front_back = ','.join([frontfile]+[backfile])
            main(front_back)

    elif listfile:
        with open(listfile,'r') as lf:
            for lf_row in lf:
                lf_row = lf_row.rstrip().split('\t')
                frontfile = lf_row[0]
                backfile = lf_row[1]
                if savefile:
                    front_back_save = ','.join([frontfile]+[backfile]+[savefile])
                    main(front_back_save)
                else:
                    front_back = ','.join([frontfile]+[backfile])
                    main(front_back)
    else:
        sys.exit(1)


