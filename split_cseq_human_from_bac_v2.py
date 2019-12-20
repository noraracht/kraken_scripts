import sys
from sys import stderr
import os
import re
import subprocess
import shlex


#run script with "python split_cseq_human_from_bac_v2.py -i NWR_NOLA.fastq -threads 24 -db /nucleus/pedigree/projects/norar/skimming/soft_install/database -confidence1 0.0 -confidence2 0.5"


if __name__ == "__main__":
    # parse user arguments
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=True, type=str, help="Input File")
    parser.add_argument('-threads', '--threads', required=False, type=int, default=1)
    parser.add_argument('-db', '--database', required=True, type=str, help="Kraken Database")
    parser.add_argument('-confidence1', '--confidence1', required=False, type=float, default=0.0, help="Confidence1")
    parser.add_argument('-confidence2', '--confidence2', required=False, type=float, default=0.5, help="Confidence2")
    args = parser.parse_args()


    if args.input == 'stdin':
        from sys import stdin;

        infile = stdin
    elif args.input.lower().endswith('.gz'):
        from gzip import open as gopen;

        infile = gopen(args.input)
    else:
        infile = open(args.input)

    # run kraken 1st time
    report_fname = "report_" + args.input
    output_fname = "output_" + args.input
    cseq_fname = "cseq_" + args.input
    ucseq_fname = "ucseq_" + args.input

    f = open(output_fname, "w")
    cmd = "kraken2 --use-names --threads {} --report {} --db {}  --confidence {} --classified-out {}  --unclassified-out {} {}".format(args.threads, report_fname, args.database, args.confidence1, cseq_fname, ucseq_fname, args.input)
    temp = subprocess.Popen(shlex.split(cmd), stdout=f)
    temp.communicate()
    f.close()


    try:
        if os.stat(cseq_fname).st_size > 0:
            pass
        else:
            print "empty file"
            sys.exit()
    except OSError:
        print ("no file")
        sys.exit()


    outfile = open("human_only.fq", 'w')
    taxid = 9606
    lines_to_print = [0, 0, 0]

    human_infile = open(cseq_fname, 'r')

    for i, line in enumerate(human_infile, 0):

        if "SRR" in line or "kraken:taxid" in line:
            if int(line.split("|")[-1]) == taxid:
                outfile.write(line)
                lines_to_print = [i+1, i+2, i+3]

        elif i in lines_to_print:
            outfile.write(line)

        else:
            pass


    try:
        if os.stat("human_only.fq").st_size > 0:
            pass
        else:
            print "no human reads found"
            sys.exit()
    except OSError:
        print ("no file")
        sys.exit()


    # run kraken 2nd time
    report_fname2 = "report2_" + args.input
    output_fname2 = "output2_" + args.input
    cseq_fname2 = "cseq2_" + args.input
    ucseq_fname2 = "ucseq2_" + args.input

    f = open(output_fname2, "w")
    cmd = "kraken2 --use-names --threads {} --report {} --db {}  --confidence {} --classified-out {}  --unclassified-out {} human_only.fq".format(args.threads, report_fname2, args.database, args.confidence2, cseq_fname2, ucseq_fname2) 
    temp = subprocess.Popen(shlex.split(cmd), stdout=f)
    temp.communicate()
    f.close()

    ucseq_final = "clean_" + args.input
    
    try:
        if os.stat(ucseq_fname2).st_size > 0:
            pass
        else:
            print("all reads passed second filter")
            os.remove(cseq_fname)
            os.remove(cseq_fname2)
            os.rename(ucseq_fname, ucseq_final)
            os.remove(ucseq_fname2)
            os.remove(output_fname)
            os.remove(output_fname2)
            os.remove("human_only.fq")
           
            archive = "report_"+ args.input + ".tar.gz"
            cmd = "tar czvf {} {} {}".format(archive, report_fname, report_fname2)
            temp = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            temp.communicate()
            f.close()
            
            os.remove(report_fname)
            os.remove(report_fname2)
            sys.exit()
    except OSError:
        print ("no file")
        sys.exit()




    # merge clean output
    f = open(ucseq_final, "w")
    cmd = "cat {} {}".format(ucseq_fname, ucseq_fname2)
    temp = subprocess.Popen(shlex.split(cmd), stdout=f)
    temp.communicate()
    f.close()


    infile.close()
    human_infile.close()
    outfile.close()

    os.remove(cseq_fname)
    os.remove(cseq_fname2)
    os.remove(ucseq_fname)
    os.remove(ucseq_fname2)
    os.remove(output_fname)
    os.remove(output_fname2)
    os.remove("human_only.fq")

    archive = "report_"+ args.input + ".tar.gz"
    #cmd = "tar czvf {} {} {} {} {} {} {} human_only.fq".format(archive, report_fname, report_fname2, cseq_fname, cseq_fname2, ucseq_fname, ucseq_fname2, output_fname, output_fname2)
    cmd = "tar czvf {} {} {}".format(archive, report_fname, report_fname2)
    temp = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    temp.communicate()
    f.close()
    os.remove(report_fname)
    os.remove(report_fname2)
