#!/usr/bin/env python
# -*- coding: utf-8 -*-

""""""

from __future__ import print_function
import argparse
import os
import subprocess
import multiprocessing
import tempfile
from collections import namedtuple
import dna_jellyfish

__author__ = "Florian Plaza Oñate"
__copyright__ = "Copyright 2019, INRA"
__maintainer__ = "Florian Plaza Oñate"
__email__ = "florian.plaza-onate@inra.fr"
__status__ = "Development"

jellyfish_bin='/export/mgps/home/fplazaonate/jellyfish-2.2.10/bin/jellyfish'

def check_program_available(program):
    try:
        with open(os.devnull, 'w') as devnull:
            subprocess.call([program], stdout=devnull, stderr=devnull)
    except OSError:
        raise RuntimeError('{0} not found or not in system path'.format(program))

def get_parameters():
    """Parse command line parameters.
    """
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('-i', dest='fastq_input_file', type=argparse.FileType('r'), required=True, 
            help='')

    parser.add_argument('-k', dest='kmer_length', type=int, choices=range(17,32,2), required=True, 
            help='')

    max_num_threads = multiprocessing.cpu_count()
    parser.add_argument('-t', dest='num_threads', type=int, choices=range(1,max_num_threads+1), default=max_num_threads, 
            help='')

    parser.add_argument('-o', dest='output_stats_file', type=argparse.FileType('w'), required=True, 
            help='')

    return parser.parse_args()

def create_jf_db(fastq_input_file, kmer_length, num_threads):
    jellyfish_db = tempfile.mkstemp(dir = '/dev/shm', suffix = '.jf')[1]
    jellyfish_count_cmd = [jellyfish_bin, 'count',
                           '-C', '-m', str(kmer_length),
                           '-s', '1G',
                           '-t', str(num_threads),
                           '-o', jellyfish_db,
                           fastq_input_file.name]

    subprocess.check_call(jellyfish_count_cmd)
    
    return jellyfish_db

def count_solid_kmers(jellyfish_db):
    jellyfish_stats_cmd = [jellyfish_bin, 'stats',
                           '-L', '2',
                           jellyfish_db]

    jellyfish_stats_proc = subprocess.Popen(jellyfish_stats_cmd, stdout=subprocess.PIPE)

    jellyfish_stats_proc.stdout.readline()
    num_distinct_kmers = jellyfish_stats_proc.stdout.readline()
    num_distinct_kmers = int(num_distinct_kmers.split()[-1])

    jellyfish_stats_proc.wait()

    return num_distinct_kmers

FastqEntry=namedtuple('FastqEntry', ['seq_id', 'seq_len', 'seq', 'qual'])
def fastq_reader(istream):
    while istream:
        seq_id = istream.next().rstrip('\n')
        seq = istream.next().rstrip('\n')
        istream.next()
        qual = istream.next().rstrip('\n')
        yield FastqEntry(seq_id, len(seq), seq, qual)

def count_mercy_kmers(fastq_input_file, jellyfish_db, kmer_length):
    num_mercy_kmers = 0
    qf = dna_jellyfish.QueryMerFile(jellyfish_db)

    for fastq_entry_id, fastq_entry in enumerate(fastq_reader(fastq_input_file), start=1):
        kmers_all_unique = True
        num_kmers = fastq_entry.seq_len - kmer_length + 1
        for i in range(0, num_kmers):
            mer = dna_jellyfish.MerDNA(fastq_entry.seq[i:(i+kmer_length)])
            mer.canonicalize()
            kmers_all_unique = bool(qf[mer] == 1)
            
            if not kmers_all_unique:
                break

        if kmers_all_unique:
            num_mercy_kmers += num_kmers

    return num_mercy_kmers

def print_output_stats(output_stats_file, num_solid_kmers, num_mercy_kmers):
    print('num_solid_kmers\t{}'.format(num_solid_kmers), file=output_stats_file)
    print('num_solid_and_mercy_kmers\t{}'.format(num_solid_kmers+num_mercy_kmers), file=output_stats_file)

def main():
    parameters = get_parameters()

    check_program_available(jellyfish_bin)

    print('Creating jellyfish database...')
    jellyfish_db = create_jf_db(parameters.fastq_input_file, parameters.kmer_length, parameters.num_threads)
    print('Done.\n')

    print('Counting solid kmers...')
    num_solid_kmers = count_solid_kmers(jellyfish_db)
    print('Done.\n')

    print('Counting mercy kmers...')
    num_mercy_kmers = count_mercy_kmers(parameters.fastq_input_file, jellyfish_db, parameters.kmer_length)
    print('Done.\n')

    print('Printing output stats...')
    print_output_stats(parameters.output_stats_file, num_solid_kmers, num_mercy_kmers)
    print('Done.\n')

    print('Cleanup...')
    os.remove(jellyfish_db)
    print('Done.\n')

if __name__ == '__main__':
    main()

