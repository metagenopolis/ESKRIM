#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""""""

import argparse
import os
import sys
import re
from collections import namedtuple
import fileinput 
import itertools
import random
import subprocess
import multiprocessing
import multiprocessing.pool
import tempfile
try:
    import dna_jellyfish
except ImportError as err:
    raise RuntimeError('Python bindings of jellyfish are not installed')

__author__ = "Florian Plaza Oñate"
__copyright__ = "Copyright 2019, INRA"
__maintainer__ = "Florian Plaza Oñate"
__email__ = "florian.plaza-onate@inra.fr"
__status__ = "Development"

def check_program_available(program):
    try:
        with open(os.devnull, 'w') as devnull:
            subprocess.call([program], stdout=devnull, stderr=devnull)
    except OSError:
        raise RuntimeError('{0} not found or not in system path'.format(program))

def num_threads_type(value):
    max_num_threads = multiprocessing.cpu_count()

    try:
        value = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError('NUM_THREADS is not an integer')

    if value <= 0:
        raise argparse.ArgumentTypeError('minimum NUM_THREADS is 1')
    elif value > max_num_threads:
        raise argparse.ArgumentTypeError('maximum NUM_THREADS is {}'.format(max_num_threads))

    return value

def get_parameters():
    """Parse command line parameters.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', dest='input_fastq_files', nargs='+', required=True,
            help='INPUT_FASTQ_FILES with reads from a single metagenomic sample (gzip and bzip2 compression accepted)')

    parser.add_argument('-s', dest='sample_name', default='NA',
            help='name of the metagenomic sample')

    parser.add_argument('-l', dest='read_length', type=int, default=80,
            help='discard reads shorter than READ_LENGTH bases and trim those exceeding this length')

    parser.add_argument('-r', dest='num_reads', type=int, default=10000000,
            help='NUM_READS to draw randomly from INPUT_FASTQ_FILES')

    parser.add_argument('-k', dest='kmer_length', type=int, choices=range(17,32,2), default=21, 
            help='length of kmers to count')

    parser.add_argument('-t', dest='num_threads', type=num_threads_type, default=4, 
            help='NUM_THREADS to launch for kmers counting')

    parser.add_argument('-o', dest='output_stats_file', type=argparse.FileType('w'), required=True, 
            help='OUTPUT_STATS_FILE with kmer richness estimates')

    return parser.parse_args()

def hook_compressed_text(filename, mode, encoding='utf8'):
    ext = os.path.splitext(filename)[1]
    if ext == '.gz':
        import gzip
        return gzip.open(filename, mode + 't', encoding=encoding)
    elif ext == '.bz2':
        import bz2
        return bz2.open(filename, mode + 't', encoding=encoding)
    else:
        return open(filename, mode, encoding=encoding)

def check_fastq_files(fastq_files):
    fastq_extensions = {'.fastq', '.fq', '.fastq.gz', '.fq.gz', '.fastq.bz2', '.fq.bz2'}
    regexp_match_fastq_extensions = '({all_extensions})$'.format(all_extensions = '|'.join(str.replace(fastq_extension, '.', '\.') for fastq_extension in fastq_extensions))
    problematic_fastq_files = [os.path.basename(fastq_file) for fastq_file in fastq_files if re.sub(regexp_match_fastq_extensions, '', fastq_file).endswith('2')]

    if problematic_fastq_files:
        print('warning: some files probably contain reverse reads ({})'.format(','.join(problematic_fastq_files)))
        print('warning: use only forward reads for accurate results\n')

FastqEntry=namedtuple('FastqEntry', ['name', 'seq', 'qual'])
def fastq_reader(istream, target_read_length):
    while istream:
        name = next(istream).rstrip('\n')
        seq = next(istream).rstrip('\n')
        next(istream)
        qual = next(istream).rstrip('\n')
        if len(seq) >= target_read_length:
            yield FastqEntry(name, seq[0:target_read_length], qual[0:target_read_length])

def fastq_formatter(fastq_entry):
    return f'@{fastq_entry.name}\n{fastq_entry.seq}\n+\n{fastq_entry.qual}\n'

def subsample_fastq_files(input_fastq_files, subsample_size, target_read_length):
    with fileinput.FileInput(input_fastq_files, openhook=hook_compressed_text) as fi:
        # Fill reservoir
        selected_reads = list(itertools.islice(fastq_reader(fi, target_read_length), subsample_size))

        if len(selected_reads) < subsample_size:
            raise RuntimeError('only {num_reads} reads of at least {read_length} bases are available in input FASTQ files.'.format(num_reads=len(selected_reads), read_length=target_read_length))

        for total_num_reads, new_read in enumerate(fastq_reader(fi, target_read_length), start=subsample_size+1):
            ind = random.randrange(0,total_num_reads)
            if ind < subsample_size:
                selected_reads[ind] = new_read

    return total_num_reads, selected_reads

def create_jf_db(reads, kmer_length, num_threads):
    jellyfish_db_path = tempfile.mkstemp(dir = '/dev/shm', suffix = '.jf')[1]
    jellyfish_count_cmd = ['jellyfish', 'count',
                           '-C', '-m', str(kmer_length),
                           '-s', '1G',
                           '-t', str(num_threads),
                           '-o', jellyfish_db_path,
                           '/dev/stdin']

    jellyfish_count_proc = subprocess.Popen(jellyfish_count_cmd, stdin=subprocess.PIPE)

    for read in reads:
        jellyfish_count_proc.stdin.write(fastq_formatter(read).encode())

    jellyfish_count_proc.stdin.close()
    jellyfish_count_proc.wait()
    
    return jellyfish_db_path

def count_distinct_kmers(jellyfish_db_path):
    jellyfish_stats_cmd = ['jellyfish', 'stats',
                           jellyfish_db_path]

    jellyfish_stats_proc = subprocess.Popen(jellyfish_stats_cmd, stdout=subprocess.PIPE)

    jellyfish_stats_proc.stdout.readline()
    num_distinct_kmers = jellyfish_stats_proc.stdout.readline()
    num_distinct_kmers = int(num_distinct_kmers.split()[-1])

    jellyfish_stats_proc.wait()

    return num_distinct_kmers

def count_solid_kmers(jellyfish_db_path):
    jellyfish_stats_cmd = ['jellyfish', 'stats',
                           '-L', '2',
                           jellyfish_db_path]

    jellyfish_stats_proc = subprocess.Popen(jellyfish_stats_cmd, stdout=subprocess.PIPE)

    jellyfish_stats_proc.stdout.readline()
    num_solid_kmers = jellyfish_stats_proc.stdout.readline()
    num_solid_kmers = int(num_solid_kmers.split()[-1])

    jellyfish_stats_proc.wait()

    return num_solid_kmers

def count_mercy_kmers_aux(params):
        reads, jellyfish_db_path, read_length, kmer_length = params
        num_mercy_kmers = 0
        jellyfish_db = dna_jellyfish.QueryMerFile(jellyfish_db_path)
        read_num_kmers = read_length - kmer_length + 1

        for read in reads:
            cur_read_kmers_all_unique = True
            for i in range(0, read_num_kmers):
                mer = dna_jellyfish.MerDNA(read.seq[i:(i+kmer_length)])
                mer.canonicalize()
                cur_read_kmers_all_unique = bool(jellyfish_db[mer] == 1)
                
                if not cur_read_kmers_all_unique:
                    break

            if cur_read_kmers_all_unique:
                num_mercy_kmers += read_num_kmers

        return num_mercy_kmers

def count_mercy_kmers(reads, jellyfish_db_path, read_length, kmer_length, num_threads, chunksize = 200000):
    chunks_parameters = ((reads[x:x+chunksize],jellyfish_db_path, read_length, kmer_length)
            for x in range(0, len(reads), chunksize))
    with multiprocessing.Pool(num_threads) as pool:
        chunks_num_mercy_kmers = pool.imap_unordered(count_mercy_kmers_aux, chunks_parameters)
        final_num_mercy_kmers = sum(chunks_num_mercy_kmers)
    
    return final_num_mercy_kmers

def main():
    parameters = get_parameters()

    check_program_available('jellyfish')
    check_fastq_files(parameters.input_fastq_files)

    print('Subsampling reads from FASTQ files...')
    total_num_reads, selected_reads = subsample_fastq_files(parameters.input_fastq_files, parameters.num_reads, parameters.read_length)
    print('Done. {num_selected_reads} reads out of {total_num_reads} selected ({proportion}%).\n'.format(
        num_selected_reads=len(selected_reads),
        total_num_reads=total_num_reads,
        proportion=round(100.0*len(selected_reads)/total_num_reads, 2)))

    print('Creating jellyfish database...')
    jellyfish_db_path = create_jf_db(selected_reads, parameters.kmer_length, parameters.num_threads)
    print('Done.\n')

    print('Counting distinct kmers...')
    num_distinct_kmers = count_distinct_kmers(jellyfish_db_path)
    m_distinct_kmers=0
    print('Done.\n')

    print('Counting solid kmers...')
    num_solid_kmers = count_solid_kmers(jellyfish_db_path)
    print('Done.\n')

    print('Counting mercy kmers...')
    num_mercy_kmers = count_mercy_kmers(selected_reads, jellyfish_db_path, parameters.read_length, parameters.kmer_length, parameters.num_threads)
    print('Done.\n')

    print('Printing output stats...')
    print('sample_name\ttotal_num_reads\tnum_selected_reads\tread_length\tkmer_length\tnum_distinct_kmers\tnum_solid_kmers\tnum_mercy_kmers', file=parameters.output_stats_file)
    print('\t'.join([
        parameters.sample_name,
        str(total_num_reads),
        str(parameters.num_reads),
        str(parameters.read_length),
        str(parameters.kmer_length),
        str(num_distinct_kmers),
        str(num_solid_kmers),
        str(num_mercy_kmers)]),
        file=parameters.output_stats_file)
    print('Done.\n')

    print('Cleanup...')
    os.remove(jellyfish_db_path)
    print('Done.')

if __name__ == '__main__':
    main()

