#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" ESKRIM: EStimate with K-mers the RIchness in a Microbiome """

import logging
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
except ImportError as import_err:
    raise RuntimeError('Python bindings of jellyfish are not installed') from import_err

__author__ = "Florian Plaza Oñate"
__copyright__ = "Copyright 2019-2023, INRAE"
__maintainer__ = "Florian Plaza Oñate"
__email__ = "florian.plaza-onate@inrae.fr"
__status__ = "Production"
__licence__ = "GNU GPLv3"
__version__ = "1.0.8"

eskrim_version = f'ESKRIM {__version__}'


def setup_logger():
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logging.basicConfig(format='%(asctime)s :: %(levelname)s :: %(message)s')


def check_program_available(program):
    try:
        subprocess.call([program], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except OSError:
        logging.error('%s not found or not in system path', program)
        sys.exit(1)


def num_threads_type(value):
    max_num_threads = multiprocessing.cpu_count()

    try:
        value = int(value)
    except ValueError as value_err:
        raise argparse.ArgumentTypeError('NUM_THREADS is not an integer') from value_err

    if value <= 0:
        raise argparse.ArgumentTypeError('minimum NUM_THREADS is 1')
    if value > max_num_threads:
        raise argparse.ArgumentTypeError(f'maximum NUM_THREADS is {max_num_threads}')

    return value


def readable_writable_dir(path):
    if not os.path.isdir(path):
        raise NotADirectoryError(path)

    if not os.access(path, os.R_OK):
        raise PermissionError(f"directory {path} is not readable")

    if not os.access(path, os.W_OK):
        raise PermissionError(f"directory {path} is not writable")

    return path


def get_parameters():
    """Parse command line parameters.
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', dest='input_fastq_files', nargs='+', required=True,
                        help=('INPUT_FASTQ_FILES with reads from a single metagenomic sample '
                              '(gzip and bzip2 compression accepted)'))

    parser.add_argument('-n', dest='sample_name',
                        default='NA',
                        help='name of the metagenomic sample')

    parser.add_argument('-l', dest='read_length',
                        type=int, default=80,
                        help=('discard reads shorter than READ_LENGTH bases '
                              'and trim those exceeding this length'))

    parser.add_argument('-r', dest='target_num_reads',
                        type=int, default=10000000,
                        help='TARGET_NUM_READS to draw randomly from INPUT_FASTQ_FILES')

    parser.add_argument('-k', dest='kmer_length',
                        type=int, choices=range(17, 32, 2), default=21,
                        help='length of kmers to count')

    parser.add_argument('-t', dest='num_threads',
                        type=num_threads_type, default=4,
                        help='NUM_THREADS to launch for kmers counting')

    parser.add_argument('-o',
                        dest='output_fastq_file', type=argparse.FileType('w'),
                        help='OUTPUT_FASTQ_FILE with the randomly selected reads')

    parser.add_argument('-s', dest='output_stats_file',
                        type=argparse.FileType('w'), required=True,
                        help='OUTPUT_STATS_FILE with kmer richness estimates')

    parser.add_argument('--tmp-dir', dest='temp_dir',
                        type=readable_writable_dir, default='/dev/shm',
                        help='Temporary directory to store the jellyfish database')

    parser.add_argument('--seed', dest='rng_seed',
                        type=int, default=0,
                        help='Seed for random number generator')

    parser.add_argument('--version', action='version', version=eskrim_version)

    return parser.parse_args()


def hook_compressed_text(filename, mode, encoding='utf8'):
    ext = os.path.splitext(filename)[1]
    if ext == '.gz':
        import gzip
        return gzip.open(filename, mode + 't', encoding=encoding)
    if ext == '.bz2':
        import bz2
        return bz2.open(filename, mode + 't', encoding=encoding)

    return open(filename, mode, encoding=encoding)


def check_fastq_files(fastq_files):
    fastq_extensions = {
        '.fastq', '.fq', '.fastq.gz', '.fq.gz', '.fastq.bz2', '.fq.bz2'}
    fastq_extensions_escape = (
        str.replace(fastq_extension, '.', r'\.') for fastq_extension in fastq_extensions)
    regexp_match_fastq_extensions = \
        '({all_extensions})$'.format(all_extensions='|'.join(fastq_extensions_escape))

    problematic_fastq_files = [os.path.basename(fastq_file) for fastq_file in fastq_files
                               if re.sub(regexp_match_fastq_extensions, '', fastq_file).endswith(
                                   ('_2', '.2', '_R2'))]
    problematic_fastq_files = ','.join(problematic_fastq_files)

    if problematic_fastq_files:
        logging.warning('Input FASTQ files probably contain reverse reads (%s)',
                        problematic_fastq_files)
        logging.warning('Use only forward reads for accurate results\n')


FastqEntry = namedtuple('FastqEntry', ['name', 'seq', 'qual'])


def fastq_reader(fastq_fi, target_read_length):
    while fastq_fi:
        try:
            name = next(fastq_fi).rstrip('\n')

            if fastq_fi.isfirstline():
                logging.info('Reading %s', fastq_fi.filename())
        except StopIteration:
            return
        try:
            seq = next(fastq_fi).rstrip('\n')
            next(fastq_fi)
            qual = next(fastq_fi).rstrip('\n')
        except StopIteration:
            logging.error('FASTQ file %s is truncated or invalid', fastq_fi.filename())
            sys.exit(1)

        fastq_reader.total_num_reads += 1

        if 'N' in seq:
            fastq_reader.num_Ns_reads_ignored += 1
        elif len(seq) < target_read_length:
            fastq_reader.num_too_short_reads_ignored += 1
        else:
            yield FastqEntry(name, seq[0:target_read_length], qual[0:target_read_length])


fastq_reader.total_num_reads = 0
fastq_reader.num_Ns_reads_ignored = 0
fastq_reader.num_too_short_reads_ignored = 0


def fastq_formatter(fastq_entry):
    return f'{fastq_entry.name}\n{fastq_entry.seq}\n+\n{fastq_entry.qual}\n'


def subsample_fastq_files(input_fastq_files, target_num_reads, target_read_length):
    with fileinput.FileInput(input_fastq_files, openhook=hook_compressed_text) as fastq_fi:
        # Fill reservoir
        selected_reads = list(
            itertools.islice(fastq_reader(fastq_fi, target_read_length), target_num_reads))

        for new_read in fastq_reader(fastq_fi, target_read_length):
            ind = random.randrange(0, fastq_reader.total_num_reads)
            if ind < target_num_reads:
                selected_reads[ind] = new_read

    if fastq_reader.total_num_reads == 0:
        logging.error('Input FASTQ files are empty')
        sys.exit(1)

    if len(selected_reads) == 0:
        logging.error('No reads of at least %d bp and no Ns are available '
                      'in input FASTQ files', target_read_length)
        sys.exit(1)

    if len(selected_reads) < target_num_reads:
        logging.warning('Only %d reads of at least %d bp and no Ns are available '
                        'in input FASTQ files', len(selected_reads), target_read_length)
        logging.warning('Selected read count (%d) is less than the target read count (%d)',
                        len(selected_reads), target_num_reads)
        logging.warning('Restart ESKRIM with a lower target read count '
                        'or discard the current sample')

    return selected_reads


def create_jf_db(reads, kmer_length, num_threads, temp_dir):
    jellyfish_db_path = tempfile.mkstemp(dir=temp_dir, suffix='.jf')[1]
    jellyfish_count_cmd = ['jellyfish', 'count',
                           '-C', '-m', str(kmer_length),
                           '-s', '1G',
                           '-t', str(num_threads),
                           '-o', jellyfish_db_path,
                           '/dev/stdin']

    with subprocess.Popen(jellyfish_count_cmd, stdin=subprocess.PIPE) as jellyfish_count_proc:
        for read in reads:
            jellyfish_count_proc.stdin.write(fastq_formatter(read).encode())

        jellyfish_count_proc.stdin.close()
        jellyfish_count_proc.wait()

    return jellyfish_db_path


def count_distinct_kmers(jellyfish_db_path):
    jellyfish_stats_cmd = ['jellyfish', 'stats',
                           jellyfish_db_path]

    with subprocess.Popen(jellyfish_stats_cmd, stdout=subprocess.PIPE) as jellyfish_stats_proc:
        jellyfish_stats_proc.stdout.readline()
        num_distinct_kmers = jellyfish_stats_proc.stdout.readline()
        jellyfish_stats_proc.wait()

    num_distinct_kmers = int(num_distinct_kmers.split()[-1])

    return num_distinct_kmers


def count_solid_kmers(jellyfish_db_path):
    jellyfish_stats_cmd = ['jellyfish', 'stats',
                           '-L', '2',
                           jellyfish_db_path]

    with subprocess.Popen(jellyfish_stats_cmd, stdout=subprocess.PIPE) as jellyfish_stats_proc:
        jellyfish_stats_proc.stdout.readline()
        num_solid_kmers = jellyfish_stats_proc.stdout.readline()
        jellyfish_stats_proc.wait()

    num_solid_kmers = int(num_solid_kmers.split()[-1])

    return num_solid_kmers


def count_mercy_kmers_aux(params):
    reads, jellyfish_db_path, read_length, kmer_length = params
    num_mercy_kmers = 0
    jellyfish_db = dna_jellyfish.QueryMerFile(jellyfish_db_path)
    read_num_kmers = read_length - kmer_length + 1

    for read in reads:
        cur_read_kmers_all_unique = True
        for i in range(0, read_num_kmers):
            mer = dna_jellyfish.MerDNA(read.seq[i:(i + kmer_length)])
            mer.canonicalize()
            cur_read_kmers_all_unique = bool(jellyfish_db[mer] == 1)

            if not cur_read_kmers_all_unique:
                break

        if cur_read_kmers_all_unique:
            num_mercy_kmers += read_num_kmers

    return num_mercy_kmers


def count_mercy_kmers(reads, jellyfish_db_path, read_length, kmer_length, num_threads):
    chunk_size = 200000
    chunks_parameters = ((reads[x:x + chunk_size], jellyfish_db_path, read_length, kmer_length)
                         for x in range(0, len(reads), chunk_size))

    with multiprocessing.Pool(num_threads) as pool:
        chunks_num_mercy_kmers = pool.imap_unordered(count_mercy_kmers_aux, chunks_parameters)
        final_num_mercy_kmers = sum(chunks_num_mercy_kmers)

    return final_num_mercy_kmers


def main():
    setup_logger()

    parameters = get_parameters()

    logging.info('%s\n', eskrim_version)

    check_program_available('jellyfish')
    check_fastq_files(parameters.input_fastq_files)
    random.seed(parameters.rng_seed)

    logging.info('Subsampling reads from FASTQ files')
    selected_reads = subsample_fastq_files(
        parameters.input_fastq_files, parameters.target_num_reads, parameters.read_length)

    prop_selected_reads = round(100.0 * len(selected_reads) / fastq_reader.total_num_reads, 2)
    logging.info('%s reads out of %d (%g%%) selected in input FASTQ files\n',
                 len(selected_reads), fastq_reader.total_num_reads, prop_selected_reads)

    if parameters.output_fastq_file:
        logging.info('Writing selected reads')
        for read in selected_reads:
            parameters.output_fastq_file.write(fastq_formatter(read))
        parameters.output_fastq_file.close()
        logging.info('Selected reads saved in %s\n', parameters.output_fastq_file.name)

    logging.info('Creating jellyfish database')
    jellyfish_db_path = create_jf_db(
        selected_reads, parameters.kmer_length, parameters.num_threads, parameters.temp_dir)
    logging.info('Jellyfish database saved in %s\n', jellyfish_db_path)

    logging.info('Counting distinct kmers')
    num_distinct_kmers = count_distinct_kmers(jellyfish_db_path)
    logging.info('%d distinct kmers found\n', num_distinct_kmers)

    logging.info('Counting solid kmers')
    num_solid_kmers = count_solid_kmers(jellyfish_db_path)
    logging.info('%d solid kmers found\n', num_solid_kmers)

    logging.info('Counting mercy kmers')
    num_mercy_kmers = count_mercy_kmers(
        selected_reads, jellyfish_db_path,
        parameters.read_length, parameters.kmer_length, parameters.num_threads)
    logging.info('%d mercy kmers found\n', num_mercy_kmers)

    print('\t'.join([
          'sample_name',
          'total_num_reads',
          'num_Ns_reads_ignored',
          'num_too_short_reads_ignored',
          'target_num_reads',
          'num_selected_reads',
          'read_length',
          'kmer_length',
          'num_distinct_kmers',
          'num_solid_kmers',
          'num_mercy_kmers']),
          file=parameters.output_stats_file)
    print('\t'.join([
          parameters.sample_name,
          str(fastq_reader.total_num_reads),
          str(fastq_reader.num_Ns_reads_ignored),
          str(fastq_reader.num_too_short_reads_ignored),
          str(parameters.target_num_reads),
          str(len(selected_reads)),
          str(parameters.read_length),
          str(parameters.kmer_length),
          str(num_distinct_kmers),
          str(num_solid_kmers),
          str(num_mercy_kmers)]),
          file=parameters.output_stats_file)
    parameters.output_stats_file.close()
    logging.info('Statistics saved in %s', parameters.output_stats_file.name)

    os.remove(jellyfish_db_path)
    logging.info('Jellyfish database removed')


if __name__ == '__main__':
    main()
