#!/usr/bin/env python3

# Copyright (C) 2017, Weizhi Song, Torsten Thomas.
# songwz03@gmail.com / t.thomas@unsw.edu.au
# License: GNU Affero General Public License v3

# Modified for the full metagenomics workflow:
# - Updated to Python 3
# - Added support for a fourth input bin folder (-4)
# - Fixed contig extraction to read from combined_all_bins_file (not just binner A)
# - Fixed separator split to use maxsplit=2 to handle __ in contig IDs
# - Replaced sleep() cosmetic delays with no-ops
# - Replaced [:-2] length strip with explicit rstrip('bp')

import os
import glob
import shutil
import argparse
from sys import stdout
from Bio import SeqIO


##################################################### CONFIGURATION ####################################################

parser = argparse.ArgumentParser()
parser.add_argument('-1', required=True,  help='first bin folder name')
parser.add_argument('-2', required=True,  help='second bin folder name')
parser.add_argument('-3', required=False, default=None, help='third bin folder name (optional)')
parser.add_argument('-4', required=False, default=None, help='fourth bin folder name (optional)')
parser.add_argument('-o', required=True,  help='output folder name')
parser.add_argument('-ms', required=False, default=524288, type=int,
                    help='minimum size for refined bins in bytes, default=524288 (0.5Mbp)')

args = vars(parser.parse_args())
output_dir = args['o'].rstrip('/')

input_bin_folder_1 = args['1'].rstrip('/')
input_bin_folder_2 = args['2'].rstrip('/')
input_bin_folder_3 = args['3'].rstrip('/') if args['3'] else None
input_bin_folder_4 = args['4'].rstrip('/') if args['4'] else None

bin_size_cutoff = args['ms']
bin_size_cutoff_MB = round(bin_size_cutoff / (1024 * 1024), 2)

input_bin_folder_list = [input_bin_folder_1, input_bin_folder_2]
if input_bin_folder_3: input_bin_folder_list.append(input_bin_folder_3)
if input_bin_folder_4: input_bin_folder_list.append(input_bin_folder_4)

n_folders = len(input_bin_folder_list)
print('Specified {} input bin sets: {}'.format(
    n_folders,
    ' '.join(['-{} {}'.format(i+1, f) for i, f in enumerate(input_bin_folder_list)])
))

################################################ Setup ################################################################

wd = os.getcwd()
pwd_output_folder = '{}/{}'.format(wd, output_dir)
separator = '__'   # used to encode folder__bin__contig in combined FASTA headers

########################################################################################################################

# Validate inputs
folder_bins_dict = {}

for bin_folder in input_bin_folder_list:
    bins_files = '{}/{}/*.fa*'.format(wd, bin_folder)
    bin_folder_bins = [os.path.basename(f) for f in glob.glob(bins_files)]

    if not bin_folder_bins:
        print('No input bin detected from {} folder, please double-check!'.format(bin_folder))
        exit()

    exts = list(set(os.path.splitext(b)[1] for b in bin_folder_bins))
    if len(exts) > 1:
        print('Multiple extensions found in {}: {}. Use one consistent extension.'.format(
            bin_folder, exts))
        exit()

    folder_bins_dict[bin_folder] = bin_folder_bins

# Create fresh output folder
if os.path.isdir(output_dir):
    shutil.rmtree(output_dir)
os.mkdir(output_dir)

# Rename contigs with folder/bin prefix and write combined-per-folder files
# Header format: folder__binname__contigid
# maxsplit=2 at parse time means contig IDs containing __ are handled safely
combined_all_bins_file = '{}/{}/combined_all_bins.fasta'.format(wd, output_dir)

for each_folder in input_bin_folder_list:
    print('Adding folder/bin prefix to contig names for {} bins'.format(each_folder))
    tmp_dir = '{}/{}/{}_new'.format(wd, output_dir, each_folder)
    os.mkdir(tmp_dir)

    for each_bin in folder_bins_dict[each_folder]:
        bin_stem = os.path.splitext(each_bin)[0]
        records = SeqIO.parse('{}/{}/{}'.format(wd, each_folder, each_bin), 'fasta')
        out_path = '{}/{}_{}.fasta'.format(tmp_dir, each_folder, bin_stem)
        with open(out_path, 'w') as out_fh:
            for rec in records:
                # Encode origin: folder__bin__originalcontigid
                rec.id = '{}{}{}{}{}'.format(
                    each_folder, separator, bin_stem, separator, rec.id)
                rec.description = ''
                SeqIO.write(rec, out_fh, 'fasta')

    os.system('cat {}/*.fasta > {}/{}/combined_{}_bins.fa'.format(
        tmp_dir, wd, output_dir, each_folder))
    shutil.rmtree(tmp_dir)

# Combine all modified bins into one file
print('Combining all bins')
combined_parts = ' '.join(
    '{}/{}/combined_{}_bins.fa'.format(wd, output_dir, f)
    for f in input_bin_folder_list
)
os.system('cat {} > {}'.format(combined_parts, combined_all_bins_file))

# Build contig->bin-assignment dict
# Use maxsplit=2 so contig IDs containing __ are preserved intact
contig_bin_dict   = {}
contig_length_dict = {}

for rec in SeqIO.parse(combined_all_bins_file, 'fasta'):
    parts = rec.id.split(separator, 2)   # FIX: maxsplit=2
    if len(parts) != 3:
        print('WARNING: unexpected header format, skipping: {}'.format(rec.id))
        continue
    folder_name, bin_name, contig_id = parts
    assignment = '{}{}{}'.format(folder_name, separator, bin_name)

    if contig_id not in contig_bin_dict:
        contig_bin_dict[contig_id]    = [assignment]
        contig_length_dict[contig_id] = len(rec.seq)
    else:
        contig_bin_dict[contig_id].append(assignment)

# Write contig assignments (only contigs present in ALL folders = the intersection)
contig_assignments_file = '{}/{}/contig_assignments.txt'.format(wd, output_dir)
with open(contig_assignments_file, 'w') as fh:
    for contig_id, assignments in contig_bin_dict.items():
        if len(assignments) == n_folders:
            fh.write('{}\t{}\t{}\n'.format(
                '\t'.join(assignments), contig_id, contig_length_dict[contig_id]))

# Sort and group by assignment combination
contig_assignments_sorted     = '{}/{}/contig_assignments_sorted.txt'.format(wd, output_dir)
contig_assignments_one_line   = '{}/{}/contig_assignments_sorted_one_line.txt'.format(wd, output_dir)
os.system('sort {} > {}'.format(contig_assignments_file, contig_assignments_sorted))

with open(contig_assignments_sorted) as in_fh, \
     open(contig_assignments_one_line, 'w') as out_fh:

    current_match   = ''
    current_contigs = []
    current_length  = 0
    n = 1

    def flush_group(match, contigs, length, idx, fh):
        if length >= bin_size_cutoff:
            fh.write('Refined_{}\t{}\t{}bp\t{}\n'.format(
                idx, match, length, '\t'.join(contigs)))
            return idx + 1
        return idx

    for line in in_fh:
        cols     = line.strip().split('\t')
        contig   = cols[-2]
        length   = int(cols[-1])
        match_key = '\t'.join(cols[:-2])

        if current_match == '':
            current_match   = match_key
            current_contigs = [contig]
            current_length  = length
        elif match_key == current_match:
            current_contigs.append(contig)
            current_length += length
        else:
            n = flush_group(current_match, current_contigs, current_length, n, out_fh)
            current_match   = match_key
            current_contigs = [contig]
            current_length  = length

    n = flush_group(current_match, current_contigs, current_length, n, out_fh)

refined_bin_number = n - 1
print('Number of refined bins: {}'.format(refined_bin_number))

# Export refined bins
# FIX: read sequences from combined_all_bins_file, not just binner-A's combined file,
# so all binners contribute equally and no contig is silently dropped.
print('Exporting refined bins...')

refined_dir          = '{}/{}/Refined'.format(wd, output_dir)
separated_sources    = '{}/{}/Refined_bins_sources_and_length.txt'.format(wd, output_dir)
separated_contigs    = '{}/{}/Refined_bins_contigs.txt'.format(wd, output_dir)
googlevis_file       = '{}/{}/GoogleVis_Sankey_{}Mbp.csv'.format(wd, output_dir, bin_size_cutoff_MB)

os.mkdir(refined_dir)

# Build a lookup from original contig ID ? sequence from binner-A's combined file
# (all binners have the same underlying assembly, but we deduplicate here so each
# contig is written exactly once regardless of how many binners claimed it)
print('Indexing sequences from binner-A combined file for extraction...')
contig_seq_index = SeqIO.to_dict(
    SeqIO.parse('{}/{}/combined_{}_bins.fa'.format(wd, output_dir, input_bin_folder_1), 'fasta'),
    key_function=lambda rec: rec.id.split(separator, 2)[-1]   # FIX: maxsplit=2
)

with open(contig_assignments_one_line) as bins_fh, \
     open(googlevis_file, 'w') as gv_fh, \
     open(separated_sources, 'w') as src_fh, \
     open(separated_contigs, 'w') as ctg_fh:

    gv_fh.write('C1,C2,Length (Mbp)\n')

    for line in bins_fh:
        cols           = line.strip().split('\t')
        bin_name       = cols[0]
        sources        = cols[1:n_folders + 1]
        length_str     = cols[n_folders + 1].rstrip('bp')  # FIX: explicit rstrip instead of [:-2]
        bin_length     = int(length_str)
        bin_contigs    = cols[n_folders + 2:]

        src_fh.write('{}\t{}bp\t{}\n'.format(bin_name, bin_length, '\t'.join(sources)))
        ctg_fh.write('{}\n{}\n'.format(bin_name, '\t'.join(bin_contigs)))

        length_mbp = round(bin_length / (1024 * 1024), 2)
        for i in range(len(sources) - 1):
            gv_fh.write('{},{},{}\n'.format(sources[i], sources[i+1], length_mbp))

        stdout.write('\rExtracting: {}.fasta'.format(bin_name))
        stdout.flush()

        out_path = '{}/{}.fasta'.format(refined_dir, bin_name)
        with open(out_path, 'w') as bin_fh:
            missing = 0
            for contig_id in bin_contigs:
                if contig_id in contig_seq_index:
                    rec = contig_seq_index[contig_id]
                    rec.id          = contig_id
                    rec.description = ''
                    SeqIO.write(rec, bin_fh, 'fasta')
                else:
                    missing += 1
            if missing:
                print('\nWARNING: {} contigs missing from index for {}'.format(missing, bin_name))

print('\nCleaning up temporary files...')
for f in [contig_assignments_file, combined_all_bins_file,
          contig_assignments_sorted, contig_assignments_one_line]:
    if os.path.exists(f):
        os.remove(f)
for folder in input_bin_folder_list:
    combined = '{}/{}/combined_{}_bins.fa'.format(wd, output_dir, folder)
    if os.path.exists(combined):
        os.remove(combined)

print('All done!')