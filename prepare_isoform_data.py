#! /usr/bin/env python


'''
Input format:
    chr1|1,10,32|5,20,50|+
    chr3|30,45|25,60|-

Externel tools:
    - bedtools
'''


import argparse
import re
import os
import os.path
import subprocess as sp
import io
import tempfile as tp
from contextlib import contextmanager
from collections import namedtuple
from functools import partial
from operator import itemgetter



class Bed:
    BED_TITLE = (
        'chr',
        'start',
        'end',
        'name',
        'score',
        'strand',
        'thickStart',
        'thickEnd',
        'itemRGB',
        'blockCount',
        'blockSizes',
        'blockStarts'
    )
    _region = namedtuple('Region', ('start', 'end'))

    def __init__(self, regions, name='.', union=False):
        self.regions = regions
        self.name = name
        self.union = union

        if not self.union:
            self.regions = self.reverse_regions_if_minus_strand(self.regions)
        else:
            self.regions = self.get_union_regions(self.regions)

        self._parse_regions(self.regions)

    @staticmethod
    def reverse_regions_if_minus_strand(regions):
        strand = regions[0][3]
        if strand == '-':
            regions = list(reversed(regions))

        return regions

    def _parse_regions(self, regions):
        self._parse_region(regions[0])
        for region in regions[1:]:
            self._add_block(region)

    def _parse_region(self, region):
        self.chrom = region[0]
        self.start = int(region[1]) - 1
        self.end = int(region[2])
        self.strand = region[3]

        self.block_count = 1
        self.block_sizes = [self.end - self.start]
        self.block_starts = [0]

    def _add_block(self, region):
        other = Bed([region])

        assert self.strand == other.strand

        self.block_count += other.block_count
        self.block_sizes += other.block_sizes

        all_starts = [
            start + self.start for start in self.block_starts
        ] + [other.start]

        if other.start < self.start:
            self.start = other.start

        self.block_starts = [start - self.start for start in all_starts]

        if other.end > self.end:
            self.end = other.end

    def get_data(self, all_fields=False):
        fields = [
            self.chrom,
            self.start,
            self.end,
            self.name,
            '.',
            self.strand
        ]

        if all_fields:
            fields += [
                self.start,
                self.end,
                0,
                self.block_count,
                self._list_to_str(self.block_sizes, sep=','),
                self._list_to_str(self.block_starts, sep=',')
            ]

        return fields

    def to_string(self, all_fields=False):
        data = self.get_data(all_fields=all_fields)
        bed_txt = self._list_to_str(data)
        return bed_txt

    @staticmethod
    def _list_to_str(list_, sep='\t', end=''):
        with io.StringIO() as tmp:
            print(*list_, sep=sep, end=end, file=tmp)
            return tmp.getvalue()

    @classmethod
    def get_union_regions(cls, regions):
        assert len(set(map(itemgetter(0), regions))) == 1, \
            "Not all regions in the same chromosome!"
        assert len(set(map(itemgetter(3), regions))) == 1, \
            "Not all regions at the same strand!"

        chr_ = regions[0][0]
        strand = regions[0][3]

        regions = sorted(
            map(
                lambda r: cls._region(r[1], r[2]),
                regions
            ),
            key=lambda r: r.start
        )

        union_regions = []
        r1 = regions[0]
        for r2 in regions[1:]:
            if r1.end < r2.start:
                union_regions.append(r1)
                r1 = r2
            else:
                if r1.end < r2.end:
                    r1 = cls._region(r1.start, r2.end)
        else:
            union_regions.append(r1)

        union_regions = tuple((chr_, start, end, strand)
                              for start, end in union_regions)
        return union_regions


def get_fasta(bed_file, genome_file, use_blocks=False, bedtools_bin='bedtools'):
    with tp.NamedTemporaryFile(dir='.') as tmp_file:
        cmd = [bedtools_bin, 'getfasta']
        cmd += ['-fi', genome_file]
        cmd += ['-bed', bed_file]
        cmd += ['-fo', tmp_file.name]
        cmd += ['-name', '-s', '-tab']

        if use_blocks:
            cmd += ['-split']

        sp.run(cmd)

        with open(tmp_file.name) as fa_in:
            for line in fa_in:
                name, fa_seq = line.rstrip('\n').split('\t')
                name = re.sub(r'\([+-]\)$', '', name)
                name = re.sub(r'::.*$', '', name)
                yield name, fa_seq


@contextmanager
def cwd(path):
    origin_pwd = os.getcwd()
    os.chdir(path)
    yield
    os.chdir(origin_pwd)


Isoform = namedtuple('Isoform', ['chr_', 'starts', 'ends', 'strand'])


def parse_isoform_data(isoform_string):
    chr_, starts, ends, strand = isoform_string.split('|')
    starts = list(map(int, starts.split(',')))
    ends = list(map(int, ends.split(',')))
    isoform_data = Isoform(chr_, starts, ends, strand)

    return isoform_data


def output_isoform_as_bed(isoform, isoform_name='.'):
    assert len(isoform.starts) == len(isoform.ends)
    assert all([i <= j for i, j in zip(isoform.starts, isoform.ends)])

    regions = [
        (isoform.chr_, start, end, isoform.strand)
        for start, end in zip(isoform.starts, isoform.ends)
    ]

    bed = Bed(regions, name=isoform_name, union=True)

    return bed


def create_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
 Input format:
    chr1|1,10,32|5,20,50|+
    chr3|30,45|25,60|-

 Externel tools:
    - bedtools"""
    )
    parser.add_argument('-g', '--ref_file', required=True)
    parser.add_argument('-i', '--circ_file', type=argparse.FileType('r'), required=True)
    parser.add_argument('-o', '--out_dir', default='isoform_data')
    parser.add_argument('--bedtools_bin', default='bedtools')

    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    ref_file = os.path.abspath(args.ref_file)

    with cwd(args.out_dir):
        with open('isoforms.bed', 'w') as bed_out, open('isoforms.length.tsv', 'w') as bed_len_out:
            for isoform_string in args.circ_file:
                isoform_string = isoform_string.rstrip('\n')
                isoform_data = parse_isoform_data(isoform_string)
                isoform_bed = output_isoform_as_bed(isoform_data, isoform_string)
                print(isoform_bed.to_string(all_fields=True), file=bed_out)

                total_len = sum(isoform_bed.block_sizes)
                print(isoform_string, total_len, sep='\t', file=bed_len_out)

        with open('isoforms.fa', 'w') as seq_out, open('isoforms.ext.fa', 'w') as seq_ext_out:
            fasta_data = get_fasta('isoforms.bed', ref_file, use_blocks=True, bedtools_bin=args.bedtools_bin)

            for fa_id, fa_seq in fasta_data:
                print(f'>{fa_id}', file=seq_out)
                print(fa_seq, file=seq_out)

                print(f'>{fa_id}', file=seq_ext_out)
                print(fa_seq + fa_seq[:30], file=seq_ext_out)
