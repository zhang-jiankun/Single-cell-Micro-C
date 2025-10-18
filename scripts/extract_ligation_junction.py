#!/usr/bin/env python

import argparse
import sys 
import pysam 

def report_ligation_junction(query_name, seg1, seg2, show_read_id=False):
    
    if not show_read_id: query_name = '.'

    if seg1.is_reverse == False:
        left_ligation_junction = [ seg1.reference_name, str(seg1.reference_end), str(seg1.reference_end + 1), query_name, ".", '-' if seg1.is_reverse else '+' ]
    else:
        left_ligation_junction = [ seg1.reference_name, str(seg1.reference_start), str(seg1.reference_start + 1), query_name, ".", '-' if seg1.is_reverse else '+' ]
    
    if seg2.is_reverse == False:
        right_ligation_junction = [ seg2.reference_name, str(seg2.reference_start), str(seg2.reference_start + 1), query_name, ".", '-' if seg2.is_reverse else '+' ]
    else:
        right_ligation_junction = [ seg2.reference_name, str(seg2.reference_end), str(seg2.reference_end + 1), query_name,  ".", '-' if seg2.is_reverse else '+' ]

    return left_ligation_junction, right_ligation_junction

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='extract ligation junction')
    parser.add_argument('sam', type=str, help='gzipped sam file')
    parser.add_argument('--readid', help='print read id', action='store_true')
    args = parser.parse_args()

    sam_file = pysam.AlignmentFile(args.sam, 'rb')
    show_read_id = args.readid
    
    read1_seg_data = dict()
    read2_seg_data = dict()
    
    for sam_read in sam_file.fetch(until_eof=True):
        
        if sam_read.is_unmapped or sam_read.is_qcfail or sam_read.is_duplicate: continue
        # supplement alignment
        if not sam_read.has_tag('SA'): continue
        
        ###
        if sam_read.mapping_quality < 20 or float(sam_read.get_tag("NM")) / (sam_read.query_alignment_end - sam_read.query_alignment_start) > 0.05: continue
        
        # separately for read1 and read2
        if sam_read.is_read1:
            if read1_seg_data.has_key(sam_read.query_name): read1_seg_data[sam_read.query_name].append(sam_read)
            else: read1_seg_data[sam_read.query_name] =  [ sam_read ]

        if sam_read.is_read2:
            if read2_seg_data.has_key(sam_read.query_name): read2_seg_data[sam_read.query_name].append(sam_read)
            else: read2_seg_data[sam_read.query_name] =  [ sam_read ]    

    for query_name in read1_seg_data.keys():

        n_seg = len(read1_seg_data[query_name])
        if n_seg != 2: continue # 1 parimary alignment, 1 supplement alignment
        seg1, seg2 = read1_seg_data[query_name]
        if seg1.is_supplementary and seg2.is_supplementary: continue
        if not seg2.is_supplementary:
            seg1, seg2 = seg2, seg1
        left_ligation_junction, right_ligation_junction = report_ligation_junction(query_name, seg1, seg2, show_read_id)

        sys.stdout.write('\t'.join(left_ligation_junction) + '\n')
        sys.stdout.write('\t'.join(right_ligation_junction) + '\n')
        
    for query_name in read2_seg_data.keys():
        
        n_seg = len(read2_seg_data[query_name])
        if n_seg != 2: continue
        seg1, seg2 = read2_seg_data[query_name]
        if seg1.is_supplementary and seg2.is_supplementary: continue
        if not seg2.is_supplementary:
            seg1, seg2 = seg2, seg1
        left_ligation_junction, right_ligation_junction = report_ligation_junction(query_name, seg1, seg2, show_read_id)

        sys.stdout.write('\t'.join(left_ligation_junction) + '\n')
        sys.stdout.write('\t'.join(right_ligation_junction) + '\n')
    