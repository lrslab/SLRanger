#!/usr/bin/env python
import re
import argparse
import pysam
import pandas as pd
# from ssw import AlignmentMgr
from pyssw.ssw_wrap import Aligner
from Bio.Seq import Seq
import random
import multiprocessing
from tqdm import tqdm
import time
from SLRanger.visualization import visualize_html

def fasta_to_dict(fasta_path):
    fasta_dict = {}
    with open(fasta_path, 'r') as f:
        current_key = None
        sequence_parts = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_key:
                    # 保存上一个序列
                    fasta_dict[current_key] = ''.join(sequence_parts)
                current_key = line[1:]  # 去掉 '>'
                sequence_parts = []
            else:
                sequence_parts.append(line)
        # 保存最后一个序列
        if current_key:
            fasta_dict[current_key] = ''.join(sequence_parts)
    return fasta_dict


def generate_mismatches(kmer, bases=['A', 'T', 'C', 'G']):
    mismatches = set()
    for i in range(len(kmer)):
        for base in bases:
            if base != kmer[i]:
                mismatched_kmer = kmer[:i] + base + kmer[i + 1:]
                mismatches.add(mismatched_kmer)
    return mismatches

def build_mismatch_index(sequences, k):
    mismatch_to_kmer = {}

    for seq_id, seq in sequences.items():
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            mismatches = generate_mismatches(kmer)
            for mismatch in mismatches:
                if mismatch not in mismatch_to_kmer:
                    mismatch_to_kmer[mismatch] = kmer

    return mismatch_to_kmer


def extract_kmers(sequences, k):
    kmer_id = {}
    # 遍历每个序列并为每个5-mer分配ID
    for key, seq in sequences.items():
        encoded_seq = []
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            encoded_seq.append(kmer)
        kmer_id[key] = encoded_seq
    return kmer_id

def find_matches(query, mismatch_to_kmer, kmer_to_refs, k):
    query_kmers = [query[i:i + k] for i in range(len(query) - k + 1)]
    matches = []

    for query_kmer in query_kmers:
        if query_kmer in mismatch_to_kmer:
            original_kmer = mismatch_to_kmer[query_kmer]
            if original_kmer in kmer_to_refs:
                match_info = kmer_to_refs[original_kmer]
                matches.extend(match_info)

    return matches

def longest_consecutive(nums):
    num_set = set(nums)
    longest_streak = 0
    max_end_value = None

    for num in num_set:
        if num - 1 not in num_set:
            current_num = num
            current_streak = 1

            while current_num + 1 in num_set:
                current_num += 1
                current_streak += 1

            if current_streak > longest_streak:
                longest_streak = current_streak
                max_end_value = current_num

    return longest_streak, max_end_value

def find_best_match(query, mismatch_to_kmer, kmer_to_refs, k):
    encoded_query = [query[i:i+k] for i in range(len(query) - k + 1)]
    ref_positions = []
    max_consecutive = 0
    best_consecutive_end = 0

    for query_kmer in encoded_query:
        if query_kmer in kmer_to_refs:
            if kmer_to_refs.count(query_kmer) > 1:
                ref_positions.extend(index for index, kmer in enumerate(kmer_to_refs) if kmer == query_kmer)
            else:
                ref_positions.append(kmer_to_refs.index(query_kmer))
        elif query_kmer in mismatch_to_kmer:
            original_kmer = mismatch_to_kmer[query_kmer]
            if original_kmer in kmer_to_refs:
                if kmer_to_refs.count(query_kmer) > 1:
                    ref_positions.extend(index for index, kmer in enumerate(kmer_to_refs) if kmer == original_kmer)
                else:
                    ref_positions.append(kmer_to_refs.index(original_kmer))
    max_intersection = len(set(ref_positions))
    if max_intersection > 0:
        max_consecutive, best_consecutive_end = longest_consecutive(ref_positions)

    if max_consecutive > len(query) - k:
    # if max_intersection > max_consecutive:
        max_intersection = max_consecutive

    return max_intersection, max_consecutive, best_consecutive_end


def soft_processed(sequence, strand, cigar):
    """
    根据 strand 和 cigar 信息提取指定长度的序列部分。

    参数:
    - sequence (str): 原始的序列字符串。这个sequence已经被reverse处理过了，与正链一致，所以取序列的时候都从开头取起
    - strand (str): 链方向，"+" 或 "-"。
    - cigar (list of tuples): CIGAR 信息，表示为列表的元组 (operation, length)。
                              例如 [(4, 5), (0, 50)] 表示 5 bp 的 soft clipping，接着是 50 bp 的匹配。
    返回:
    - str: 根据条件截取后的序列。再多加两个碱基在后面
    """
    if not cigar:
        return ""  # Return empty string if no CIGAR info

    if strand == "+":
        if cigar[0][0] == 4:
            soft_clip_length = cigar[0][1]
            # 返回从序列开头截取 soft_clip_length 长度的序列部分
            return sequence[:soft_clip_length + 2]  #
    elif strand == "-":
        if cigar[-1][0] == 4:
            soft_clip_length = cigar[-1][1]
            # 返回从序列开头截取 soft_clip_length 长度的序列部分
            return sequence[:soft_clip_length + 2]  #
    else:
        # 不符合条件时，返回原始序列
        return ""


def ssw_wrapper(seq1, seq2, match=1, mismatch=1, gap_open=1, gap_extend=1):
    """
    parameter are write inside the function
    seq1 is ref and seq2 is query
    # todo : leave a api to change matrix
    """

    ref_seq = str(seq1)
    read_seq = str(seq2)
    # reduce the gap open score from 3 to 1 for nanopore reads
    aligner = Aligner(ref_seq,
                      match, mismatch, gap_open, gap_extend,
                      report_cigar=True, )

    aln = aligner.align(read_seq)  # min_score=20, min_len=10)

    return aln

def consensus(ref, query, cigar_string, ref_shift=0, query_shift = 0):
    """
    to get a concensus sequence as the longest cons of both sequence,
    and fill all disagreement with N
    for example,
    query:AAATA-TAGAA
    ref:  AAACACTA-AA
    cons: AAANANCANAA

    use query to get the con
    then get two chain file, 0 based
    """

    # pattern = re.compile('([0-9]*)([MIDNSHP=X])')
    pattern = re.compile(r"([0-9]+)b?'?([MIDNSHP=X])")
    query_out = []
    ref_out = []
    cons = []

    query_pos = query_shift
    ref_pos = ref_shift

    cigar_pattern = pattern.findall(cigar_string)
    for i, (length, code) in enumerate(cigar_pattern):
        if i == 0 or i == len(cigar_pattern) - 1:
            if code == "S":
                continue

        length = int(length)

        if code == "M":
            for i in range(length):
                q_s = query[query_pos].upper()
                try:
                    r_s = ref[ref_pos].upper()
                except IndexError as e:
                    print
                    "indexerror", ref_pos

                query_out.append(q_s)
                ref_out.append(r_s)

                cons.append(q_s) if q_s == r_s else cons.append("m")

                ref_pos += 1
                query_pos += 1

        elif code in "D":
            for i in range(length):
                r_s = ref[ref_pos]

                ref_out.append(r_s)
                query_out.append("-")
                cons.append("d")

                ref_pos += 1

        elif code in "IHS":
            for i in range(length):
                q_s = query[query_pos]

                query_out.append(q_s)
                ref_out.append("-")
                cons.append("i")
                query_pos += 1

    return "".join(cons)   # "".join(cons), "".join(query_out), "".join(ref_out)

def score_calculate(sw_score, max_intersection, max_consecutive, ref_end, ref_length, seq_start, seq_length, soft_length):
    # factor = get_value(best_end, ref_length - k)

    score_seq = (1 - seq_start/soft_length) # * (seq_length/soft_length)  # 在soft clipping中的位置，越靠头部越好
    # score_ref = 1 - (ref_length - ref_end)/ref_length
    score_ref = 1 - (ref_length - ref_end)/25   #(SL1_len - k)*2
    score_kmer = 0.5 * (max_intersection + max_consecutive - 1)
    score = score_seq * score_ref * (sw_score + score_kmer)

    final_score = max(score, 0)
    return final_score

def length_index(SL, ref_seq, kmer, mismatch_to_kmer, k):
    dict = {}

    for i in range(k, len(ref_seq) + 1):
        ref_seq_s = ref_seq[-i:]
        sw_aln = ssw_wrapper(ref_seq, ref_seq_s)  # score>10,
        sw_score = sw_aln.score
        # ref_start = sw_aln.ref_begin
        # ref_end = sw_aln.ref_end
        read_start = sw_aln.query_begin
        read_end = sw_aln.query_end + 1
        # sw_cigar = sw_aln.cigar_string

        corrected_sequence_sw = ref_seq_s
        # last_three_chars = corrected_sequence_sw[-3:]
        seq_s_length = len(corrected_sequence_sw)
        ref_length = len(ref_seq)
        seq_length = read_end - read_start
        kmer_to_refs = kmer[SL]
        max_intersection, max_consecutive, best_end = find_best_match(corrected_sequence_sw, mismatch_to_kmer,
                                                                      kmer_to_refs, k)
        final_score = score_calculate(sw_score, max_intersection, max_consecutive, ref_length, ref_length,
                                      read_start, seq_s_length, len(ref_seq_s))
        dict[i] = final_score

    return dict

def final_score_process(final_score, ref_length, seq_s_length, length_score):
    if seq_s_length <= ref_length:
        SL_score = final_score * (final_score / length_score[seq_s_length])
        # final_score_normalized = 100 * (final_score / length_score[seq_s_length])
    else:
        SL_score = final_score * (final_score / (seq_s_length / ref_length * length_score[ref_length]))
        # final_score_normalized = 100 * (final_score / (seq_s_length / ref_length * length_score[ref_length]))
    return SL_score # final_score_normalized

def random_score(random_sequences_dict, random_kmer, random_mismatch_to_kmer, length_score, corrected_sequence, k):
    random_sw_score_max = 0
    random_final_score_max = 0
    random_SL_score_max = 0
    for key, random_seq in random_sequences_dict.items():
        random_length_score = length_score
        random_aln_sw = ssw_wrapper(random_seq, corrected_sequence)
        random_sw_score = random_aln_sw.score
        random_read_start = random_aln_sw.query_begin
        random_read_end = random_aln_sw.query_end + 1
        corrected_sequence_sw = corrected_sequence[random_read_start:random_read_end]
        seq_s_length = len(corrected_sequence_sw)

        if seq_s_length >= 5:
            random_ref_end = random_aln_sw.ref_end + 1
            random_ref_length = len(random_seq)
            kmer_to_refs = random_kmer[key]
            max_intersection, max_consecutive, best_end = find_best_match(corrected_sequence_sw,
                                                                          random_mismatch_to_kmer,
                                                                          kmer_to_refs, k)
            random_final_score = score_calculate(random_sw_score, max_intersection, max_consecutive,
                                                 random_ref_end, random_ref_length, random_read_start,
                                                 seq_s_length, len(corrected_sequence))
            random_SL_score = final_score_process(random_final_score, random_ref_length,
                                                  seq_s_length, random_length_score)
            if random_SL_score > random_SL_score_max:
                random_sw_score_max = random_sw_score
                random_final_score_max = random_final_score
                random_SL_score_max = random_SL_score

    return random_sw_score_max, random_final_score_max, random_SL_score_max

def update(*a):
    if a[0] is not None:
        outfile.write(a[0])
    pbar.update(1)

def calculation_per_process(item,sl_dict,length_scores,random_sequences_dict ,random_kmer,random_mismatch_to_kmer,k,kmer,mismatch_to_kmer):
    query_name = item[0]
    query_seq = item[1]
    strand = item[2]

    corrected_sequence = soft_processed(query_seq, strand, item[3])
    aligned_len = item[4]

    if corrected_sequence is not None:
        soft_length = len(corrected_sequence)
        mes = query_name + '\t' + strand + '\t' + str(soft_length) + '\t' + str(aligned_len) + '\t'  # with soft_processed
    else:
        soft_length = 0
        mes = query_name + '\t' + strand + '\t' + 'NA' + '\t' + str(aligned_len) + '\t'  # with soft_processed

    if soft_length < k:  # 太短了就不要了
        SL_sw_dict = {'read_end': 'NA', 'query_length': 'NA',
                      'consensus': 'NA', 'random_sw_score': 'NA',
                      'random_final_score': 'NA', 'random_SL_score': 'NA',
                      'sw_score': 'NA', 'final_score': 'NA',
                      'SL_score': 'NA'}
        mes = mes + '\t'.join([str(value) for value in SL_sw_dict.values()]) + '\t' + 'random' + "\n"
        return mes

    ### use SL1 seq for SW check
    SL_sw_dict = {}
    for SL, SEQ in sl_dict.items():
        length_score = length_scores[SL]
        # corrected_sequence = 'CAAG'
        sw_aln = ssw_wrapper(SEQ, corrected_sequence)  # score>10,
        sw_score = sw_aln.score
        ref_start = sw_aln.ref_begin
        ref_end = sw_aln.ref_end + 1
        read_start = sw_aln.query_begin
        read_end = sw_aln.query_end + 1
        sw_cigar = sw_aln.cigar_string

        corrected_sequence_sw = corrected_sequence[read_start:read_end]

        seq_s_length = len(corrected_sequence_sw)
        random_sw_score_max, random_final_score_max, random_SL_score_max = random_score(random_sequences_dict,
                                                                                        random_kmer,
                                                                                        random_mismatch_to_kmer,
                                                                                        length_scores['SL1'],
                                                                                        corrected_sequence, k)
        if seq_s_length >= k:
            sw_ref = SEQ[ref_start:ref_end]

            # print("\n".join(consensus(corrected_sequence_sw, sw_ref, sw_cigar, 0)))
            # last_three_chars = corrected_sequence_sw[-3:]
            con_seq = consensus(sw_ref, corrected_sequence_sw, sw_cigar, 0)
            ref_length = len(SEQ)
            kmer_to_refs = kmer[SL]
            max_intersection, max_consecutive, best_end = find_best_match(corrected_sequence_sw, mismatch_to_kmer,
                                                                          kmer_to_refs, k)
            final_score = score_calculate(sw_score, max_intersection, max_consecutive, ref_end, ref_length,
                                          read_start, seq_s_length,
                                          soft_length)  # sw_score, max_intersection, max_consecutive, ref_end, ref_length, seq_start, seq_length, soft_length
            SL_score = final_score_process(final_score, ref_length, seq_s_length, length_score)

            SL_sw_dict[SL] = {'read_end': read_end,
                              'query_length': seq_s_length,
                              'consensus': con_seq,
                              'random_sw_score': random_sw_score_max,
                              'random_final_score': round(random_final_score_max, 2),
                              'random_SL_score': round(random_SL_score_max, 2),
                              'sw_score': sw_score,
                              'final_score': round(final_score, 2),
                              'SL_score': round(SL_score, 2)
                              }
        else:
            SL_sw_dict[SL] = {'read_end': read_end,
                              'query_length': None,
                              'consensus': None,
                              'random_sw_score': 0,
                              'random_final_score': 0,
                              'random_SL_score': 0,
                              'sw_score': 0,
                              'final_score': 0,
                              'SL_score': 0
                              }

    SL_sw_df = pd.DataFrame.from_dict(SL_sw_dict, orient='index')
    SL_sw_df_s = SL_sw_df[SL_sw_df['SL_score'] > SL_sw_df['random_SL_score']]
    sl_max = SL_sw_df_s['SL_score'].max()
    # filtered_df = SL_sw_df_s[SL_sw_df_s['SL_score'] == sl_max]
    filtered_df = SL_sw_df_s[SL_sw_df_s['SL_score'] == sl_max]

    if len(filtered_df) > 1:
        mes = mes + '\t'.join([str(value) for value in filtered_df.iloc[0]]) + '\t' + filtered_df.index[
            0] + '_unknown' + "\n"
    elif len(filtered_df) == 1:
        mes = mes + '\t'.join([str(value) for value in filtered_df.iloc[0]]) + '\t' + filtered_df.index[0] + "\n"
    else:
        filtered_df_s = SL_sw_df[SL_sw_df['SL_score'] == SL_sw_df['SL_score'].max()]
        mes = mes + '\t'.join([str(value) for value in filtered_df_s.iloc[0]]) + '\t' + 'random' + "\n"
    return mes


def main(args):
    global outfile, pbar
    """
    NM comparison between SL1 and SL2
    read reads in bam
    write out a dataframe that including:
        query name, 22nt sequence, SW score, SL1 score, SL2 score, SL1 cigar, SL2 cigar ,SL type
        query name, selected 22nt sequence with soft clipping...
    """
    sl_dict = fasta_to_dict(args.ref)

    # 生成10个长度为SL mean长度的碱基的随机序列
    ref_lengths = [len(key) for key in sl_dict.values()]
    random_seq_len = round(sum(ref_lengths) / len(ref_lengths))
    random.seed(826)
    random_sequences_dict = {}
    for i in range(10):
        random_sequence = ''.join([random.choice('AGTC') for _ in range(random_seq_len)])
        random_sequences_dict[i] = random_sequence
    min_len = 7
    k = 5
    kmer = extract_kmers(sl_dict, k)
    mismatch_to_kmer = build_mismatch_index(sl_dict, k)
    random_kmer = extract_kmers(random_sequences_dict, k)
    random_mismatch_to_kmer = build_mismatch_index(random_sequences_dict, k)

    length_scores = {}
    for SL, SEQ in sl_dict.items():
        length_score = length_index(SL, SEQ, kmer, mismatch_to_kmer, k)
        length_scores[SL] = length_score

    bam_file = pysam.AlignmentFile(args.bam, 'rb')
    timestamp = int(time.time())
    tmp_output_name = f"tmp_{timestamp}.csv"
    outfile = open(tmp_output_name, "w")
    outfile.write("query_name\tstrand\tsoft_length\taligned_length\tread_end\tquery_length\tconsensus\trandom_sw_score\trandom_final_score\trandom_SL_score\tsw_score\tfinal_score\tSL_score\tSL_type\n")

    # 迭代每个read
    print('Loading the BAM file')
    bam_list=[]
    for read in bam_file.fetch():
        if read.is_supplementary or read.is_secondary:
            continue
        ### sequence extract
        query_name = read.query_name     # 获取query的名称
        full_query_sequence = read.query_sequence  # 获取query的序列
        if read.is_reverse:
            query_seq = str(Seq(full_query_sequence).reverse_complement())   # 序列映射到负链，需要反向互补处理
            strand = '-'
        else:
            query_seq = full_query_sequence    # 序列映射到正链，直接使用
            strand = '+'
        bam_list.append([query_name,query_seq,strand,read.cigartuples,read.query_alignment_length])
    pbar = tqdm(total=len(bam_list), position=0, leave=True)
    with multiprocessing.Pool(processes=args.cpu) as pool:
        for item in bam_list:
            pool.apply_async(calculation_per_process, args=(item,sl_dict,length_scores,random_sequences_dict,random_kmer,random_mismatch_to_kmer,k,kmer,mismatch_to_kmer), callback=update)
        pool.close()
        pool.join()
    pbar.close()
    outfile.close()
    df = pd.read_csv(tmp_output_name,sep='\t')
    df.sort_values(by=['query_name'], inplace=True)
    df.to_csv(args.output, index=False,sep='\t')
    if args.visualization:
        visualize_html(args.output, args.cutoff)
    print('Finished')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="help to know spliced leader and distinguish SL1 and SL2")
    parser.add_argument("-r", "--ref", type=str, required=True, help="SL reference")
    parser.add_argument("-b", "--bam", type=str, required=True, help="input the bam file")
    parser.add_argument("-o", "--output", type=str, default="SLRanger.txt",
                        help="output file")
    parser.add_argument( "--visualization", action='store_true', help='Turn on the visualization mode')
    parser.add_argument("-t", "--cpu", type=int,
                        default=4, help="CPU number")
    parser.add_argument("-c", "--cutoff", type=float, default=4, help="cutoff of high confident SL sequence")
    args = parser.parse_args()
    main(args)
