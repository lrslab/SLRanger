#!/usr/bin/env python
import re
import argparse
from pathlib import Path

import pandas as pd

try:
    from SLRanger.run_ex_function import run_track_cluster
except ImportError:
    run_track_cluster = None

# 解析GFF文件并构建DataFrame
def parse_gff(gff_file):
    genes = []
    with open(gff_file, 'r') as file:
        for line in file:
            if line.startswith("#") or line.strip() == '':
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            if parts[2] == 'gene':
                attr_field = parts[8]
                gene_id = None
                for attr in attr_field.split(';'):
                    if attr.startswith('ID='):
                        gene_id = attr.split('ID=')[1].split(';')[0]
                if gene_id is None:
                    continue
                genes.append({
                    'gene': gene_id,
                    'chromosome': parts[0],
                    'start': int(parts[3]),
                    'end': int(parts[4]),
                    'strand': parts[6]
                })
    return pd.DataFrame(
        genes,
        columns=['gene', 'chromosome', 'start', 'end', 'strand'],
    )

def parse_cds_gene(gff_file):
    """
    解析 GFF 文件，查找哪些 gene 具有 CDS。
    """
    gene_map = {}  # 存储 transcript -> gene 的映射
    cds_parents = set()  # 存储有 CDS 的 transcript ID
    genes_with_cds = set()  # 存储最终有 CDS 的 gene ID

    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue  # 跳过注释行

            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue  # 确保是完整的 GFF 行

            feature_type = fields[2]
            attributes = fields[8]

            if feature_type in ['transcript', 'mRNA']:
                # 提取 transcript ID 和其对应的 gene ID
                transcript_match = re.search(r'ID=([^;]+)', attributes)
                parent_match = re.search(r'Parent=([^;]+)', attributes)
                if transcript_match and parent_match:
                    transcript_id = transcript_match.group(1)
                    gene_id = parent_match.group(1)
                    gene_map[transcript_id] = gene_id

            elif feature_type == "CDS":
                # 提取 CDS 的 Parent（即 transcript）
                cds_match = re.search(r'Parent=([^;]+)', attributes)
                if cds_match:
                    transcript_id = cds_match.group(1).split(',')[0]
                    cds_parents.add(transcript_id)

    # 通过 transcript 找 gene
    for transcript_id in cds_parents:
        if transcript_id in gene_map:
            genes_with_cds.add(gene_map[transcript_id])

    # 转换为 DataFrame 输出
    df = pd.DataFrame({"gene": list(genes_with_cds)})
    return df

# 排序基因并计算基因之间的距离
def sort_and_calc_distance(df):
    output_columns = list(df.columns) + ['rank', 'intergenic_distance']
    output_columns = list(dict.fromkeys(output_columns))
    if df.empty:
        return pd.DataFrame(columns=output_columns)

    df_pos_list = []
    df_neg_list = []

    for chrom in df['chromosome'].unique():
        # 正链
        df_pos_chrom = df[(df['chromosome'] == chrom) & (df['strand'] == '+')].sort_values('end').reset_index(
            drop=True)
        df_pos_chrom['rank'] = range(1, len(df_pos_chrom) + 1)
        df_pos_chrom['intergenic_distance'] = df_pos_chrom['start'] - df_pos_chrom['end'].shift(1)
        df_pos_list.append(df_pos_chrom)

        # 负链
        df_neg_chrom = df[(df['chromosome'] == chrom) & (df['strand'] == '-')].sort_values('end',
                                                                                           ascending=False).reset_index(
            drop=True)
        df_neg_chrom['intergenic_distance'] = df_neg_chrom['start'].shift(1) - df_neg_chrom['end']
        df_neg_chrom['rank'] = range(1, len(df_neg_chrom) + 1)
        df_neg_list.append(df_neg_chrom)

    df_pos_final = pd.concat(df_pos_list).reset_index(drop=True)
    df_neg_final = pd.concat(df_neg_list).reset_index(drop=True)

    df_pos_final['rank'] = df_pos_final.groupby('chromosome').cumcount() + 1
    df_neg_final['rank'] = df_neg_final.groupby('chromosome').cumcount() + 1
    df_pos_final['strand'] = '+'
    df_neg_final['strand'] = '-'

    df_final = pd.concat([df_pos_final, df_neg_final])
    # df_final_s = df_final[['gene', 'strand', 'chromosome', 'rank', 'intergenic_distance']]
    return df_final

def sw_ratio(df, cols):
    if df.empty:
        return pd.DataFrame(columns=[*cols, 'ratio'])
    df_long = pd.melt(df, value_vars=cols, var_name='group', value_name='score')
    df_counts = df_long.groupby(['score', 'group']).size().reset_index(name='count')
    # len_dict={
    #     'random':df_long[df_long['group']=='random'].shape[0],
    #     'sw':df_long[df_long['group']=='sw'].shape[0]
    # }
    # df_counts['count'] = df_counts.apply(lambda x: len_dict[x['group']]-x['count'], axis=1)
    df_wide = df_counts.pivot(index='score', columns='group', values='count').fillna(0)
    df_wide = df_wide.reindex(columns=cols, fill_value=0)

    df_wide['ratio'] = df_wide[cols[1]] / df_wide[cols[0]]
    return df_wide

def cutoff(data, cf):
    if data.empty:
        return None
    df_wide_sw = sw_ratio(data, ['random', 'sw'])
    df_wide_sw.reset_index(inplace=True)
    df_wide_sw_s = df_wide_sw[df_wide_sw['score'] > 3] ### 有些值太低会有问题
    # sw_sum = df_wide_sw['sw'][df_wide_sw['ratio'] > 5].sum()
    candidates = df_wide_sw_s.loc[df_wide_sw_s['ratio'] > cf, 'score']
    if candidates.empty:
        return None
    return candidates.min()

def parse_sl_map(value):
    if value is None:
        return set()
    return {item.strip() for item in str(value).split(',') if item.strip()}

def standardize_sl_type(sl_type, sl1_refs, sl2_refs, legacy_mapping=False):
    sl_type = str(sl_type).strip()
    if sl_type in {'random', 'SL_unknown', 'SL1_unknown'}:
        return None
    if sl_type in sl1_refs:
        return 'SL1'
    if sl_type in sl2_refs:
        return 'SL2'
    if legacy_mapping:
        return 'SL2'
    return None

def sl_process(path, cf, sl1_refs, sl2_refs, legacy_mapping=False):
    output_columns = ['query_name', 'SL_type', 'SL']
    try:
        sl = pd.read_csv(path, sep='\t')
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=output_columns)

    required_columns = {'query_name', 'SL_type', 'SL_score', 'random_SL_score'}
    missing_columns = sorted(required_columns - set(sl.columns))
    if missing_columns:
        raise ValueError(
            'SL input is missing required column(s): ' + ', '.join(missing_columns)
        )

    sl = sl.copy()
    sl['SL_score'] = pd.to_numeric(sl['SL_score'], errors='coerce')
    sl['random_SL_score'] = pd.to_numeric(sl['random_SL_score'], errors='coerce')
    sl = sl.dropna(subset=list(required_columns))
    if sl.empty:
        return pd.DataFrame(columns=output_columns)

    sl['query_name'] = sl['query_name'].astype(str)
    sl['random'] = (sl['random_SL_score'] * 2).round() / 2
    sl['sw'] = (sl['SL_score'] * 2).round() / 2
    cutoff_value = cutoff(sl, cf)
    if cutoff_value is None or pd.isna(cutoff_value):
        return pd.DataFrame(columns=output_columns)

    # cutoff is calculated from the rounded "sw" score distribution, so the
    # filter must use the same score domain.  >= also retains the boundary bin.
    sl_s = sl[sl['sw'] >= cutoff_value]
    sl_s = sl_s.copy()
    sl_s['SL'] = sl_s['SL_type'].apply(
        lambda x: standardize_sl_type(
            x, sl1_refs, sl2_refs, legacy_mapping=legacy_mapping
        )
    )
    sl_s = sl_s.dropna(subset=['SL'])
    return sl_s[output_columns]

def fusion_expand(df, genes_dict):
    result_df = pd.DataFrame(columns=df.columns)
    df_with_or = df[df['gene'].str.contains(';', regex=True)]
    df_without_or = df[~df['gene'].str.contains(';', regex=True)]
    # 处理每一行
    for index, row in df_with_or.iterrows():
        gene = row['gene']
        split_genes = gene.split(";")
        temp_dict = {'gene': [], 'strand': [], 'chromosome': [],  'start': [], 'end': []}
        for gene in split_genes:
            if gene in genes_dict:
                temp_dict['gene'].append(gene)
                temp_dict['strand'].append(genes_dict[gene]['strand'])
                temp_dict['chromosome'].append(genes_dict[gene]['chromosome'])
                temp_dict['start'].append(genes_dict[gene]['start'])
                temp_dict['end'].append(genes_dict[gene]['end'])
        new_row = row.copy()
        if len(temp_dict['gene']) > 1:
            all_plus = all(info == '+' for info in temp_dict['strand'])
            all_minus = all(info == '-' for info in temp_dict['strand'])

            if all_plus:
                min_start = min(temp_dict['start'])
                min_start_index = temp_dict['start'].index(min_start)  # x[1] 是 start
                min_end = temp_dict['end'][min_start_index]
                max_start = min_start + (min_end - min_start)/2
                for i, start in enumerate(temp_dict['start']):
                    if min_start <= start <= max_start:  # 与最小 start 差值不超过 300bp
                        new_row['gene'] = temp_dict['gene'][i]
                        result_df = pd.concat([result_df, pd.DataFrame([new_row])], ignore_index=True)
            elif all_minus:
                max_end = max(temp_dict['end'])
                max_end_index = temp_dict['end'].index(max_end)  # x[1] 是 start
                max_start = temp_dict['start'][max_end_index]
                min_end = max_end - (max_end - max_start)/2
                for i, end in enumerate(temp_dict['end']):
                    if min_end <= end <= max_end:  # 与最小 start 差值不超过 300bp
                        new_row['gene'] = temp_dict['gene'][i]
                        result_df = pd.concat([result_df, pd.DataFrame([new_row])], ignore_index=True)
        elif len(temp_dict['gene']) == 1:
            new_row['gene'] = temp_dict['gene'][0]
            result_df = pd.concat([result_df, pd.DataFrame([new_row])], ignore_index=True)

    result_df = pd.concat([df_without_or, result_df], ignore_index=True)
    return result_df, df_with_or

def fusion_to_ref(df, genes_df):
    if df.empty:
        return []

    df = df.copy()
    df['group_id'] = [f'FS{str(i + 1).zfill(4)}' for i in range(len(df))]

    # 第二步：按"||"分割并展开
    df_expanded = df.assign(gene=df['gene'].str.split(';')).explode('gene')
    filtered_df = df_expanded[df_expanded['gene'].isin(genes_df['gene'])]
    value_counts = filtered_df['group_id'].value_counts()
    # 筛选出出现次数大于1的值
    duplicated_values = value_counts[value_counts > 1].index
    # 保留这些值的行
    filtered_df_s = filtered_df[filtered_df['group_id'].isin(duplicated_values)]
    if filtered_df_s.empty:
        return []
    filtered_df_s = pd.merge(filtered_df_s, genes_df, on='gene', how='left')
    # 第三步：为每个基因添加带序号的编码
    filtered_df_s = filtered_df_s.sort_values(['group_id', 'rank'])
    filtered_df_s['gene_id'] = filtered_df_s.groupby('group_id').cumcount() + 1

    filtered_df_s['gene_id'] = filtered_df_s['group_id'] + '_' + filtered_df_s['gene_id'].astype(str)
    filtered_df_ss = filtered_df_s[filtered_df_s['count']>3]
    result = [tuple(group['gene']) for _, group in filtered_df_ss.groupby('group_id')]
    # result = filtered_df_ss.groupby('group_id')['gene'].apply(list).to_dict()

    return result

def reshape(df_input):
    gene_data = {}
    for index, row in df_input.iterrows():
        gene = row['gene']  # 获取基因ID
        sl_type = row['SL']  # 获取SL类型
        count = int(row['count'])  # 获取数量

        # 如果基因还不在字典中，初始化
        if gene not in gene_data:
            gene_data[gene] = {'SL1': 0, 'SL2': 0}

        # 更新相应的SL计数
        gene_data[gene][sl_type] = count
    df = pd.DataFrame.from_dict(gene_data, orient='index').reset_index()
    df.columns = ['gene', 'SL1', 'SL2']
    return df

def operon_ref_process(path):
    operon = pd.read_csv(path, sep='\t', header=None)
    operon_s = operon[[8]]
    operon_s.columns = ['attributes']
    expanded_rows = []
    gene_add_row = []
    for row in operon_s['attributes']:
        parts = row.split(';')
        name = parts[0].split('=')[1]
        genes = parts[1].split('=')[1]
        gene_list = genes.split(',')
        gene_all = []
        for gene in gene_list:
            gene2 = 'Gene:' + gene
            expanded_rows.append({'operon': name, 'gene': gene2})
            gene_all.append(gene2)
        gene_add_row.append({'operon': name, 'gene': ",".join(gene_all)})
    gene_add_df = pd.DataFrame(gene_add_row)
    gene_add_df['gene'] = gene_add_df['gene'].apply(lambda x: ','.join(sorted(x.split(','))))

    return gene_add_df, pd.DataFrame(expanded_rows)

def extract_operon_names(df, count_fusion, median_value_sl2):
    operons = []
    i = 0
    n = df.shape[0]
    while i < n:
        if df.loc[i, 'type2'] == 'SL2':
            operon_genes = []
            sum_counts = 0
            if i > 0:
                if df.loc[i - 1, 'type'] == 'SL1' and df.loc[i - 1, 'SL1'] > 1:
                    operon_genes.append(df.loc[i - 1, 'gene'])
                elif (pd.isna(df.loc[i - 1, 'type'])
                      and df.loc[i, 'type'] == 'SL2'
                      and df.loc[i, 'SL2'] >= median_value_sl2):
                    operon_genes.append(df.loc[i - 1, 'gene'])
                elif (any((df.loc[i - 1, 'gene'] in tup and df.loc[i, 'gene'] in tup) for tup in count_fusion)):
                    operon_genes.append(df.loc[i - 1, 'gene'])
            operon_genes.append(df.loc[i, 'gene'])
            sum_counts += df.loc[i, 'sum_count']
            i += 1

            while i < n and df.loc[i, 'type2'] == 'SL2':
                operon_genes.append(df.loc[i, 'gene'])
                sum_counts += df.loc[i, 'sum_count']
                i += 1
            # 将得到的operon_genes加入operons列表
            if len(operon_genes) > 1 and sum_counts >= 3:
                operons.append(operon_genes)
            elif len(operon_genes) == 1 and df.loc[i-1, 'type'] == 'SL2' and sum_counts >= median_value_sl2:
                operons.append(operon_genes)
        else:
            i += 1

    return operons

def group_genes_into_operons(df, count_fusion, distance, median_value_sl2):
    operons = []
    current_operon = []

    for idx, row in df.iterrows():
        if pd.isna(row['intergenic_distance']) or row['intergenic_distance'] > distance:
            # 当遇到距离>5000的基因时，判断并保存之前的operon
            if current_operon:
                operons.append(current_operon)
                current_operon = []
        current_operon.append(row)
    # 保存最后一个operon
    if current_operon:
        operons.append(current_operon)

    # 对operon内的基因进行SL判断并重组
    operon_list = []
    for operon in operons:
        operon_df = pd.DataFrame(operon)

        sl2_genes = operon_df[operon_df['type2'] == 'SL2']['gene'].tolist()

        if sl2_genes and len(operon_df) > 1:
            operon_df = operon_df.reset_index(drop=True)
            # value = operon_df['rank'].iloc[0]
            # chr= operon_df['chromosome'].iloc[0]
            operons_t = extract_operon_names(operon_df, count_fusion, median_value_sl2)

            if operons_t:
                operon_list.extend(operons_t)

    dup_operon = list(set(tuple(sublist) for sublist in operon_list))
    return dup_operon


def merge_single_gene_sublists(gene_list, gene_df):
    modified_list = list(set(gene_list))
    # 找到所有只有一个基因的子列表
    single_gene_sublists = [sublist for sublist in gene_list if len(sublist) == 1]
    multiple_gene_sublists = [sublist for sublist in gene_list if len(sublist) > 1]
    gene_to_sublists = {}
    for sublist in multiple_gene_sublists:
        for gene in sublist:
            gene_to_sublists[gene] = sublist  # 每个基因都映射到它所在的子列表

    # 遍历单基因子列表
    for sublist in single_gene_sublists:
        gene = sublist[0]
        gene_row = gene_df[gene_df['gene'] == gene]
        chromosome, strand, rank = gene_row['chromosome'].tolist()[0], gene_row['strand'].tolist()[0], gene_row['rank'].tolist()[0]
        rank_minus_1 = rank - 1
        rank_plus_1 = rank + 1
        filtered_minus = gene_df[(gene_df['chromosome'] == chromosome) &
                                  (gene_df['strand'] == strand) &
                                  (gene_df['rank'] == rank_minus_1)]
        rank_minus_gene = filtered_minus['gene'].tolist()[0] if not filtered_minus.empty else None

        filtered_plus = gene_df[(gene_df['chromosome'] == chromosome) &
                                (gene_df['strand'] == strand) &
                                (gene_df['rank'] == rank_plus_1)]
        rank_plus_gene = filtered_plus['gene'].tolist()[0] if not filtered_plus.empty else None
        target_sublist = []
        if rank_minus_gene in gene_to_sublists:
            target_sublist = list(gene_to_sublists[rank_minus_gene])
            target_sublist.append(gene)
        elif rank_plus_gene in gene_to_sublists:
            target_sublist = list(gene_to_sublists[rank_plus_gene])
            target_sublist.insert(0, gene)

        if target_sublist:
            modified_list.append(target_sublist)
    # 去掉是另一个子列表子集的列表
    filtered_sublists = []
    for sublist in modified_list:
        if not any(set(sublist).issubset(set(other)) and sublist != other for other in modified_list):
            filtered_sublists.append(sublist)

    # 仅返回被修改过的子列表 + 没有被合并的单基因子列表
    return filtered_sublists

def generate_operon_gff(gene_combinations, gene_dict):
    """
    gene_combinations: list of lists, 每个子list包含一组gene ID
    gene_dict: dict, key是gene ID, value是包含start和end的dict/tuple
    chromosome: 染色体编号，默认"I"
    """
    operon_gff_dict = {'chromosome': [], 'source': [], 'type':[], 'start': [], 'end': [], 'score': [],
                       'strand': [], 'phase': [],'gene': []}
    gene_gff_dict = {'chromosome': [], 'source': [], 'type':[], 'start': [], 'end': [], 'score': [],
                       'strand': [], 'phase': [],'gene': []}
    # 处理每个gene组合
    for idx, genes in enumerate(gene_combinations):
        if not genes:  # 跳过空组合
            continue
        formatted_idx = f"LRS{idx + 1:04d}"   # idx + 1 因为要从 0001 开始，04d 确保 4 位补零
        # 获取所有gene的start和end
        starts = []
        ends = []
        chromosomes = []
        strands = []  # 检查方向是否一致

        for gene in genes:
            if gene in gene_dict:
                chromosome = gene_dict[gene]['chromosome']
                start = gene_dict[gene]['start']
                end = gene_dict[gene]['end']
                strand = gene_dict[gene]['strand']

                gene_gff_dict['chromosome'].append(chromosome)
                gene_gff_dict['strand'].append(strand)
                gene_gff_dict['start'].append(start)
                gene_gff_dict['end'].append(end)
                gene_entry = f"ID={gene};Parent={formatted_idx}"
                gene_gff_dict['gene'].append(gene_entry)

                starts.append(start)
                ends.append(end)
                chromosomes.append(chromosome)
                strands.append(strand)

        if not starts:  # 如果没有有效基因，跳过
            continue

        # 计算operon的范围
        chromosome = set(chromosomes).pop()
        operon_start = min(starts)
        operon_end = max(ends)
        strand = set(strands).pop()  # 默认正链，可根据需要修改逻辑
        genes_str = ",".join(genes)
        gene_entry = f"ID={formatted_idx};genes={genes_str}"

        operon_gff_dict['chromosome'].append(chromosome)
        operon_gff_dict['strand'].append(strand)
        operon_gff_dict['start'].append(operon_start)
        operon_gff_dict['end'].append(operon_end)
        operon_gff_dict['gene'].append(gene_entry)

    operon_gff_dict['source'] = ['SLRanger'] * len(operon_gff_dict['gene'])
    operon_gff_dict['type'] = ['operon'] * len(operon_gff_dict['gene'])
    operon_gff_dict['score'] = ['.'] * len(operon_gff_dict['gene'])
    operon_gff_dict['phase'] = ['.'] * len(operon_gff_dict['gene'])

    gene_gff_dict['source'] = ['SLRanger'] * len(gene_gff_dict['gene'])
    gene_gff_dict['type'] = ['gene'] * len(gene_gff_dict['gene'])
    gene_gff_dict['score'] = ['.'] * len(gene_gff_dict['gene'])
    gene_gff_dict['phase'] = ['.'] * len(gene_gff_dict['gene'])

    operon_gff_df = pd.DataFrame(operon_gff_dict)
    gene_gff_df = pd.DataFrame(gene_gff_dict)
    df = pd.concat([operon_gff_df, gene_gff_df])
    df_sorted = df.sort_values(by=['chromosome', 'start'])

    return df_sorted

def count_process(df, df_pos, genes_dict):
    empty_counts = pd.DataFrame(
        columns=['gene', 'SL1', 'SL2', 'sum_count', 'sl2_ratio', 'type', 'type2']
    )
    df = df.dropna(subset=['gene', 'SL'])
    if df.empty:
        return empty_counts, []

    counts = df.value_counts().reset_index()
    counts_expand, counts_fusion = fusion_expand(counts, genes_dict)
    fusion_ref = fusion_to_ref(counts_fusion, df_pos)
    if counts_expand.empty:
        return empty_counts, fusion_ref

    counts_s = counts_expand.groupby(['gene', 'SL'])['count'].sum().reset_index()
    # has_semicolon = counts_filtered[counts_filtered['sort_fusion'].str.contains(';', na=False)]
    counts_re = reshape(counts_s)
    counts_re['sum_count'] = counts_re['SL2'] + counts_re['SL1']
    counts_re['sl2_ratio'] = counts_re['SL2'] / counts_re['sum_count']

    counts_re['type'] = counts_re.apply(lambda row: 'SL2' if row['sl2_ratio'] > 0.5 else 'SL1', axis=1)
    counts_re['type2'] = counts_re['type']
    # 满足条件的行赋值为 'SL2'
    counts_re.loc[(counts_re['SL2'] > 5) & (counts_re['sl2_ratio'] >= 0.25), 'type2'] = 'SL2'

    return counts_re, fusion_ref

def expand_gene_associations(df):
    """Return one row per unique read/gene/SL association.

    A semicolon-delimited mapping such as ``geneA;geneB`` contributes once to
    each gene.  This table reports gene-read associations, so its total can be
    larger than the number of unique reads when a read maps to multiple genes.
    """
    columns = ['query_name', 'gene', 'SL']
    if df.empty:
        return pd.DataFrame(columns=columns)

    associations = df[columns].dropna(subset=columns).copy()
    associations['gene'] = associations['gene'].astype(str).str.split(';')
    associations = associations.explode('gene')
    associations['gene'] = associations['gene'].str.strip()
    associations = associations[associations['gene'] != '']
    return associations.drop_duplicates(columns).reset_index(drop=True)


def build_gene_sl_table(df, df_pos):
    associations = expand_gene_associations(df)
    if associations.empty:
        counts_re = pd.DataFrame(columns=['gene', 'SL1', 'SL2'])
    else:
        counts_s = (
            associations.groupby(['gene', 'SL'])['query_name']
            .nunique()
            .reset_index(name='count')
        )
        counts_re = reshape(counts_s)

    annotation_columns = ['gene', 'chromosome', 'strand', 'rank']
    annotation = df_pos[
        [column for column in annotation_columns if column in df_pos.columns]
    ].drop_duplicates('gene')
    counts_re = pd.merge(counts_re, annotation, how='left', on='gene')
    for col in ['SL1', 'SL2']:
        if col not in counts_re.columns:
            counts_re[col] = 0
        counts_re[col] = counts_re[col].fillna(0).astype(int)

    counts_re['sum_count'] = counts_re['SL1'] + counts_re['SL2']
    counts_re['sl2_ratio'] = counts_re.apply(
        lambda r: r['SL2'] / r['sum_count'] if r['sum_count'] > 0 else 0.0, axis=1
    )
    counts_re['type'] = counts_re.apply(
        lambda r: 'SL2' if r['sl2_ratio'] > 0.5 else 'SL1', axis=1
    )
    counts_re['type2'] = counts_re['type']
    counts_re.loc[
        (counts_re['SL2'] > 5) & (counts_re['sl2_ratio'] >= 0.25), 'type2'
    ] = 'SL2'
    counts_re = counts_re[counts_re['sum_count'] > 0].copy()

    gene_sl_cols = [
        'gene', 'chromosome', 'strand', 'rank',
        'SL1', 'SL2', 'sum_count', 'sl2_ratio', 'type', 'type2',
    ]
    return counts_re[[c for c in gene_sl_cols if c in counts_re.columns]]


def resolve_mapping(args, gff_file):
    mapping_path = getattr(args, 'mapping', None)
    if mapping_path:
        return mapping_path

    bam_path = getattr(args, 'bam', None)
    if not bam_path:
        raise ValueError('Either a BAM file or a read-to-gene mapping file is required.')
    if run_track_cluster is None:
        raise RuntimeError(
            'BAM input requires the installed SLRanger package. '
            'Alternatively, provide a mapping file with -m/--mapping.'
        )
    return run_track_cluster(gff_file, bam_path)


def read_mapping(path):
    try:
        mapping = pd.read_csv(
            path,
            sep='\t',
            header=None,
            usecols=[0, 1],
            names=['query_name', 'gene'],
        )
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=['query_name', 'gene'])

    mapping = mapping.dropna(subset=['query_name', 'gene'])
    mapping = mapping[
        ~(
            (mapping['query_name'].astype(str) == 'query_name')
            & (mapping['gene'].astype(str) == 'gene')
        )
    ]
    mapping['query_name'] = mapping['query_name'].astype(str)
    mapping['gene'] = mapping['gene'].astype(str)
    # Keep the original mapping rows for backward-compatible operon counts.
    # Per-gene reporting performs its own read/gene/SL deduplication.
    return mapping.reset_index(drop=True)


def report_skipped_operon_output(output):
    output_path = Path(output)
    if output_path.exists():
        print(
            'Warning: operon prediction was skipped, so the existing output file '
            + str(output_path)
            + ' was not updated.'
        )
    else:
        print('Operon output was not generated: ' + str(output_path))


def main(args):
    gff_file = getattr(args, 'gff', None) or getattr(args, 'refer', None)
    if not gff_file:
        raise ValueError('A GFF annotation file is required.')

    df_genes = parse_gff(gff_file)
    df_genes_filter = parse_cds_gene(gff_file)
    if df_genes_filter.empty:
        df_genes_with_cds = df_genes.copy()
    else:
        df_genes_with_cds = pd.merge(
            df_genes, df_genes_filter, on='gene', how='right'
        )
    df_pos = sort_and_calc_distance(df_genes_with_cds)
    df_pos_dict = (
        df_genes_with_cds.drop_duplicates('gene').set_index('gene').to_dict('index')
    )
    mapping_path = resolve_mapping(args, gff_file)
    map_gene = read_mapping(mapping_path)

    sl1_map_value = getattr(args, 'sl1_map', getattr(args, 'SL1_map', None))
    sl2_map_value = getattr(args, 'sl2_map', getattr(args, 'SL2_map', None))
    legacy_mapping = sl1_map_value is None and sl2_map_value is None
    sl1_refs = {'SL1'} if legacy_mapping else parse_sl_map(sl1_map_value)
    sl2_refs = set() if legacy_mapping else parse_sl_map(sl2_map_value)
    sl_ss = sl_process(
        args.input,
        args.cutoff,
        sl1_refs,
        sl2_refs,
        legacy_mapping=legacy_mapping,
    )
    sl_ss_gene = pd.merge(sl_ss, map_gene, how='left', on='query_name')
    sl_ss_gene_for_count = sl_ss_gene.dropna(subset=['gene', 'SL'])

    gene_sl_table = getattr(args, 'gene_sl_table', None)
    if gene_sl_table is None:
        output_path = Path(args.output)
        gene_sl_table = str(
            output_path.with_name(output_path.stem + '_gene_sl1_sl2.tsv')
        )
        args.gene_sl_table = gene_sl_table

    gene_sl_df = build_gene_sl_table(
        sl_ss_gene_for_count[['query_name', 'gene', 'SL']], df_pos
    )
    gene_sl_path = Path(gene_sl_table)
    gene_sl_path.parent.mkdir(parents=True, exist_ok=True)
    gene_sl_df.to_csv(gene_sl_path, sep='\t', index=False)
    print('Per-gene SL1/SL2 counts written to ' + str(gene_sl_path))

    detected_types = set(sl_ss['SL'])
    if not detected_types:
        print(
            'No high-confidence SL1 or SL2 reads were detected. '
            'The per-gene SL table was still generated; operon prediction was skipped.'
        )
        report_skipped_operon_output(args.output)
        return 0
    if detected_types == {'SL1'}:
        print(
            'Only SL1 reads were detected. Operon prediction requires SL2 reads '
            '(normally together with SL1), so only the per-gene SL table was generated.'
        )
        report_skipped_operon_output(args.output)
        return 0
    if sl_ss_gene_for_count.empty or 'SL2' not in set(sl_ss_gene_for_count['SL']):
        print(
            'SL2 reads were detected, but none could be assigned to a gene. '
            'The per-gene SL table was generated; operon prediction was skipped.'
        )
        report_skipped_operon_output(args.output)
        return 0

    counts_re, count_fusion = count_process(sl_ss_gene_for_count[['gene', 'SL']], df_pos, df_pos_dict)
    #median_value_all = counts_re['sum_count'].median()
    median_value_sl2 = counts_re[counts_re['type2'] == 'SL2']['SL2'].median()
    counts_re = pd.merge(counts_re, df_pos, how='right', on='gene')
    counts_re['sum_count'] = counts_re['sum_count'].fillna(0)

    operon_result = group_genes_into_operons(counts_re, count_fusion, args.distance, median_value_sl2)
    updated_gene_list = merge_single_gene_sublists(operon_result, df_pos)
    # operon_combination = pd.DataFrame([','.join(sublist) for sublist in updated_gene_list])
    operon_combination_gff = generate_operon_gff(updated_gene_list, df_pos_dict)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    operon_combination_gff.to_csv(output_path, sep='\t', index=False, header=False)

    print('Operon detected to ' + str(output_path))
    return 0


def build_parser():
    parser = argparse.ArgumentParser(
        description="help to know spliced leader and distinguish SL1 and SL2")
    parser.add_argument(
        "-g", "--gff", "-r", "--refer", dest="gff", required=True,
        help="GFF annotation file",
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("-b", "--bam", type=str, help="BAM file")
    input_group.add_argument(
        "-m", "--mapping", type=str, help="existing read-to-gene mapping file"
    )
    parser.add_argument("-i", "--input", type=str, required=True, help="input the SL detection file")
    parser.add_argument("-o", "--output", type=str, default="SLRanger.gff", help="output operon detection file")
    parser.add_argument("--gene-sl-table", type=str, default=None, help="per-gene SL1/SL2 count table")
    parser.add_argument(
        "--sl1-map", "--SL1_map", "-SL1_map", dest="sl1_map", default=None,
        help="comma-separated SL_type values treated as SL1",
    )
    parser.add_argument(
        "--sl2-map", "--SL2_map", "-SL2_map", dest="sl2_map", default=None,
        help=(
            "comma-separated SL_type values treated as SL2; when neither map "
            "option is supplied, legacy SLRanger classification is used"
        ),
    )
    parser.add_argument("-d", "--distance", type=int, default=5000, help="promoter scope")
    parser.add_argument("-c", "--cutoff", type=float, default=4, help="cutoff of high confident SL sequence")
    return parser


if __name__ == '__main__':
    parser = build_parser()
    raise SystemExit(main(parser.parse_args()))
