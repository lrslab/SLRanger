#!/usr/bin/env python
import re
import argparse
import pandas as pd
from SLRanger.run_ex_function import run_track_cluster
import warnings

warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)  # 如果使用plotnine可能需要
# 解析GFF文件并构建DataFrame
def parse_gff(gff_file):
    genes = []
    with open(gff_file, 'r') as file:
        for line in file:
            if line.startswith("#") or line.strip() == '':
                continue
            parts = line.strip().split('\t')
            if parts[2] == 'gene':
                attr_field = parts[8]
                gene_id = None
                for attr in attr_field.split(';'):
                    if attr.startswith('ID='):
                        gene_id = attr.split('ID=')[1].split(';')[0]
                genes.append({
                    'gene': gene_id,
                    'chromosome': parts[0],
                    'start': int(parts[3]),
                    'end': int(parts[4]),
                    'strand': parts[6]
                })
    return pd.DataFrame(genes)

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
    df_final_s = df_final[['gene', 'strand', 'chromosome', 'rank', 'intergenic_distance']]
    return df_final_s

def sw_ratio(df, cols):
    df_long = pd.melt(df, value_vars=cols, var_name='group', value_name='score')
    df_counts = df_long.groupby(['score', 'group']).size().reset_index(name='count')
    # len_dict={
    #     'random':df_long[df_long['group']=='random'].shape[0],
    #     'sw':df_long[df_long['group']=='sw'].shape[0]
    # }
    # df_counts['count'] = df_counts.apply(lambda x: len_dict[x['group']]-x['count'], axis=1)
    df_wide = df_counts.pivot(index='score', columns='group', values='count').fillna(0)

    df_wide['ratio'] = df_wide[cols[1]] / df_wide[cols[0]]
    return df_wide

def cutoff(data, cf):
    df_wide_sw = sw_ratio(data, ['random', 'sw'])
    df_wide_sw.reset_index(inplace=True)
    # sw_sum = df_wide_sw['sw'][df_wide_sw['ratio'] > cf].sum()
    sw_min = df_wide_sw['score'][df_wide_sw['ratio'] > cf].min()
    return sw_min

def sl_process(path, cf):
    sl = pd.read_csv(path, sep='\t')
    sl = sl.dropna()
    sl['SL_score'] = sl['SL_score'].astype(float)
    sl['random'] = (sl['random_SL_score'] * 2).round() / 2
    sl['sw'] = (sl['SL_score'] * 2).round() / 2
    sl_s = sl[sl['SL_score'] > cutoff(sl, cf)]
    sl_s = sl_s[(sl_s['SL_type'] != 'random') & (sl_s['SL_type'] != 'SL1_unknown')]
    sl_s['SL'] = sl_s['SL_type'].apply(lambda x: 'SL1' if x == 'SL1' else 'SL2')
    return sl_s[['query_name', 'SL']]

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
                min_start_gene = temp_dict['gene'][min_start_index]
                min_end = temp_dict['end'][min_start_index]
                max_start = min_start + (min_end - min_start)/2
                new_row['gene'] = min_start_gene
                result_df = pd.concat([result_df, pd.DataFrame([new_row])], ignore_index=True)
                for i, start in enumerate(temp_dict['start']):
                    if min_start <= start <= max_start:  # 与最小 start 差值不超过 300bp
                        new_row['gene'] = temp_dict['gene'][i]
                        result_df = pd.concat([result_df, pd.DataFrame([new_row])], ignore_index=True)
            elif all_minus:
                max_end = max(temp_dict['end'])
                max_end_index = temp_dict['end'].index(max_end)  # x[1] 是 start
                max_end_gene = temp_dict['gene'][max_end_index]
                max_start = temp_dict['start'][max_end_index]
                min_end = max_end - (max_end - max_start)/2
                new_row['gene'] = max_end_gene
                result_df = pd.concat([result_df, pd.DataFrame([new_row])], ignore_index=True)
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
    df['group_id'] = [f'FS{str(i + 1).zfill(4)}' for i in range(len(df))]

    # 第二步：按"||"分割并展开
    df_expanded = df.assign(gene=df['gene'].str.split(';')).explode('gene')
    filtered_df = df_expanded[df_expanded['gene'].isin(genes_df['gene'])]
    value_counts = filtered_df['group_id'].value_counts()
    # 筛选出出现次数大于1的值
    duplicated_values = value_counts[value_counts > 1].index
    # 保留这些值的行
    filtered_df_s = filtered_df[filtered_df['group_id'].isin(duplicated_values)]
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
                elif (pd.isna(df.loc[i - 1, 'type']) and df.loc[i, 'type'] == 'SL2'
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
    counts = df.value_counts().reset_index()
    counts_expand, counts_fusion = fusion_expand(counts, genes_dict)
    fusion_ref = fusion_to_ref(counts_fusion, df_pos)
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

def main(args):
    print("Detection start ...")
    gff_file = args.gff
    df_genes = parse_gff(gff_file)
    df_genes_filter = parse_cds_gene(gff_file)
    df_genes_with_cds = pd.merge(df_genes, df_genes_filter, on='gene', how='right')
    df_pos = sort_and_calc_distance(df_genes_with_cds)
    df_pos_dict = df_genes_with_cds.set_index('gene').to_dict('index')
    args.mapping=run_track_cluster(args.gff,args.bam)
    map_gene = pd.read_csv(args.mapping, sep='\t', names=['query_name', 'gene'])
    sl_ss = sl_process(args.input, args.cutoff)
    sl_ss_gene = pd.merge(sl_ss, map_gene, how='left', on='query_name')
    counts_re, count_fusion = count_process(sl_ss_gene[['gene', 'SL']], df_pos, df_pos_dict)
    # median_value_all = counts_re['sum_count'].median()
    median_value_sl2 = counts_re[counts_re['type2'] == 'SL2']['SL2'].median()
    counts_re = pd.merge(counts_re, df_pos, how='right', on='gene')
    counts_re['sum_count'] = counts_re['sum_count'].fillna(0)

    operon_result = group_genes_into_operons(counts_re, count_fusion, args.distance, median_value_sl2)
    updated_gene_list = merge_single_gene_sublists(operon_result, df_pos)

    # operon_combination = pd.DataFrame([','.join(sublist) for sublist in updated_gene_list])
    operon_combination_gff = generate_operon_gff(updated_gene_list, df_pos_dict)
    operon_combination_gff.to_csv(args.output, sep='\t', index=False, header=False)

    print('Operon detected to ' + args.output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="help to know spliced leader and distinguish SL1 and SL2")
    parser.add_argument("-g", "--gff", type=str, required=True, help="GFF annotation file")
    parser.add_argument("-b", "--bam", type=str, required=True, help="bam file")
    parser.add_argument("-i", "--input", type=str, required=True, help="input the SL detection file")
    parser.add_argument("-o", "--output", type=str,  default="SLRanger.gff",help="output operon detection file")
    parser.add_argument("-d", "--distance", type=int, default=5000, help="promoter scope")
    parser.add_argument("-c", "--cutoff", type=float, default=4, help="cutoff of high confident SL sequence")
    args = parser.parse_args()
    main(args)
