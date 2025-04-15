#!/usr/bin/env python
from trackcluster.utils import myexe, is_bin_in, get_file_prefix,del_files,get_file_location
from trackcluster.tracklist import read_bigg, write_bigg, list_to_dic
from trackcluster.pre import tracklist_add_gene,get_gendic_bedinter,group_bigg_by_gene
import os
import argparse

def wrapper_bedtools_intersect2_select(bedfile1,bedfile2,outfile,fraction_bed1, fraction_bed2):
    """
    Using two bedfile to get the intsersection of pairs
    use simular fraction as substract instead of 0.2, the bigg can contain large intron in first intersect

    # -f 0.01= 1% from bigg, -F 0.05, means 5% from gene annotation

    :param bigg_one: the read
    :param bigg_two: the ref
    :param: fraction, may need to inlcude some read track with very long
    :return:
    """
    # generate the bedfile, -r is reciprocal overlap, -s is strandedness/keep the same strand

    # sometime the name check for an intact bed and an non-intact one can result in a blank file
    # add -nonamecheck
    cmd="bedtools intersect -nonamecheck -wa -wb  -s -f {f1} -F {f2} -a {bedfile1} -b {bedfile2}>{out}".format(
        bedfile1=bedfile1, bedfile2=bedfile2, out=outfile, f1=fraction_bed1, f2=fraction_bed2)

    _=myexe(cmd)

    return outfile

def flow_add_gene(wkdir, prefix, bigg_gff_file, bigg_nano_file):
    os.chdir(wkdir)
    # make sure one read one track
    bigg_raw=read_bigg(bigg_nano_file)
    bigg_dedup=list(list_to_dic(bigg_raw).values())
    print("raw bigg number: {}; after dedup:{}".format(len(bigg_raw), len(bigg_dedup)))
    outbed=prefix+"_dedup.bed"
    write_bigg(bigg_dedup, outbed)

    ### get two parts, the gene part and the novel part,

    # the gene part
    out_single = prefix + "_single_inter.bed"
    # write the outfile to disk
    wrapper_bedtools_intersect2_select(outbed, bigg_gff_file, outfile=out_single,
                                       fraction_bed1=0.01, fraction_bed2=0.05)
    out_fusion = prefix + "_fusion_inter.bed"
    # write the outfile to disk
    wrapper_bedtools_intersect2_select(outbed, bigg_gff_file, outfile=out_fusion,
                                       fraction_bed1=0.33, fraction_bed2=0.33)

    read_gene_single = get_gendic_bedinter(out_single)
    read_gene_fusion = get_gendic_bedinter(out_single)
    print("read number in genes:", len(read_gene_single))

    bigg_single = tracklist_add_gene(bigg_dedup, read_gene_single)
    bigg_fusion = tracklist_add_gene(bigg_dedup, read_gene_fusion)

    # cleanup
    del_files([outbed, out_single, out_fusion])

    return bigg_single, bigg_fusion

def flow_fusion_annotation(wkdir,isoform_bed, prefix=None):
    if prefix is None:
        prefix=get_file_prefix(isoform_bed, sep=".")

    os.chdir(wkdir)
    bigg_isoform=read_bigg(isoform_bed)

    # make the fusion class separately
    fusion_d = {}
    for bigg in bigg_isoform:
        genename_l=bigg.geneName.split("||")
        if len(genename_l)>1:
            fusion_d[bigg.name]=genename_l

    # IO for fusion
    with open(prefix+"_fusion.txt", "w") as fw:
        for k, v in fusion_d.items():
            v_str = ";".join(v)
            fw.write(k + "\t" + v_str + "\n")

    return prefix+"_fusion.txt"

def flow_single_annotation(wkdir,isoform_bed, prefix=None):
    if prefix is None:
        prefix=get_file_prefix(isoform_bed, sep=".")

    os.chdir(wkdir)
    bigg_isoform=read_bigg(isoform_bed)

    # make the fusion class separately
    with open(prefix + "_single.txt", "w") as fw:
        for bigg in bigg_isoform:
            if bigg.geneName != 'none':
                genename_l=bigg.geneName.split("||")
                if len(genename_l) == 1:
                    fw.write(bigg.name + "\t" + genename_l[0] + "\n")

    return prefix + "_single.txt"

def concatenate_and_remove(file1, file2, output_file):
    # 打开输出文件，以写入模式
    with open(output_file, 'w') as outfile:
        # 读取并写入第一个文件的内容
        with open(file1, 'r') as f1:
            outfile.write(f1.read())
        # 读取并写入第二个文件的内容
        with open(file2, 'r') as f2:
            outfile.write(f2.read())

if __name__ == '__main__':
# def addgene(self):
    parser = argparse.ArgumentParser(
        description="Used to add gene annotation for read bigg tracks, useful in some analysis, "
                    "the process is included in cluster runs. the new bigg file will be prefix_gene.bed"
    )
    parser.add_argument("-d", "--folder", default=os.getcwd(),
                        help="the folder contains all the seperated tracks in different locus/genes, default is the current dir")
    parser.add_argument("-s", "--sample",
                        help="the bigg format of the read track, with the key of GeneName")
    parser.add_argument("-r", "--reference",
                        help="the bigg format of the reference annotation track")

    args = parser.parse_args()
    args.prefix = get_file_prefix(args.sample, sep=".")

    os.chdir(args.folder)
    bigg_single, bigg_fusion = flow_add_gene(wkdir=args.folder,
                                             prefix=args.prefix,
                                             bigg_gff_file=args.reference,
                                             bigg_nano_file=args.sample,
                                             )
    out_single = args.prefix + "_single_gene.bed"
    out_fusion = args.prefix + "_fusion_gene.bed"
    write_bigg(bigg_single, out_single)
    write_bigg(bigg_fusion, out_fusion)

    fusion_f = flow_fusion_annotation(wkdir=args.folder,
                                      isoform_bed=out_fusion,
                                      prefix=args.prefix)

    single_f = flow_single_annotation(wkdir=args.folder,
                                      isoform_bed=out_single,
                                      prefix=args.prefix)
    output_file = args.prefix + "_gene.bed"
    concatenate_and_remove(single_f, fusion_f, output_file)
    del_files([out_single, out_fusion, single_f, fusion_f])
