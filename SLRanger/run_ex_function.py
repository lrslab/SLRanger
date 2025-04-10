import os
import time
from collections import OrderedDict
import subprocess
import numpy as np
import pandas as pd
from tqdm import tqdm
import multiprocessing
import warnings
from datetime import datetime
import shutil


def identify_file_path(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError("File do not exist! Please check your path : " + file_path)

def run_cmd(cmd):
    try:
        # 执行命令并捕获输出
        print(cmd)
        subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT, universal_newlines=True)
        # 处理命令输出
        # print(output)
    except subprocess.CalledProcessError as e:
        # 处理命令执行错误
        print("Command execution failed:", cmd)
        print("Reason:", e.output)
        raise RuntimeError('There are some errors in the cmd as below, please check your env ')

def run_track_cluster(gff_file,bam_file):
    identify_file_path(gff_file)
    identify_file_path(bam_file)

    # 获取当前时间戳（格式：YYYYMMDD_HHMMSS）
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # 创建文件夹名称（SLRange_view_中间加时间戳）
    folder_name = f"SLRange_temp_{timestamp}/"
    # 创建文件夹
    os.makedirs(folder_name, exist_ok=True)
    try:
        gff2bigg_cmd = 'gff2bigg.py -i '+gff_file+' -o '+folder_name+'ref.bed'
        run_cmd(gff2bigg_cmd)
        bam2bigg_cmd = 'bam2bigg.py -b '+bam_file+' -o '+folder_name+'read.bed'
        run_cmd(bam2bigg_cmd)
        bedtools_cmd = 'bedtools sort -i '+folder_name+'read.bed'+' >'+folder_name+'read_sort.bed'
        run_cmd(bedtools_cmd)
        track_cmd = "add_gene.py -r "+folder_name+'ref.bed'+" -s "+folder_name+"read_sort.bed"
        run_cmd(track_cmd)
    except Exception as e:
        print(e)
        # shutil.rmtree(folder_name)
        raise RuntimeError('There are some errors in the cmd as below, please check your env ')
    # shutil.rmtree(folder_name)
    return "read_sort_gene.bed"
