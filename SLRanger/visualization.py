import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from plotnine import *
from plotnine.ggplot import PlotnineWarning
import re
import markdown
from pathlib import Path
from datetime import datetime
import seaborn as sns
import matplotlib.colors as mcolors
import warnings
# 隐藏特定的警告
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)  # 如果使用plotnine可能需要
warnings.filterwarnings('ignore', category=PlotnineWarning)


# Function to calculate ratio
def sw_ratio(df, cols):
    df_long = pd.melt(df, value_vars=cols, var_name='group', value_name='score')
    df_counts = df_long.groupby(['score', 'group']).size().reset_index(name='count')
    # len_dict={
    #     'random_seq':df_long[df_long['group']=='random_seq'].shape[0],
    #     'SL_reference':df_long[df_long['group']=='SL_reference'].shape[0]
    # }
    # df_counts['count'] = df_counts.apply(lambda x: len_dict[x['group']]-x['count'], axis=1)
    df_wide = df_counts.pivot(index='score', columns='group', values='count').fillna(0)

    df_wide['ratio'] = df_wide[cols[1]] / df_wide[cols[0]]
    return df_wide

def plot_cumulative_line(data, output_name, cf):
    df_wide_sw = sw_ratio(data, ['random_seq', 'SL_reference'])
    df_wide_sw.reset_index(inplace=True)
    # sw_sum = df_wide_sw['SL_reference'][df_wide_sw['ratio'] > 5].sum()
    sw_min = df_wide_sw['score'][df_wide_sw['ratio'] > cf].min()

    # Set ratios > 50 to infinity
    # df_wide_sw.loc[df_wide_sw['ratio'] > 50, 'ratio'] = np.inf
    # df_wide_sw['SL_reference'] = df_wide_sw['SL_reference'].sum() - df_wide_sw['SL_reference']
    # df_wide_sw['random_seq'] = df_wide_sw['random_seq'].sum() - df_wide_sw['random_seq']
    # df_wide_sw['cumulative_sw'] = df_wide_sw['SL_reference'].cumsum()
    # df_wide_sw['cumulative_random'] = df_wide_sw['random_seq'].cumsum()
    df_wide_sw['SL_reference'] = df_wide_sw['SL_reference'][::-1].cumsum()[::-1]
    df_wide_sw['random_seq'] = df_wide_sw['random_seq'][::-1].cumsum()[::-1]
    plot_df = pd.melt(df_wide_sw,id_vars=['score'],value_vars=['SL_reference', 'random_seq'])

    # Create cumulative line plot
    tick_step = 5
    x_max = int(np.floor(plot_df['score'].max()))
    # 生成整数刻度（不超过 x_max）
    base_breaks = list(range(0, x_max + 1, tick_step))
    x_breaks = sorted(set(base_breaks + [sw_min]))
    p_cumulative = (
            ggplot(plot_df, aes(x='score', y='value', color='group'))
            + geom_line()
            + geom_point(size=1)
            + geom_vline(xintercept=sw_min, linetype='dashed', color='black')
            + scale_x_continuous(breaks=x_breaks)
            + labs(title='Cumulative Counts', x='Score', y='Cumulative Count')
            + theme_bw()
            + theme(
        plot_title=element_text(ha='center'),
        legend_position='right',
        axis_text_y=element_text(margin={'r': 10}),
    )
    )
    # Save the plot
    p_cumulative.save(output_name, width=5, height=3.5,format='png',dpi=300)

    return  sw_min

def plot_aligned_length(df,folder_name):

    def extract_number(s,unique_number):
        match = re.search(r'SL(\d+)', s)  # 使用正则表达式提取数字
        if s=='SL2_unknown':
            return unique_number - 1
        elif s =='SL1_unknown':
            return unique_number - 2
        if match:
            return int(match.group(1))  # 返回数字部分并转换为整数
        return 0  # 如果没有匹配到，返回 0（可以根据需要调整）
    # unknown will be merged together and SL1 separated
    # df['SL_type'].apply(lambda x: if 'unknown'(x))

    result_list=['SL1','SL2','SL3','SL4','SL5','SL6','SL7','SL8','SL9','SL10','SL11','SL12','SL13','SL1_unknown','SL2_unknown']
    colors_13 = [
        "#1F77B4",  # SL1 - 蓝色
        "#FF7F0E",  # SL2 - 橙色
        "#2CA02C",  # SL3 - 绿色
        "#D62728",  # SL4 - 红色
        "#9467BD",  # SL5 - 紫色
        "#8C564B",  # SL6 - 棕色
        "#E377C2",  # SL7 - 粉红
        "#00008B",  # SL8 - 深紫
        "#BCBD22",  # SL9 - 橄榄黄
        "#17BECF",  # SL10 - 青色
        "#AEC7E8",  # SL11 - 浅蓝
        "#FFBB78",  # SL12 - 浅橙
        "#98DF8A"  # SL13 - 浅绿
    ]

    # 创建完整的颜色映射字典
    color_mapping = {
        f"SL{i + 1}": colors_13[i] for i in range(13)
    }
    color_mapping['SL1_unknown']='#231815' # black
    color_mapping['SL2_unknown'] = 'darkgrey'
    df['SL_type'] = df['SL_type'].apply(lambda x: 'SL2_unknown' if 'unknown' in x and x != 'SL1_unknown' else x)
    data = df['SL_type'].value_counts().reset_index()
    category = pd.api.types.CategoricalDtype(categories=result_list, ordered=True)
    data['SL_type'] = data['SL_type'].astype(category)
    data = data.sort_values('SL_type')
    def autopct_func(pct, allvals):
        absolute = round(pct / 100. * sum(allvals))
        return f"{pct:.1f}%\n({absolute:d})"

    labels = data['SL_type'].values
    sizes = data['count'].values
    # 生成15个清新风格的颜色
    # 使用Seaborn的"husl"色系生成14个颜色，然后添加指定的灰色
    base_colors = sns.color_palette("pastel", len(labels)-1)  # 生成14个清新颜色
    # 将RGB颜色转换为HEX
    base_hex_colors = [mcolors.to_hex(color) for color in base_colors]
    # 添加最后一个指定的灰色

    plt.figure(figsize=(5, 5))
    plt.pie(sizes, labels=labels,
            colors=[color_mapping[label] for label in labels],
            pctdistance=0.65,  # 百分比标签距离中心的距离
            autopct=lambda pct: autopct_func(pct, sizes),
            startangle=90,
            textprops={'fontsize': 8},  # 调整字体大小
            labeldistance=1.1,  # 标签距离中心的距离
            rotatelabels=True)  # 旋转标签以避免重叠
    # 确保饼图为圆形
    plt.axis('equal')
    # 保存图形
    plt.savefig(folder_name+'type_pie.png', dpi=300, bbox_inches='tight')

    plot = ggplot(df, aes(x='query_length')) \
           + theme_bw()\
           + geom_histogram(binwidth=0.5) \
           + theme(
        figure_size=(8, 4),
        axis_text=element_text(size=12,),
        axis_title=element_text(size=12,),
        panel_grid_minor=element_blank(),
        title=element_text(size=12,),
        strip_background=element_rect(alpha=0),
        strip_text=element_text(size=12,),
        legend_position='bottom',
        legend_text=element_text(size=8,),
    )
    plot.save(folder_name+'query_length.png',dpi=300,format='png')

    Q1 = np.percentile(df['aligned_length'], 25)
    Q2 = np.percentile(df['aligned_length'], 50)
    Q3 = np.percentile(df['aligned_length'], 75)
    Q_list = [Q1, Q2, Q3]
    Q_axis_list = [[0, 0],
                   [1000, 1000],
                   [2000, 2000],
                   [3000, 3000],
                   ['Q1', Q_list[0]],
                   ['Q2', Q_list[1]],
                   ['Q3', Q_list[2]]]
    Q_axis_list = pd.DataFrame(Q_axis_list)
    df = df[df['aligned_length'] < 3500]
    Q_axis_list.sort_values(by=[1], inplace=True)
    plot = ggplot(df, aes(x='aligned_length')) \
           + theme_bw() \
           + geom_density(fill='#D7D7D7') \
           + scale_x_continuous(breaks=Q_axis_list[1].values.tolist(), labels=Q_axis_list[0].values.tolist()) \
           + theme(
        figure_size=(8, 4),
        axis_text=element_text(size=12, ),
        axis_title=element_text(size=12, ),
        panel_grid_minor=element_blank(),
        title=element_text(size=12, ),
        strip_background=element_rect(alpha=0),
        strip_text=element_text(size=12, ),
        legend_position='bottom',
        legend_text=element_text(size=8, ),
    )
    for item in Q_list:
        plot = plot + geom_vline(xintercept=item, color='black', alpha=0.5, linetype='dashdot')
    plot.save(folder_name + 'aligned_length.png', dpi=300, format='png')

    return data

def create_image_gallery_md_html(df, image_paths, output_md_path, output_html_path):
    # Markdown模板 - 添加表格部分
    md_template = """
# Data Visualization Results

This document contains the visualization results from the SLRanger.

## Data Summary Table

{data_table}

## Visualization Gallery

{image_sections}

---

*Generated on {current_date}. If you think this tool is pretty GOOD, don’t forget to give our [Git](https://github.com/lrslab/SLRanger) a star! *
    """

    # HTML样式模板 - 添加表格样式
    html_template = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Data Visualization Results</title>
    <style>
        body {{
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            line-height: 1.6;
        }}
        h1 {{
            color: #333;
            text-align: center;
            border-bottom: 2px solid #eee;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #444;
            margin-top: 30px;
        }}
        img {{
            max-width: 100%;
            height: auto;
            display: block;
            margin: 20px auto;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            border-radius: 4px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}
        th, td {{
            padding: 8px;
            text-align: left;
            border: 1px solid #ddd;
        }}
        th {{
            background-color: #f5f5f5;
            color: #333;
        }}
        tr:nth-child(even) {{
            background-color: #fafafa;
        }}
        .footer {{
            text-align: center;
            color: #666;
            font-size: 0.9em;
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid #eee;
        }}
    </style>
</head>
<body>
    {html_content}
</body>
</html>
"""

    # 将DataFrame转换为Markdown表格
    data_table = df.to_markdown(index=False)

    # 创建每个图片的Markdown片段
    image_sections = ""
    image_captions = [
        "Cumulative Counts (SW). The **dashed black line** represents the recommended **cutoff** value.",
        "Cumulative Counts (SL). The **dashed black line** represents the recommended **cutoff** value.",
        "Query Length Distribution",
        "Aligned Length Distribution",
        "SL Type Distribution"
    ]

    for img_path, caption in zip(image_paths, image_captions):
        image_sections += f"""
### {caption}

![{caption}]({img_path})

---
        """

    # 获取当前日期
    current_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # 生成完整的Markdown内容
    md_content = md_template.format(
        data_table=data_table,
        image_sections=image_sections,
        current_date=current_date
    )

    # 写入Markdown文件
    with open(output_md_path, 'w', encoding='utf-8') as f:
        f.write(md_content)

    # 将Markdown转换为HTML
    html_content = markdown.markdown(
        md_content,
        extensions=['markdown.extensions.extra', 'markdown.extensions.toc', 'markdown.extensions.tables']
    )

    # 将HTML内容嵌入到模板中
    final_html = html_template.format(html_content=html_content)

    # 写入HTML文件
    with open(output_html_path, 'w', encoding='utf-8') as f:
        f.write(final_html)

# 主程序

def visualize_html(output_file, cf):
    import os
    from datetime import datetime

    # 获取当前时间戳（格式：YYYYMMDD_HHMMSS）
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")


    path = output_file
    out_put_name = path.split('.')[0]

    # 创建文件夹名称（SLRange_view_中间加时间戳）
    folder_name = f"{out_put_name}_{timestamp}/"
    # 创建文件夹
    os.makedirs(folder_name, exist_ok=True)
    # Read data
    df = pd.read_csv(path, sep='\t')
    reads_all = len(df)
    df = df.dropna()
    reads_na = len(df)
    # print(reads_na / reads_all)

    # Data processing
    df['SL_score'] = df['SL_score'].astype(float)

    # Calculate AUC
    # SW processing
    data = df.copy()
    data['random_seq'] = data['random_sw_score'].round(0)
    data['SL_reference'] = data['sw_score'].round(0)
    sw_min = plot_cumulative_line(data, folder_name+'cumulative_int_sw.png', cf)

    # SL processing
    data = df.copy()
    data['random_seq'] = (data['random_SL_score'] * 2).round() / 2
    data['SL_reference'] = (data['SL_score'] * 2).round() / 2
    sl_min = plot_cumulative_line(data, folder_name+'cumulative_int_sl.png', cf)

    df = df[df['SL_type'] != 'random']
    reads_sw_solid = len(df[df['sw_score'] >= sw_min])
    potential_read = len(df)
    df = df[df['SL_score'] > sl_min]
    reads_sl_solid = len(df)
    # 图片文件路径（根据你的代码生成的图片名称）
    df['query_length'] = df['query_length'].astype(int).astype(str).astype(int)
    type_table = plot_aligned_length(df,folder_name)

    output_table = pd.DataFrame({
        'SL_type': ['Total', 'Candidate', 'Potential SL', 'SW Solid SL', 'SLRanger Solid SL'],
        'count': [reads_all, reads_na, potential_read, reads_sw_solid, reads_sl_solid]
    })

    output_table = pd.concat([output_table, type_table])
    output_table.columns = ['Variable', 'Read Count']
    output_table['Proportion to total reads (%)'] = output_table['Read Count'] / reads_all
    output_table['Proportion to Candidate reads (%)'] = output_table['Read Count'] / reads_na
    output_table['Proportion to Potential reads (%)'] = output_table['Read Count'] / potential_read
    output_table['Proportion to SLRanger Solid SL reads (%)'] = output_table['Read Count'] / reads_sl_solid

    proportion_cols = [col for col in output_table.columns if 'Proportion' in col]
    for col in proportion_cols:
        output_table[col] = (output_table[col] * 100).round(2)
        output_table.loc[output_table[col] >= 100, col] = '/'

    image_paths = [
        "cumulative_int.png",
        "cumulative_int_sl.png",
        "query_length.png",
        "aligned_length.png",
        "type_pie.png"
    ]

    # 输出文件路径
    output_md = folder_name+"visualization_results.md"
    output_html = folder_name+"visualization_results.html"
    output_table_path = folder_name+"summary_table.csv"
    # 创建Markdown和HTML文件
    output_table.to_csv(output_table_path, index=False)
    create_image_gallery_md_html(output_table,image_paths, output_md, output_html)
    print(f"Markdown file has been created at: {output_md}")
    print(f"HTML file has been created at: {output_html}")