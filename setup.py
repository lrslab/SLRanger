from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="SLRanger",
    version="0.0.0.9",
    author="GUO Zhihao",
    author_email="qhuozhihao@icloud.com",
    description='An integrated approach for spliced leader detection and operon prediction in eukaryotes using long RNA reads',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lrslab/SLRanger",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7.0,<=3.11.7',
    install_requires=[
        'bio>=1.5.0',
        'numpy>=1.23.0,<2.0.0',
        'pandas>=2.1.0',
        'plotnine>=0.12.4',
        'tqdm>=4.0.0',
        "pysam>=0.21.0",
        "biopython>=1.80",
        'pyssw==0.1.7',
        'trackcluster==0.1.7',
        "scikit-learn>=1.0.2",
        'tabulate>=0.8.0',
        'Markdown>=3.5',
        'seaborn>0.12.0'
    ],
    scripts=['SLRanger/SL_detect.py','SLRanger/operon_predict.py','SLRanger/add_gene.py']
)