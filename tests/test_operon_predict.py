import argparse
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from SLRanger import operon_predict as OPERON


def write_minimal_case(temp, sl_type, mapped_gene='geneB'):
    gff = temp / 'genes.gff3'
    mapping = temp / 'mapping.tsv'
    sl_input = temp / 'sl.tsv'
    gff.write_text(
        'chr1\ttest\tgene\t1\t100\t.\t+\t.\tID=geneA\n'
        'chr1\ttest\tmRNA\t1\t100\t.\t+\t.\tID=txA;Parent=geneA\n'
        'chr1\ttest\tCDS\t1\t100\t.\t+\t0\tParent=txA\n'
        'chr1\ttest\tgene\t201\t300\t.\t+\t.\tID=geneB\n'
        'chr1\ttest\tmRNA\t201\t300\t.\t+\t.\tID=txB;Parent=geneB\n'
        'chr1\ttest\tCDS\t201\t300\t.\t+\t0\tParent=txB\n'
    )
    mapping.write_text('read1\t' + mapped_gene + '\n')
    pd.DataFrame(
        {
            'query_name': ['read1'],
            'random_SL_score': [1.0],
            'SL_score': [10.0],
            'SL_type': [sl_type],
        }
    ).to_csv(sl_input, sep='\t', index=False)
    return gff, mapping, sl_input


class OperonPredictTests(unittest.TestCase):
    def test_empty_fusion_returns_empty_list(self):
        fusion = pd.DataFrame(columns=['gene', 'SL', 'count'])
        genes = pd.DataFrame(columns=['gene', 'rank'])
        self.assertEqual(OPERON.fusion_to_ref(fusion, genes), [])

    def test_single_read_is_retained_at_cutoff_boundary(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            sl_path = Path(temp_dir) / 'sl.tsv'
            pd.DataFrame(
                {
                    'query_name': ['read1'],
                    'random_SL_score': [1.0],
                    'SL_score': [10.0],
                    'SL_type': ['SL2'],
                }
            ).to_csv(sl_path, sep='\t', index=False)
            result = OPERON.sl_process(
                sl_path, 4, {'SL1'}, {'SL2'}, legacy_mapping=False
            )
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]['SL'], 'SL2')

    def test_header_only_input_returns_fixed_empty_table(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            sl_path = Path(temp_dir) / 'sl.tsv'
            pd.DataFrame(
                columns=['query_name', 'random_SL_score', 'SL_score', 'SL_type']
            ).to_csv(sl_path, sep='\t', index=False)
            result = OPERON.sl_process(
                sl_path, 4, {'SL1'}, set(), legacy_mapping=True
            )
        self.assertEqual(list(result.columns), ['query_name', 'SL_type', 'SL'])
        self.assertTrue(result.empty)

    def test_legacy_mapping_keeps_sl2_variants(self):
        self.assertEqual(
            OPERON.standardize_sl_type(
                'SL3_unknown', {'SL1'}, set(), legacy_mapping=True
            ),
            'SL2',
        )

    def test_gene_table_expands_fusion_and_deduplicates_reads(self):
        reads = pd.DataFrame(
            {
                'query_name': ['read1', 'read1', 'read2'],
                'gene': ['geneA;geneB', 'geneA;geneB', 'geneA'],
                'SL': ['SL1', 'SL1', 'SL2'],
            }
        )
        annotation = pd.DataFrame(
            {
                'gene': ['geneA'],
                'chromosome': ['chr1'],
                'strand': ['+'],
                'rank': [1],
            }
        )
        result = OPERON.build_gene_sl_table(reads, annotation).set_index('gene')
        self.assertEqual(result.loc['geneA', 'SL1'], 1)
        self.assertEqual(result.loc['geneA', 'SL2'], 1)
        self.assertEqual(result.loc['geneB', 'SL1'], 1)
        self.assertEqual(result.loc['geneB', 'sum_count'], 1)

    def test_operon_fusion_expansion_does_not_duplicate_anchor(self):
        counts = pd.DataFrame(
            {'gene': ['geneA;geneB'], 'SL': ['SL2'], 'count': [1]}
        )
        genes = {
            'geneA': {
                'strand': '+', 'chromosome': 'chr1', 'start': 1, 'end': 100,
            },
            'geneB': {
                'strand': '+', 'chromosome': 'chr1', 'start': 50, 'end': 80,
            },
        }
        expanded, _ = OPERON.fusion_expand(counts, genes)
        self.assertEqual(expanded['gene'].value_counts().to_dict(), {
            'geneA': 1,
            'geneB': 1,
        })

    def test_sl2_only_without_fusion_completes(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            gff, mapping, sl_input = write_minimal_case(temp, 'SL2')
            output = temp / 'operon.gff'
            gene_table = temp / 'gene_sl.tsv'

            args = argparse.Namespace(
                gff=str(gff),
                bam=None,
                mapping=str(mapping),
                input=str(sl_input),
                output=str(output),
                gene_sl_table=str(gene_table),
                sl1_map='SL1',
                sl2_map='SL2',
                distance=5000,
                cutoff=4.0,
            )
            return_code = OPERON.main(args)
            counts = pd.read_csv(gene_table, sep='\t')

            self.assertEqual(return_code, 0)
            self.assertTrue(output.exists())
            self.assertEqual(len(counts), 1)
            self.assertEqual(counts.iloc[0]['gene'], 'geneB')
            self.assertEqual(counts.iloc[0]['SL2'], 1)

    def test_sl1_only_writes_counts_and_skips_operon_output(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            gff, mapping, sl_input = write_minimal_case(temp, 'SL1')
            output = temp / 'operon.gff'
            gene_table = temp / 'gene_sl.tsv'
            args = argparse.Namespace(
                gff=str(gff),
                bam=None,
                mapping=str(mapping),
                input=str(sl_input),
                output=str(output),
                gene_sl_table=str(gene_table),
                sl1_map='SL1',
                sl2_map='SL2',
                distance=5000,
                cutoff=4.0,
            )

            return_code = OPERON.main(args)
            counts = pd.read_csv(gene_table, sep='\t')

            self.assertEqual(return_code, 0)
            self.assertFalse(output.exists())
            self.assertEqual(counts.iloc[0]['SL1'], 1)
            self.assertEqual(counts.iloc[0]['SL2'], 0)

    def test_parser_preserves_original_bam_cli(self):
        args = OPERON.build_parser().parse_args(
            ['-g', 'genes.gff', '-b', 'reads.bam', '-i', 'SLRanger.txt']
        )
        self.assertEqual(args.gff, 'genes.gff')
        self.assertEqual(args.bam, 'reads.bam')
        self.assertIsNone(args.mapping)
        self.assertEqual(args.output, 'SLRanger.gff')


if __name__ == '__main__':
    unittest.main()
