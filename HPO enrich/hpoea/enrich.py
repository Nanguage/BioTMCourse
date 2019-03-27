from hpoea.utils.log import get_logger
log = get_logger(__name__)


class GSEA(object):
    """Class for Gene Set Enrichment Analyze. 

    Args
    ----
    gaf : str, optional
        Path GAF(Gene Associative File).
        This file provide the link between genes and HPO term.
        (see: https://hpo.jax.org/app/download/annotation)
        If not specified, will download it automatically.
    """
    def __init__(self, gaf=None):
        from scipy import stats
        if not gaf_file:  # download HPO GAF
            from hpoea.utils.download import get_hpo_gaf
            gaf = get_hpo_gaf()
        from hpoea.utils.parse import parse_hpo_gaf
        self.gaf = parse_hpo_gaf(gaf)

    def enrich(self, gene_list, gene_id_source="entrez_id"):
        """Perform gene set enrichment analyze.

        Args
        ----
        gene_list : iterable
            Input gene list, it's the subset of genes in the GAF.
        gene_id_source : {"entrez_id", "entrez_symbol", "ensembl"}, optional
            The source of input genes, default in "entrez_id",
            same to the ID source used in GAF.
            If specified other formats, will convert to "entrez_id".

        Return
        ------
        enrichment_table : pandas.DataFrame
        """
        import pandas as pd
        gene_ids = self._to_entrez_id(gene_list, gene_id_source)
        possi_terms = self._get_all_possible_terms(gene_ids)

        pvals_uncorr = []
        counts = []
        for term_id in possi_terms:  # calculate the p-value of each possible terms
            counts = self._get_counts(term_id)
            pval = _calc_pvalue(*counts)
            pvals_uncorr.append(pval)
            counts.append(counts)
        study_count, n_study, population_count, n_population = zip(*counts)
        
        enrichment_table = pd.DataFrame({
            'HPO_term_ID': possi_terms,
            'study_count': study_count,
            'n_study': n_study,
            'population_count': population_count,
            'n_population': n_population,
            'pvalue': pvals_uncorr,
        })
        self.enrichment_table = enrichment_table

    def multiple_test_corretion(self, method='fdr_bh'):
        """Perform multiple test correction, adjust the p-value.

        Args
        ----
        method : method, optional
            Methods used for correction. default 'fdr_bh'.
            All supported method list see `statsmodels.stats.multitest.multipletests`:
                http://www.statsmodels.org/devel/generated/statsmodels.stats.multitest.multipletests.html#statsmodels.stats.multitest.multipletests
        """
        from statsmodels.stats.multitest import multipletests
        df = self.enrichment_table
        pvals = df['pvalue']
        padjs = multipletests(pvals)
        df['padj'] = padjs

    def _calc_pvalue(self, study_count, study_n, pop_count, pop_n):
        """Calculate uncorrected p-value.

        See Also
        --------
        scipy.stats.fisher_exact :
            http://docs.scipy.org/doc/scipy-0.17.0/reference/generated/scipy.stats.fisher_exact.html
        """
        avar = study_count
        bvar = study_n - study_count
        cvar = pop_count - study_count
        dvar = pop_n - pop_count - bvar
        assert cvar >= 0, "(pop_count - study_count) must large than zero."
        oddsratio, p_val = stats.fisher_exact([[avar, bvar], [cvar, dvar]])
        return p_val

    def _get_counts(self, term_id):
        """Get the counts value used for fisher exact test."""
        pop_cnt = self.gaf[self.gaf.HPO_Term_ID == term_id].shape[0]
        pop_n = self.gaf.shape[0]
        study_cnt = self._study[self._study.HPO_Term_ID == term_id].shape[0]
        study_n = self._study.shape[0]
        counts = (study_cnt, study_n, pop_cnt, pop_n)
        return counts

    def _get_all_possible_terms(self, gene_list):
        study = self.gaf.entrez_gene_id.isin(gene_list)
        self._study = study
        possible_hpo_term_ids = list(study.HPO_Term_ID.drop_duplicates())
        return possible_hpo_term_ids
    
    def _to_entrez_id(self, gene_list, source):
        gene_list = list(gene_list)
        if source == "entrez_id":
            return gene_list
        elif source == "entrez_symbol":
            from hpoea.utils.idconvert import EntrezConvert
            self._entrez_converter = EntrezConvert(self.gaf)
            ids = self._entrez_converter.symbol2id(gene_list)
            return ids
        elif source == "ensembl":
            from hpoea.utils.idconvert import EntrezEnsemblConvert
            self._entrez_ensembl_converter = EntrezEnsemblConvert()
            ids = self._entrez_ensembl_converter.ensembl2entrez(gene_list)
            return ids 
        else:
            raise NotImplemented("gene id source only support: entrez_id | entrez_symbol | ensembl")
