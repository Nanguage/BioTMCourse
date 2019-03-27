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
        gene_ids = self._to_entrez_id(gene_list, gene_id_source)
        possi_items = self._get_all_possible_items(gene_ids)

    def _get_all_possible_items(self, gene_list):
        pass
    
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
