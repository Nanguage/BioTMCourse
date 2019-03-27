from hpoea.utils.log import get_logger

log = get_logger(__name__)


class EntrezEnsemblConvert(object):
    """Conversion between entrez_id and ensembl_id.
    """

    def __init__(self):
        import mygene
        self.mg = mygene.MyGeneInfo()

    def ensembl2entrez(self, gene_list):
        log.info("Perform gene ID conversion: from ensembl id to entrez id.")
        resp = self.mg.querymany(gene_list, scopes='ensembl.gene', fields='entrezgene')

        self._conversion = {}
        res = []
        for r_ in resp:
            entrez = r_.get('entrezgene')
            if entrez:
                self._conversion[r_['query']] = int(entrez)
                res.append(int(entrez))
            else:
                self._conversion[r_['query']] = None

        if len(res) < len(gene_list):
            failed = [k for k,v in self._conversion.items() if v is None]
            log.warning("{} gene conversion failed, failed genes:".format(len(failed)))
            log.warning(",".join(failed))
        return res



class EntrezConvert(object):
    """Conversion between entrez_id and entrez_name
    """
    def __init__(self, gaf_df):
        self.gaf = gaf_df

    def id2symbol(self, ids):
        df_ = self.gaf[['entrez_gene_id', 'entrez_gene_symbol']]
        df_ = df_.drop_duplicates()
        df_.index = df_.entrez_gene_id
        symbols = list( df_[ids].entrez_gene_symbol )
        return symbols

    def symbol2id(self, names):
        df_ = self.gaf[['entrez_gene_id', 'entrez_gene_symbol']]
        df_ = df_.drop_duplicates()
        df_.index = df_.entrez_gene_
        ids = list( df_[ids].entrez_gene_id )
        return ids