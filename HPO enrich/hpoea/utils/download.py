import os
import pathlib

import requests
from tqdm import tqdm

from hpoea.utils.log import get_logger

log = get_logger(__name__)

DATA_DIR = os.path.join(os.path.expanduser("~"), ".hpoea")

def make_data_dir(path):
    if not os.path.exists(path):
        os.mkdir(path)
    return path

def download_file(url, path):
    """Download file and save it."""
    log.info("Download file from: {}".format(url))
    response = requests.get(url, stream=True)
    with open(path, "wb") as handle:
        for data in tqdm(response.iter_content()):
            handle.write(data)
    log.info("Data download done, file saved to: {}".format(path))
    return path

HPO_GAF_URL = "http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastSuccessfulBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt"

def get_hpo_gaf(url=HPO_GAF_URL, file_name="gaf.txt", file_dir=DATA_DIR):
    """Download the HPO GAF(Gene Associative File)
    It provide the link between genes and HPO term.

    More detail see:
        https://hpo.jax.org/app/download/annotation
    """
    make_data_dir(file_dir)
    path = os.path.join(file_dir, file_name)
    if os.path.exists(path):
        return path
    else:
        return download_file(url, path)
