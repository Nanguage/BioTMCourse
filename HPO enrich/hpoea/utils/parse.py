

def hpo_table_header(path):
    with open(path) as f:
        line = f.readline().strip()
    columns = line.replace("#Format: ", "").split("<tab>")
    columns = [i.replace("-", "_") for i in columns]
    return columns


def parse_hpo_gaf(path):
    """Parse the HPO GAF,
    return a pandas dataframe.
    """
    import pandas as pd
    # parse header
    columns = hpo_table_header(path)
    df = pd.read_table(path, skiprows=1, header=None)
    df.columns = columns
    return df


def parse_hpo_obo(path):
    """Parse HPO obo file.
    return a network MultiDiGraph object
    """
    import obonet
    graph = obonet.read(path)
    return graph
