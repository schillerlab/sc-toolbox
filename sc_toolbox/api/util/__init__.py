from typing import List

from anndata import AnnData


def timestamp():
    """
    Custom timestamp in common EU format.

    Returns:
        datetime timestamp
    """
    import datetime
    import time

    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime("%d-%m-%Y %H:%M:%S")

    return st


def binarize_score(adata, score_label: str, threshold: float):
    """
    Binarizes a provided key of an AnnData object by labeling values over a threshold as 'positive' or 'negative'.

    Args:
        adata: AnnData object to perform the binarization for
        score_label: Label in adata.obs which will be used for the thresholding
        threshold: The threshold to evaluate for

    Returns:
        List of 'positive' for all scores > threshold and 'negative' else
    """
    result = ["positive" if score > threshold else "negative" for score in adata.obs[score_label]]

    return result


def read_concatenate_to_adata(file_paths: List[str]) -> AnnData:
    """
    Parses a list of file paths and concatenates them memory efficiently into a single AnnData object.

    Args:
        file_paths: List of file paths

    Returns:
        Single AnnData object containing all concatenated files
    """
    import anndata as ad
    import gc
    import scipy.sparse

    for index, file_name in enumerate(file_paths):
        if index == 1:
            adata_concatenated = ad.read(file_name)
            adata_concatenated.X = scipy.sparse.csr_matrix(adata_concatenated.X)
        else:
            adata = ad.read(file_name)
            adata.X = scipy.sparse.csr_matrix(adata.X)
            gc.collect()
            adata_concatenated = adata_concatenated.concatenate(adata)

    return adata_concatenated
