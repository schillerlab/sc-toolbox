def timestamp():
    """
    Custom timestamp in common EU format.
    """
    import datetime
    import time

    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime("%d-%m-%Y %H:%M:%S")

    return st


def binarize_score(adata, score_label: str, threshold: float):
    """

    Args:
        adata: AnnData object to perform the binarization for
        score_label: Label in adata.obs which will be used for the thresholding
        threshold: The threshold to evaluate for

    Returns:
        List of 'positive' for all scores > threshold and 'negative' else

    """
    result = ["positive" if score > threshold else "negative" for score in adata.obs[score_label]]

    return result
