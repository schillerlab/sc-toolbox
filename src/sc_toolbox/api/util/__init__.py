def timestamp():
    """
    Custom timestamp in common EU format.
    """
    import datetime
    import time

    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime("%d-%m-%Y %H:%M:%S")

    return st
