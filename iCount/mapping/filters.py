""".. Line to protect from pydocstyle D205, D400.

Filter hits
-----------

"""


def remove_duplicates(hits):
    """
    Remove duplicate hits that map to same position and have same randomer.

    Records must grouped by contig and sorted by increasing position on
    contig.
    """
    return


def group_by_start(hits):
    """
    Group hits by start.

    Records must grouped by contig and sorted by increasing position on
    contig.
    """
    return


def group_by_end(hits):
    """
    Group hits by their end position.

    Records must grouped by contig and sorted by increasing position on
    contig.
    """
    return


def remove_wrong_assignments(hits_list):
    """
    Remove wrong_assignments.

    Remove low-frequency hits mapped to same position as frequent hits
    from other experiments.

    Records must grouped by contig and sorted by increasing position on
    contig.
    """
    return
