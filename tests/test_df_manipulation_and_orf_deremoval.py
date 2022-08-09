


'''
def test_count_additions_orf_unremover(duplicate_orf_map:List[List[str]]) -> int:
    """
    Test the number of additions to the rxn dataframe that were made by the orf unremover
    :param duplicate_orf_map: duplicate orf map
    :return: number of additions to the rxn dataframe
    """

    return sum([len(rxn)-1
                for rxn
                in duplicate_orf_map])'''