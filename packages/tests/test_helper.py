import packages.vaccine.helper as hp
import pandas as pd
import numpy as np

def test_get_n_length_seq():
    assert hp.get_n_length_seq('12345678', 10) == '12345678JJ'
    assert hp.get_n_length_seq('123456789', 9) == '123456789'
    assert hp.get_n_length_seq('1234567890', 9) == '123456789'


def test_split_spike_sequence():
    test_array = ['M', 'F', 'V', 'F', 'L', 'V', 'L', 'L', 'P', 'L', 'V', 'S']
    results = []
    for each in hp.split_spike_sequence(test_array, 9):
        results.append(each)

    assert len(results) == 4
    first = results[0]
    second = results[1]
    assert first[0] == 0
    assert first[1] == 8
    assert first[2] == ['M', 'F', 'V', 'F', 'L', 'V', 'L', 'L', 'P']
    assert second[0] == 1
    assert second[1] == 9
    assert second[2] == ['F', 'V', 'F', 'L', 'V', 'L', 'L', 'P', 'L']


def test_encode_seq():
    test_array = ['M', 'F', 'V', 'F']
    result = hp.encode_seq(test_array)
    assert result.tolist() == [77, 70, 86, 70]
    assert chr(result[0]) == 'M'


def test_split_into_features():
    test_array = ['M', 'F', 'V', 'F']
    result = hp.split_into_features(hp.encode_seq(test_array))

    assert result.shape[0] == 1
    assert result.shape[1] == 4


def test_get_pIC50_from_data():
    # test dataset
    c1 = ['HLA-A*01:02', 'HLA-A*02:01', 'HLA-A*01:01']
    c2 = ['CVADYSVLY', 'CVADYSVLY', 'CVADYSVLY']
    c3 = [7, 7.8, 7.874116247000001]
    test_dataframe = pd.DataFrame(list(zip(c1, c2, c3)), columns=['HLA Allele', 'seq_standard', 'pIC50'])

    # match exists
    result = hp.get_pIC50_from_data(test_dataframe, 'HLA-A*01:01', 'CVADYSVLY')
    assert result == 7.874116247000001

    # match doesn't exist
    result = hp.get_pIC50_from_data(test_dataframe, 'HLA-A*02:01', 'CVADYSVLY')
    if not result:
        assert False

def test_get_sequence_string():
    assert hp.get_sequence_string(np.array(['M', 'F', 'V'])) == 'MFV'

def test_calculate_metrics():
    test_data = [('HLA-A*01:01', 'NCVADYSVL', 0, 8, 4.691983176635914),
                 ('HLA-A*01:01', 'CVADYSVLY', 1, 9, 7.874116247000001),
                 ('HLA-A*02:01', 'NCVADYSVL', 0, 8, 5.292583412663456),
                 ('HLA-A*02:01', 'CVADYSVLY', 1, 9, 5.3770656045923655)]

    result = hp.calculate_metrics('II', test_data)
    assert result[0] == 'II'
    assert result[1] == 4.691983176635914
    assert result[2] == 7.874116247000001
    assert result[3] == 5.808937110222935
    assert result[4] == 5.334824508627911
    assert result[5] == 1.2212375992301565
