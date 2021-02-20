from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import RandomForestRegressor
from numpy import mean
from numpy import std
import numpy as np
import pandas as pd

def get_n_length_seq(seq ,n):
    '''
    Returns a string of the first n characters of the given sequence.

    If given sequence has less than n characters, appends J until
    the desired length has been reached, then returns that string.

    Note that J does not code for a natural amino acid.
    '''
    seq_length = len(seq)
    if seq_length < n:
        # append dummy variable J to make it length n
        addition = (n - seq_length) * 'J'
        result = seq + addition
    else:
        # retain the first 9 characters of sequence
        result = seq[:n]
    return result


def get_allele_data(df,allele):
    '''
    Returns a DataFrame with data only for the given allele.
    '''
    return df[df['HLA Allele'] == allele]


def split_feature_and_target(df):
    '''
    Returns a tuple of (X,y), where X is a DataFrame
    of the features of the dataset and y is a Series
    of the target.

    Features consist only of the sequence columns
    with individual encoded letters.

    Target is the pIC50 column.
    '''
    drop_list = ['species', 'HLA Allele', 'peptide_length', 'sequence', 'seq_standard', 'IC50 (nM concentration)',
                 'pIC50']
    X = df.drop(drop_list, axis=1)
    y = df['pIC50']
    return X,y


def visualize_data(df, allele,a1,a2,axs):
    '''
    Transforms high dimensional data into two dimensions then plots it.
    '''
    X,y = split_feature_and_target(get_allele_data(df,allele))
    rows = X.shape[0]

    if rows > 0:
        # use PCA to visualize data
        pca = PCA(n_components=2)
        pca.fit(X)
        X_pca = pca.transform(X)
        axs[a1,a2].scatter(X_pca[:,0],X_pca[:,1])


def plot_separate_alleles(mhc_df, class_list):
    '''
    Plots each allele in the class list in a separate plot.
    '''

    fig, axs = plt.subplots(3, 9, figsize=(20, 10))

    a1 = 0
    a2 = 0
    for allele in class_list:
        if a2 > 8:
            a1 += 1
            a2 = 0

        visualize_data(mhc_df, allele, a1, a2, axs)

        a2 += 1


def score_regression_model(df, allele):
    '''
    Trains a linear regression model on data for
    a given allele, using 5-fold cross validation.

    Returns the mean and standard deviation of the
    cross validation scores as a tuple (mean,std_dev).

    Returns (0.0,0.0) if no data exists for the given
    allele.
    '''
    df_allele = get_allele_data(df, allele)
    X, y = split_feature_and_target(df_allele)
    results = (0.0, 0.0)

    rows = X.shape[0]
    if rows > 0:  # data exists for allele

        # setup crossfold validation
        cv = KFold(n_splits=5, random_state=1, shuffle=True)

        # learn model
        regr = RandomForestRegressor(max_depth=2, random_state=0)
        scores = cross_val_score(regr, X, y, scoring='r2', cv=cv)

        results = (mean(scores), std(scores))

    return results


def get_spike_sequence(path_to_fasta, filename):
    '''
    Reads file in FASTA format and returns
    the sequence as a string.
    '''

    with open(path_to_fasta) as f:
        sequence = ''.join(f.read().split('\n')[1:])
        return sequence


def get_spike_array(sequence):
    '''
    Returns a numpy character array of the given sequence.
    '''
    return np.array([char for char in sequence])


def split_spike_sequence(seq, k):
    '''
    Returns a generator for k-mers.

    Each iteration of the generator returns a tuple containing
    the next k-mer of the sequence in following form:
    (start index, end index, k-mer)
    '''

    start = 0
    end = start + (k - 1)
    last_position = len(seq) - 1

    while end <= last_position:
        yield (start, end, seq[start:end + 1])

        # get next sequence
        start = start + 1
        end = start + (k - 1)


def encode_seq(sequence):
    '''
    Converts numpy character array to numpy integer array
    where each character is the ASCII encoded value.
    '''
    encoder = lambda t: ord(t)
    vfunc = np.vectorize(encoder)
    return vfunc(sequence)


def split_into_features(seq):
    '''
    Takes an encoded sequence and returns a
    DataFrame with each element as its own
    column.
    '''
    return pd.DataFrame(seq).T


def get_pIC50_from_data(df, allele, sequence):
    '''
    Checks given DataFrame for a row with the given allele
    and sequence. Returns pIC50 value if a matching row is
    found, returns None otherwise.
    '''

    # find all matches in DataFrame for given allele and sequence
    matches = df[(df['HLA Allele'] == allele) & (df['seq_standard'] == sequence)]

    if matches.shape[0] > 0:

        # take the first match and extract the pIC50 value from it
        pIC50 = matches['pIC50'].iloc[0]

    else:

        pIC50 = None

    return pIC50


def get_sequence_string(seq):
    '''
    Returns a string of the given sequence.
    '''
    return ''.join(seq.tolist())


def predict_pIC50_values_for_allele(pIC50_df, model, allele, spike_array, k):
    '''
    Given a regression model, allele, spike sequence array,
    and k-mer length, returns a list of pIC50 predictions in
    the form of a tuple (allele,peptide,start,end,pIC50).

    If an exact sequence match is found in the pIC50 data itself,
    uses that pIC50 value instead of predicting a value.
    '''
    results = []

    # executes once for each k-mer in the sequence
    for start, end, arr in split_spike_sequence(spike_array, k):

        # determine if the k-mer exists in the data
        peptide = get_sequence_string(arr)
        pIC50 = get_pIC50_from_data(pIC50_df, allele, peptide)

        if not pIC50:
            # predict the pIC50 value
            encoded = encode_seq(arr)
            features = split_into_features(encoded)
            pIC50 = model.predict(features)[0]

        # save result
        results.append((allele, peptide, start, end, pIC50))

    return results


def train_regresson_model(df, allele):
    '''
    Fits a regression model to the data for a
    given allele and returns the model.

    Returns None if no data exists for the given
    allele.
    '''
    df_allele = get_allele_data(df, allele)
    X, y = split_feature_and_target(df_allele)
    model = None

    rows = X.shape[0]
    if rows > 0:  # data exists for allele

        model = RandomForestRegressor(max_depth=2, random_state=0)
        model.fit(X, y)

    return model


def predict_pIC50_values(allele_list, df, seq, k):
    '''
    For each allele in the allele list, trains a regression
    model using the given DataFrame, then uses that model to
    predict pIC50 for each k-mer in the given protein sequence.

    Returns a list of all results.
    '''

    results = []

    for allele in allele_list:

        print('.', end='')  # visualizing progress

        # train model for this allele
        model = train_regresson_model(df, allele)

        # use model as long as it is not None
        if model:
            # predict pIC50 for each k-mer in the sequence
            result = predict_pIC50_values_for_allele(df, model, allele, seq, k)
            results = results + result

    print()  # newline for visualization
    return results


def calculate_metrics(mhc_class, predictions):
    '''
    Calculates min, max, mean, median, and standard deviation
    of all results for a given MHC class.

    Returns the metrics in a tuple with the form:
    (min, max, mean, median, std_dev)
    '''

    # create numpy array of the pIC50 values
    a = np.array([t[4] for t in predictions])

    return [mhc_class, np.min(a), np.max(a), np.mean(a), np.median(a), np.std(a)]
