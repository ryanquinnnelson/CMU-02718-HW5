import pandas as pd


def get_list_above_p_min(predictions, p_min):
    """
    Returns a list of predictions which have a pIC50 value above p_min.
    :param predictions:
    :param p_min:
    :return:
    """

    results = []

    for prediction in predictions:
        pIC50 = prediction[4]
        if pIC50 >= p_min:
            results.append(prediction)

    return results


def get_peptide_coverage_dict(predictions_list):
    """
    Builds and returns a dictionary in which keys are every
    peptide in the prediction list and the value for each
    key is a set of alleles which are covered by that peptide.
    :param predictions_list:
    :return:
    """

    # maintain a set of peptides which cover each allele
    coverage_dict = dict()

    for allele, peptide, start, end, pIC50 in predictions_list:

        # define key to avoid collisions
        key = peptide + '_' + str(start) + '_' + str(end)

        if key in coverage_dict:
            # append allele to the list
            coverage_dict[key].add(allele)
        else:
            # add set to dictionary
            new_set = set()
            new_set.add(allele)
            coverage_dict[key] = new_set

    return coverage_dict


def build_row_df(coverage_dict, key, alleles):
    """
    Builds and returns a DataFrame consisting of a single row
    with the following columns:

    - peptide
    - start position
    - end position
    - allele flag for each allele in alleles
    - count

    Allele flag is 1 if allele binds with this peptide,
    0 otherwise.

    Count is a sum of all allele flags for this peptide
    :param coverage_dict:
    :param key:
    :param alleles:
    :return:
    """

    peptide, start, end = key.split('_')

    new_dict = dict()
    new_dict['peptide'] = [peptide]
    new_dict['start'] = [int(start)]
    new_dict['end'] = [int(end)]

    # create row for DataFrame of all alleles
    counter = 0
    for allele in alleles:

        if allele in coverage_dict[key]:
            new_dict[allele] = [1]
            counter += 1
        else:
            new_dict[allele] = [0]

    new_dict['count'] = [counter]
    return pd.DataFrame.from_dict(new_dict)


def build_peptide_coverage_df(coverage_dict, alleles):
    """
    Builds and returns a DataFrame which contains a row for
    each peptide in the coverage dictionary with columns
    for each allele in alleles.
    :param coverage_dict:
    :param alleles:
    :return:
    """

    counter = 0
    df = None
    for key in coverage_dict:

        # build row
        df_row = build_row_df(coverage_dict, key, alleles)

        if counter == 0:
            df = df_row
        else:
            df = df.append(df_row)

        counter += 1
    return df


def get_peptide_coverage_df(predictions_list, p_min, alleles):
    """
    Selects from the given predictions list those
    peptides which have pIC50 at or above p_min, then
    builds a DataFrame representing the allele coverage
    for each of the selected peptides.
    :param predictions_list:
    :param p_min:
    :param alleles:
    :return:
    """
    above = get_list_above_p_min(predictions_list, p_min)
    if len(above) > 0:
        p_dict = get_peptide_coverage_dict(above)
        df = build_peptide_coverage_df(p_dict, alleles)
    else:
        df = pd.DataFrame()  # empty DataFrame
    return df


def get_allele_least_coverage(df_coverage, uncovered):
    """
    Given the peptide coverage DataFrame, returns the
    allele with the minimum number of peptides that cover it
    which has a coverage greater than zero.
    :param df_coverage:
    :param uncovered:
    :return:
    """
    # sum the coverage for each allele column
    # then sort the values so the smallest count is at the top
    counted_and_sorted_series = df_coverage.drop(['peptide', 'start', 'end', 'count'], axis=1).sum(axis=0).sort_values()
    df = pd.DataFrame(counted_and_sorted_series, columns=['total']).reset_index()
    df = df.rename(columns={'index': 'HLA Allele'})

    # drop all alleles with a count of zero
    df = df[df['total'] > 0]

    # retain only those alleles which are in uncovered
    df = df[df['HLA Allele'].isin(uncovered)]

    if df.shape[0] > 0:

        # there is an allele which can still be covered
        allele = df.iloc[0]['HLA Allele']
        total = df.iloc[0]['total']
    else:
        allele = 'NA'
        total = -1

    return (allele, total)


def find_first_peptide_cover(df_coverage, allele):
    """
    Returns a dictionary of the first peptide in the
    DataFrame of available peptides which covers the given allele, or an
    empty dictionary if no cover is found.
    :param df_coverage:
    :param allele:
    :return:
    """

    # find all peptides which covers allele
    df = df_coverage[df_coverage[allele] == 1]

    if df.shape[0] > 0:

        # select the peptide with the highest coverage over all alleles
        df_peptide = df.sort_values(by=['count'], ascending=False).head(1)

        # return dictionary of the peptide row
        return df_peptide.to_dict()
    else:

        # return an empty dictionary
        return dict()


def update_selected_peptides(peptide_dict, selected):
    """
    Adds peptide string to selected list.
    :param peptide_dict:
    :param selected:
    :return:
    """
    # add to list of selected
    peptide = peptide_dict['peptide'][0]
    selected.append(peptide)


def update_uncovered_alleles(peptide_dict, uncovered, num_covered):
    """
    Removes each allele covered by given peptide from uncovered list.
    :param peptide_dict:
    :param uncovered:
    :param num_covered:
    :return:
    """
    n = num_covered

    for key in peptide_dict:

        # for each allele flag in the peptide dictionary
        if key in uncovered:

            # mark off the alleles covered by this peptide if they have not been covered
            flag = peptide_dict[key][0]
            if flag == 1:
                # peptide covers this allele
                uncovered.remove(key)
                n += 1
    return n


def get_nonoverlapping_peptides(peptide_dict, df):
    """
    Removes any peptides from the given DataFrame
    if their sequence overlaps with that of the selected
    peptide, then returns a new DataFrame of the
    remaining peptides which are available.
    :param peptide_dict:
    :param df:
    :return:
    """

    # remove any peptides from the df which overlap with the selected peptide
    start = peptide_dict['start'][0]
    end = peptide_dict['end'][0]
    removed1 = df.loc[(df['start'] >= start) & (df['start'] <= end)]
    removed2 = df.loc[(df['end'] >= start) & (df['end'] <= end)]
    retained = df.loc[(df['end'] < start) | (df['start'] > end)]

    return (removed1, removed2, retained)


def select_peptides(df_coverage, alleles, design_limit):
    """
    Selects peptides based on predetermined criteria.
    :param df_coverage:
    :param alleles:
    :param design_limit:
    :return:
    """

    # list of selected peptides
    selected = []

    # list of alleles not yet covered
    uncovered = alleles.copy()
    num_covered = 0

    # select the allele with the lowest nonzero coverage
    allele, total = get_allele_least_coverage(df_coverage, uncovered)
    selected_peptide = find_first_peptide_cover(df_coverage, allele)

    # update data structures based on selection
    update_selected_peptides(selected_peptide, selected)
    num_covered = update_uncovered_alleles(selected_peptide, uncovered, num_covered)
    r1, r2, retained = get_nonoverlapping_peptides(selected_peptide, df_coverage)

    # as long as
    # 1 - the number of peptides selected has not exceeded the design limit
    # 2 - there are alleles which have not yet been covered
    # 3 - there are peptides left to check
    while len(selected) < design_limit and len(uncovered) != 0 and retained.shape[0] > 0:

        # continue searching
        allele, total = get_allele_least_coverage(retained, uncovered)

        if total == -1:

            # no alleles remain which can be covered
            break

        else:

            # an allele can still be covered
            selected_peptide = find_first_peptide_cover(retained, allele)

            # update data structures based on selection
            update_selected_peptides(selected_peptide, selected)
            num_covered = update_uncovered_alleles(selected_peptide, uncovered, num_covered)
            r1, r2, retained = get_nonoverlapping_peptides(selected_peptide, retained)

    # end while - search is over
    return (num_covered, selected)


def design_coverage(allele_list, mhc_results, design_limits, p_mins):
    """

    :param allele_list:
    :param mhc_results:
    :param design_limits:
    :param p_mins:
    :return:
    """
    designs = []
    coverages = []

    for d in design_limits:
        for p in p_mins:

            coverage_df = get_peptide_coverage_df(mhc_results, p, allele_list)
            if coverage_df.shape[0] > 0:

                # there is at least one peptide which has pIC50 >= p_min
                num_selected, peptides = select_peptides(coverage_df, allele_list, d)

            else:

                # no peptides have pIC50 >= p_min
                num_selected, peptides = (0, [])

            designs.append((d, p, num_selected, peptides))
            coverages.append(coverage_df)

    return (designs, coverages)


def get_design_coverage(design, coverage_df):
    """
    Returns a DataFrame containing only the rows from the
    coverage DataFrame which are peptides in the design.
    :param design:
    :param coverage_df:
    :return:
    """
    return coverage_df[coverage_df['peptide'].isin(design)]


def get_allele_distribution_df(design, coverage_df):
    """
    Returns a DataFrame of alleles and coverage
    counts for each allele.
    :param design:
    :param coverage_df:
    :return:
    """

    # select only rows for peptides in design
    df = get_design_coverage(design, coverage_df)

    # remove non-allele columns
    df = df.drop(['peptide', 'start', 'end', 'count'], axis=1)

    # sum each allele column then convert Series back to DataFrame
    df = df.sum(axis=0).to_frame().reset_index()

    # rename DataFrame columns
    df = df.rename(columns={'index': 'HLA Allele', 0: 'Count'})

    return df


def select_allele_with_min_coverage(distribution_df, alleles):
    """
    From alleles with nonzero coverage, selects one which
    has the minimum coverage count and which is also in the
    given list of alleles. Returns None if no minimum
    allele is found.
    :param distribution_df:
    :param alleles:
    :return:
    """
    first = None

    # select alleles with nonzero counts
    df = distribution_df[(distribution_df['Count'] > 0) & distribution_df['HLA Allele'].isin(alleles)]

    if df.shape[0] > 0:
        # sort by Count
        df = df.sort_values(by=['Count'], ascending=True)

        # select the first record
        first = df['HLA Allele'].iloc[0]

    return first


def get_available_coverage(design, coverage_df):
    """
    Returns a DataFrame containing only the rows from the
    coverage DataFrame which are not peptides in the design.
    :param design:
    :param coverage_df:
    :return:
    """
    return coverage_df[~coverage_df['peptide'].isin(design)]


def update_allele_distribution(distribution_df, peptide_dict):
    """
    Updates the coverage count for each allele covered by the
    given peptide and returns the updated DataFrame.
    :param distribution_df:
    :param peptide_dict:
    :return:
    """

    # get list of alleles in distribution
    allele_list = distribution_df['HLA Allele'].tolist()

    for key in peptide_dict:

        # for each allele flag in the peptide dictionary
        if key in allele_list:

            # mark off the alleles covered by this peptide if they have not been covered
            flag = peptide_dict[key][0]
            if flag == 1:
                # peptide covers this allele
                prev = distribution_df.loc[distribution_df['HLA Allele'] == key].iloc[0]['Count']
                distribution_df.loc[distribution_df['HLA Allele'] == key, 'Count'] = prev + 1


def increase_design_fairness(design, coverage_df, design_limit):
    """
    While the design has less peptides than the design_limit,
    adds additional peptides which cover the least-covered
    alleles.

    Returns the supplemented design as a list of peptides.
    :param design:
    :param coverage_df:
    :param design_limit:
    :return:
    """
    updated_design = design.copy()

    # calculate allele distribution
    dist = get_allele_distribution_df(design, coverage_df)

    # get list of alleles with nonzero coverage
    allele_list = dist[dist['Count'] > 0]['HLA Allele'].tolist()

    # get list of peptides not currently in design
    available = coverage_df[~coverage_df['peptide'].isin(design)]

    while len(updated_design) < design_limit and available.shape[0] > 0 and len(allele_list) > 0:

        # select allele with least nonzero coverage
        allele = select_allele_with_min_coverage(dist, allele_list)

        # remove allele from list of consideration
        allele_list.remove(allele)

        # select peptide to cover it from available choices
        selected_peptide = find_first_peptide_cover(available, allele)

        if selected_peptide:
            # update allele distribution
            update_allele_distribution(dist, selected_peptide)

            # update available peptides
            r1, r2, available = get_nonoverlapping_peptides(selected_peptide, available)

            # update design
            update_selected_peptides(selected_peptide, updated_design)

    return updated_design


def execute_design_supplement(mhc_class, designs, coverages):
    """

    :param mhc_class:
    :param designs:
    :param coverages:
    :return:
    """
    designs_after = []
    for i, design in enumerate(designs):
        coverage = coverages[i]

        design_limit, p_min, num_selected, design_before = design
        print('before', design_limit, p_min, num_selected, design_before)
        design_after = increase_design_fairness(design_before, coverage, design_limit)
        print('after ', design_limit, p_min, num_selected, design_after)
        designs_after.append((design_limit, p_min, num_selected, design_after))

    return designs_after


def get_uncovered_alleles(design, coverage_df, allele_list):
    """
    Returns a list of alleles which are not covered by the given design.
    :param design:
    :param coverage_df:
    :param allele_list:
    :return:
    """

    # limit coverage to only peptides in the design
    df = coverage_df[coverage_df['peptide'].isin(design)]

    # sum up counts for alleles
    dist = get_allele_distribution_df(design, df)

    # get a list of alleles which are not covered
    return dist[dist['Count'] < 1]['HLA Allele'].tolist()


def select_best_10(supple_10_9, orig_10_6, supple_20_6, coverage_supple_20_6):
    """
    Selects best peptides from the given designs using predefined
    criteria and returns a list of the 10 best peptides.

    Criteria:
     - Select all peptides from the supplemented design for
       design_limit = 10 and p_min = 9.0.

     - Select all peptides from the minimum coverage design for
       design_limit = 10 and p_min = 6.0.

     - With any remaining space, from the peptides in the
       supplemented design for design_limit = 20 and p_min = 6.0
       which have not been selected already, chose those
       which have the highest individual cover over alleles.
    :param supple_10_9:
    :param orig_10_6:
    :param supple_20_6:
    :param coverage_supple_20_6:
    :return:
    """

    best = []

    # select all from 9.0 supplemented cover
    best += supple_10_9

    # Add all from 6.0 non-supplemented cover to ensure coverage of all alleles.
    best += orig_10_6

    # Remove duplicates
    best = set(best)

    # Determine which peptides have not been selected from supplemented 6.0 cover.
    supple_20_6_set = set(supple_20_6)
    remaining = supple_20_6_set - best

    # Determine how many additional peptides are needed to reach 10 total
    num_needed = 10 - len(best)

    # Fill remaining space with unselected peptides from the
    # supplemented 6.0 cover which cover the most alleles.
    selected_df = coverage_supple_20_6[coverage_supple_20_6['peptide'] \
        .isin(remaining)] \
        .sort_values(by=['count'], ascending=False) \
        .head(num_needed)
    selected = selected_df['peptide'].tolist()
    return list(best) + selected


def build_1_5_list(start_num, mhc_class, chosen, predictions):
    """

    :param start_num:
    :param mhc_class:
    :param chosen:
    :param predictions:
    :return:
    """
    temp = []

    for allele, peptide, start, end, pIC50 in predictions:
        if peptide in chosen:
            # add tuple to temporary list
            temp.append((peptide, allele, mhc_class, pIC50))

    # sort temp by peptide
    temp = sorted(temp, key=lambda x: x[0])

    # set number for each peptide
    curr_num = start_num
    temp2 = []
    visited = []

    for peptide, allele, mhc_class, pIC50 in temp:
        if peptide not in visited:
            visited.append(peptide)
            curr_num += 1

        temp2.append((curr_num, peptide, allele, mhc_class, pIC50))

    return temp2
