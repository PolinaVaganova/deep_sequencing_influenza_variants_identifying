import pandas as pd
import glob
import os
from pathlib import Path


def get_statistics(path_to_input_df):
    df = pd.read_csv(path_to_input_df, sep='\s+', header=None)

    freq = df[3].str.strip('%').astype(float)

    return [freq.mean(), freq.std()]


def collect_statistics(path_to_control_dfs):
    total_statistics = pd.DataFrame()

    for path_to_df in path_to_control_dfs:
        statistics = get_statistics(path_to_df)
        statistics_df = pd.DataFrame(statistics, index=['mean', 'std'],
                                     columns=[f"{os.path.basename(path_to_df)[:-4]}"])
        total_statistics = pd.concat([total_statistics, statistics_df], axis=1)

    total_statistics.to_csv(os.path.join(path_to_wd, 'controls_statistics.csv'), sep='\t')

    return (total_statistics)


def main():
    path_to_roommate_df = os.path.join(Path.cwd(), 'reference.roommate.variants.csv')
    paths_to_control_dfs = glob.glob(os.path.join(Path.cwd(), 'reference.control_*.csv'))

    roommate_df = pd.read_csv(path_to_roommate_df, sep='\s+', header=None)
    roommate_df[3] = roommate_df[3].str.strip('%').astype(float)

    total_boolean_df = pd.DataFrame()

    for ind in range(total_statistics.shape[1]):
        mean = total_statistics.loc['mean', total_statistics.columns[ind]]
        std = total_statistics.loc['std', total_statistics.columns[ind]]
        boolean_df = (roommate_df[3] > mean + 3 * std) | (roommate_df[3] < mean - 3 * std)
        total_boolean_df = pd.concat([total_boolean_df, boolean_df], axis=1)

    total_boolean_df['result'] = total_boolean_df.all(axis=1)

    filtered_roomate_variants = roommate_df[total_boolean_df['result']]

    filtered_roomate_variants.rename(columns={0: 'ref',
                                              1: 'pos',
                                              2: 'alt',
                                              3: 'freq'},
                                     inplace=True)

    filtered_roomate_variants.to_csv(os.path.join(path_to_wd, 'filtered_roommate_variants.csv'), sep='\t')


if __name__ == "__main__":
    main()