import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import re
import numpy as np

MAX_K = 20

def span_length_k_keys(code, k):
    pattern = r'pre_(\d+)_end_(\d+)'
    match = re.match(pattern, code)
    if match:
        num1, num2 = map(int, match.groups())  # Extract and convert to integers
        return num2 == num1 + k  # Check if second number is first + k
    return False

def is_span_length_less_euqal_than_k(code, k):
    pattern = r'pre_(\d+)_end_(\d+)'
    match = re.match(pattern, code)
    if match:
        num1, num2 = map(int, match.groups())  # Extract and convert to integers
        return num2 - num1 <= k  # Check if second number is first + k
    return False

def reduce_span_entropy_to_average(data_df):
    
    result_df = data_df.copy()  # Start with a copy of the original DataFrame
    print(f"delete spans that length > {MAX_K}")
    result_df = result_df[result_df['Key'].apply(lambda x: is_span_length_less_euqal_than_k(x, MAX_K))]
    print(f"new df length = {len(result_df)}" )

    reduce_records_plan_to_append = []
    # Loop over k and generate the new rows
    for k in tqdm(range(2, MAX_K)):
        # Filter rows where the 'Key' matches the span_length_k_keys condition
        filtered_result_df = result_df[result_df['Key'].apply(lambda x: span_length_k_keys(x, k))]
        
        # Compute the average for each 'id'
        averages = (
            filtered_result_df
            .groupby('id')['Value']
            .mean()
            .reset_index()
            .rename(columns={'Value': 'average_value'})
        )
        
        # Create the new span_k rows
        span_k_rows = averages.copy()
        span_k_rows['Key'] = f'span_{k}'
        span_k_rows['Value'] = span_k_rows['average_value']
        span_k_rows = span_k_rows[['id', 'Key', 'Value']]  # Keep necessary columns
        
        # Concatenate the new rows into the final DataFrame
        reduce_records_plan_to_append.append(span_k_rows.copy())

        result_df = result_df[~result_df['Key'].apply(lambda x: span_length_k_keys(x, k))]
        print(f"    - new DF length = {len(result_df)}" )
    result_df = pd.concat([result_df] + reduce_records_plan_to_append, ignore_index=True)
    return result_df


def load_dataset_to_data_frame(file_path):
    print(f"load data from: {file_path}")
    return pd.read_csv(file_path, sep=',')
    
    
def one_way_anova(key_name):
    from scipy.stats import f_oneway

    seizure_data = seizures[seizures.Key == key_name]['Value'].dropna()
    preepileptic_data = preepileptic[preepileptic.Key == key_name]['Value'].dropna()
    normal_data = normal[normal.Key == key_name]['Value'].dropna()
    return f_oneway(seizure_data, preepileptic_data, normal_data)


def pairwise_comparisons(key_name):
    from scipy.stats import mannwhitneyu
    from statsmodels.stats.multitest import multipletests

    # Filter the data based on the key_name
    seizure_data = seizures[seizures.Key == key_name]['Value'].dropna()
    preepileptic_data = preepileptic[preepileptic.Key == key_name]['Value'].dropna()
    normal_data = normal[normal.Key == key_name]['Value'].dropna()

    # Pairwise comparisons (Mann-Whitney U test)
    groups = {
        'Seizure vs Preepileptic': (seizure_data, preepileptic_data),
        'Seizure vs Normal': (seizure_data, normal_data),
        'Preepileptic vs Normal': (preepileptic_data, normal_data)
    }

    results = []
    for comparison, (group1, group2) in groups.items():
        stat, p_value = mannwhitneyu(group1, group2)
        results.append((comparison, stat, p_value))

    # Adjust p-values using Bonferroni correction
    p_values = [result[2] for result in results]
    adjusted_p_values = multipletests(p_values, method='bonferroni')[1]

    return adjusted_p_values

def pairwise_comparisons_KL(seizures_df, preepileptic_df, normal_df, key_name):
    import numpy as np
    from scipy.stats import mannwhitneyu
    from scipy.stats import entropy
    from statsmodels.stats.multitest import multipletests
    from sklearn.neighbors import KernelDensity


    # Filter the data based on the key_name
    seizure_data = seizures_df[seizures_df.Key == key_name]['Value'].dropna()
    preepileptic_data = preepileptic_df[preepileptic_df.Key == key_name]['Value'].dropna()
    normal_data = normal_df[normal_df.Key == key_name]['Value'].dropna()

    # Pairwise comparisons (Mann-Whitney U test)
    groups = {
        'Seizure vs Preepileptic': (seizure_data, preepileptic_data),
        'Seizure vs Normal': (seizure_data, normal_data),
        'Preepileptic vs Normal': (preepileptic_data, normal_data)
    }


    results = []
    for comparison, (group1, group2) in groups.items():
        # Step 1: Estimate PDFs using Kernel Density Estimation (KDE)
        kde1 = KernelDensity(kernel='gaussian', bandwidth=0.2).fit(np.array(group1).reshape(-1, 1))
        kde2 = KernelDensity(kernel='gaussian', bandwidth=0.2).fit(np.array(group2).reshape(-1, 1))

        # Create a range of points to evaluate the PDFs
        x = np.linspace(min(group1.min(), group2.min()), max(group1.max(), group2.max()), 1000)

        # Evaluate PDFs at points in range
        p = np.exp(kde1.score_samples(x[:, None]))  # PDF of sequence1
        q = np.exp(kde2.score_samples(x[:, None]))  # PDF of sequence2

        # Normalize PDFs to ensure they sum to 1
        p /= np.sum(p)
        q /= np.sum(q)

        # Step 2: Compute KL divergence
        kl_divergence = entropy(p, q)  # Computes D_KL(P || Q)

        results.append((comparison, kl_divergence))

    # Print results with adjusted p-values
    
    return results



def evaluate_mutual_information(seizure_reduced, preepileptic_reduced, normal_reduced, key_name):
    from sklearn.feature_selection import mutual_info_classif


    seizure_data = seizure_reduced[seizure_reduced.Key == key_name]['Value'].dropna()
    preepileptic_data = preepileptic_reduced[preepileptic_reduced.Key == key_name]['Value'].dropna()
    normal_data = normal_reduced[normal_reduced.Key == key_name]['Value'].dropna()
    groups = {
        'Seizure vs Normal': (seizure_data, normal_data),
        'Preepileptic vs Normal': (preepileptic_data, normal_data)
    }

    results = []
   
    for comparison, (group1, group2) in groups.items():
        # Convert groups to arrays
        array_group1 = np.array(group1)
        array_group2 = np.array(group2)

        # Create the feature and target arrays
        features = np.concatenate([array_group1, array_group2]).reshape(-1, 1)
        predictions = np.concatenate([
            np.ones(len(group1)),  # Label group1 as 1
            np.zeros(len(group2)) # Label group2 as 0
        ])

        # Compute mutual information
        mi = mutual_info_classif(features, predictions, discrete_features=False)

        # Store results
        results.append((comparison, mi[0]))

    return results


seizures = load_dataset_to_data_frame('/data1/seizures_data.csv' )
preepileptic = load_dataset_to_data_frame('/data1/preepileptic_data.csv' )
normal = load_dataset_to_data_frame('/data1/normal_data.1.csv' )

## important: remove code inner block in the future. 
seizures = seizures[seizures['id'] < 1000]
preepileptic = preepileptic[preepileptic['id'] < 1000]
normal = normal[normal['id'] < 1000]
###


print(f'1. Reduce span entropy records to avergae')
seizure_reduced = reduce_span_entropy_to_average(seizures)
normal_reduced = reduce_span_entropy_to_average(normal)
preepileptic_reduced = reduce_span_entropy_to_average(preepileptic)

print(f'2. Evaluate mutual information and K-L divergency.')
results = {}
profiled_items = set(
    seizure_reduced['Key'].unique()
).intersection(
    preepileptic_reduced['Key'].unique(),
    normal_reduced['Key'].unique()
)

for profiled_item in tqdm(profiled_items):
    print(f'evalue Key = {profiled_item}')
    results[profiled_item] = {}
    mi_comparision = evaluate_mutual_information(seizure_reduced, preepileptic_reduced, normal_reduced, profiled_item)
    results[profiled_item]["mi_comparision"] = mi_comparision
    kl_comparision = pairwise_comparisons_KL(seizure_reduced, preepileptic_reduced, normal_reduced, profiled_item)
    results[profiled_item]["kl_comparision"] = kl_comparision

np.save("feature_evaluation_results.dict", results, allow_pickle=True)