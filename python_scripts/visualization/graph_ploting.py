
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
    
def plot_distribution_of_record_items(key_name):
    
    plt.figure(figsize=(10, 6))

    # Create a seaborn distribution plot
    #sns.histplot(seizures[seizures.Key == 'derivation_entropy']['Value'], kde=True, bins=30)  # Adjust bins as needed
    sns.kdeplot(seizures[seizures.Key == key_name]['Value'], color="red", label="KDE")
    sns.kdeplot(preepileptic[preepileptic.Key == key_name]['Value'], color="blue", label="KDE")
    sns.kdeplot(normal[normal.Key == key_name]['Value'], color="green", label="KDE")

    # Customize the plot
    plt.title(f'Distribution of {key_name}')
    plt.xlabel(key_name)
    plt.ylabel('Frequency')
    plt.legend(["Seizure", "Pre-epileptic (length = 5 minutes)", "Normal"])
    
def box_plot_record_items(key_name):
    
    plt.figure(figsize=(10, 6))

    seizure_data = seizures[seizures.Key == key_name]['Value']
    preepileptic_data = preepileptic[preepileptic.Key == key_name]['Value']
    normal_data = normal[normal.Key == key_name]['Value']
    data = pd.DataFrame({
        'Distribution': ['Seizure'] * len(seizure_data) + ['Pre-epileptic (length = 5 minutes)']*len(preepileptic_data) + ['Normal']*len(normal_data),
        'Value': list(seizure_data) + list(preepileptic_data) + list(normal_data)
    })
    sns.boxplot(x='Distribution', y='Value', data=data)

    plt.title(f'Box plot of {key_name}')
    plt.xlabel(key_name)
    plt.ylabel('Value')
    #plt.legend(["Seizure", "Pre-epileptic (length = 5 minutes)", "Normal"])

def plot_all_graph(key_name):
    plot_distribution_of_record_items(key_name)
    box_plot_record_items(key_name)