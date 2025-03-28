import argparse
import pandas as pd
from sklearn.metrics import confusion_matrix, precision_score, recall_score, f1_score
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import re
from cyvcf2 import VCF # type: ignore

def split_snarl(input_str):
    # Split the string and filter empty elements, then convert to integers
    return [str(num) for num in re.split(r'[><]', input_str) if num]

def parse_snarl_path_file_dict(path_file:str) -> dict :
    
    snarl_paths = {}
    df = pd.read_csv(path_file, sep='\t', dtype=str)
    df['paths'] = df['paths'].str.split(',')
    for paths, pos in zip(df['paths'], df['pos']):
        for path in paths :
            list_node = split_snarl(path)
            # all case : SNP, INS, DEL, COMPLEX
            node = int(list_node[1])
            snarl_paths[node] = pos

    return snarl_paths

# Function to process the frequency file and get result list with differences
def process_file(freq_file, threshold=0.2):
    df = pd.read_csv(freq_file, sep='\t')
    freq_path_list = []
    true_labels = []
    list_diff = []
 
    # Iterate in pairs and calculate the differences
    for i in range(0, len(df) - 1, 2):
        row_1 = df.iloc[i]
        row_2 = df.iloc[i + 1]
        diff = abs(row_1['freq'] - row_2['freq'])

        # If group freq diff > threshold, then true label is 0, else 1
        if diff > threshold:
            true_labels.append(0)
        else :
            true_labels.append(1)

        freq_path_list.append(f"{int(row_1['start_node'])}_{int(row_1['next_node'])}")
        list_diff.append(float(diff))

    return freq_path_list, true_labels, list_diff

def check_valid_snarl(start_node_1, next_node_1, start_node_2, next_node_2, snarl_list):
    """
    Check if a valid snarl exists in the snarl_list.

    A valid snarl is defined as:
    - At least one path contains both start_node_1 and next_node_1.
    - At least one path contains both start_node_2 and next_node_2.
    - start_node_1 and start_node_2 must be the same str.
    - One path can not containt the 4 node.
    """

    if start_node_1 != start_node_2:  # First condition: start nodes must match
        raise ValueError("start_node_1 != start_node_2")

    contains_first_pair = False
    contains_second_pair = False

    for path in snarl_list:
        # Extract numbers from the path as integers
        decomposed_path = {int(num) for num in re.findall(r'\d+', path)}
        
        # Check if the current path contains both nodes of the first pair
        if int(start_node_1) in decomposed_path and int(next_node_1) in decomposed_path:
            contains_first_pair = True

        # Check if the current path contains both nodes of the second pair
        if int(start_node_2) in decomposed_path and int(next_node_2) in decomposed_path:
            contains_second_pair = True

        # If both conditions are satisfied, we can exit early
        if contains_first_pair and contains_second_pair:
            return True

    # Return true only if both pairs are satisfied
    return contains_first_pair and contains_second_pair

def match_snarl(freq_path_list, true_labels, list_diff, p_value_file, paths_file):

    p_value_df = pd.read_csv(p_value_file, sep='\t')
    paths_df = pd.read_csv(paths_file, sep='\t')['paths']
    split = p_value_df['SNARL'].str.split('_')

    # To store predicted labels
    predicted_labels_10_2 = []
    predicted_labels_10_5 = []
    predicted_labels_10_8 = []
    cleaned_true_labels = []
    clean_list_diff = []
    pvalue = []
    num_sample =[]

    for idx in range(0, len(freq_path_list) - 1, 2):  # Step by 2 to process pairs

        start_node_1, next_node_1 = map(int, freq_path_list[idx].split('_'))
        start_node_2, next_node_2 = map(int, freq_path_list[idx+1].split('_'))

        # We want to know if the snarl is in the range/containt of the snarl in the p_value file
        matched_row = p_value_df[(split.str[1].astype(int) <= start_node_1) & (split.str[0].astype(int) >= next_node_1) |
                                 (split.str[0].astype(int) <= start_node_1) & (split.str[1].astype(int) >= next_node_1)]

        # Case where the snarl is found 
        if not matched_row.empty:
            indices = matched_row.index
            split_paths = [paths_df[idx] for idx in indices]

            # Check if at least one path in the snarl contains the start node followed by the next node
            for idx_paths, list_path in enumerate(split_paths):

                if check_valid_snarl(start_node_1, next_node_1, start_node_2, next_node_2, list_path.split(',')) : 
                    matched_p_value = matched_row.loc[indices[idx_paths]]
                    if type_ == 'binary':  
                        p_value = matched_p_value['P_FISHER']
                    elif type_ == 'quantitative':
                        p_value = matched_p_value['P']
                    else :
                        raise ValueError("type_ must be binary or quantitative")

                    try:
                        p_value = float(p_value)  # Try converting p_value to a float
                        if p_value != p_value:  # Check if it's NaN
                            raise ValueError("p_value is NaN")
                    except (ValueError, TypeError):  # Catch invalid float or type errors
                        # Continue or handle the case where p_value is not a valid float
                        continue

                    predicted_labels_10_2.append(0 if p_value < 0.01 else 1)
                    predicted_labels_10_5.append(0 if p_value < 0.00001 else 1)
                    predicted_labels_10_8.append(0 if p_value < 0.00000001 else 1)
                    cleaned_true_labels.append(true_labels[idx])
                    clean_list_diff.append(list_diff[idx])
                    pvalue.append(p_value)
                    try :
                        allele_num = matched_row.loc[indices[idx_paths]]['ALLELE_NUM']
                    except :
                        allele_num = 200
                    num_sample.append(allele_num)

    return predicted_labels_10_2, predicted_labels_10_5, predicted_labels_10_8, cleaned_true_labels, clean_list_diff, pvalue, num_sample

def conf_mat_maker(p_val, predicted_labels, true_labels, output):
        
    # Calculate confusion matrix for p-value < p_val
    print(f"\nMetrics for p-value < {p_val}:")

    # Inverse because I want the X axis to be the truth labels and Y axis to be the predicted labels
    cm = confusion_matrix(predicted_labels, true_labels)
    print(f"Confusion Matrix for p-value < {p_val}:\n{cm}")
    prec = precision_score(predicted_labels, true_labels)
    recall = recall_score(predicted_labels, true_labels)
    f1 = f1_score(predicted_labels, true_labels)
    print(f"Precision: {prec:.3f}")
    print(f"Recall: {recall:.3f}")
    print(f"F1 Score: {f1:.3f}")

    # Plot confusion matrix for p-value < p_val
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', cbar=False,
                xticklabels=['Positive', 'Negative'], 
                yticklabels=['Positive', 'Negative'],
                annot_kws={"size": 30})
    plt.xticks(fontsize=16)  
    plt.yticks(fontsize=16)  
    plt.title(f'Confusion Matrix for p-value < {p_val}', fontsize=18)  # Increase title font size
    plt.xlabel('Truth Labels', fontsize=20)  # Increase x-label font size
    plt.ylabel('Predicted Labels', fontsize=20)  # Increase y-label font size
    plt.savefig(output + f'_{p_val}.png', format='png', dpi=300)

def print_confusion_matrix(predicted_labels_10_2, predicted_labels_10_5, predicted_labels_10_8, true_labels, output):
    
    p_val_10_2 = 0.01
    p_val_10_5 = 0.00001
    p_val_10_8 = 0.00000001

    conf_mat_maker(p_val_10_2, predicted_labels_10_2, true_labels, output)
    conf_mat_maker(p_val_10_5, predicted_labels_10_5, true_labels, output)
    conf_mat_maker(p_val_10_8, predicted_labels_10_8, true_labels, output)

def p_value_distribution(test_predicted_labels, cleaned_true_labels, list_diff, p_value, num_sample, output):
    
    false_positive_indices = [
        i for i, (pred, true) in enumerate(zip(test_predicted_labels, cleaned_true_labels)) 
        if pred == 0 and true == 1]

    print("len(false_positive_indices) : " , len(false_positive_indices))

    true_positive_indices = [
        i for i, (pred, true) in enumerate(zip(test_predicted_labels, cleaned_true_labels)) 
        if pred == 0 and true == 0]
    
    print("len(true_positive_indices) : " , len(true_positive_indices))

    false_negative_indices = [
        i for i, (pred, true) in enumerate(zip(test_predicted_labels, cleaned_true_labels)) 
        if pred == 1 and true == 0]

    print("len(false_negative_indices) : " , len(false_negative_indices))

    diff_false_positive = [list_diff[i] for i in false_positive_indices]
    pvalue_false_positive = [p_value[i] for i in false_positive_indices]
    minsample_false_positive = [num_sample[i] for i in false_positive_indices]

    diff_true_positives = [list_diff[i] for i in true_positive_indices]
    pvalue_true_positives = [p_value[i] for i in true_positive_indices]
    minsample_true_positives = [num_sample[i] for i in true_positive_indices]

    diff_false_negative = [list_diff[i] for i in false_negative_indices]
    pvalue_false_negative = [p_value[i] for i in false_negative_indices]
    minsample_false_negative = [num_sample[i] for i in false_negative_indices]

    # Create a DataFrame for easy plotting
    data = {
        'P-Value': pvalue_false_positive + pvalue_true_positives + pvalue_false_negative,
        'Difference': diff_false_positive + diff_true_positives + diff_false_negative,
        'Min Sample': minsample_false_positive + minsample_true_positives + minsample_false_negative,
        'Type': ['False Positives'] * len(pvalue_false_positive) + ['True Positives'] * len(pvalue_true_positives) + ['False Negatives'] * len(pvalue_false_negative)
    }

    df = pd.DataFrame(data)

    # Create the interactive scatter plot
    fig = px.scatter(
        df, 
        x='P-Value', 
        y='Difference', 
        size='Min Sample', 
        color='Type',
        hover_name=df.index,
        hover_data={
            "P-Value": True,  # Include P-Value in hover box
            "Difference": True,  # Include Difference in hover box
            "Min Sample": True,  # Optionally hide Min Sample (size is already shown)
        },
        title="Distribution of P-Values for False Positives and True Positives",
        labels={"P-Value": "P-Value", "Difference": "Simulated Effect (Difference in Probabilities)"},
        size_max=20
    )

    fig.update_layout(
        xaxis_title="P-Value",
        yaxis_title="Simulated Effect (Difference in Probabilities)",
        legend_title="Type",
        template="plotly_white",
        xaxis_title_font=dict(size=30),  # Increase font size for x-axis title
        yaxis_title_font=dict(size=30),  # Increase font size for y-axis title
        xaxis=dict(tickfont=dict(size=25)),  # Increase font size for x-axis ticks
        yaxis=dict(tickfont=dict(size=25)),  # Increase font size for y-axis ticks
    )

    # Show the interactive plot
    fig.show()

    # Optionally, save the plot as an HTML file
    fig.write_html(f'{output}_pvalue_interactive.html')

def plot_diff_distribution(test_predicted_labels:list, cleaned_true_labels:list, clean_list_diff:list, output:str, pvalue:list):

    # True label 0 : difference freq significative
    # True label 1 : No difference freq significative
    # Predict label 0 : pvalue significative
    # Predict label 1 : No pvalue significative

    # ----------------------------- True Positive -----------------------------
    true_positive_indices= [
        i for i, (pred, true) in enumerate(zip(test_predicted_labels, cleaned_true_labels)) 
        if pred == 0 and true == 0]

    true_positive_diffs = [clean_list_diff[i]*100 for i in true_positive_indices]

    plt.figure(figsize=(10, 6))
    sns.histplot(true_positive_diffs, bins=20, kde=True, color='blue')
    plt.title("Distribution of Differences for True Positive", fontsize=16)
    plt.xlabel("Difference (%)", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.grid(False)
    plt.savefig(output + pvalue + '_distribution_true_positive.png', format='png', dpi=300)

    # ----------------------------- False Positive -----------------------------
    false_positive_indices = [
        i for i, (pred, true) in enumerate(zip(test_predicted_labels, cleaned_true_labels)) 
        if pred == 0 and true == 1]

    false_positive_diffs = [clean_list_diff[i]*100 for i in false_positive_indices]

    plt.figure(figsize=(10, 6))
    sns.histplot(false_positive_diffs, bins=20, kde=True, color='blue')
    plt.title("Distribution of Differences for False Positive", fontsize=16)
    plt.xlabel("Difference (%)", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.grid(False)
    plt.savefig(output + pvalue + '_distribution_false_positive.png', format='png', dpi=300)

    # ----------------------------- False Negative -----------------------------
    false_negative_indices = [
        i for i, (pred, true) in enumerate(zip(test_predicted_labels, cleaned_true_labels)) 
        if pred == 1 and true == 0]

    false_negative_diffs = [clean_list_diff[i]*100 for i in false_negative_indices]

    plt.figure(figsize=(10, 6))
    sns.histplot(false_negative_diffs, bins=20, kde=True, color='blue')
    plt.title("Distribution of Differences for False Negative", fontsize=16)
    plt.xlabel("Difference (%)", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.grid(False)
    plt.savefig(output + pvalue + '_distribution_false_negative.png', format='png', dpi=300)

    # ----------------------------- True Negative -----------------------------
    true_negative_indices = [
        i for i, (pred, true) in enumerate(zip(test_predicted_labels, cleaned_true_labels)) 
        if pred == 1 and true == 1]

    true_negative_diffs = [clean_list_diff[i]*100 for i in true_negative_indices]

    plt.figure(figsize=(10, 6))
    sns.histplot(true_negative_diffs, bins=20, kde=True, color='blue')
    plt.title("Distribution of Differences for True Negatives", fontsize=16)
    plt.xlabel("Difference (%)", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.grid(False)
    plt.savefig(output + pvalue + '_distribution_true_negatives.png', format='png', dpi=300)
    
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Parse a file with start_node, next_node, group, and freq columns.")
    parser.add_argument("--freq", help="Path to allele frequence file")
    parser.add_argument("--p_value", help="Path to p_value gwas output file")
    parser.add_argument("--paths", help="Path to snarl paths list file")
    parser.add_argument("-t", "--threshold", type=float, required=False, help="Threshold to define the truth label")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-b", "--binary", action='store_true', help="binary test")
    group.add_argument("-q", "--quantitative", action='store_true', help="quantitative test")

    args = parser.parse_args()

    if args.binary:
        output = "tests/binary_tests_output/"
        output_diff = f"{output}/binary_test"
        type_ = 'binary'

    elif args.quantitative:
        output = "tests/quantitative_tests_output/"
        output_diff = f"{output}/quantitative_test"
        type_ = 'quantitative'

    # THRESHOLD_FREQ : Threshold to define the truth label by comparing the frequence difference between both group
    THRESHOLD_FREQ = 0.2

    # Define Truth label from freq file
    freq_test_path_list, test_true_labels, test_list_diff = process_file(args.freq, THRESHOLD_FREQ)
    assert len(freq_test_path_list) == len(test_true_labels) == len(test_list_diff)

    test_predicted_labels_10_2, test_predicted_labels_10_5, test_predicted_labels_10_8, cleaned_true_labels, clean_list_diff, _, _ = match_snarl(freq_test_path_list, test_true_labels, test_list_diff, args.p_value, args.paths)
    
    # Plot confusion matrix
    print_confusion_matrix(test_predicted_labels_10_2, test_predicted_labels_10_5, test_predicted_labels_10_8, cleaned_true_labels, f"{output}/confusion_matrix_{THRESHOLD_FREQ}")

    if args.binary or args.quantitative :
        # THRESHOLD_FREQ = 0.0 : Case where just a difference between both group snarl is considered like Truth label
        THRESHOLD_FREQ = 0.0
        freq_test_path_list, test_true_labels, test_list_diff = process_file(args.freq, THRESHOLD_FREQ)
        test_predicted_labels_10_2, test_predicted_labels_10_5, test_predicted_labels_10_8, cleaned_true_labels, clean_list_diff, pvalue, num_sample = match_snarl(freq_test_path_list, test_true_labels, test_list_diff, args.p_value, args.paths)
        print_confusion_matrix(test_predicted_labels_10_2, test_predicted_labels_10_5, test_predicted_labels_10_8, cleaned_true_labels, f"{output}/confusion_matrix_{THRESHOLD_FREQ}")
        assert len(cleaned_true_labels) == len(clean_list_diff)

        print("Pourcentage of paths tested : ", (len(pvalue)/len(freq_test_path_list))*100*2) # *2 because we jump 2 per 2 the paths
        # Plot distribution of p-values for false negatives and true positives
        plot_diff_distribution(test_predicted_labels_10_2, cleaned_true_labels, clean_list_diff, output_diff, "10^-2")
        p_value_distribution(test_predicted_labels_10_2, cleaned_true_labels, clean_list_diff, pvalue, num_sample, output_diff)
    
    """
    python3 tests/verify_truth.py --freq tests/simulation/quantitative/pg.snarls.freq.tsv \
    --p_value output/cpp.quantitative_analysis.tsv --paths tests/simulation/quantitative/snarl_paths.tsv -q

    python3 tests/verify_truth.py --freq tests/simulation/binary/pg.snarls.freq.tsv \
    --p_value output/binary_analysis.tsv --paths tests/simulation/binary/snarl_paths.tsv -b
    """
