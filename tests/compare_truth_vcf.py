from cyvcf2 import VCF
import argparse
from sklearn.metrics import confusion_matrix, classification_report, precision_score, recall_score, f1_score
import pandas as pd
import plotly.express as px
import matplotlib.pyplot as plt
import seaborn as sns

def parse_vcf(vcf_file, key):
    """Parses the VCF file and extracts the required INFO field values using cyvcf2."""
    data = {}
    
    # Open the VCF file using cyvcf2
    vcf = VCF(vcf_file)

    # Iterate over each variant in the VCF file
    for record in vcf:
        pos = record.POS  # Position
        list_list_genotypes = record.genotypes  # Genotypes for each sample
        
        # If the value is a list, convert it to float if necessary
        if key == 'FR':
            info_value = record.INFO.get(key).split(',')
            info_value = [float(x) if isinstance(x, (int, float, str)) else x for x in info_value]
            iterator = -2
        elif key == 'PV':
            info_value = [record.INFO.get(key)]
            iterator = -1
        
        # Calculate the number of alleles from the genotypes
        num_alleles = 0
        for list_genotypes in list_list_genotypes:
            for genotype in list_genotypes[:iterator]:
                num_alleles += 1 if int(genotype) != -1 else 0

        # Store the extracted data (INFO field value and number of alleles)
        data[pos] = (info_value, num_alleles)
    
    return data

def evaluate_conditions(fr_data, pv_data):
    """Evaluates the conditions and prepares data for the confusion matrix."""
    y_true = []  # Actual values (0 or 1 based on PV condition)
    y_pred = []  # Predicted values (0 if FR condition is not met, 1 if met)
    num_sample = []
    pvalues = []
    list_diff = []

    for key in fr_data:
        if key in pv_data:
            fr_values = fr_data[key][0]  # List of values expected
            pv_value = pv_data[key][0][0]  # Single value expected
            num_sample.append(pv_data[key][1])

            fr_condition = any(f == 0.0 for f in fr_values)
            list_diff.append(max(fr_values[0], fr_values[1]))
            pv_condition = pv_value > 0.01
            pvalues.append(pv_value)

            y_pred.append(int(fr_condition))
            y_true.append(int(pv_condition))
    
    return y_true, y_pred, num_sample, list_diff, pvalues

def p_value_distribution(y_true, y_pred, num_sample, list_diff, p_value, output="plink"):
    
    false_positive_indices = [
        i for i, (true, pred) in enumerate(zip(y_true, y_pred)) if pred == 1 and true == 0]

    print("len(false_positive_indices) : " , len(false_positive_indices))

    true_positive_indices = [
        i for i, (true, pred) in enumerate(zip(y_true, y_pred)) if pred == 0 and true == 0]
    
    print("len(true_positive_indices) : " , len(true_positive_indices))

    false_negative_indices = [
       i for i, (true, pred) in enumerate(zip(y_true, y_pred)) if pred == 0 and true == 1]

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
        labels={"P-Value": "P-Value", "Difference": "Simulated Effect (Difference in Probabilities)"}
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

def conf_mat_maker(cm, output="plink"):

    # Plot confusion matrix for p-value < p_val
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', cbar=False,
                xticklabels=['Positive', 'Negative'], 
                yticklabels=['Positive', 'Negative'],
                annot_kws={"size": 30})
    plt.xticks(fontsize=16)  
    plt.yticks(fontsize=16)  
    plt.title(f'Confusion Matrix for p-value < 0.01', fontsize=18)  # Increase title font size
    plt.xlabel('Truth Labels', fontsize=20)  # Increase x-label font size
    plt.ylabel('Predicted Labels', fontsize=20)  # Increase y-label font size
    plt.savefig(output + f'_0.01.png', format='png', dpi=300)

def main(vcf1, vcf2):
    fr_data = parse_vcf(vcf1, 'FR')
    pv_data = parse_vcf(vcf2, 'PV')
    
    y_true, y_pred, num_sample, list_diff, pvalues = evaluate_conditions(fr_data, pv_data)
    
    cm = confusion_matrix(y_true, y_pred)
    precision = precision_score(y_true, y_pred)
    recall = recall_score(y_true, y_pred)
    f1 = f1_score(y_true, y_pred)

    print("Confusion Matrix:")
    conf_mat_maker(cm)
    print("\nClassification Report:")
    print(classification_report(y_true, y_pred))
    print(f"\nPrecision: {precision:.3f}")
    print(f"Recall: {recall:.3f}")
    print(f"F1-score: {f1:.3f}")
    p_value_distribution(y_true, y_pred, num_sample, list_diff, pvalues)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process two VCF files and compute a confusion matrix.")
    parser.add_argument("vcf_truth", help="VCF file containing FR values")
    parser.add_argument("vcf_plink", help="VCF file containing PV values")
    args = parser.parse_args()
    
    main(args.vcf_truth, args.vcf_plink)

# python3 tests/compare_truth_vcf.py tests/quantitative.truth.deconstruct.norm.vcf tests/quantitative.plink_vcf.vcf
