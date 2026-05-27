import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

def create_dot_plot(data, title, output_file=None):
    """Generate a dot plot with an additional shape factor."""
    plt.figure(figsize=(16, 9))
    scatter = sns.scatterplot(
        data=data,
        x='Amino Acid Position',
        y='Percentile Rank',
        hue='MHC Allele',
        style='Host',  # Use the new shape factor
        palette='viridis',
        s=100,
        marker='o'
    )
    plt.title(title, fontsize=14)
    plt.xlabel('Amino acid Position', fontsize=12)
    plt.ylabel('Percentile Rank', fontsize=12)
    
    # Retrieve handles and labels from the legend
    handles, labels = scatter.get_legend_handles_labels()
    plt.legend(handles=handles, labels=labels, 
               title='Legend', bbox_to_anchor=(1, 1), loc='upper left', ncol=2)
    
    plt.grid(True)
    plt.tight_layout()
    
    # Save the plot or show it
    if output_file:
        plt.savefig(output_file)
        plt.close()
    else:
        plt.show()

def process_sheets(input_file):
    """Process each sheet in the Excel file and generate plots ."""
    excel_data = pd.ExcelFile(input_file)
    columns_of_interest = ['allele', 'start', 'consensus_percentile_rank', 'shape_column']
    
    for sheet in excel_data.sheet_names:
        sheet_data = excel_data.parse(sheet)
        sheet_data.rename(columns=lambda x: x.strip(), inplace=True)
        
        if set(columns_of_interest).issubset(sheet_data.columns):
            cleaned_data = sheet_data[columns_of_interest].dropna()
            cleaned_data.rename(columns={
                'allele': 'MHC Allele',
                'start': 'Amino Acid Position',
                'consensus_percentile_rank': 'Percentile Rank',
                'shape_column': 'Host'  # Rename Column E to 'Shape'
            }, inplace=True)
            
            # Generate the dot plot 
            output_file = f"{sheet}_dot_plot.png"
            create_dot_plot(cleaned_data, f'{sheet}: Dot Plot ', output_file)
            print(f"Saved plot to {output_file}")
        else:
            print(f"Sheet {sheet} does not contain all required columns.")

def main():
    parser = argparse.ArgumentParser(description="Generate dot plots from an Excel file.")
    parser.add_argument("input_file", help="Path to the input Excel file.")
    args = parser.parse_args()

    process_sheets(args.input_file)

if __name__ == "__main__":
    main()
