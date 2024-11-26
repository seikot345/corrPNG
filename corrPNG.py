import argparse
import csv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as st
import seaborn as sns
import sys
import os

# Calculate correlation coefficient
def read_file(file_path):
    # tsv or csv
    _, ext = os.path.splitext(file_path)
    delimiter = '\t' if ext == '.tsv' or '.Rtab' else ','

    with open(file_path, newline='', encoding='utf-8') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=delimiter)
        header = next(csvreader)
        return [row for row in csvreader]

def pearson(input_file, pheno, output_file, threshold):
    print(f"Executing spearman with option: -i {input_file} -p {pheno} -o {output_file} -t {threshold}")
    output_file = output_file + '.txt'
    with open(output_file, 'w') as file:
        pheno_data = read_file(pheno)
        gene_data = read_file(input_file)
        threshold = float(threshold)
        threshold2 = -1*threshold

        # write header
        file.write("pheno\tgene\tcoef\tpvalue\n")

        for row1 in pheno_data:
            for row2 in gene_data:
                X = np.array(row1[1:], dtype=np.float64)
                Y = np.array(row2[1:], dtype=np.float64)

                # caluculate pearson correlation coefficient
                pec, pep = st.pearsonr(X, Y)

                # write result if over threshold
                if pec >= threshold or pec <= threshold2:
                    file.write(f"{row1[0]}\t{row2[0]}\t{pec}\t{pep}\n")

def spearman(input_file, pheno, output_file, threshold):
    print(f"Executing spearman with option: -i {input_file} -p {pheno} -o {output_file} -t {threshold}")
    output_file = output_file + '.txt'
    with open(output_file, 'w') as file:
        pheno_data = read_file(pheno)
        gene_data = read_file(input_file)
        threshold = float(threshold)
        threshold2 = -1*threshold

        # write header
        file.write("pheno\tgene\trho\tpvalue\n")

        for row1 in pheno_data:
            for row2 in gene_data:
                X = np.array(row1[1:], dtype=np.float64)
                Y = np.array(row2[1:], dtype=np.float64)

                # caluculate spearman rank correlation coefficient
                spr, spp = st.spearmanr(X, Y)

                # write result if over threshold
                if spr >= threshold or spr <= threshold2:
                    file.write(f"{row1[0]}\t{row2[0]}\t{spr}\t{spp}\n")

def kendall(input_file, pheno, output_file, threshold):
    print(f"Executing kendall with option: -i {input_file} -p {pheno} -o {output_file} -t {threshold}")
    output_file = output_file + '.txt'
    with open(output_file, 'w') as file:
        pheno_data = read_file(pheno)
        gene_data = read_file(input_file)
        threshold = float(threshold)
        threshold2 = -1*threshold

        # write header
        file.write("pheno\tgene\ttau\tpvalue\n")

        for row1 in pheno_data:
            for row2 in gene_data:
                X = np.array(row1[1:], dtype=np.float64)
                Y = np.array(row2[1:], dtype=np.float64)

                # caluculate spearman rank correlation coefficient
                ket, kep = st.kendalltau(X, Y)

                # write result if over threshold
                if ket >= threshold or ket <= threshold2:
                    file.write(f"{row1[0]}\t{row2[0]}\t{ket}\t{kep}\n")

# Visualization
def determine_delimiter(file_path):
    with open(file_path, 'r') as f:
        first_line = f.readline()

    if first_line.count('\t') > first_line.count(','):
        return '\t'
    else:
        return ','

def plot(input_file, pheno, output_file, gene_id, pheno_id, r, jx, jy):
    print(f"generating the correlation diagram with option: -i {input_file} -p {pheno} -o {output_file} -g {gene_id} -n {pheno_id}")

    sep_g = determine_delimiter(input_file)
    sep_p = determine_delimiter(pheno)

    inputA = pd.read_csv(input_file, sep=sep_g, index_col=0)
    inputB = pd.read_csv(pheno, sep=sep_p, index_col=0)

    if gene_id not in inputA.index:
        raise ValueError(f"{gene_id} not found in {input_file}")
    if pheno_id not in inputB.index:
        raise ValueError(f"{pheno_id} not found in {pheno}")

    x_values = inputA.loc[gene_id]
    y_values = inputB.loc[pheno_id]

    # Standardization of sample names
    x_values.index = x_values.index.str.lower().str.replace(" ", "")
    y_values.index = y_values.index.str.lower().str.replace(" ", "")

    # Only common samples are taken
    common_samples = x_values.index.intersection(y_values.index)
    x_values = x_values[common_samples].astype(float)
    y_values = y_values[common_samples].astype(float)

    if x_values.empty or y_values.empty:
        raise ValueError(f"No common samples found between {input_file} and {pheno}.")

    plt.figure(figsize=(8, 6))

    # option for jitter plot
    jit_x = float(x_values.max())/50
    jit_y = float(y_values.max())/50
    jittered_x = x_values + np.random.normal(0, jit_x, size=len(x_values))
    jittered_y = y_values + np.random.normal(0, jit_y, size=len(y_values))

    if r == True:
        if jx == True:
            if jy == True:
                sns.regplot(x=jittered_x, y=jittered_y, scatter_kws={'s': 10}, line_kws={'color': 'red'})
            else:
                sns.regplot(x=jittered_x, y=y_values, scatter_kws={'s': 10}, line_kws={'color': 'red'})
        elif jy == True:
            sns.regplot(x=x_values, y=jittered_y, scatter_kws={'s': 10}, line_kws={'color': 'red'})
        else:
            sns.regplot(x=x_values, y=y_values, scatter_kws={'s': 10}, line_kws={'color': 'red'})
    elif jx == True:
        if jy == True:
            sns.regplot(x=jittered_x, y=jittered_y, scatter_kws={'s': 10}, line_kws={'color': 'red'}, fit_reg=False)
        else:
            sns.regplot(x=jittered_x, y=y_values, scatter_kws={'s': 10}, line_kws={'color': 'red'}, fit_reg=False)
    else:
        if jy == True:
            sns.regplot(x=x_values, y=jittered_y, scatter_kws={'s': 10}, line_kws={'color': 'red'}, fit_reg=False)
        else:
            sns.regplot(x=x_values, y=y_values, scatter_kws={'s': 10}, line_kws={'color': 'red'}, fit_reg=False)

    plt.xlabel(f"Number of homologs for {gene_id}")
    plt.ylabel(f"Value of {pheno_id}")
    plt.title(f"Correlation diagram between {gene_id} and {pheno_id}")

    # Set x-axis to integer values
    max_x = int(np.ceil(x_values.max()))  # Adjust the maximum value of the x-axis to the maximum value of the data
    plt.xticks(range(0, max_x + 1))  # Set the scale to an integer between 0 and the maximum value

    output_file = output_file + '.png'
    plt.savefig(output_file, format='png')
    plt.close()

# data formatting
def sort_column(input_file, output_file, samplelist):
    print(f"Executing sort_column with option: -i {input_file} -s {samplelist} -o {output_file}")

    sep = determine_delimiter(input_file)

    # load sort order
    with open(samplelist, 'r') as f:
        order = [line.strip() for line in f.readlines()]
    order.insert(0, 'Gene')

    data = pd.read_csv(input_file, sep=sep)

    # Sorts columns in a specified order
    ordered_data = data[order]

    # Save the sorted data
    if output_file.endswith('.tsv'):
        output_sep = '\t'
    elif output_file.endswith('.csv'):
        output_sep = ','
    else:
        output_sep = sep
    ordered_data.to_csv(output_file, sep=output_sep, index=False)

def transpose_table(input_file, output_file):
    print(f"Executing transpose_table with option: -i {input_file} -o {output_file}")

    sep = determine_delimiter(input_file)

    data = pd.read_csv(input_file, sep=sep)

    # Extract the original header
    original_header = data.columns.tolist()

    transposed_data = data.transpose()

    # Add the original header row as the first column
    transposed_data.insert(0, 'Original Header', original_header)

    if output_file.endswith('.tsv'):
        output_sep = '\t'
    elif output_file.endswith('.csv'):
        output_sep = ','
    else:
        output_sep = sep

    transposed_data.to_csv(output_file, sep=output_sep, index=False, header=False)

class CustomHelpFormatter(argparse.RawTextHelpFormatter):
    def add_arguments(self, actions):
        # Override the `add_arguments` method to display nothing
        pass

def main():
    parser = argparse.ArgumentParser(add_help=False)
    subparsers = parser.add_subparsers(dest="command")

    parser_1 = subparsers.add_parser("pearson", help="pearson command")
    parser_1.add_argument("-i", "--input", required=True, metavar="STR", help="Name of input tsv file")
    parser_1.add_argument("-p", "--pheno", required=True, metavar="STR", help="Name of input phenotype file")
    parser_1.add_argument("-o", "--output", default="output", metavar="STR", help="Name of output file")
    parser_1.add_argument("-t", "--threshold", default="0.7", metavar="FLOAT", help="threshold")

    parser_2 = subparsers.add_parser("spearman", help="spearman command")
    parser_2.add_argument("-i", "--input", required=True, metavar="STR", help="Name of input tsv file")
    parser_2.add_argument("-p", "--pheno", required=True, metavar="STR", help="Name of input phenotype file")
    parser_2.add_argument("-o", "--output", default="output", metavar="STR", help="Name of output file")
    parser_2.add_argument("-t", "--threshold", default="0.7", metavar="FLOAT", help="threshold")

    parser_3 = subparsers.add_parser("kendall", help="kendall command")
    parser_3.add_argument("-i", "--input", required=True, metavar="STR", help="Name of input tsv file")
    parser_3.add_argument("-p", "--pheno", required=True, metavar="STR", help="Name of input phenotype file")
    parser_3.add_argument("-o", "--output", default="output", metavar="STR", help="Name of output file")
    parser_3.add_argument("-t", "--threshold", default="0.7", metavar="FLOAT", help="threshold")

    parser_4 = subparsers.add_parser("plot", help="plot command")
    parser_4.add_argument("-i", "--input", required=True, metavar="STR", help="Name of input tsv file")
    parser_4.add_argument("-p", "--pheno", required=True, metavar="STR", help="Name of input phenotype file")
    parser_4.add_argument("-o", "--output", default="output", metavar="STR", help="Name of output file")
    parser_4.add_argument("-g", "--gene_id", required=True, metavar="STR", help="Gene name")
    parser_4.add_argument("-n", "--pheno_id", required=True, metavar="STR", help="Phenotype name")
    parser_4.add_argument("-r", action="store_true", help="Plot the regression line")
    parser_4.add_argument("-jx", action="store_true", help="Jitter plot for x-axis")
    parser_4.add_argument("-jy", action="store_true", help="Jitter plot for y-axis")

    parser_5 = subparsers.add_parser("sort_column", help="sort_column command")
    parser_5.add_argument("-i", "--input", required=True, metavar="STR", help="Name of input tsv or csv file")
    parser_5.add_argument("-o", "--output", default="output.sorted.tsv", metavar="STR", help="Name of output file")
    parser_5.add_argument("-s", "--samplelist", required=True, metavar="STR", help="Name of input sample list file")

    parser_6 = subparsers.add_parser("transpose_table", help="transpose_table command")
    parser_6.add_argument("-i", "--input", required=True, metavar="STR", help="Name of input tsv or csv file")
    parser_6.add_argument("-o", "--output", default="output.sorted.tsv", metavar="STR", help="Name of output file")

    if len(sys.argv) < 2 or sys.argv[1] in ["-h", "--help"]:
        print("A tool for calculating the correlation coefficient between Phenotype data and the Number of Genes.\n\n\
Version: 1.0\n\n\
Usage:\n\
    {} [command] option\n\n\
Command list:\n\
    pearson          Calculate pearson correlation coefficient\n\
    spearman         Calculate spearman rank correlation coefficient\n\
    kendall          Calculate kendall rank correlation coefficient\n\
    plot             Visualize correlation diagram between one gene and one phenotype\n\
    sort_column      Change the column notation of tsv file\n\
    transpose_table  Transpose the table\n".format(sys.argv[0]))
        sys.exit(0)

    args = parser.parse_args()

    if args.command == "pearson":
        pearson(args.input, args.pheno, args.output, args.threshold)
    elif args.command == "spearman":
        spearman(args.input, args.pheno, args.output, args.threshold)
    elif args.command == "kendall":
        kendall(args.input, args.pheno, args.output, args.threshold)
    elif args.command == "plot":
        if args.r:
            if args.jx:
                if args.jy:
                    plot(args.input, args.pheno, args.output, args.gene_id, args.pheno_id, r=True, jx=True, jy=True)
                else:
                    plot(args.input, args.pheno, args.output, args.gene_id, args.pheno_id, r=True, jx=True, jy=False)
            elif args.jy:
                plot(args.input, args.pheno, args.output, args.gene_id, args.pheno_id, r=True, jx=False, jy=True)
            else:
                plot(args.input, args.pheno, args.output, args.gene_id, args.pheno_id, r=True, jx=False, jy=False)
        elif args.jx:
            if args.jy:
                plot(args.input, args.pheno, args.output, args.gene_id, args.pheno_id, r=False, jx=True, jy=True)
            else:
                plot(args.input, args.pheno, args.output, args.gene_id, args.pheno_id, r=False, jx=True, jy=False)
        else:
            if args.jy:
                plot(args.input, args.pheno, args.output, args.gene_id, args.pheno_id, r=False, jx=False, jy=True)
            else:
                plot(args.input, args.pheno, args.output, args.gene_id, args.pheno_id, r=False, jx=False, jy=False)
    elif args.command == "sort_column":
        sort_column(args.input, args.output, args.samplelist)
    elif args.command == "transpose_table":
        transpose_table(args.input, args.output)
    else:
        print("Error: {} [command] \n\n    [pearson] [spearman] [kendall] [plot] [sort_column] [transpose_table]\n".format(sys.argv[0]))
        sys.exit(1)

if __name__ == "__main__":
    main()
