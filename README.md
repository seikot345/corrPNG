# corrPNG : A tool for calculating the correlation coefficient between <ins>P</ins>henotype data and the <ins>N</ins>umber of <ins>G</ins>enes. 

```
#Calculate pearson correlation coefficient 
python corrPNG.py pearson -i genePAV.Rtab -p pheno -o output -t 0.7  

#Calculate spearman rank correlation coefficient 
python corrPNG.py spearman -i genePAV.Rtab -p pheno -o output -t 0.7  

#Calculate kendall rank correlation coefficient 
python corrPNG.py kendall -i genePAV.Rtab -p pheno -o output -t 0.7  

#Visualization correlation diagram
python corrPNG.py plot -i genePAV.Rtab -p pheno.tsv -o output -g gene_id -n phenotype_id -r

#Change the column notation of tsv file 
python corrPNG.py sort_column -i genePAV.Rtab -o output_sort -s samplelist

#Change row and column
python corrPNG.py transpose_table -i input.tsv -o output.tsv
```

# Requirements
The following software must be installed on your machine:  
Python : tested with version 3.8

# introduction 
corrPNG.py calculate correlation coefficient between the phenotype data and copy number variation of gene including the gene presence/absence. Run the following command in [pangene][pangene] and create the tab-delimited file. 
```
k8 pangene.js gfa2matrix graph.gfa > genePAV.Rtab
k8 pangene.js gfa2matrix -c graph.gfa > geneCNV.Rtab
```

# Calculate correlation coefficient
Use obtained tsv file for corrPNG.py as input. Phenotype data is averrable csv and tsv file. Set the data as that row is phenotype name and column is sample name. corrPNG.py has three commands for calculate correlation coefficient:
**pearson** for pearson’s coefficient (r),
**spiaman** for spiarman rho (ρ), and 
**kendoll** for kendoll tau (τ).
And output the file coefficient and p value over the threshold. Threshold should be one value (set 0.7 in default) and it works as positive and negative (e.g. 0.7 and -0.7).  
<img width="750" alt="fig1" src=https://github.com/user-attachments/assets/66988a6a-fb45-4055-9d14-8c1e00702b36>

# Visualization
**plot** generate correlation diagram (only scatter plot) between one gene and one phenotype. Set gene name as -g and phenotype name as -n. If you add -r option, write regression line. If you add -jx and -jy option, generate jitter plot.  
<img width="750" alt="fig2" src=https://github.com/user-attachments/assets/7ccaed08-e3ca-4130-8deb-a0f20c756b0b>

# Data formatting
**sort_column** is command allows you to arbitrarily change the column notation of genePAV.tsv obtained by gff2matrix in [pangene][pangene].  
**transpose_table** is change row and column of csv or tsv file.
Please include the extension for the options of these two command.  
<img width="500" alt="fig3" src=https://github.com/user-attachments/assets/813bab6d-44c7-4570-8829-8f943de699cb>
<img width="400" alt="fig4" src=https://github.com/user-attachments/assets/2faeaca0-8769-4fbf-a204-9f800d43c899>


# Citation
xxx

[pangene]: https://github.com/lh3/pangene
