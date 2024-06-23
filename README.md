# MNMO
Discover driver genes from a multi-omics data based multi-layer network model
## Installation

You can install the development version of MNMO like so:

```{r, eval=FALSE}
library(psych)
install.packages("/path_to_package/MNMO_xxx.tar.gz", repos = NULL)
library(MNMO)
```

## Example

This is a basic example which shows you how to get a list of potential driver genes for breast cancer:

```{r, eval=FALSE}
library(psych)
install.packages("/path_to_package/MNMO_xxx.tar.gz", repos = NULL)
library(MNMO)
```
Now, you can load data about breast cancer:
```{r, eval=FALSE, message=FALSE, warning = FALSE}
data(network)
data(protein_info)
data(mir_gene_network)
data(brca_normal_s)
data(brca_cancer_s)
data(brca_normal_si)
data(brca_cancer_si)
data(brca_G_mut)
```

Of course, you can also read your own data by specifying the path:
```{r, eval=FALSE, message=FALSE, warning = FALSE}
brca_normal_s <- read.table('normal_sample.txt',sep='\t',TRUE)
brca_cancer_s <- read.table('cancer_sample.txt',sep='\t',TRUE)
...
```
Then we use the following function to extract differentially expressed genes and miRNAs:

```{r, eval=FALSE, message=FALSE, warning = FALSE}
library(MNMO)
degs_05 <- cal_degs(brca_normal_s,brca_cancer_s,0.5)
mirna_value <- cal_degs(brca_normal_si,brca_cancer_si,0.5)
```
Now, we used multi-omics data to extract multi-layer networks and calculate the control ability scores of genes:

```{r, eval=FALSE, message=FALSE, warning = FALSE}
library(MNMO)
construct_net(brca_G_mut,network,protein_info,degs_05)
```
After run the construct_net function,  multiple output files will get, including the second layer network and the control ability score of the gene.
Read the obtained files into memory:

```{r, eval=FALSE, message=FALSE, warning = FALSE}
layer2_gene_list <- read.csv("layer2_gene_list.csv")
layer2_gene_list <- layer2_gene_list[,1]
layer2_network <- read.csv("layer2_network.csv")
gene_expscore_2 <-read.csv("gene_expscore_2.csv")
gene_expscore_3 <-read.csv("gene_expscore_3.csv")
```

Finally, the following function is used to calculate the final score of the gene and rank the potential driver genes:
```{r, eval=FALSE, message=FALSE, warning = FALSE}
library(MNMO)
brca_list <- MNMO(layer2_gene_list,layer2_network,brca_cancer_s,brca_G_mut,
                  mirna_value,mir_gene_network,gene_expscore_2,gene_expscore_3)
```
