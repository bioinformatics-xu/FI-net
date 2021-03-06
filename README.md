# FI-net (R project)
FI-net: identification of cancer driver genes by using functional impact prediction neural network

Developer: Xiaolu Xu lu891111@mail.dlut.edu.cn

Faculty of Electronic Information and Electrical Engineering 

Dalian University of Technology, China
# Instructions to FI-net
## Dependencies
Some R package should be imported to apply FI-net, including:
- data.table
- plyr
- stringr
- h2o
- fitdistrplus
## Identify driver genes based on FI-net
Identify driver genes based no an artificial neural network model and hierarchical clustering algorithm using ./FInet.R. 
### 1. Run FI-net
Run ./FInet.R to identify driver genes, and the gene list and the corresponding qValue will be returned.
### 2. Input
- input.matrix: The input matrix contains the gene feature and the functional impact score.
- N: The pre-defined expected number of genes in each cluster.
### 3. Output
- driver gene list.
## Folder Structure
```
FI-net
|__ README
|__ FInet.R
|__ input.matrix.csv
```
## Usage
1. Make sure that the dependent R package (data.table, plyr, stringr, h2o, and fitdistrplus) are successfully installed.

2. Run the following R program:
```
library(data.table)
library(stringr)
library(plyr)
library(h2o)
source("FInet.R")

input.matrix <- read.csv("input.matrix.csv")
result <- FInet(input.matrix,3000)
driver.genes <- result$gene[which(result$qValue <= 0.05)]
```
