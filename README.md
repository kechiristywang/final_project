# Differential Gene Expression in TCGA within Stage A Chronic Lymphocytic Leukemia comparing Progression and Stable using DeSEQ2

## Download Data

* I used the data from https://dcc.icgc.org/search. Examining clinical data, there are 139 tumor samples, and 37 are defined by me as progression patients and 102 are identified as stable patients. The classification standard is the category data in the clinical data. All donor data is contained in the same file. Keep donor and seq files in all files for subsequent analysis. After the download was completed I removed extraneous information from the donor file, such as age at last follow-up. 

* The raw data I used and the data I got are in the final_project_of_Ke_Wang/scripts/Data/, but one of the raw files cannot be uploaded to Github. I used the Google Drive to share my files. [Here is the link to Google Drive.](https://drive.google.com/drive/folders/1vO7PULg_fK82M5A3ZRaQFXb1fVS_e0zy?usp=sharing)

![image](https://user-images.githubusercontent.com/89613437/144685483-feda143f-6846-4e0f-b2a7-ffbab55cc322.png)

## Step 1: Change the downloaded data format to a format suitable for DESesq2 analysis. Change gene_id to symbols.

Packages that need to be installed in advance：

install packages ggplot2 in R
```
install.packages("ggplot2")
```
install packages DESeq2, org.Hs.eg.db, apeglm, IHW, and vsn in R
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("apeglm")
BiocManager::install("IHW")
BiocManager::install("vsn")
```
install packages ashr in R
```
install.packages("devtools")
install_github("stephens999/ashr")
```
install packages pheatmap in R
```
install.packages("pheatmap")
```




(1) load Htseq-counts data in R, and keep the required data. Change data format (length data to width data).

```
rm(list = ls())
pasAnno <- read.table(file = 'seq.tsv', sep = '\t', header = TRUE) 
pasAnno1<-subset(pasAnno,select = c(1,8,10))
require(tidyverse)
dat <-  pasAnno1 %>% 
  group_by(icgc_donor_id) %>%
  mutate(index = row_number()) %>%
  pivot_wider(names_from = icgc_donor_id, 
              values_from = raw_read_count) %>%
  dplyr::select(-index)
pasAnno2<-dat
head(pasAnno)
head(pasAnno1)
head(pasAnno2)
```

* This is the raw RNA-seq data from the ICGC Data Portal.

![image](https://user-images.githubusercontent.com/89613437/144563579-f829d0bb-960c-430f-b660-5e2e143bb50a.png)

* Eliminate unwanted genes and keep Htseq-counts.

![image](https://user-images.githubusercontent.com/89613437/144564320-b2553c60-775e-4c1d-810e-a3bd55cb8970.png)

* Long format data is not conducive to subsequent analysis; turn long format data into wide data format.

![image](https://user-images.githubusercontent.com/89613437/144564556-b11f5427-c5ee-4541-b5ef-cc6e04999e33.png)

(2) Delete invalid data (Delete all NULL rows, incomplete data, and all zero rows.)：

```
pasAnno3<- filter(pasAnno2, pasAnno2$DO51966 !='NULL') 
pasAnno3<- dplyr::select(pasAnno3,-c(DO51973,DO6940,DO6635,DO6630,DO6625,DO6620,DO6730,DO6408,DO6406,DO6400,DO6454,DO6450,DO6452,DO6448,DO6444,DO6438,DO6430,DO6426,DO6422,DO6424,DO6498,DO6494,DO6488,DO6484,DO6486,DO6480,DO6482,DO6476,DO6478,DO6472,DO6474,DO6466,DO6464,DO6460,DO6531,DO6534,DO6525,DO6519,DO6516,DO6510,DO6507,DO6504,DO6501,DO6565,DO6561,DO6558,DO6555,DO6552,DO6546,DO6543,DO6585,DO6615,DO6374,DO6360,DO6356,DO6351,DO6396,DO6394,DO6390,DO7108,DO7136,DO7112,DO7048,DO7084,DO6760,DO6418))
pasAnno4<-as.data.frame(pasAnno3)
row.names(pasAnno4)<-pasAnno4[,1]
pasAnno4<-pasAnno4[,-1]
pasAnno4 = pasAnno4[rowSums(pasAnno4) !=0,]
head(pasAnno4)
```

* The changes I made in this step are not displayed in the picture, but the columns changed.

![image](https://user-images.githubusercontent.com/89613437/144565828-b1c4906d-3519-446e-a354-012f3ea11854.png)

(3) Change gene_id to symbol

```
pasAnno5<- cbind(rownames(pasAnno4), data.frame(pasAnno4, row.names=NULL))
colnames(pasAnno5)[1] = 'gene_id'
pasAnno5<-separate(pasAnno5,col=gene_id,into=c("ensembl_id","b"),sep="[.]" )
pasAnno5<- dplyr::select(pasAnno5,-c(b))
library(org.Hs.eg.db)
#ls("package:org.Hs.eg.db")
g2s<-toTable(org.Hs.egSYMBOL)
g2e<-toTable(org.Hs.egENSEMBL)
b<-merge(pasAnno5,g2e,by='ensembl_id',all.x=T)
d<-merge(b,g2s,by='gene_id',all.x=T)
d<-d[,c(142,2:141,1)]
detach("package:org.Hs.eg.db");
pasAnno6<- dplyr::select(d,-c(ensembl_id,gene_id))
require(tidyverse)
pasAnno6<- filter(pasAnno6, pasAnno6$symbol !='na') %>% 
dplyr::distinct(symbol, .keep_all = TRUE)
head(pasAnno6)
```

* Symbols were used to represent genes.

![image](https://user-images.githubusercontent.com/89613437/144566316-6ddb97a8-63b1-4a13-83c5-68897bc7773d.png)

(4) Store the processed file for subsequent steps.

```
write.csv(pasAnno6, file = "counts.csv", row.names = F, quote = F)
```

## Step 2: Using data in step 1 to analyse difference between two donor type (progression and stable). Draw MA plot, plot counts, variance plot, heatmap, sample-to-sample heatmap, and PCA plot.

(1) Load data (donor information and Htseq-counts from last step and ICGC) in R for the DESeq2.

```
rm(list = ls())
counts0<-read.csv("counts.csv",row.names = 1)
donor0<-read.csv("donor.csv",row.names = 1)
colnames(counts0) == donor0$icgc_donor_id
```

* The donor information in the two files corresponds.

![image](https://user-images.githubusercontent.com/89613437/144567249-e563d105-ac56-4fe2-a289-7f5c28ef74ab.png)

(2)load data in DESeq2 (donor type, donor sex, and donor age at diagnosis).

```
library("DESeq2")
dds<-DESeqDataSetFromMatrix(countData = counts0,
                            colData = donor0,
                            design = ~donor_sex + donor_type + donor_age_at_diagnosis)
```

(3)Pre-filtering: keep the data have more then 10 values in a raw.

```
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
```

(4)Note on factor levels (donor type).

```
dds$donor_type <- factor(dds$donor_type, levels = c("progression","stable"))
```

(5) Differential expression analysis and save data

```
dds<-DESeq(dds)
res <- results(dds, alpha=0.05)
write.csv(as.data.frame(res), 
          file="donor_type_progression_stable.csv")
head(res)
```

* res shows in picture. P value is also in the picture.

![image](https://user-images.githubusercontent.com/89613437/144568250-7d767259-b1fc-49d2-a361-e65ee47335c1.png)

(6)Log fold change shrinkage for visualization and ranking by LFC estimates.

```
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="donor_type_stable_vs_progression", type="apeglm")
sum(resLFC$padj < 0.05, na.rm=TRUE)
sum(resLFC$padj < 0.01, na.rm=TRUE)
```

* There were 1364 P values of genes less than 0.05. There were 452 P values of genes less than 0.01.

(7) draw MS plot with apeglm, normal, and ashr.

```
plotMA(res, ylim=c(-0.3,0.3))
resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(0.1,1e6); ylim <- c(-2,2)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
```

* The MA plot by res is not good.

![image](https://user-images.githubusercontent.com/89613437/144570471-a89d08c2-db02-4586-b356-11c6ae8a3a28.png)

* Using apeglm, normal, and ashr to imporve MA plot. I picked apeglm because it is beautiful and reasonable. It removes the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.

![image](https://user-images.githubusercontent.com/89613437/144570871-49dcce72-6cf0-48ec-adef-d087efc324d6.png)

(8)draw MA-plot with apeglm.

```
plotMA(resLFC, ylim=c(-2,2))
```

* The non-significantly different genes (grey) and the significantly different genes (blue) can be seen on the MA plot. They are symmetrically distributed on both sides of the straight line. The blue genes on the straight line are UP group, and the blue genes under the straight line are DOWN group.

![image](https://user-images.githubusercontent.com/89613437/144573283-f4c1954c-8901-435f-96ed-aca07b57c327.png)

(9) Draw Plot counts with dds. Here I specify the gene which had the smallest p value from the results table created above. You can select the gene to plot by row name or by numeric index.

```
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="donor_type", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=donor_type, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(50,500,5000))
```

* There was no significant difference in counts between the two groups of donors for the gene with the smallest P value.

![image](https://user-images.githubusercontent.com/89613437/144573696-c6ebda97-0c95-497c-99e0-cce972a6d446.png)

(10) Extracting transformed values.

```
vsd <- vst(dds, blind=FALSE)
ntd <- normTransform(dds)
```

(11) Effects of transformations on the variance.

```
library("vsn")
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
```

*  There is a significant change in the variance calculated by ntd.

![image](https://user-images.githubusercontent.com/89613437/144575237-6b002deb-9804-44ab-a275-a39fa44bfc03.png)

*  The variance calculated by vsd is stable. I picked vsd for subsequent steps.

![image](https://user-images.githubusercontent.com/89613437/144575656-7483438a-3292-4f1f-9848-1eadb1652a91.png)

(12) Draw heatmap of the count matrix. The heatmap shows the results of 139 donors and 20 genes. If you need to view specific genes, you can change [1:20].

```
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("donor_type","donor_sex","donor_age_at_diagnosis")])
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df,show_colnames = FALSE)

```

*  The obvious difference in expression can be seen in the picture, but there is no obvious correlation with my variables. If you need to see specific genes, you can modify 1 and 20.

![image](https://user-images.githubusercontent.com/89613437/144576534-ac09014a-e880-4653-a86d-8334e2446d93.png)

(13) Draw heatmap of the sample-to-sample distances.

```
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```

*  The plot shows the differences between the samples, and the PAC plots will be used subsequently to determine whether the differences are due to donor type.

![image](https://user-images.githubusercontent.com/89613437/144578586-6469a5b1-39b5-4b1d-bd2d-856df4f37f9e.png)

(14) Principal component plot of the samples (PCA plot)

```
pcaData <- plotPCA(vsd, intgroup=c("donor_type", "donor_sex"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=donor_sex, shape=donor_type)) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
  
  pcaData <- plotPCA(vsd, intgroup=c("donor_type", "donor_age_at_diagnosis"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=donor_age_at_diagnosis, shape=donor_type)) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
```

*  PCA plot 1 (donor type and donor sex): The donor type in the plot did not show significant differences, but were divided into two groups based on gender.

![image](https://user-images.githubusercontent.com/89613437/144579083-becab4d4-8d9f-4fd1-8c46-835a518be5cc.png)

*  PCA plot 2 (donor type and donor age at diagnosis): Although the data were divided into two groups, there was no relationship with donor type and donor age at diagnosis.

![image](https://user-images.githubusercontent.com/89613437/144579192-d7094ad3-34e3-469a-b380-653d112c183d.png)

**(15) Conclusion 1: The differences can be analyzed according to Htseq-counts, but the differences are not related to the variables I want to explore. I found that the expression differences were related to gender. In the following I will start a new analysis of stage A chronic lymphocytic leukemia according to gender.**

## Step 3: Using data in step 1 to analyse difference between two donor sex (male and female). Draw MA plot, plot counts, variance plot, heatmap, sample-to-sample heatmap, and PCA plot.

(1)load data in R for the DESeq2.

```
rm(list = ls())
counts0<-read.csv("counts.csv",row.names = 1)
donor0<-read.csv("donor.csv",row.names = 1)
```

(2) Donor information is categorized by gender, and Htseq-counts are corresponded to donor information.

```
donor0<- donor0[order(donor0$donor_sex),]
counts1<- t(counts0)
counts1<- cbind(rownames(counts1), data.frame(counts1, row.names=NULL))
counts1<- dplyr::rename(counts1, "icgc_donor_id" = "rownames(counts1)")
counts2<- dplyr::select(donor0,c(icgc_donor_id))
counts2<- dplyr::inner_join(counts2,counts1,by='icgc_donor_id')
row.names(counts2) <- counts2[, 1]
counts2<- counts2[, -1]
counts2<- t(counts2)
colnames(counts2) == donor0$icgc_donor_id
```

* The donor information in the two files corresponds.

![image](https://user-images.githubusercontent.com/89613437/144682587-d7f8e19d-ab89-463f-8f44-b4cd64933985.png)

(3)load data in DESeq2 (donor sex).

```
library("DESeq2")
dds<-DESeqDataSetFromMatrix(countData = counts2,
                            colData = donor0,
                            design = ~donor_sex)
```

(4)Pre-filtering: keep the data have more then 10 values in a raw.

```
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
```

(5) Note on factor levels.

```
dds$donor_type <- factor(dds$donor_type, levels = c("male","female"))
```

(6) Differential expression analysis and save data

```{r}
dds<-DESeq(dds)
res <- results(dds, alpha=0.05)
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), 
          file="donor_sex_male_vs_female_p_order.csv")
summary(res)
sum(res$padj < 0.05, na.rm=TRUE)
```

* There were 24,300 genes analyzed. 117 genes had p-values less than 0.05.

![image](https://user-images.githubusercontent.com/89613437/144683189-f3536db1-73bc-40fe-8cb9-a3f9cca6d15e.png)

(7)Log fold change shrinkage for visualization and ranking by LFC estimates.

```
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="donor_sex_male_vs_female", type="apeglm")
summary(resLFC)
sum(resLFC$padj < 0.05, na.rm=TRUE)
```

* After LFC estimates， the results do not change.

![image](https://user-images.githubusercontent.com/89613437/144683321-4cd637fd-77e9-42aa-940f-f28bd7bcce90.png)

(8) draw MS plot with apeglm, normal, and ashr.

```
plotMA(res,  xlim=c(0.01,1e5), ylim=c(-10,15))
#resNorm <- lfcShrink(dds, coef=2, type="normal")
#resAsh <- lfcShrink(dds, coef=2, type="ashr")
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(0.1,1e6); ylim <- c(-1.5,1.5)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
```

* The MA plot by res is not perfect.

![image](https://user-images.githubusercontent.com/89613437/144683514-10a34260-bce7-4471-8159-da8dc869b680.png)

* Using apeglm, normal, and ashr to imporve MA plot.

![image](https://user-images.githubusercontent.com/89613437/144683564-9e6d0f9e-d5f0-4bf9-951b-408a140b1820.png)

(9)draw MA-plot with ashr because it retains more functions and is suitable for handling large amounts of data.

```{r}
plotMA(resAsh,xlim=c(0.1,1e6), ylim=c(-6,6))
```

* This picture is not perfect, but it shows the difference genes clearly.

![image](https://user-images.githubusercontent.com/89613437/144683611-e04308de-13ff-4f81-9652-fe24b04884d5.png)

(10) Draw Plot counts with dds. Here I specify the gene which had the smallest p value from the results table created above. You can select the gene to plot by row name or by numeric index.

```
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="donor_sex", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=donor_sex, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(50,500,10000))
```

* It is clear from the figure that there are expression differences between the two groups. The male expression is significantly higher than the female. Some women do not even express this gene. Five males expressed significantly less than others because of individual differences, but all males expressed this gene.

![image](https://user-images.githubusercontent.com/89613437/144683772-dd01b292-2f00-4a61-bf9e-0a8a5c1f1b45.png)

(11)Extracting transformed values

```
vsd <- vst(dds, blind=FALSE)
ntd <- normTransform(dds)
```

(12) Effects of transformations on the variance

```
library("vsn")
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
```

* There is a significant change in the variance calculated by ntd.

![image](https://user-images.githubusercontent.com/89613437/144684065-ef52b24b-23fd-4587-9e05-431151db8fd4.png)

* The variance calculated by vsd is stable. I picked vsd for subsequent steps.

![image](https://user-images.githubusercontent.com/89613437/144684114-03ea4a7e-ec3f-4b70-8023-02a158b7c401.png)

(13) Heatmap of the count matrix. The heatmap shows the results of 139 donors and 20 genes. If you need to view specific genes, you can change [626:645].

```
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[626:645]
df <- dplyr::select(donor0,c(icgc_donor_id,donor_sex))
rownames(df) = df[,1]
df<- dplyr::select(df,-c(icgc_donor_id))
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=T,
         cluster_cols=FALSE, annotation_col=df, show_colnames = FALSE)

```

* The color of the DDX3Y gene and the sex of the donor in the figure have a corresponding relationship. There are 116 other similar genes. Other genes can be observed in the heatmap by modifying the code（626:645）.

![image](https://user-images.githubusercontent.com/89613437/144684190-5915f911-d68a-42a3-942a-a22eb087867b.png)

(14) Heatmap of the sample-to-sample distances

```
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```

* The plot shows the differences between the samples. This difference is caused by gender.

![image](https://user-images.githubusercontent.com/89613437/144684622-0fc10989-c31e-4f25-b228-e4bae4253a95.png)

(15) Principal component plot of the samples (PCA plot)

```
pcaData <- plotPCA(vsd, intgroup=c("donor_sex"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=donor_sex)) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
```

* The data were clearly divided into two groups. One male was assigned to the female group, which may be due to his individuality. In the plot counts, there is a male and a female with similar data. This outlier is caused by the same person as the outlier in plot counts. The other five outliers are also caused by the same factors.

![image](https://user-images.githubusercontent.com/89613437/144684700-056e55ae-4e74-493b-8dba-a89a250a2897.png)

**(16) Conclusion 2: There were 117 genes that produced significant differences between the two groups of men and women. Which of these genes are associated with chronic lymphocytic leukemia requires further analysis.**

## Evaluation genes

(1) I will use the GSEA website to evaluate genes (https://www.gsea-msigdb.org/gsea/index.jsp). Genes with a p-value of 0 were  analyzed using the investigate gene sites. There were seven genes judged to have a p-value equal to 0: TXLNGY, PRKY, TTTY15, ZFY, KDM5D, USP9Y, DDX3Y, and EIF1AY. The results show that four of these genes and three gene sets overlap. There are two down-regulated genes in PTEN_DN.V2_DN [143]; two up-regulated genes in RELA_DN.V1_UP [149]; and two down-regulated genes in PGF_UP.V1_DN [190]. I will continue to find what DDX3Y, ZFY, USP9Y, and KDM5D affect on human.

![image](https://user-images.githubusercontent.com/89613437/144692782-47e4de55-c73e-4685-820f-a05622653cce.png)

![image](https://user-images.githubusercontent.com/89613437/144692793-c101356a-fbf7-43cf-aadb-27d7acf729b3.png)

(2) Search Gene Sets: PTEN_DN.V2_DN [143], RELA_DN.V1_UP [149], and PGF_UP.V1_DN [190] are not related to blood cancers. Therefore, there was no data to determine in the GSEA study whether the genes I screened for were associated with chronic lymphocytic leukemia.

## Known Issues

(1) The number of donors is relatively small.

(2) The sample-to-sample heat map does not reflect much information.

(3) I only selected the genes with the smallest p-values for analysis, and there are still some eligible genes that can be analyzed.

## Conclusions

According to Inherited predisposition to chronic lymphocytic leukemia, CXCR4, SMAD7, and DAPK may cause chronic lymphocytic leukemia. Therefore, the genes I picked are not in the article. After I search the foctions of the genes, I believe that these genes are related to differences in expression between  male and female, and they are not related to chronic lymphocytic leukemia. This guess is based on the fact that these genes are all highly expressed in males, but expressed in low or no amounts in females. On the GeneCards, DDX3Y, ZFY, USP9Y, and KDM5D are all related to the Y chromosome. Therefore, my htseq-counts data were divided in two groups by gender not chronic lymphocytic leukemia.



