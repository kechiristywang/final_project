

## Step 1: Change the downloaded data format to a format suitable for DESesq2 analysis. Change gene_id to symbols.

(1) load Htseq-counts data in R, and keep the required data. Change data format (length data to width data).

* This is the raw RNA-seq data from the ICGC Data Portal.

![image](https://user-images.githubusercontent.com/89613437/144563579-f829d0bb-960c-430f-b660-5e2e143bb50a.png)

* Eliminate unwanted genes and keep Htseq-counts.

![image](https://user-images.githubusercontent.com/89613437/144564320-b2553c60-775e-4c1d-810e-a3bd55cb8970.png)

* Long format data is not conducive to subsequent analysis; turn long format data into wide data format.

![image](https://user-images.githubusercontent.com/89613437/144564556-b11f5427-c5ee-4541-b5ef-cc6e04999e33.png)

(2) Delete invalid data (Delete all NULL rows, incomplete data, and all zero rows.)：

* The changes I made in this step are not displayed in the picture, but the columns changed.

![image](https://user-images.githubusercontent.com/89613437/144565828-b1c4906d-3519-446e-a354-012f3ea11854.png)

(3) Change gene_id to symbol

* Symbols were used to represent genes.

![image](https://user-images.githubusercontent.com/89613437/144566316-6ddb97a8-63b1-4a13-83c5-68897bc7773d.png)

(4) Store the processed file for subsequent steps.

## Step 2: Using data in step 1 to analyse difference between two donor type (progression and stable). Draw MA plot, plot counts, variance plot, heatmap, sample-to-sample heatmap, and PCA plot.

(1) Load data (donor information and Htseq-counts from last step and ICGC) in R for the DESeq2.

* The donor information in the two files corresponds.

![image](https://user-images.githubusercontent.com/89613437/144567249-e563d105-ac56-4fe2-a289-7f5c28ef74ab.png)

(2)load data in DESeq2 (donor type, donor sex, and donor age at diagnosis).

(3)Pre-filtering: keep the data have more then 10 values in a raw.

(4)Note on factor levels (donor type).

(5) Differential expression analysis and save data

* res shows in picture. P value is also in the picture.

![image](https://user-images.githubusercontent.com/89613437/144568250-7d767259-b1fc-49d2-a361-e65ee47335c1.png)

(6)Log fold change shrinkage for visualization and ranking by LFC estimates.

* There were 1364 P values of genes less than 0.05. There were 452 P values of genes less than 0.01.

(7) draw MS plot with apeglm, normal, and ashr.

* The MA plot by res is not good.

![image](https://user-images.githubusercontent.com/89613437/144570471-a89d08c2-db02-4586-b356-11c6ae8a3a28.png)

* Using apeglm, normal, and ashr to imporve MA plot. I picked apeglm because it is beautiful and reasonable. It removes the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.

![image](https://user-images.githubusercontent.com/89613437/144570871-49dcce72-6cf0-48ec-adef-d087efc324d6.png)

(8)draw MA-plot with apeglm.

* The non-significantly different genes (grey) and the significantly different genes (blue) can be seen on the MA plot. They are symmetrically distributed on both sides of the straight line. The blue genes on the straight line are UP group, and the blue genes under the straight line are DOWN group.

![image](https://user-images.githubusercontent.com/89613437/144573283-f4c1954c-8901-435f-96ed-aca07b57c327.png)

(9) Draw Plot counts with dds. Here I specify the gene which had the smallest p value from the results table created above. You can select the gene to plot by row name or by numeric index.

* There was no significant difference in counts between the two groups of donors for the gene with the smallest P value.

![image](https://user-images.githubusercontent.com/89613437/144573696-c6ebda97-0c95-497c-99e0-cce972a6d446.png)

(10) Extracting transformed values.

(11) Effects of transformations on the variance.

*  There is a significant change in the variance calculated by ntd.

![image](https://user-images.githubusercontent.com/89613437/144575237-6b002deb-9804-44ab-a275-a39fa44bfc03.png)

*  The variance calculated by vsd is stable. I picked vsd for subsequent steps.

![image](https://user-images.githubusercontent.com/89613437/144575656-7483438a-3292-4f1f-9848-1eadb1652a91.png)

(12) Draw heatmap of the count matrix. The heatmap shows the results of 139 donors and 20 genes. If you need to view specific genes, you can change [1:20].

*  The obvious difference in expression can be seen in the picture, but there is no obvious correlation with my variables. If you need to see specific genes, you can modify 1 and 20.

![image](https://user-images.githubusercontent.com/89613437/144576534-ac09014a-e880-4653-a86d-8334e2446d93.png)

(13) Draw heatmap of the sample-to-sample distances.

*  The plot shows the differences between the samples, and the PAC plots will be used subsequently to determine whether the differences are due to donor type.

![image](https://user-images.githubusercontent.com/89613437/144578586-6469a5b1-39b5-4b1d-bd2d-856df4f37f9e.png)

(14) Principal component plot of the samples (PCA plot)

*  PCA plot 1 (donor type and donor sex): The donor type in the plot did not show significant differences, but were divided into two groups based on gender.

![image](https://user-images.githubusercontent.com/89613437/144579083-becab4d4-8d9f-4fd1-8c46-835a518be5cc.png)

*  PCA plot 2 (donor type and donor age at diagnosis): Although the data were divided into two groups, there was no relationship with donor type and donor age at diagnosis.

![image](https://user-images.githubusercontent.com/89613437/144579192-d7094ad3-34e3-469a-b380-653d112c183d.png)

(15) **Conclusion: The differences can be analyzed according to Htseq-counts, but the differences are not related to the variables I want to explore. I found that the expression differences were related to gender. In the following I will start a new analysis of stage A chronic lymphocytic leukemia according to gender.**














