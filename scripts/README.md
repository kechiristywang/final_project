

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

* Using apeglm, normal, and ashr to imporve MA plot. I picked apeglm because it removes the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.
![image](https://user-images.githubusercontent.com/89613437/144570871-49dcce72-6cf0-48ec-adef-d087efc324d6.png)













