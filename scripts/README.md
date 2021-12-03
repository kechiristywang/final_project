

## Step 1:
Change the downloaded data format to a format suitable for DESesq2 analysis. Change gene_id to symbols.

(1) load Htseq-counts data in R, and keep the required data. Change data format (length data to width data).

* This is the raw RNA-seq data from the ICGC Data Portal.
![image](https://user-images.githubusercontent.com/89613437/144563579-f829d0bb-960c-430f-b660-5e2e143bb50a.png)
* Eliminate unwanted genes and keep Htseq-counts.
![image](https://user-images.githubusercontent.com/89613437/144564320-b2553c60-775e-4c1d-810e-a3bd55cb8970.png)
* Long format data is not conducive to subsequent analysis; turn long format data into wide data format.
![image](https://user-images.githubusercontent.com/89613437/144564556-b11f5427-c5ee-4541-b5ef-cc6e04999e33.png)










