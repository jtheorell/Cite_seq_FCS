# Cite seq FCS
### View 82 human surface markers simultaneously in flow cytometry software. Data from [Kotliarov et al, Nat. Med. 2020](https://www.nature.com/articles/s41591-020-0769-8).

Here, a set of flow cytometry standard (FCS) files are presented, that have been extracted from the 10X dataset of healthy individuals undergoing a vaccination trial, presented in: Kotliarov, Y., Sparks, R. et al. Broad immune activation underlies shared set point signatures for vaccine responsiveness in healthy individuals and disease activity in patients with lupus. Nat. Med. [DOI](https://doi.org/10.1038/s41591-020-0769-8) (2020) . 

The raw data can be retrieved [here](https://nih.figshare.com/articles/dataset/CITE-seq_protein-mRNA_single_cell_data_from_high_and_low_vaccine_responders_to_reproduce_Figs_4-6_and_associated_Extended_Data_Figs_/11349761?backTo=/collections/).

The preprocessing steps are shown in the file "Sce2fcs_processing.R". The preprocessing steps closely follow Robert Amezquita and Aaron Luns work ["Orchestrating Single-Cell Analysis with Bioconductor"](https://bioconductor.org/books/release/OSCA/) . In this case, it is especially [chapter 20](https://bioconductor.org/books/release/OSCA/integrating-with-protein-abundance.html) that is used.

To look at the data in the flow cytometry software, download the files (one per donor, name format: X_yyy, where X = batch and xxx = individual) and add them like fcs files are normally added to your software of choice. I have found arcsinh transformation to work well, but have tweaked the settings differently for each marker. 

NB! Marker CD133, which is a stem cell marker, appears to be expressed in a subset of B cells. When going through the transcriptome data of the CD133 "positive" cells, there is not a single one that expresses a CD133 transcript. This is also true for another, independent dataset with B cells. Therefore, this is for now considered a technical artifact. 
