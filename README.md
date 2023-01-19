# Mapping single cell data from cellenONE, to proteomics data


This code allows you to map single cell data from the cellenONE to your proteomics data when preparing samples via [nPOP](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02817-5). Auto maps your cell types to the wells of the 384 plate you inject from.

The example code here is for some 2-plex plexDIA sets. It also inculdes some quick data quality QCs that are relevant when using plexDIA. Minor changes to the script need to be made for mapping to 3plex.

Code for the downstream QCs will not apply for the 14plex TMT data here but the function for mapping cellenONE data will still be useful.
