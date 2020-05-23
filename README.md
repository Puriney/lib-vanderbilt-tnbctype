# Gene signatures for TNBC subtypes

Genes are downloaded from the TNBCtype paper "Identification of human triple-negative breast cancer subtypes and preclinical models for selection of targeted therapies" (https://www.jci.org/articles/view/45014#top).


# Result

Use the two files:
- vanderbilt_tnbc_subtypes_gene_signatures_plus.gmt
- vanderbilt_tnbc_subtypes_gene_signatures_minus.gmt

The 'plus' and 'minus' postfix are used rather than 'up'/'dn' because of the [vision package](https://yoseflab.github.io/VISION/articles/Signatures.html#create-your-own--gmt-files) so that the signature autocorrelation analysis is ready to do.

The original gene symbols are converted the official ones using limma::alias2Symbol. Keep the original ones if no mappings.
