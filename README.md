# SPACED
SPACED (Spatial Classification of mRNA Expression Data) analysis is a computational method to infer
spatial organization of celltypes identified in scRNA-seq data.
This method requires the following parameters as input:

                 1.obj -> seurat object
                 2.degs -> differential expressed genes
                 3.ident_use -> cell identity class
                 4.celltype ->  cell identity

SPACED exhibited greater performance with the top10 most cluster specific genes 
SPACED is restricted to the brain regions which have clear anatomic structures and distinguishable borders
