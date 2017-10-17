# LazyCircos
Create the configuration files for generating circos from genelists


Loads an Excel workbook or openDocumentSpreadsheet containing lists of genes, finds overlaps between these lists, and generates the configuration files for generating circos plot.

This still requires the perl implementation of circos to be installed, and only creates the files for it.


With an excel workbook comprising of a list of genes on a sheet per condition (e.g. up and down regulated genes in treated vs untreated, would be 4 sheets in the workbook), the configuration files for the circos can be derived in as few as 4 lines


install_github("Ben-Cossins/LazyCircos")
library(LazyCircos)
read_genelist("yourfilehere")
initialise_circos(Comparison = x, Condition = y)
