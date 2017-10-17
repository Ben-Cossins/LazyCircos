#' read_genelist
#' 
#' Reads in an XLSX file or ODS file containing multiple sheets each with a single list of genes and creates a list of dataframes
#' @param file Location of the file to be read in
#' @param type Type of file to be read in. Default to "xlsx" (other option is "ods")
#' @param header Specifies whether sheets contain headers or not. Defaults to TRUE
#' @import openxlsx
#' @import readODS
#' @examples  
#' read_genelist(file = "test.xlsx", type = "xlsx", header = TRUE)
#' read_genelist(file = "test.ods", type = "ods", header = FALSE)
#' @export


read_genelist = function(file = "", type =  "xlsx", header = TRUE){
    a2 = list()
    if (toupper(type) == "XLSX"){
        a = loadWorkbook(file)
        a1 = sheets(a)
        for (i in 1:length(a1)){
            dat = readWorkbook(file, sheet = i, colNames = header)
            a2[[a1[i]]] = dat
        }
    } else if (toupper(type) == "ODS"){
        a1 = ods_sheets(file)
        for (i in 1:length(a1)){
            dat = suppressMessages(read_ods(file, i, col_names =  header))
            a2[[a1[i]]] = dat
        }
    } else {
        print("type of file unknown")
        stop()
    }
    
    duplicatedgenes = list()
    for (i in 1:length(a2)){
        duplicatedgenes[[i]] = a2[[i]][duplicated(a2[[i]][,1]),1]
    }
    if (sum(lengths(duplicatedgenes)) > 0){
        cat("\n\nInputted lists contain duplicates, whilst these will be plotted on the chromosomes, they will not contribute to the ribbons linking them \n\n")
        for (i in 1:length(duplicatedgenes)){
            if (length(duplicatedgenes[[i]] > 0 )){
                cat(paste("Duplicated genes on sheet: ", a1[i], "\n", sep = "" ))
                cat(duplicatedgenes[[i]],"\n")
                cat("\n")
            }
        }
    }
    
    
    sheet_list <<- gsub("[[:space:]]", ".", a1)
    gene_list <<-a2
    
    
    .chr_length = integer()
    for (i in 1:length(sheet_list)){
        .chr_length[i] = length(gene_list[[i]][,1])
    }
    .chr_length <<- .chr_length
}
