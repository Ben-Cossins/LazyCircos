#' create_karyotype
#' 
#' Generate karyotype file. run after initialising circos as it is dependent on variables declared within that function
#' The karyotype.txt file will be saved within the current working directory in a folder called data
#' @export 




create_karyotype = function(){
    m = 0
    j = 1
    k = 1
    x = .Condition +1
    sink("data/karyotype.txt")
    for (i in 1:(((.Condition + 1) * .Comparison) + 1)){
        e = (((.Condition + 1) * .Comparison) + 1)
        if (((i+.Condition) %% x) == 0){
            if (i == 1 || i == e){
                cat(paste("chr", "-", paste("Blk.", j, sep=""), i, m, m + (.gap/2), "white\n",  sep = "\t"))
                m = m + (.gap/2)
                j = j + 1
            } else {
                cat(paste("chr", "-", paste("Blk.", j, sep=""), i, m, m + .gap, "white\n",  sep = "\t"))
                m = m + .gap
                j = j + 1
            }
        } else {
            cat(paste("chr", "-", sheet_list[[k]], i, m, m + nrow(gene_list[[k]]), "black\n",  sep = "\t"))
            m = m + nrow(gene_list[[k]])
            k = k + 1
        }
    }
    sink(NULL)
}






