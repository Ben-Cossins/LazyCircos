#' reverse_lists
#' 
#' split the lists of conditions and reverse the order of the second comparison.
#' @export

reverse_lists = function(){
    if (.Comparison == 2){
        order_sheet_list = list()
        order_sheet_list[[1]] = sheet_list[1:.Condition]
        order_sheet_list[[2]] = rev(sheet_list[.Condition+1:.Condition])
    
        order_gene_list = list()
        i = 1
        x = length(gene_list)
        for (k in 1:.Condition){
            order_gene_list[[i]] = gene_list[[k]]
            i = i+1
        }
        for (k in x:(1+.Condition)){
            order_gene_list[[i]] = gene_list[[k]]
            .chr_length[i] <<- nrow(gene_list[[k]])
            i = i+1
        }
        .sheet_list <<- sheet_list
        sheet_list <<- append(order_sheet_list[[1]], order_sheet_list[[2]])
        .gene_list <<- gene_list
        gene_list <<- order_gene_list
        
        
    } else {print("More than 2 comparisons, shouldn't reverse conditions.")}
}
