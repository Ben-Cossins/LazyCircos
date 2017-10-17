#' create_tile
#' 
#' Creates the files for each condition
#' The karyotype.txt file will be saved within the current working directory in a folder called data
#' @export


create_tile = function(){
    .karyotype <<- read.csv("data/karyotype.txt", sep = "\t", header = FALSE)
    i = 2
    j = (2*.Condition)+2
        for (k in 1:.Condition){
        file = paste("data/Tile",k,".txt", sep = "")
        sink(file = file)
        if (.Comparison == 2){
            cat(paste(.karyotype[i,3], .karyotype[i,5], .karyotype[i,6], "\n", sep = "\t"))
            cat(paste(.karyotype[j,3], .karyotype[j,5], .karyotype[j,6], "\n", sep = "\t"))
            #print(paste(i, j))
            i = i + 1
            j = j - 1
        } else {
            for (i in 1:(((.Condition + 1)* .Comparison)+1)){
                if (((i + .Condition)- k) %% (.Condition + 1) == 0){
                    cat(paste(.karyotype[i,3], .karyotype[i,5], .karyotype[i,6], "\n", sep = "\t"))
                }
            }
        }
        sink(NULL)
    }
}

