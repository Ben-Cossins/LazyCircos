#' create_2x3_links
#' 
#' Calculate ribbons for 2 comparisons with 3 conditions 
#' @export


create_2x3_links = function(){
    a.1.4 = intersect(gene_list[[1]][,1], gene_list[[4]][,1])
    a.1.5 = intersect(gene_list[[1]][,1], gene_list[[5]][,1])
    a.1.6 = intersect(gene_list[[1]][,1], gene_list[[6]][,1])
    
    a.2.4 = intersect(gene_list[[2]][,1], gene_list[[4]][,1])
    a.2.5 = intersect(gene_list[[2]][,1], gene_list[[5]][,1])
    a.2.6 = intersect(gene_list[[2]][,1], gene_list[[6]][,1])
    
    a.3.4 = intersect(gene_list[[3]][,1], gene_list[[4]][,1])
    a.3.5 = intersect(gene_list[[3]][,1], gene_list[[5]][,1])
    a.3.6 = intersect(gene_list[[3]][,1], gene_list[[6]][,1])
    
    a.1.u = setdiff(gene_list[[1]][,1], c(a.1.4, a.1.5, a.1.6))
    a.2.u = setdiff(gene_list[[2]][,1], c(a.2.4, a.2.5, a.2.6))
    a.3.u = setdiff(gene_list[[3]][,1], c(a.3.4, a.3.5, a.3.6))
    a.4.u = setdiff(gene_list[[4]][,1], c(a.1.4, a.2.4, a.3.4))
    a.5.u = setdiff(gene_list[[5]][,1], c(a.1.5, a.2.5, a.3.5))
    a.6.u = setdiff(gene_list[[6]][,1], c(a.1.6, a.2.6, a.3.6))
    
    ribbonlist = list(list(a.1.6, a.1.5, a.1.4, a.1.u),
                      list(a.2.6, a.2.5, a.2.4, a.2.u),
                      list(a.3.6, a.3.5, a.3.4, a.3.u),
                      list(a.3.4, a.2.4, a.1.4, a.4.u),
                      list(a.3.5, a.2.5, a.1.5, a.5.u),
                      list(a.3.6, a.2.6, a.1.6, a.6.u))
    
    
    a.1.cg = sum(lengths(ribbonlist[[1]][1:3]))
    a.2.cg = sum(lengths(ribbonlist[[2]][1:3]))
    a.3.cg = sum(lengths(ribbonlist[[3]][1:3]))
    a.4.cg = sum(lengths(ribbonlist[[4]][1:3]))
    a.5.cg = sum(lengths(ribbonlist[[5]][1:3]))
    a.6.cg = sum(lengths(ribbonlist[[6]][1:3]))
    
    
    
    r.t = function(x){return(as.numeric(lengths(x) != 0))} # calulates if there are any links present requiring ribbons
    
    a.1.r = sum(r.t(ribbonlist[[1]])[1:3])
    a.2.r = sum(r.t(ribbonlist[[2]])[1:3])
    a.3.r = sum(r.t(ribbonlist[[3]])[1:3])
    a.4.r = sum(r.t(ribbonlist[[4]])[1:3])
    a.5.r = sum(r.t(ribbonlist[[5]])[1:3])
    a.6.r = sum(r.t(ribbonlist[[6]])[1:3])
    
    spacing.table = data.frame(Chrsize = .chr_length, 
                               crossgenes = c(a.1.cg, a.2.cg, a.3.cg, a.4.cg, a.5.cg, a.6.cg),
                               ribbons = c(a.1.r, a.2.r, a.3.r, a.4.r, a.5.r, a.6.r))
    spacing.table$ribbon.gap = floor((spacing.table$Chrsize - spacing.table$crossgenes) 
                                     / (spacing.table$ribbons))
    
    NoGapKaryotype = .karyotype[c(2,3,4,6,7,8),]
    
    chr.coordinates = list()
    for(i in 1:(.Comparison * .Condition)){
        dat = data.frame(chr = character(), start = integer(), end = integer(), ribbon = integer(), 
                         stringsAsFactors=FALSE)
        start = NoGapKaryotype[i,5]
        for (k in 1:.Condition){
            dat[k,1] = sheet_list[i]
            dat[k,2] = start
            dat[k,3] = start + length(ribbonlist[[i]][[k]])
            dat[k,4] = r.t(ribbonlist[[i]])[k]
            start = dat[k,3]
            if (dat[k,4] == 1){start = start + spacing.table[i,4]}
            #cat(paste(dat[k,1], dat[k,2], dat[k,3], dat[k,4], (dat[k,3]-dat[k,2]), "\n", sep="\t"))
        }
        chr.coordinates[[i]] = dat
    }
    
    
    sink(file = "data/Ribbon1.txt")
    if(chr.coordinates[[1]][1,4] == 1){
        cat(paste(chr.coordinates[[1]][1,1], "\t", chr.coordinates[[1]][1,2], "\t", chr.coordinates[[1]][1,3], "\t",
                  chr.coordinates[[6]][3,1], "\t", chr.coordinates[[6]][3,2], "\t", chr.coordinates[[6]][3,3], "\n",
                  sep = ""))}
    if(chr.coordinates[[1]][2,4] == 1){
        cat(paste(chr.coordinates[[1]][2,1], "\t", chr.coordinates[[1]][2,2], "\t", chr.coordinates[[1]][2,3], "\t",
                  chr.coordinates[[5]][3,1], "\t", chr.coordinates[[5]][3,2], "\t", chr.coordinates[[5]][3,3], "\n",
                  sep = ""))}
    if(chr.coordinates[[1]][3,4] == 1){
        cat(paste(chr.coordinates[[1]][3,1], "\t", chr.coordinates[[1]][3,2], "\t", chr.coordinates[[1]][3,3], "\t",
                  chr.coordinates[[4]][3,1], "\t", chr.coordinates[[4]][3,2], "\t", chr.coordinates[[4]][3,3], "\n",
                  sep = ""))}
    
    if(chr.coordinates[[2]][1,4] == 1){
        cat(paste(chr.coordinates[[2]][1,1], "\t", chr.coordinates[[2]][1,2], "\t", chr.coordinates[[2]][1,3], "\t",
                  chr.coordinates[[6]][2,1], "\t", chr.coordinates[[6]][2,2], "\t", chr.coordinates[[6]][2,3], "\n",
                  sep = ""))}
    if(chr.coordinates[[2]][2,4] == 1){
        cat(paste(chr.coordinates[[2]][2,1], "\t", chr.coordinates[[2]][2,2], "\t", chr.coordinates[[2]][2,3], "\t",
                  chr.coordinates[[5]][2,1], "\t", chr.coordinates[[5]][2,2], "\t", chr.coordinates[[5]][2,3], "\n",
                  sep = ""))}
    if(chr.coordinates[[2]][3,4] == 1){
        cat(paste(chr.coordinates[[2]][3,1], "\t", chr.coordinates[[2]][3,2], "\t", chr.coordinates[[2]][3,3], "\t",
                  chr.coordinates[[4]][2,1], "\t", chr.coordinates[[4]][2,2], "\t", chr.coordinates[[4]][2,3], "\n",
                  sep = ""))}
    
    if(chr.coordinates[[3]][1,4] == 1){
        cat(paste(chr.coordinates[[3]][1,1], "\t", chr.coordinates[[3]][1,2], "\t", chr.coordinates[[3]][1,3], "\t",
                  chr.coordinates[[6]][1,1], "\t", chr.coordinates[[6]][1,2], "\t", chr.coordinates[[6]][1,3], "\n",
                  sep = ""))}
    if(chr.coordinates[[3]][1,4] == 1){
        cat(paste(chr.coordinates[[3]][2,1], "\t", chr.coordinates[[3]][2,2], "\t", chr.coordinates[[3]][2,3], "\t",
                  chr.coordinates[[5]][1,1], "\t", chr.coordinates[[5]][1,2], "\t", chr.coordinates[[5]][1,3], "\n",
                  sep = ""))}
    if(chr.coordinates[[3]][1,4] == 1){
        cat(paste(chr.coordinates[[3]][3,1], "\t", chr.coordinates[[3]][3,2], "\t", chr.coordinates[[3]][3,3], "\t",
                  chr.coordinates[[4]][1,1], "\t", chr.coordinates[[4]][1,2], "\t", chr.coordinates[[4]][1,3], "\n",
                  sep = ""))}
    sink(NULL)
    
    
    
    ribbonsizes = data.frame(Chr1 = character(), Chr2 = character(), size = integer(), 
                             stringsAsFactors=FALSE)
    tmp = read.csv("data/Ribbon1.txt", header = FALSE, sep = "\t") 
    for (i in 1:nrow(tmp)){
        ribbonsizes[i,1] = as.character(tmp[i,1])
        ribbonsizes[i,2] = as.character(tmp[i,4])
        ribbonsizes[i,3] = tmp[i,3] - tmp[i,2]
    }
    wb = createWorkbook()
    addWorksheet(wb, "Ribbon1", gridLines = TRUE)
    writeData(wb, sheet = 1, ribbonsizes, rowNames = FALSE)
    saveWorkbook(wb, "RibbonSizes.xlsx", overwrite = TRUE)
}
