#' create_3x2_links
#' 
#' Calculate ribbons for 3 comparisons with 2 conditions 
#' @export


create_3x2_links = function(){
    a.1.3 = intersect(gene_list[[1]][,1], gene_list[[3]][,1])
    a.1.4 = intersect(gene_list[[1]][,1], gene_list[[4]][,1])
    a.1.5 = intersect(gene_list[[1]][,1], gene_list[[5]][,1])
    a.1.6 = intersect(gene_list[[1]][,1], gene_list[[6]][,1])
    
    a.2.3 = intersect(gene_list[[2]][,1], gene_list[[3]][,1])
    a.2.4 = intersect(gene_list[[2]][,1], gene_list[[4]][,1])
    a.2.5 = intersect(gene_list[[2]][,1], gene_list[[5]][,1])
    a.2.6 = intersect(gene_list[[2]][,1], gene_list[[6]][,1])
    
    a.3.5 = intersect(gene_list[[3]][,1], gene_list[[5]][,1])
    a.3.6 = intersect(gene_list[[3]][,1], gene_list[[6]][,1])
    a.4.5 = intersect(gene_list[[4]][,1], gene_list[[5]][,1])
    a.4.6 = intersect(gene_list[[4]][,1], gene_list[[6]][,1])
    
    b.1.3.5 = intersect(a.1.3, a.1.5)
    b.1.3.6 = intersect(a.1.3, a.1.6)
    b.1.4.5 = intersect(a.1.4, a.1.5)
    b.1.4.6 = intersect(a.1.4, a.1.6)
    
    b.2.3.5 = intersect(a.2.3, a.2.5)
    b.2.3.6 = intersect(a.2.3, a.2.6)
    b.2.4.5 = intersect(a.2.4, a.2.5)
    b.2.4.6 = intersect(a.2.4, a.2.6)
    
    a.1.3 = setdiff(a.1.3, c(b.1.3.5, b.1.3.6))
    a.1.4 = setdiff(a.1.4, c(b.1.4.5, b.1.4.6))
    a.1.5 = setdiff(a.1.5, c(b.1.3.5, b.1.4.5))
    a.1.6 = setdiff(a.1.6, c(b.1.3.5, b.1.4.6))
    
    a.2.3 = setdiff(a.2.3, c(b.2.3.5, b.2.3.6))
    a.2.4 = setdiff(a.2.4, c(b.2.4.5, b.2.4.6))
    a.2.5 = setdiff(a.2.5, c(b.2.3.5, b.2.4.5))
    a.2.6 = setdiff(a.2.6, c(b.2.3.6, b.2.4.6))
    
    a.3.5 = setdiff(a.3.5, c(b.1.3.5, b.2.3.5))
    a.3.6 = setdiff(a.3.6, c(b.1.3.6, b.2.3.6))
    a.4.5 = setdiff(a.4.5, c(b.1.4.5, b.2.4.5))
    a.4.6 = setdiff(a.4.6, c(b.1.4.6, b.2.4.6))
    
    a.1.u = setdiff(gene_list[[1]][,1], c(a.1.3,a.1.4,a.1.5,a.1.6,b.1.3.5,b.1.3.6,b.1.4.5,b.1.4.6))
    a.2.u = setdiff(gene_list[[2]][,1], c(a.2.3,a.2.4,a.2.5,a.2.6,b.2.3.5,b.2.3.6,b.2.4.5,b.2.4.6))
    a.3.u = setdiff(gene_list[[3]][,1], c(a.1.3,a.2.3,a.3.5,a.3.6,b.1.3.5,b.1.3.6,b.2.3.5,b.2.3.6))
    a.4.u = setdiff(gene_list[[4]][,1], c(a.1.4,a.2.4,a.4.5,a.4.6,b.1.4.5,b.1.4.6,b.2.4.5,b.2.4.6))
    a.5.u = setdiff(gene_list[[5]][,1], c(a.1.5,a.2.5,a.3.5,a.4.5,b.1.3.5,b.1.4.5,b.2.3.5,b.2.4.5))
    a.6.u = setdiff(gene_list[[6]][,1], c(a.1.6,a.2.6,a.3.6,a.4.6,b.1.3.6,b.1.4.6,b.2.3.6,b.2.4.6))
    

    ribbonlist = list(list(a.1.6, b.1.4.6, b.1.3.6, a.1.5, a.1.4, b.1.4.5, b.1.3.5, a.1.3, a.1.u),
                      list(a.2.6, b.2.4.6, b.2.4.5, a.2.5, a.2.4, b.2.3.6, b.2.3.5, a.2.3, a.2.u),
                      list(a.2.3, b.2.3.5, b.2.3.6, a.1.3, a.3.6, b.1.3.6, b.1.3.5, a.3.5, a.3.u),
                      list(a.2.4, b.2.4.6, b.1.4.6, a.1.4, a.4.6, b.1.4.5, b.2.4.5, a.4.5, a.4.u),
                      list(a.4.5, b.2.4.5, b.1.4.5, a.3.5, a.2.5, b.2.3.5, b.1.3.5, a.1.5, a.5.u),
                      list(a.4.6, b.2.4.6, b.2.3.6, a.3.6, a.2.6, b.1.3.6, b.1.4.6, a.1.6, a.6.u))
    
        
    # calculate number of ribbons, total number of genes linked, and size of spacers between ribbons
    a.1.cg = sum(lengths(ribbonlist[[1]][1:8]))
    a.2.cg = sum(lengths(ribbonlist[[2]][1:8]))
    a.3.cg = sum(lengths(ribbonlist[[3]][1:8]))
    a.4.cg = sum(lengths(ribbonlist[[4]][1:8]))
    a.5.cg = sum(lengths(ribbonlist[[5]][1:8]))
    a.6.cg = sum(lengths(ribbonlist[[6]][1:8]))
    
    r.t = function(x){return(as.numeric(lengths(x) != 0))} # calulates if there are any links present requiring ribbons
    
    a.1.r = sum(r.t(ribbonlist[[1]])[1:8])
    a.2.r = sum(r.t(ribbonlist[[2]])[1:8])
    a.3.r = sum(r.t(ribbonlist[[3]])[1:8])
    a.4.r = sum(r.t(ribbonlist[[4]])[1:8])
    a.5.r = sum(r.t(ribbonlist[[5]])[1:8])
    a.6.r = sum(r.t(ribbonlist[[6]])[1:8])

    
    

    
    spacing.table = data.frame(Chrsize = .chr_length, 
                               crossgenes = c(a.1.cg, a.2.cg, a.3.cg, a.4.cg, a.5.cg, a.6.cg),
                               ribbons = c(a.1.r, a.2.r, a.3.r, a.4.r, a.5.r, a.6.r))
    spacing.table$ribbon.gap = floor((spacing.table$Chrsize - spacing.table$crossgenes) 
                                     / (spacing.table$ribbons))
    
    NoGapKaryotype = .karyotype[c(2,3,5,6,8,9),]
    
    chr.coordinates = list()
    for(i in 1:(.Comparison * .Condition)){
        dat = data.frame(chr = character(), start = integer(), end = integer(), ribbon = integer(), stringsAsFactors=FALSE)
        start = NoGapKaryotype[i,5]
        for (k in 1:8){
            dat[k,1] = sheet_list[i]
            dat[k,2] = start
            dat[k,3] = start + length(ribbonlist[[i]][[k]])
            dat[k,4] = r.t(ribbonlist[[i]])[k]
            start = dat[k,3]
            if (dat[k,4] == 1){start = start + spacing.table[i,4]}
        }
        chr.coordinates[[i]] = dat
    }
    

    sink(file = "data/Ribbon1.txt")
    if(chr.coordinates[[1]][1,4] == 1){
        cat(paste(chr.coordinates[[1]][1,1], chr.coordinates[[1]][1,2], chr.coordinates[[1]][1,3],
                  chr.coordinates[[6]][8,1], chr.coordinates[[6]][8,2], chr.coordinates[[6]][8,3],
                  "\n", sep = "\t"))}
    if(chr.coordinates[[1]][4,4] == 1){
        cat(paste(chr.coordinates[[1]][4,1], chr.coordinates[[1]][4,2], chr.coordinates[[1]][4,3],
                  chr.coordinates[[5]][8,1], chr.coordinates[[5]][8,2], chr.coordinates[[5]][8,3],
                  "\n", sep = "\t"))}
    if(chr.coordinates[[1]][5,4] == 1){
        cat(paste(chr.coordinates[[1]][5,1], chr.coordinates[[1]][5,2], chr.coordinates[[1]][5,3],
                  chr.coordinates[[4]][4,1], chr.coordinates[[4]][4,2], chr.coordinates[[4]][4,3],
                  "\n", sep = "\t"))}
    if(chr.coordinates[[1]][8,4] == 1){
        cat(paste(chr.coordinates[[1]][8,1], chr.coordinates[[1]][8,2], chr.coordinates[[1]][8,3],
                  chr.coordinates[[3]][4,1], chr.coordinates[[3]][4,2], chr.coordinates[[3]][4,3],
                  "\n", sep = "\t"))}
    
    if(chr.coordinates[[2]][1,4] == 1){
        cat(paste(chr.coordinates[[2]][1,1], chr.coordinates[[2]][1,2], chr.coordinates[[2]][1,3],
                  chr.coordinates[[6]][5,1], chr.coordinates[[6]][5,2], chr.coordinates[[6]][5,3],
                  "\n", sep = "\t"))}
    if(chr.coordinates[[2]][4,4] == 1){
        cat(paste(chr.coordinates[[2]][4,1], chr.coordinates[[2]][4,2], chr.coordinates[[2]][4,3],
                  chr.coordinates[[5]][5,1], chr.coordinates[[5]][5,2], chr.coordinates[[5]][5,3],
                  "\n", sep = "\t"))}
    if(chr.coordinates[[2]][5,4] == 1){
        cat(paste(chr.coordinates[[2]][5,1], chr.coordinates[[2]][5,2], chr.coordinates[[2]][5,3],
                  chr.coordinates[[4]][1,1], chr.coordinates[[4]][1,2], chr.coordinates[[4]][1,3],
                  "\n", sep = "\t"))}
    if(chr.coordinates[[2]][8,4] == 1){
        cat(paste(chr.coordinates[[2]][8,1], chr.coordinates[[2]][8,2], chr.coordinates[[2]][8,3],
                  chr.coordinates[[3]][1,1], chr.coordinates[[3]][1,2], chr.coordinates[[3]][1,3],
                  "\n", sep = "\t"))}
    
    if(chr.coordinates[[3]][5,4] == 1){
        cat(paste(chr.coordinates[[3]][5,1], chr.coordinates[[3]][5,2], chr.coordinates[[3]][5,3],
                  chr.coordinates[[6]][4,1], chr.coordinates[[6]][4,2], chr.coordinates[[6]][4,3],
                  "\n", sep = "\t"))}
    if(chr.coordinates[[3]][8,4] == 1){
        cat(paste(chr.coordinates[[3]][8,1], chr.coordinates[[3]][8,2], chr.coordinates[[3]][8,3],
                  chr.coordinates[[5]][4,1], chr.coordinates[[5]][4,2], chr.coordinates[[5]][4,3],
                  "\n", sep = "\t"))}
    if(chr.coordinates[[4]][5,4] == 1){
        cat(paste(chr.coordinates[[4]][5,1], chr.coordinates[[4]][5,2], chr.coordinates[[4]][5,3],
                  chr.coordinates[[6]][1,1], chr.coordinates[[6]][1,2], chr.coordinates[[6]][1,3],
                  "\n", sep = "\t"))}
    if(chr.coordinates[[4]][8,4] == 1){
        cat(paste(chr.coordinates[[4]][8,1], chr.coordinates[[4]][8,2], chr.coordinates[[4]][8,3],
                  chr.coordinates[[5]][1,1], chr.coordinates[[5]][1,2], chr.coordinates[[5]][1,3],
                  "\n", sep = "\t"))}
    sink(NULL)
    
    
    
    sink(file = "data/Ribbon2.txt")
    if(chr.coordinates[[1]][2,4] == 1){
        cat(paste(chr.coordinates[[1]][2,1], chr.coordinates[[1]][2,2], chr.coordinates[[1]][2,3],
                  chr.coordinates[[4]][3,1], chr.coordinates[[4]][3,2], chr.coordinates[[4]][3,3],
                  "\n", sep = "\t"))
        cat(paste(chr.coordinates[[1]][2,1], chr.coordinates[[1]][2,2], chr.coordinates[[1]][2,3],
                  chr.coordinates[[6]][7,1], chr.coordinates[[6]][7,2], chr.coordinates[[6]][7,3],
                  "\n", sep = "\t"))
        cat(paste(chr.coordinates[[4]][3,1], chr.coordinates[[4]][3,2], chr.coordinates[[4]][3,3],
                  chr.coordinates[[6]][7,1], chr.coordinates[[6]][7,2], chr.coordinates[[6]][7,3],
                  "\n", sep = "\t"))
        }
    
    if(chr.coordinates[[1]][3,4] == 1){
        cat(paste(chr.coordinates[[1]][3,1], chr.coordinates[[1]][3,2], chr.coordinates[[1]][3,3],
                  chr.coordinates[[3]][6,1], chr.coordinates[[3]][6,2], chr.coordinates[[3]][6,3],
                  "\n", sep = "\t"))
        cat(paste(chr.coordinates[[1]][3,1], chr.coordinates[[1]][3,2], chr.coordinates[[1]][3,3],
                  chr.coordinates[[6]][6,1], chr.coordinates[[6]][6,2], chr.coordinates[[6]][6,3],
                  "\n", sep = "\t"))
        cat(paste(chr.coordinates[[3]][6,1], chr.coordinates[[3]][6,2], chr.coordinates[[3]][6,3],
                  chr.coordinates[[6]][6,1], chr.coordinates[[6]][6,2], chr.coordinates[[6]][6,3],
                  "\n", sep = "\t"))
        }
    
    if(chr.coordinates[[1]][6,4] == 1){
        cat(paste(chr.coordinates[[1]][6,1], chr.coordinates[[1]][6,2], chr.coordinates[[1]][6,3],
                  chr.coordinates[[4]][6,1], chr.coordinates[[4]][6,2], chr.coordinates[[4]][6,3],
                  "\n", sep = "\t"))
        cat(paste(chr.coordinates[[1]][6,1], chr.coordinates[[1]][6,2], chr.coordinates[[1]][6,3],
                  chr.coordinates[[5]][3,1], chr.coordinates[[5]][3,2], chr.coordinates[[5]][3,3],
                  "\n", sep = "\t"))
        cat(paste(chr.coordinates[[4]][6,1], chr.coordinates[[4]][6,2], chr.coordinates[[4]][6,3],
                  chr.coordinates[[5]][3,1], chr.coordinates[[5]][3,2], chr.coordinates[[5]][3,3],
                  "\n", sep = "\t"))
        }
    
    if(chr.coordinates[[1]][7,4] == 1){
        cat(paste(chr.coordinates[[1]][7,1], chr.coordinates[[1]][7,2], chr.coordinates[[1]][7,3],
                  chr.coordinates[[3]][7,1], chr.coordinates[[3]][7,2], chr.coordinates[[3]][7,3],
                  "\n", sep = "\t"))
        cat(paste(chr.coordinates[[1]][7,1], chr.coordinates[[1]][7,2], chr.coordinates[[1]][7,3],
                  chr.coordinates[[5]][7,1], chr.coordinates[[5]][7,2], chr.coordinates[[5]][7,3],
                  "\n", sep = "\t"))
        cat(paste(chr.coordinates[[3]][7,1], chr.coordinates[[5]][7,2], chr.coordinates[[5]][7,3],
                  chr.coordinates[[5]][7,1], chr.coordinates[[3]][7,2], chr.coordinates[[3]][7,3],
                  "\n", sep = "\t"))
        }
    
    if(chr.coordinates[[2]][2,4] == 1){
        cat(paste(chr.coordinates[[2]][2,1], chr.coordinates[[2]][2,2], chr.coordinates[[2]][2,3],
                  chr.coordinates[[4]][2,1], chr.coordinates[[4]][2,2], chr.coordinates[[4]][2,3],
                  "\n", sep = "\t"))
        cat(paste(chr.coordinates[[2]][2,1], chr.coordinates[[2]][2,2], chr.coordinates[[2]][2,3],
                  chr.coordinates[[6]][2,1], chr.coordinates[[6]][2,2], chr.coordinates[[6]][2,3],
                  "\n", sep = "\t"))
        cat(paste(chr.coordinates[[4]][2,1], chr.coordinates[[4]][2,2], chr.coordinates[[4]][2,3],
                  chr.coordinates[[6]][2,1], chr.coordinates[[6]][2,2], chr.coordinates[[6]][2,3],
                  "\n", sep = "\t"))
        }
    
    if(chr.coordinates[[2]][3,4] == 1){
        cat(paste(chr.coordinates[[2]][3,1], chr.coordinates[[2]][3,2], chr.coordinates[[2]][3,3],
                  chr.coordinates[[4]][7,1], chr.coordinates[[4]][7,2], chr.coordinates[[4]][7,3],
                  "\n", sep = "\t"))
        cat(paste(chr.coordinates[[2]][3,1], chr.coordinates[[2]][3,2], chr.coordinates[[2]][3,3],
                  chr.coordinates[[5]][2,1], chr.coordinates[[5]][2,2], chr.coordinates[[5]][2,3],
                  "\n", sep = "\t"))
        cat(paste(chr.coordinates[[5]][2,1], chr.coordinates[[5]][2,2], chr.coordinates[[5]][2,3],
                  chr.coordinates[[4]][7,1], chr.coordinates[[4]][7,2], chr.coordinates[[4]][7,3],
                  "\n", sep = "\t"))
        }
    
    if(chr.coordinates[[2]][6,4] == 1){
        cat(paste(chr.coordinates[[2]][6,1], chr.coordinates[[2]][6,2], chr.coordinates[[2]][6,3],
                  chr.coordinates[[3]][3,1], chr.coordinates[[3]][3,2], chr.coordinates[[3]][3,3],
                  "\n", sep = "\t"))
        cat(paste(chr.coordinates[[2]][6,1], chr.coordinates[[2]][6,2], chr.coordinates[[2]][6,3],
                  chr.coordinates[[6]][3,1], chr.coordinates[[6]][3,2], chr.coordinates[[6]][3,3],
                  "\n", sep = "\t"))
        cat(paste(chr.coordinates[[3]][3,1], chr.coordinates[[3]][3,2], chr.coordinates[[3]][3,3],
                  chr.coordinates[[6]][3,1], chr.coordinates[[6]][3,2], chr.coordinates[[6]][3,3],
                  "\n", sep = "\t"))
        }
    
    if(chr.coordinates[[2]][7,4] == 1){
        cat(paste(chr.coordinates[[2]][7,1], chr.coordinates[[2]][7,2], chr.coordinates[[2]][7,3],
                  chr.coordinates[[3]][2,1], chr.coordinates[[3]][2,2], chr.coordinates[[3]][2,3],
                  "\n", sep = "\t"))
        cat(paste(chr.coordinates[[2]][7,1], chr.coordinates[[2]][7,2], chr.coordinates[[2]][7,3],
                  chr.coordinates[[5]][6,1], chr.coordinates[[5]][6,2], chr.coordinates[[5]][6,3],
                  "\n", sep = "\t"))
        cat(paste(chr.coordinates[[3]][2,1], chr.coordinates[[3]][2,2], chr.coordinates[[3]][2,3],
                  chr.coordinates[[5]][6,1], chr.coordinates[[5]][6,2], chr.coordinates[[5]][6,3],
                  "\n", sep = "\t"))
        }
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
    
    ribbonsizes2 = data.frame(Chr1 = character(), Chr2 = character(), Chr3 = character(),
                              size = integer(), stringsAsFactors=FALSE)
    tmp = read.csv("data/Ribbon2.txt", header = FALSE, sep = "\t") 
    for (i in 1:(nrow(tmp) / 3)){
        ribbonsizes2[i,1] = as.character(tmp[i,1])
        ribbonsizes2[i,2] = as.character(tmp[i,4])
        ribbonsizes2[i,2] = as.character(tmp[(i+1),4])
        ribbonsizes2[i,3] = tmp[i,3] - tmp[i,2]
    }    
    addWorksheet(wb, "Ribbon2", gridLines = TRUE)
    writeData(wb, sheet = 1, ribbonsizes2, rowNames = FALSE)
    
    saveWorkbook(wb, "RibbonSizes.xlsx", overwrite = TRUE)
    
}
