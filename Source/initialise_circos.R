#' Initialise_Circos
#' 
#' Initial parameters for creating files.
#' @param Comparison Number of comparisons being made, This is the experimental design (e.g. treated vs control). Maximum of 3 for grouped 
#' @param Condition Number of conditions within each comparison e.g. Up- & Down- regulated genes.
#' @param Spacing Do you wish to alter the size of the condition segment of the chromosome. If so, takes the largest comparison segment and makes all segment this size. For 2 conditions both are joined together, for 3 way comparisons will centre the chromosome. Default is FALSE
#' @param Gap Size of the gap between chromosomes, reccomended to be even and approximately 10% of largest segment. Default is 10
#' @examples 
#' initialise_circos(Comparison = 3, Condition = 2, Spacing = FALSE, Gap = 10)
#' initialise_circos(Comparison = 3, Condition = 2)
#' @export

initialise_circos = function(Comparison = integer(), Condition = integer(), Spacing = FALSE, Gap = 10){
    if (!dir.exists("data")){
        dir.create("data")
    } 
    if (length(sheet_list) != (Comparison * Condition)){
        cat("The number of defined Comparisons & Conditions does not match the number of sheets")
        cat(paste("Number of sheets:", length(sheet_list), "Number of Conditions", (Comparison*Condition)))
    } else {
        cat(paste("\nThe number of sheets in the workbook (", length(sheet_list), ") matches what is expepected using the defined comparisons (",Comparison, ") & Conditions (", Condition, ")\n\n", sep = ""))
        .Comparison <<- Comparison
        .Condition <<- Condition
        .gap <<- Gap
        .Spacing <<- Spacing
    }
    if (Comparison == 2){reverse_lists()}
    create_karyotype()
    create_tile()
    if (Comparison == 2 && Condition == 1){create_2x1_links()}
    if (Comparison == 2 && Condition == 2){create_2x2_links()}
    if (Comparison == 2 && Condition == 3){create_2x3_links()}
    if (Comparison == 3 && Condition == 1){create_3x1_links()}
    if (Comparison == 3 && Condition == 2){create_3x2_links()}
    circos_conf()
    create_runfile()
}
