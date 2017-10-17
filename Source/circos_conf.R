#' Circos.conf
#' 
#' Create the circos.conf file required for generating the plots. This is a generic function that creates file for all circos comparisons, as the perl script doesn't care if the files are empty, so long as they exist.
#' @param TileColours Pass a list of colours defining the tile colour of the conditions within the comparison. Defaults to c("blue", "yellow", "red").
#' @param TileThickness Defines the thickness of the tiles
#' @param RibbonColours Pass a list of colours for ribbons connecting genes, Defaults to c("grey_a1", "dred_a1").
#' 



circos_conf = function(TileColours = c("blue", "yellow", "red"), TileThickness = 50, RibbonColours = c("grey_a1", "dred_a1")){
    head.conf = paste("# circos.conf

karyotype = data/karyotype.txt
chromosome_units = 1
chromosomes_display_default = yes
                    
<ideogram>
    <spacing>
        default = 0.000r
    </spacing>
    radius    = 0.9r
    thickness = 20p
    fill      = yes
</ideogram>
                    
<image>
    # Included from Circos distribution.
    <<include ../etc/image.conf>>
</image>
                    
# RGB/HSV color definitions, color lists, location of fonts, fill patterns.
# Included from Circos distribution.
<<include etc/colors_fonts_patterns.conf>>
                    
# Debugging, I/O an dother system parameters
# Included from Circos distribution.
<<include etc/housekeeping.conf>>\n\n")

    plot = list()
    for (i in 1:.Condition ){
        plot[i] = paste("    <plot>
        file = data/Tile",i,".txt
        color = ", TileColours[i],"
        thickness = ",TileThickness,"p
        r1 = 0.95r
        r0 = 0.95r
    </plot>\n", sep = "")
    }   
    plots = paste(plot, collapse = '')
    plots.conf = paste("<plots>\n    type = tile\n", plots, "</plots>\n", sep = "")                   
   
    
    link = list()
    for (i in 1:(.Comparison - 1)){
        link[i] = paste("    <link>
        file = data/Ribbon",i,".txt
        ribbon = yes
        flat = yes
        color = ", RibbonColours[i],"
        radius = 0.92r
    </link>\n", sep = "")
    }   
    links = paste(link, collapse = '')
    links.conf = paste("<links>\n", links, "</links>\n", sep = "")   
            


    conf.file <- paste(head.conf, plots.conf, links.conf)

    write.table(conf.file, "circos.conf", quote = FALSE, col.names = FALSE, row.names = FALSE)
}
