#' create_runfile
#' 
#' Creates the runfile for perl script
#' @export

create_runfile = function(){
    sink("run")
    cat("#!/bin/bash
        
../bin/circos -conf circos.conf -debug_group summary,timer > run.out")
    sink(NULL)
    Sys.chmod("run", mode = "0775")
}
