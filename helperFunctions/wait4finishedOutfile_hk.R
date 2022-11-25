wait4finishedOutfile <- function (outfile_fn, waitingtime_seconds=600, erkennungsmarke = "TOTAL TIME :") {
  
  #outfile_fn = "07_a_coloc_get_eQTL_data.R.out" 
  
  check = F
  while(check ==F) {
    outfile = try(suppressWarnings(readLines(outfile_fn)))
    if(class(outfile)=="try-error") stop("File \n'", outfile_fn, "'\n does not exist or is not readable...stopping....")
    check = any(grepl(erkennungsmarke, tail(outfile, 10)))
    print(tail(outfile, 6))
    if(check==T) message("Skript producing file '", outfile_fn, "' seems to be SUCCESFULLY FINISHED...") else {
      message("Skript producing file '", outfile_fn, "' seems NOT TO BE FINISHED at ", Sys.time(), " ...waiting...")
      Sys.sleep(time = waitingtime_seconds )
      
    }
  }
}
