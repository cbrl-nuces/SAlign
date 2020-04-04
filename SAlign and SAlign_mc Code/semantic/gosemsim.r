# CHECK THE FILE PATHS BEFORE RUNNING THIS SCRIPT

library(GOSemSim)

terms <- c('MF', 'BP')

sp1 <- net1
sp2 <- net2

for(term_i in 1:2) # MF and BP
{ 
  term <- terms[term_i]
  goterm <- tolower(term)
  align_file <- paste(c("1-net1.txt_net2.txt.alignment"),collapse="") # INPUT FILE
  alignment_file <- scan(align_file, what="character", sep='\n',quiet=TRUE)
  len <- length(alignment_file)
  d <- godata('org.Hs.eg.db', ont=term, computeIC=FALSE)
  fileConn <- paste(c("1-net1_net2_", goterm, ".txt") ,collapse="") # OUTPUT FILE

  path1 <- paste(c('go-ids/', sp1) ,collapse="")
  path2 <- paste(c('go-ids/', sp2) ,collapse="")

  match1 <- paste(c(path1, '/', goterm, '-go.txt'), collapse='')
  match2 <- paste(c(path2, '/', goterm, '-go.txt'), collapse='')
  main_stream <- ''
  
  for(i in (1:len))
  {
    local_stream <- ''
    p<- strsplit(alignment_file[i], ' ')
    p1 <- p[[1]][1]
    p2 <- p[[1]][2]
    line1 <- grep(p1, readLines(match1), value = TRUE)
    line2 <- grep(p2, readLines(match2), value = TRUE)
    val1 <- unlist(strsplit(line1, '\t'))
    val2 <- unlist(strsplit(line2, '\t'))
    if(is.null(val1[2]))
    {
      print(c(paste(p1, "\t", sp1, "\t", "skipped")))
      next
    }
    if(is.null(val2[2])) {
      print(c(paste(p2, "\t", sp2, "\t", "skipped")))
      next
    }
    go1 <- strsplit(val1[2], ',')
    go2 <- strsplit(val2[2], ',')
    sm <- 0
    sm <- mgoSim(go1[[1]], go2[[1]], semData=d, measure="Wang", combine="BMA")
    local_stream <- c(paste(p1, "\t", p2, "\t", sm, sep = ""))
    main_stream <- c(trimws(paste(main_stream,'\n', local_stream, sep = "")))
  }
  writeLines(main_stream,fileConn) 
}
print(c(paste("DONE WITH SEMANTIC CALCULATIONS.")))

