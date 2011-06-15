write.bedGraph <- function (x, raw = TRUE, vals = TRUE, directory = ".") 
{
  for (i in 1:ncol(x)) {
    if (class(x) == "ACMECalcSet") {
      if (vals) {
        filename <- file.path(directory, sprintf("%s_thresh%3.2f.bed", 
                                                 sampleNames(x)[i], threshold(x)))
        cat(filename, "\n")
        f = file(filename,'w')
        writeLines(sprintf('track type=bedGraph name=%s description=%s',
                             basename(filename),basename(filename)),con=f)
        write.table(data.frame(chromosome = chromosome(x), 
                               start = start(x), end = end(x),
                               -log10(vals(x)[, i] + min(vals(x)[vals(x)[,i] > 0, i]))),
                    file = f, sep = "\t", 
                    col.names = FALSE, row.names = FALSE, quote = FALSE)
        close(f)
      }
    }
    if (raw) {
      filename <- file.path(directory, sprintf("%s_raw.bed", 
                                               sampleNames(x)[i]))
      cat(filename, "\n")
      f = file(filename,'w')
      writeLines(sprintf('track type=bedGraph name=%s description=%s',
                         filename,filename),con=f)
      write.table(data.frame(chromosome = chromosome(x), 
                             start = start(x), end=end(x),
                             exprs(x)[, i]),
                  file = f, 
                  sep = "\t", col.names = FALSE, row.names = FALSE, 
                  quote = FALSE)
      close(f)
    }
  }
}
