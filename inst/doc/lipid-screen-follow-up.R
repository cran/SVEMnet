## ----worked-results, echo=FALSE, results='asis'-------------------------------
snapshot <- system.file("examples", "lipid-screen-follow-up.md", package = "SVEMnet",
                        mustWork = TRUE)
cat(readLines(snapshot, encoding = "UTF-8"), sep = "\n")

