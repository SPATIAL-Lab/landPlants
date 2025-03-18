# Read base model and write to tempdir
basemod = readLines("code/models/forwardFranksMulti.R")

writeLines(basemode, file.path(tempdir(), "fullFranks.txt"))

# Find and remove likelihood expressions for Pl and l, write
rl = c(grep("GCLab[i, 1]*", basemod), grep("GCWab[i, 1]*", basemod))
Donly = basemod[-rl]

writeLines(Donly, file.path(tempdir(), "DonlyFranks.txt"))

# Find and remove likelihood expression for D, write
rl = grep("Dab[i, 1]*", Donly)
d13Conly = Donly[-rl]

writeLines(d13Conly, file.path(tempdir(), "d13ConlyFranks.txt"))
