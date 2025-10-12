library(devtools)
devtools::document()
devtools::load_all()
conflicts <- intersect(ls(), ls("package:varmasim"))
if (length(conflicts) > 0) {
  cat("Conflicts detected:\n")
  print(conflicts)
  cat("You may need to remove these objects:\n")
  print(conflicts)
}
cat("Package 'varmasim' loaded successfully for development.\n")
