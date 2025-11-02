library(devtools)
devtools::document()
devtools::load_all()
conflicts <- intersect(ls(), ls("package:varmapack"))
if (length(conflicts) > 0) {
  cat("Conflicts detected:\n")
  print(conflicts)
  cat("You may need to remove these objects:\n")
  print(conflicts)
}
cat("Package 'varmapack' loaded successfully for development.\n")
