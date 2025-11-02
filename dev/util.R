vigview <- function() {
  vignette <- 'varmapack-intro.Rmd'
  render(here('vignettes', vignette, '.Rmd'), output_dir=here('doc'))
  browseURL(here('doc', vignette))
}
