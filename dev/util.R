vigview <- function() {
  vignette <- 'varmasim-intro.Rmd'
  render(here('vignettes', vignette, '.Rmd'), output_dir=here('doc'))
  browseURL(here('doc', vignette))
}
