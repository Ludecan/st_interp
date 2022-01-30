#install.packages(c("devtools", "formatR"))
#devtools::install_github("ontherunvaro/formatFileAddin")

options(
  formatR.comment=TRUE,
  formatR.blank=TRUE,
  formatR.arrow=FALSE,
  formatR.brace.newline=FALSE,
  formatR.indent=4,
  formatR.wrap=FALSE,
  formatR.width=160,
  formatR.args.newline=TRUE
)

getOption('formatR.wrap')
