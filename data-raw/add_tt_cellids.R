### Adding string vector of cell line ids for tumour type (haem vs solid) addition.
load("C:/Users/hfy163/OneDrive - Queen Mary, University of London/Documents/ReMEA_dev_Phase2/Linear_model_scripting/haem_cells.rda")
load("C:/Users/hfy163/OneDrive - Queen Mary, University of London/Documents/ReMEA_dev_Phase2/Linear_model_scripting/solid_cells.rda")

require("devtools")
usethis::use_data(haem.cells, overwrite = TRUE)
usethis::use_data(solid.cells, overwrite = TRUE)

devtools::build()
