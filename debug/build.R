# Steps to create/check/install package from directory
rm(list = ls())
bb = strsplit(getwd(),"/")[[1]]
pack_dir = paste(bb[-length(bb)],collapse = "/")
pack = strsplit(pack_dir,"/")[[1]]
pack = pack[length(pack)]
pack

if( pack %in% installed.packages()[,1] ){
	remove.packages(pack)
	q("no")
}

Rcpp::compileAttributes(pkgdir = pack_dir)
devtools::document(pkg = pack_dir)
usethis::use_gpl3_license()
devtools::check(pkg = pack_dir,
	manual = TRUE,cran = TRUE,
	error_on = "note")
devtools::install(pack_dir)

