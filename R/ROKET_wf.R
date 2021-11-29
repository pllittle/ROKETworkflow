# ----------
# ROKET workflow functions
# ----------

# ----------
# Minor functions
# ----------
setdirs = function(git_dir,work_dir,dataset,verbose = TRUE){
	
	# Make directories
	if( verbose ) cat(sprintf("%s: Make directories ...\n",date()))
	proj_dir 	= file.path(work_dir,"kSMASH"); smart_mkdir(proj_dir)
	ref_dir 	= file.path(proj_dir,"REF"); 		smart_mkdir(ref_dir)
	panc_dir 	= file.path(proj_dir,"PanCan"); smart_mkdir(panc_dir)
	type_dir 	= file.path(proj_dir,dataset); 	smart_mkdir(type_dir)
	spms_dir 	= file.path(type_dir,"SPMs"); 	smart_mkdir(spms_dir)
	out_dir 	= file.path(type_dir,"output"); smart_mkdir(out_dir)
	ithOT_dir	= file.path(type_dir,"ITH_OT"); smart_mkdir(ithOT_dir)
	reg_dir		= file.path(out_dir,"regress"); smart_mkdir(reg_dir)
	fin_dir		= file.path(proj_dir,"final");	smart_mkdir(fin_dir)
	
	# Output
	if( verbose ) cat(sprintf("%s: Output ...\n",date()))
	list(git_dir = git_dir,proj_dir = proj_dir,ref_dir = ref_dir,
		panc_dir = panc_dir,type_dir = type_dir,spms_dir = spms_dir,
		dataset = dataset,ithOT_dir = ithOT_dir,out_dir = out_dir,
		reg_dir = reg_dir,fin_dir = fin_dir)
	
}
make_sampID = function(origIDs){
	cat(sprintf("%s: Making universal ID ...\n",date()))
	sapply(origIDs,function(xx) 
		paste(strsplit(xx,"-")[[1]][1:3],collapse="-"),USE.NAMES=FALSE)
}
smart_DTR = function(...){
	data.table::fread(...,
		data.table = FALSE)
}
make_DP_VAF = function(DATA){
	req_names = c("t_alt_count","t_ref_count")
	smart_reqNames(DATA = DATA,REQ = req_names)
	
	DATA$tDP = DATA$t_alt_count + DATA$t_ref_count
	DATA$tVAF = DATA$t_alt_count / DATA$tDP
	DATA$tVAF[DATA$tDP == 0] = 0
	DATA
}
name_TUMORTYPES = function(){
	out = c()
	out = rbind(out,smart_df(V1 = "BLCA",V2 = "Bladder Urothelial Carcinoma"))
	out = rbind(out,smart_df(V1 = "BRCA",V2 = "Breast invasive carcinoma"))
	out = rbind(out,smart_df(V1 = "CESC",V2 = "Cervical squamous cell carcinoma and endocervical adenocarcinoma"))
	out = rbind(out,smart_df(V1 = "COAD",V2 = "Colon and Rectum adenocarcinoma"))
	out = rbind(out,smart_df(V1 = "ESCA",V2 = "Esophageal carcinoma"))
	out = rbind(out,smart_df(V1 = "GBM",V2 = "Glioblastoma multiforme"))
	out = rbind(out,smart_df(V1 = "HNSC",V2 = "Head and Neck squamous cell carcinoma"))
	out = rbind(out,smart_df(V1 = "KIRC",V2 = "Kidney renal clear cell carcinoma"))
	out = rbind(out,smart_df(V1 = "LGG",V2 = "Brain Lower Grade Glioma"))
	out = rbind(out,smart_df(V1 = "LIHC",V2 = "Liver hepatocellular carcinoma"))
	out = rbind(out,smart_df(V1 = "LUAD",V2 = "Lung adenocarcinoma"))
	out = rbind(out,smart_df(V1 = "LUSC",V2 = "Lung squamous cell carcinoma"))
	out = rbind(out,smart_df(V1 = "PAAD",V2 = "Pancreatic adenocarcinoma"))
	out = rbind(out,smart_df(V1 = "SARC",V2 = "Sarcoma"))
	out = rbind(out,smart_df(V1 = "SKCM",V2 = "Skin Cutaneous Melanoma"))
	out = rbind(out,smart_df(V1 = "STAD",V2 = "Stomach adenocarcinoma"))
	out = rbind(out,smart_df(V1 = "UCEC",V2 = "Uterine Corpus Endometrial Carcinoma"))
	
	return(out)
}


# ----------
# Check/Download files
# ----------






#' @importFrom smartr smart_reqNames
NULL

# Steps to create/check/install package from directory
# bb = strsplit(getwd(),"/")[[1]]; pack_dir = paste(bb[-length(bb)],collapse = "/")
# pack = strsplit(pack_dir,"/")[[1]]; pack = pack[length(pack)]
# if( pack %in% installed.packages()[,1] ){ remove.packages(pack); q("no")}
# Rcpp::compileAttributes(pkgdir = pack_dir)
# devtools::document(pkg = pack_dir); usethis::use_gpl3_license()
# devtools::check(pkg = pack_dir,manual = TRUE,cran = TRUE,error_on = "note")
# devtools::install(pack_dir)

###

