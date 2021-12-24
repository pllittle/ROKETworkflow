# ----------
# ROKET workflow functions
# ----------

# ----------
# Minor functions
# ----------

#' @title setdirs
#' @param git_dir A string containing the full
#'	path to the parent directory containing
#'	the cloned ROKETworkflow repository
#' @param work_dir A string containing the full
#'	path to the working directory for downloading,
#'	processing, and analyzing real data
#' @param dataset A string for the abbreviated
#'	cancer type such as BRCA, LGG, STAD, etc.
#' @param verbose Boolean set to TRUE to display
#'	verbose output
#' @export
setdirs = function(git_dir,work_dir,dataset,verbose = TRUE){
	
	# Make directories
	if( verbose ) cat(sprintf("%s: Make directories ...\n",date()))
	if( !dir.exists(work_dir) ) stop(sprintf("Make working directory = %s",work_dir))
	proj_dir 	= file.path(work_dir,"kSMASH"); smart_mkdir(proj_dir)
	ref_dir 	= file.path(proj_dir,"REF"); 		smart_mkdir(ref_dir)
	panc_dir 	= file.path(proj_dir,"PanCan"); smart_mkdir(panc_dir)
	type_dir 	= file.path(proj_dir,dataset); 	smart_mkdir(type_dir)
	spms_dir 	= file.path(type_dir,"SPMs"); 	smart_mkdir(spms_dir)
	out_dir 	= file.path(type_dir,"output"); smart_mkdir(out_dir)
	ithOT_dir	= file.path(type_dir,"OT"); smart_mkdir(ithOT_dir)
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
# Real Data: Check/Download files
# ----------
get_REF = function(my_dirs){
	
	# GDC's hg38 GTF file
	gtf_fn = file.path(my_dirs$ref_dir,"gencode.gene.info.v22.tsv")
	if( !file.exists(gtf_fn) ){
		url_link = "https://api.gdc.cancer.gov/data/b011ee3e-14d8-4a97-aed4-e0b10f6bbe82"
		stop(sprintf("Download %s and rename %s",url_link,gtf_fn))
	}
	
	# Cytolytic activity
	url_link = "https://www.cell.com/cms/10.1016"
	url_link = sprintf("%s/j.cell.2014.12.033/attachment",url_link)
	url_link = sprintf("%s/fdbc71ba-1786-469f-9a7b-6c05e2ab76c4/mmc1.xlsx")
	cyt_fn = file.path(my_dirs$panc_dir,"mmc1 CYT.xlsx")
	if( !file.exists(cyt_fn) ){
		stop(sprintf("Download %s and rename %s",url_link,cyt_fn))
	}
	
}
down_PanCan = function(my_dirs,getCNA = FALSE,
	samp = NULL,override = FALSE,show = TRUE){
	
	get_REF(my_dirs = my_dirs)
	dataset = my_dirs$dataset
	
	# Source: https://gdc.cancer.gov/about-data/publications/pancanatlas
	
	# Download absolute purity/ploidy from TCGA PanCan
	if( show ) cat(sprintf("%s: Get Absolute purity/ploidy ...\n",date()))
	absPP_fn = file.path(my_dirs$panc_dir,"TCGA_mastercalls.abs_tables_JSedit.fixed.txt")
	if( !file.exists(absPP_fn) ){
		url_link = "http://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5"
		stop(sprintf("Download %s and rename %s",url_link,absPP_fn))
	}
	
	# Download absolute copy number from TCGA PanCan
	if( show ) cat(sprintf("%s: Get Absolute copy number ...\n",date()))
	absCN_fn = file.path(my_dirs$panc_dir,"TCGA_mastercalls.abs_segtabs.fixed.txt")
	absCN_fn2 = sprintf("%s.gz",absCN_fn)
	if( !file.exists(absCN_fn2) ){
		url_link = "http://api.gdc.cancer.gov/data/0f4f5701-7b61-41ae-bda9-2805d1ca9781"
		stop(sprintf("Download %s, unzip and rename %s",url_link,absCN_fn))
	}
	absCN_fn = absCN_fn2
	
	# Download TCGA PanCan clinical outcomes
	if( show ) cat(sprintf("%s: Get TCGA PanCan clinical outcomes ...\n",date()))
	clin_fn = file.path(my_dirs$panc_dir,"TCGA-CDR-SupplementalTableS1.xlsx")
	if( !file.exists(clin_fn) ){
		url_link = "https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81"
		stop(sprintf("Download %s and rename %s",url_link,clin_fn))
	}
	
	# Import/Polish clin data
	if( show ) cat(sprintf("%s: Import/Prep clinical data ...\n",date()))
	clin = suppressWarnings(readxl::read_excel(clin_fn,sheet = "TCGA-CDR"))
	clin = as.data.frame(clin,stringsAsFactors = FALSE)
	clin = clin[,names(clin) != "...1"]
	clin = name_change(clin,"bcr_patient_barcode","ID")
	for(vv in names(clin)){
		idx = which(clin[,vv] %in% c("[Not Available]","[Not Applicable]"))
		if( length(idx) > 0 ) clin[idx,vv] = NA
	}
	# dim(clin); clin[1:3,]
	
	# Subset tumor type (clin, copy number, and SPMs)
	if( show ) cat(sprintf("%s: Subset tumortype = %s ...\n",date(),dataset))
	smart_table(clin$type)
	if( dataset == "COAD" ){
		clin = clin[which(clin$type %in% c("COAD","READ")),]
	} else if( dataset %in% c("BLCA","BRCA","CESC","ESCA","GBM","HNSC","KIRC",
			"LAML","LGG","LIHC","LUAD","LUSC","OV","PAAD","SARC","SKCM",
			"STAD","UCEC") ){
		clin = clin[which(clin$type == dataset),]
	} else {
		stop("Add code in down_PanCan()")
	}
	# dim(clin); clin[1:3,]
	
	if( show ) cat(sprintf("%s: Import ABSOLUTE purity/ploidy data ...\n",date()))
	absPP = smart_DTR(absPP_fn,sep = "\t",header = TRUE)
	absPP$ID = make_sampID(origIDs = absPP$array)
	absPP = absPP[which(absPP$ID %in% clin$ID),]
	
	# remove samples with multiple CN
	if( show ) cat(sprintf("%s: Remove sample duplicates with lower purity ...\n",date()))
	tab = table(absPP$ID)
	tab = tab[tab > 1]
	for(ID in names(tab)){
		idx = which(absPP$ID == ID)
		# remove duplicate with lower purity
		tmp_df = absPP[idx,]
		tmp_df = tmp_df[order(tmp_df$purity),]
		rm_samp = tmp_df$sample[1]
		absPP = absPP[which(absPP$sample != rm_samp),]
	}
	# dim(absPP); absPP[1:3,]
	dat = smart_merge(clin,absPP,all.x = TRUE); rm(absPP)
	
	# Get follow up CLIN pancancer data
	fu_fn = file.path(my_dirs$panc_dir,"clinical_PANCAN_patient_with_followup.tsv")
	fu = smart_DTR(fu_fn,header = TRUE,sep = "\t",data.table = FALSE)
	fu = name_change(fu,"bcr_patient_barcode","ID")
	int_subj = intersect(fu$ID,dat$ID); length(int_subj)
	fu = fu[fu$ID %in% int_subj,] # get tumor-type specific subjects
	fu = apply(fu,2,function(xx){
		xx[xx %in% c("","[Not Applicable]","[Not Available]",
			"[Unknown]","[Not Reported]","[Not Evaluated]","[Discrepancy]")] = NA
		xx
		})
	fu = smart_df(fu)
	num_uniq_value = apply(fu,2,function(xx) length(unique(xx))) # remove columns with no variation
	table(num_uniq_value)
	fu = fu[,num_uniq_value > 1]
	int_vars = intersect(names(fu),names(dat)); length(int_vars)
	int_vars = int_vars[int_vars != "ID"]
	fu = fu[,!(names(fu) %in% int_vars)]
	if( "height" %in% names(fu) ) fu$height = as.numeric(fu$height)
	if( "weight" %in% names(fu) ) fu$weight = as.numeric(fu$weight)
	if( all(c("height","weight") %in% names(fu)) ){
		fu$BMI = fu$weight / (fu$height / 100)^2
	}
	dim(fu); fu[1:3,]
	dim(dat)
	dat = smart_merge(dat,fu,all = TRUE)
	dim(dat)
	
	if( getCNA ){
		if( show ) cat(sprintf("%s: Import ABSOLUTE CN data ...\n",date()))
		absCN = smart_DTR(absCN_fn,sep = "\t",header = TRUE)
		absCN = absCN[which(absCN$Sample %in% dat$array),]
		absCN$ID = make_sampID(origIDs = absCN$Sample)
		absCN$Chr = paste0("chr",absCN$Chromosome)
		absCN = absCN[which(!is.na(absCN$Chromosome)),]
		absCN$tCN = absCN$Modal_HSCN_1 + absCN$Modal_HSCN_2
		# dim(absCN); absCN[1:3,]
	} else {
		absCN = NULL
	}
	
	# Other PanCan info like sequencing center and WGA status
	if( show ) cat(sprintf("%s: Get WGA and seq center data ...\n",date()))
	clin2_fn = file.path(my_dirs$panc_dir,"12864_2017_3770_MOESM2_ESM.xlsx")
	if( !file.exists(clin2_fn) ){
		url_link = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5467262/bin/12864_2017_3770_MOESM2_ESM.xlsx"
		stop(sprintf("Download %s and rename %s",url_link,clin2_fn))
	}
	
	clin_ex_fn = file.path(my_dirs$panc_dir,"clin_ex.rds")
	if( !file.exists(clin_ex_fn) ){
		clin2 = read_excel(clin2_fn,sheet = 1)
		clin2 = as.data.frame(clin2,stringsAsFactors = FALSE)
		saveRDS(clin2,clin_ex_fn)
	}
	clin2 = readRDS(clin_ex_fn)
	
	if( dataset == "COAD" ){
		clin2 = clin2[which(clin2$CANCER %in% c("COAD","READ")),]
	} else if( dataset %in% c("BLCA","BRCA","CESC","ESCA","GBM","HNSC","KIRC",
			"LAML","LGG","LIHC","LUAD","LUSC","OV","PAAD","SARC","SKCM",
			"STAD","UCEC") ){
		clin2 = clin2[which(clin2$CANCER == dataset),]
	} else {
		stop("Add code for this dataset")
	}
	
	clin2$ID = make_sampID(origIDs = clin2$TCGA_ID)
	clin2 = smart_rmcols(OBJ = clin2,rm_names = c("TCGA_ID","VCF_ID"))
	dim(clin2); clin2[1:3,]
	
	if( show ) cat(sprintf("%s: Merge clinical and sequencing data ...\n",date()))
	dat = smart_merge(dat,clin2[,c("ID","SAMPLE_TYPE","WGA_STATUS",
		"CENTER","SEQUENCER","KIT","BWA","X20","YEAR_PUBLISHED")],all = TRUE)
	
	if( dataset == "BRCA" ){
		url_link = "https://www.cell.com/cms/10.1016/j.ccell.2018.03.014"
		url_link = sprintf("%s/attachment/dd9a9dca-90c0-42f4-8d9d-304ca0cffd7e/mmc4.xlsx",url_link)
		pam50_fn = file.path(my_dirs$panc_dir,
			"1-s2.0-S1535610818301193-mmc4_Berger2018_gynec_breast.xlsx")
		if( !file.exists(pam50_fn) ){
			stop(sprintf("Download %s and rename %s",url_link,pam50_fn))
		}
		pam50 = suppressWarnings(read_excel(pam50_fn,skip = 1))
		pam50 = as.data.frame(pam50,stringsAsFactors = FALSE)
		pam50 = name_change(pam50,"Sample.ID","ID")
		smart_table(pam50$Tumor.Type,pam50$vital_status)
		pam50 = pam50[which(pam50$Tumor.Type == "BRCA"),]
		pam50 = name_change(pam50,"BRCA_Subtype_PAM50","PAM50")
		pam50$PAM50[pam50$PAM50 == "NA"] = NA
		smart_table(pam50$PAM50,pam50$vital_status)
		pam50 = pam50[,c("ID","PAM50")]
		dim(pam50); pam50[1:3,]
		dat = smart_merge(dat,pam50,all = TRUE)
	}
	
	if( !is.null(samp) ){ # Subset sample's data
		cat(sprintf("%s: Subset sample data ...\n",date()))
		absCN = absCN[which(absCN$ID == samp),]
		dat = dat[which(dat$ID == samp),]
	}
	
	list(absCN = absCN,dat = dat)
}
down_SPMs = function(my_dirs){
	dataset = my_dirs$dataset
	
	if( dataset == "BLCA" ){
		spms_fn = file.path(my_dirs$spms_dir,
			"TCGA.BLCA.mutect.0e239d8f-47b0-4e47-9716-e9ecc87605b9.DR-10.0.somatic.maf.gz")
	} else if( dataset == "BRCA" ){
		spms_fn = file.path(my_dirs$spms_dir,
			"TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf.gz")
	} else if( dataset == "CESC" ){
		spms_fn = file.path(my_dirs$spms_dir,
			"TCGA.CESC.mutect.5ffa70b1-61b4-43d1-b10a-eda412187c17.DR-10.0.somatic.maf.gz")
	} else if( dataset == "ESCA" ){
		spms_fn = file.path(my_dirs$spms_dir,
			"TCGA.ESCA.mutect.7f8e1e7c-621c-4dfd-8fad-af07c739dbfc.DR-10.0.somatic.maf.gz")
	} else if( dataset == "GBM" ){
		spms_fn = file.path(my_dirs$spms_dir,
			"TCGA.GBM.mutect.da904cd3-79d7-4ae3-b6c0-e7127998b3e6.DR-10.0.somatic.maf.gz")
	} else if( dataset == "HNSC" ){
		spms_fn = file.path(my_dirs$spms_dir,
			"TCGA.HNSC.mutect.1aa33f25-3893-4f37-a6a4-361c9785d07e.DR-10.0.somatic.maf.gz")
	} else if( dataset == "KIRC" ){
		spms_fn = file.path(my_dirs$spms_dir,
			"TCGA.KIRC.mutect.2a8f2c83-8b5e-4987-8dbf-01f7ee24dc26.DR-10.0.somatic.maf.gz")
	} else if( dataset == "LAML" ){
		spms_fn = file.path(my_dirs$spms_dir,
			"TCGA.LAML.mutect.27f42413-6d8f-401f-9d07-d019def8939e.DR-10.0.somatic.maf.gz")
	} else if( dataset == "LGG" ){
		spms_fn = file.path(my_dirs$spms_dir,
			"TCGA.LGG.mutect.1e0694ca-fcde-41d3-9ae3-47cfaf527f25.DR-10.0.somatic.maf.gz")
	} else if( dataset == "LIHC" ){
		spms_fn = file.path(my_dirs$spms_dir,
			"TCGA.LIHC.mutect.a630f0a0-39b3-4aab-8181-89c1dde8d3e2.DR-10.0.somatic.maf.gz")
	} else if( dataset == "LUAD" ){
		spms_fn = file.path(my_dirs$spms_dir,
			"TCGA.LUAD.mutect.0458c57f-316c-4a7c-9294-ccd11c97c2f9.DR-10.0.somatic.maf.gz")
	} else if( dataset == "LUSC" ){
		spms_fn = file.path(my_dirs$spms_dir,
			"TCGA.LUSC.mutect.95258183-63ea-4c97-ae29-1bae9ed06334.DR-10.0.somatic.maf.gz")
	} else if( dataset == "OV" ){
		spms_fn = file.path(my_dirs$spms_dir,
			"TCGA.OV.mutect.b22b85eb-2ca8-4c9f-a1cd-b77caab999bd.DR-10.0.somatic.maf.gz")
	} else if( dataset == "PAAD" ){
		spms_fn = file.path(my_dirs$spms_dir,
			"TCGA.PAAD.mutect.fea333b5-78e0-43c8-bf76-4c78dd3fac92.DR-10.0.somatic.maf.gz")
	} else if( dataset == "SARC" ){
		spms_fn = file.path(my_dirs$spms_dir,
			"TCGA.SARC.mutect.cc207fe8-ee0a-4b65-82cb-c8197d264126.DR-10.0.somatic.maf.gz")
	} else if( dataset == "SKCM" ){
		spms_fn = file.path(my_dirs$spms_dir,
			"TCGA.SKCM.mutect.4b7a5729-b83e-4837-9b61-a6002dce1c0a.DR-10.0.somatic.maf.gz")
	} else if( dataset == "STAD" ){
		spms_fn = file.path(my_dirs$spms_dir,
			"TCGA.STAD.mutect.c06465a3-50e7-46f7-b2dd-7bd654ca206b.DR-10.0.somatic.maf.gz")
	} else if( dataset == "THCA" ){
		spms_fn = file.path(my_dirs$spms_dir,
			"TCGA.THCA.mutect.13999735-2e70-439f-a6d9-45d831ba1a1a.DR-10.0.somatic.maf.gz")
	} else if( dataset == "UCEC" ){
		spms_fn = file.path(my_dirs$spms_dir,
			"TCGA.UCEC.mutect.d3fa70be-520a-420e-bb6d-651aeee5cb50.DR-10.0.somatic.maf.gz")
	} else {
		stop("No code for that tumor type")
	}
	
	if( dataset != "COAD" && !file.exists(spms_fn) ){
		stop(sprintf("Download %s SPMs and rename %s",
			dataset,spms_fn))
	} else if( dataset == "COAD" ){
		coad_spms_fn = file.path(my_dirs$spms_dir,
			"TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.gz")
		if( !file.exists(coad_spms_fn) )
			stop(sprintf("Download COAD SPMs and rename %s",coad_spms_fn))
		
		read_spms_fn = file.path(my_dirs$spms_dir,
			"TCGA.READ.mutect.faa5f62a-2731-4867-a264-0e85b7074e87.DR-10.0.somatic.maf.gz")
		if( !file.exists(read_spms_fn) )
			stop(sprintf("Download READ SPMs and rename %s",read_spms_fn))
		spms_fn = c(coad_spms_fn,read_spms_fn)
	
	}
	
	spms_fn
}
get_PATH = function(my_dirs,gmt_fn = NULL){
	
	if( is.null(gmt_fn) ){
		gmt_fn = file.path(my_dirs$git_dir,
			"ROKETworkflow/inst/extdata",
			"c2.cp.reactome.v7.2.symbols.gmt")
	}
	
	gmt = list()
	tmp_gmt = readLines(gmt_fn)
	tmp_gmt = t(sapply(tmp_gmt,function(xx){
		yy = strsplit(xx,"\t")[[1]][-2]
		c(yy[1],paste(yy[-1],collapse=" "))
		},USE.NAMES = FALSE))
	dim(tmp_gmt); tmp_gmt[1:4,]
	
	# make gene by pathway matrix
	nPATH = nrow(tmp_gmt)
	uGENE = unique(strsplit(paste(tmp_gmt[,2],collapse = " ")," ")[[1]])
	nGENE = length(uGENE); nGENE
	pp = matrix(0,nGENE,nPATH)
	pp = smart_names(pp,ROW = uGENE,COL = tmp_gmt[,1])
	pp[1:4,1:2]
	for(ii in seq(nPATH)){
		# ii = 1
		smart_progress(ii = ii,nn = nPATH,iter = 2e1,iter2 = 1e3)
		tmp_genes = strsplit(tmp_gmt[ii,2]," ")[[1]]
		pp[tmp_genes,ii] = 1
		rm(tmp_genes)
	}
	pp = pp[rowSums(pp) > 0,]
	
	return(pp)
}
get_GTF = function(my_dirs){
	
	gtf_fn = file.path(my_dirs$ref_dir,"gencode.gene.info.v22.tsv")
	gtf = data.table::fread(gtf_fn,sep = "\t",header = TRUE)
	table(gtf$gene_type)
	gtf = gtf[gtf$gene_type == "protein_coding" 
		& gtf$seqname %in% paste0("chr",1:22),]
	gtf = name_change(gtf,"gene_id","ENSG")
	gtf = name_change(gtf,"gene_name","GENE")
	gtf = name_change(gtf,"seqname","Chr")
	gtf$ENSG0 = sapply(gtf$ENSG,function(xx) 
		strsplit(xx,"[.]")[[1]][1],USE.NAMES = FALSE)
	dim(gtf); gtf[1:3,]
	
	return(gtf)
	
}


# ----------
# Real Data: Pre-Processing functions
# ----------
keep_GDCFILTER = function(DATA,VAR = "GDC_FILTER",rmWGA = TRUE){
	cat(sprintf("%s: Filter ndp ...\n",date()))
	DATA = DATA[!grepl("ndp",DATA[[VAR]]),,drop=FALSE]
	cat(sprintf("%s: Filter NonExonic ...\n",date()))
	DATA = DATA[!grepl("NonExonic",DATA[[VAR]]),,drop=FALSE]
	cat(sprintf("%s: Filter exac ...\n",date()))
	DATA = DATA[!grepl("common_in_exac",DATA[[VAR]]),,drop=FALSE]
	cat(sprintf("%s: Filter bitgt ...\n",date()))
	DATA = DATA[!grepl("bitgt",DATA[[VAR]]),,drop=FALSE]
	cat(sprintf("%s: Filter gdc_pon ...\n",date()))
	DATA = DATA[!grepl("gdc_pon",DATA[[VAR]]),,drop=FALSE]
	if( rmWGA ){
		cat(sprintf("%s: Filter wga ...\n",date()))
		DATA = DATA[!grepl("wga",DATA[[VAR]]),,drop=FALSE]
	}
	DATA
}
keep_varClass = function(DATA,VAR = "Variant_Classification"){
	cat(sprintf("%s: Filter specific variant classes ...\n",date()))
	keep_classes = c("Frame_Shift_Del","Frame_Shift_Ins",
		"In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation",
		"Splice_Region","Splice_Site","Translation_Start_Site")
	DATA[which(DATA[[VAR]] %in% keep_classes),,drop=FALSE]
}
remove_lowImpact = function(DATA,VAR = "Consequence"){
	cat(sprintf("%s: Filter synonymous ...\n",date()))
	DATA = DATA[!grepl("synonymous",DATA[[VAR]]),,drop=FALSE]
	cat(sprintf("%s: Filter 3_prime_UTR ...\n",date()))
	DATA = DATA[!grepl("3_prime_UTR",DATA[[VAR]]),,drop=FALSE]
	cat(sprintf("%s: Filter 5_prime_UTR ...\n",date()))
	DATA = DATA[!grepl("5_prime_UTR",DATA[[VAR]]),,drop=FALSE]
	cat(sprintf("%s: Filter NMD ...\n",date()))
	DATA = DATA[!grepl("NMD_",DATA[[VAR]]),,drop=FALSE]
	cat(sprintf("%s: Filter intron ...\n",date()))
	DATA = DATA[!grepl("intron_",DATA[[VAR]]),,drop=FALSE]
	cat(sprintf("%s: Filter non_coding_transcript ...\n",date()))
	DATA = DATA[!grepl("non_coding_transcript",DATA[[VAR]]),,drop=FALSE]
	cat(sprintf("%s: Filter stop_retain ...\n",date()))
	DATA = DATA[!grepl("stop_retain",DATA[[VAR]]),,drop=FALSE]
	DATA
}
prep_SPMs = function(my_dirs,samp = NULL,rmWGA = TRUE){
	rm_names = c("all_effects","src_vcf_id","tumor_bam_uuid",
		"normal_bam_uuid","case_id","Matched_Norm_Sample_Barcode",
		"Tumor_Sample_UUID","Matched_Norm_Sample_UUID")
	
	# Get SPMs
	dataset = my_dirs$dataset
	raw_spms_fn = file.path(my_dirs$spms_dir,"raw_spms.rds")
	if( !file.exists(raw_spms_fn) ){
		spms_fn = down_SPMs(my_dirs = my_dirs)
		spms = c()
		for(fn in spms_fn){
			tmp_spms = smart_DTR(fn,header = TRUE,sep = "\t")
			tmp_spms = smart_rmcols(OBJ = tmp_spms,rm_names = rm_names)
			tmp_spms = make_DP_VAF(DATA = tmp_spms)
			spms = rbind(spms,tmp_spms)
		}
		spms$ID = make_sampID(origIDs = spms$Tumor_Sample_Barcode)
		saveRDS(spms,raw_spms_fn)
	}
	spms = readRDS(raw_spms_fn)
	
	# Filter VCs
	spms = keep_GDCFILTER(DATA = spms,rmWGA = rmWGA)
	spms = keep_varClass(DATA = spms)
	spms = remove_lowImpact(DATA = spms)
	
	if( !is.null(samp) ){ # Subset sample data
		cat(sprintf("%s: Subset sample's SPMs ...\n",date()))
		spms = spms[which(spms$ID == samp),,drop=FALSE]
	}
	
	spms
}
make_SPM_mat = function(my_dirs,rmWGA = TRUE){
	
	dataset = my_dirs$dataset
	mat_fn = file.path(my_dirs$spms_dir,sprintf("%s_mat.rds",dataset))
	if( file.exists(mat_fn) ){
		mut_spm = readRDS(mat_fn)
		return(mut_spm)
	}
	
	# Process gtf: Convert all genes to Hugo symbols
	gtf = get_GTF(my_dirs = my_dirs)
	dim(gtf); gtf[1:3,]
	
	# Get SPMs
	spms = prep_SPMs(my_dirs = my_dirs,rmWGA = rmWGA)
	spms = spms[which(spms$t_alt_count > 0),]
	length(unique(spms$Gene))
	length(unique(spms$Hugo_Symbol))
	
	# Make SPMs MUT matrix
	cat(sprintf("%s: Make SPM mutation matrix ...\n",date()))
	spms_ids = unique(spms$ID)
	nsamp_spms = length(spms_ids); nsamp_spms
	mut_spm = matrix(0,nrow(gtf),nsamp_spms)
	mut_spm = smart_names(mut_spm,ROW = gtf$ENSG0,COL = spms_ids)
	# mut_spm[1:5,1:4]
	cnt = 1
	for(samp in spms_ids){
		# samp = spms_ids[1]; samp
		smart_progress(ii = cnt,nn = nsamp_spms,iter = 5,iter2 = 1e2)
		idx = which(spms$ID == samp)
		ensgs = unique(spms$Gene[idx])
		ensgs = intersect(ensgs,rownames(mut_spm))
		# ensgs[1:10]
		mut_spm[ensgs,samp] = 1
		cnt = cnt + 1
		rm(idx,ensgs)
	}
	
	mut_spm = mut_spm[rowSums(mut_spm) >= 1,]
	gtf = gtf[which(gtf$ENSG0 %in% rownames(mut_spm)),]
	gtf = gtf[match(rownames(mut_spm),gtf$ENSG0),]
	rownames(mut_spm) = gtf$GENE
	dim(mut_spm); mut_spm[1:10,1:5]
	saveRDS(mut_spm,mat_fn)
	
	return(mut_spm)
}


# ----------
# Real Data: Gene Similarity functions
# ----------
get_geneSim_path = function(gene1,gene2,path,NN){
	
	if( gene1 == gene2 ){
		return(0)
	}
	
	# NN = num balls = num genes
	idx1 = path[gene1,] == 1
	idx2 = path[gene2,] == 1
	
	# KK = num white balls = num genes with same pathways as gene1 and gene2
	idx_path = which(idx1 & idx2); # idx_path
	num_path = length(idx_path); num_path
	
	if( num_path > 0 ){
		idx_gene = rowSums(path[,idx_path,drop = FALSE]) == num_path
		KK = sum(idx_gene); KK
	} else {
		KK = 0
	}
	
	if( KK < 2 ){
		# stop(sprintf("check %s and %s",gene1,gene2))
		return(1)
	}
	
	# Calculate prob of drawing two genes that are on the same pathway
	prob = dhyper(
		x = 2,       # num white balls drawn
		m = KK,      # num white balls in urn
		n = NN - KK, # num black balls in urn
		k = 2        # num balls drawn
		)
	prob
	
}
get_full_geneSim_path = function(use_genes,path){
	
	use_genes2 = intersect(use_genes,rownames(path))
	ngenes = length(use_genes2); ngenes
	mat = matrix(NA,ngenes,ngenes)
	mat = smart_names(mat,ROW = use_genes2,COL = use_genes2)
	NN = nrow(path)
	
	cnt = 1; tot = ngenes^2
	for(ii in seq(ngenes)){
	for(jj in seq(ngenes)){
		smart_progress(ii = cnt,nn = tot,iter = 1e2,iter2 = 4e3)
		
		if( ii <= jj ){
			# ii = 5; jj = ii
			tmp_out = tryCatch(get_geneSim_path(
				gene1 = use_genes2[ii],
				gene2 = use_genes2[jj],
				path = path,NN = NN),
				warning = function(ww){
					cmd = sprintf("check %s and %s:",
						use_genes2[ii],use_genes2[jj])
					stop(cmd)
					})
			
			mat[ii,jj] = tmp_out
			mat[jj,ii] = mat[ii,jj]
			
		}
		cnt = cnt + 1
	}}
	
	return(mat)
}
get_geneSim_ME = function(gene1,gene2,mSPM,log10_TMB){
	if(FALSE){
		gene1 = "DMD"
		gene2 = "ABCA13"
		# mSPM
		
	}
	
	tab = table(mSPM[gene1,],mSPM[gene2,])
	names(attr(tab,"dimnames")) = c(gene1,gene2); tab
	marg_out = fisher.test(tab,alternative = "less")
	MARG = list(OR = as.numeric(marg_out$estimate),
		PVAL = marg_out$p.value)
	
	log10_TMB2 = (log10_TMB)^2
	glm_out1 = tryCatch(glm(mSPM[gene2,] ~ mSPM[gene1,]
		+ log10_TMB + log10_TMB2,
		family = "binomial"),warning = function(ww){NULL})
	# summary(glm_out)
	
	glm_out2 = tryCatch(glm(mSPM[gene1,] ~ mSPM[gene2,]
		+ log10_TMB + log10_TMB2,
		family = "binomial"),warning = function(ww){NULL})
	# summary(glm_out2)
	
	if( is.null(glm_out1) | is.null(glm_out2) | MARG$OR == 0 ){
		COND = MARG
	} else {
		PVAL_1 = pnorm(summary(glm_out1)$coefficients[2,3])
		PVAL_2 = pnorm(summary(glm_out2)$coefficients[2,3])
		PVAL = sqrt(PVAL_1 * PVAL_2)
		COND = list(PVAL = PVAL)
	}
	
	list(MARG = MARG,COND = COND)
}
get_full_geneSim_ME = function(use_genes,mSPM,
	strat = "all",TYPE = "marg"){
	
	if(FALSE){
		min_gene_mut = 45
		use_genes = names(tb1)[which(tb1 >= min_gene_mut)]
		length(use_genes)
		
		use_genes = use_genes
		mSPM = matSPM
		strat = c("all","nHM","HM")[1]
		TYPE = c("marg","cond")[2]
		
	}
	
	DAT = smart_df(ID = colnames(mSPM),
		TMB = colSums(mSPM))
	DAT$TMB_bin = ifelse(log10(DAT$TMB) > 2.5,1,0)
	
	if( strat %in% c("HM","nHM") ){
		if( strat == "nHM" ) 	DAT = DAT[which(DAT$TMB_bin == 0),]
		if( strat == "HM" ) 	DAT = DAT[which(DAT$TMB_bin == 1),]
	}
	mSPM = mSPM[,DAT$ID]
	TMB = DAT$TMB
	log10_TMB = log10(TMB)
	
	dim(mSPM)
	int_genes = intersect(use_genes,rownames(mSPM))
	ngenes = length(int_genes); ngenes
	mat = matrix(NA,ngenes,ngenes)
	mat = smart_names(mat,int_genes,int_genes)
	
	cnt = 1; tot = ngenes^2
	for(ii in seq(ngenes)){
	for(jj in seq(ngenes)){
		# ii = 1; jj = 2
		smart_progress(ii = cnt,nn = tot,iter = 5e2,iter2 = 1e4)
		
		if( ii < jj ){
			
			# int_genes[ii]; int_genes[jj]
			
			ana_out = tryCatch(get_geneSim_ME(
				gene1 = int_genes[ii],
				gene2 = int_genes[jj],
				mSPM = mSPM,log10_TMB = log10_TMB),
				warning = function(ww){NULL})
			
			if( is.null(ana_out) ){
				stop("Check gene pair")
			}
			
			mat[ii,jj] = ifelse(TYPE == "marg",
				ana_out$MARG$PVAL,ana_out$COND$PVAL)
			mat[jj,ii] = mat[ii,jj]
			
		}
		
		cnt = cnt + 1
		
	}}
	diag(mat) = 0
	# dim(mat); mat[1:5,1:5]
	# summary(diag(mat))
	
	sim_mat = mat
	if(FALSE){ # Make PVAL matrix symmetric
		sim_mat[sim_mat > -1] = NA
		for(ii in seq(ngenes)){
		for(jj in seq(ngenes)){
			if( ii == jj ){
				sim_mat[ii,ii] = mat[ii,ii]
			} else if( ii < jj ){
				# Average the logs() <=> tilde(p_ij) = sqrt(p_ij * p_ji)
				sim_mat[ii,jj] = sqrt(mat[ii,jj] * mat[jj,ii])
				sim_mat[jj,ii] = sim_mat[ii,jj]
			}
		}}
	}
	
	sgenes = intersect(rownames(mat),
		c("TP53","KRAS","BRAF","PIK3CA","TTN"))
	
	# Show top gene pairs
  tmp_mat = sim_mat
  tmp_mat[lower.tri(tmp_mat,diag = TRUE)] = 1
  xx = tmp_mat
  xx = xx[upper.tri(xx)]
	top_num = min(c(20,length(xx))); top_num
  xx2 = unique(sort(xx))[1:top_num]; # xx2
  tmp_idx = which(tmp_mat < xx2[top_num],arr.ind = TRUE); # tmp_idx
	res = c()
	for(jj in seq(nrow(tmp_idx))){
		# jj = 1
		tmp_mat[tmp_idx[jj,1],tmp_idx[jj,2],drop = FALSE]
		res = rbind(res,smart_df(G1 = rownames(tmp_mat)[tmp_idx[jj,1]],
			G2 = rownames(tmp_mat)[tmp_idx[jj,2]],
			ME = tmp_mat[tmp_idx[jj,1],tmp_idx[jj,2]]))
	}
	rm(tmp_mat)
	res = res[order(res$ME),]
	rownames(res) = NULL
	# print(res)
	
	# Calculate geneSim
	gsim = -log(sim_mat)
	upnorm = quantile(gsim[upper.tri(gsim)],0.99); upnorm
	# upnorm = max(gsim[upper.tri(gsim)])
	gsim = gsim / upnorm
	# gsim[sgenes,sgenes]
	# any(gsim > 1)
	gsim[gsim > 1] = 1
	# diag(gsim) = 1 
	
	list(PVAL = mat,SIM = sim_mat,RES = res,GS = gsim)
}
make_tb1 = function(my_dirs,rmWGA = TRUE,show = TRUE){
	
	dataset = my_dirs$dataset
	
	# Save gene mutation frequency
	tb1_fn = file.path(my_dirs$spms_dir,sprintf("%s_tb1.rds",dataset))
	if( !file.exists(tb1_fn) ){
		
		# Get SPM data
		print(dataset)
		print(rmWGA)
		spms = prep_SPMs(my_dirs = my_dirs,rmWGA = rmWGA)
		print(spms[1:10,c("Tumor_Sample_Barcode","Hugo_Symbol")])
		u_spms = unique(spms[,c("Hugo_Symbol","Tumor_Sample_Barcode")])
		print(dim(u_spms))
		
		tb1 = table(u_spms$Hugo_Symbol)
		print(sort(tb1,decreasing = TRUE)[1:30])
		rm(spms,u_spms)
		saveRDS(tb1,tb1_fn)
		
	}
	tb1 = readRDS(tb1_fn)
	
	for(freq in c(1,2,5,10,15,20,30,40,50,60)){
		# freq = 1
		if( !show ) next
		cat(sprintf("\nGenes mutated at least %s times ...\n",freq))
		print(table(tb1 >= freq))
	}
	
	return(tb1)
}

#' @title run_full_geneSim
#' @param my_dirs A list generated by \code{setdirs()}
#'	of directories specific to a cancer type
#' @param min_gene_mut A non-negative integer to filter
#'	genes, based on mutation frequency, to include in distance 
#'	calculations
#' @param rmWGA Boolean set to TRUE to exclude WGA samples
#' @export
run_full_geneSim = function(my_dirs,min_gene_mut,rmWGA = TRUE){
	if(FALSE){
		dataset = "BLCA"
		min_gene_mut = 5
		rmWGA = TRUE
		
	}
	
	dataset = my_dirs$dataset
	tb1 = make_tb1(my_dirs = my_dirs,rmWGA = rmWGA)
	
	# Make SPM matrix
	matSPM = make_SPM_mat(my_dirs = my_dirs,rmWGA = rmWGA)
	
	# Set reference genes
	fav_genes = c("TP53","ATM","KRAS","BRAF",
		"NRAS","ERBB2","APC","SOX9","TTN")
	
	# Output GO-based of all genes
	all_GOsim_fn = file.path(my_dirs$spms_dir,"gsim_GO_all.rds")
	if( !file.exists(all_GOsim_fn) ){
		
		# Get hsGO
		cat(sprintf("%s: Get hsGO ...\n",date()))
		hsGO = godata("org.Hs.eg.db",
			keytype = "SYMBOL",ont = "BP",computeIC = FALSE)
		
		cat(sprintf("%s: Get GO-based gene_sim ...\n",date()))
		gene_sim_GO = mgeneSim(genes = names(tb1),
			semData = hsGO,measure = "Wang",verbose = TRUE,
			combine = "BMA")
		saveRDS(gene_sim_GO,all_GOsim_fn)
	}
	gene_sim_GO = readRDS(all_GOsim_fn)
	
	# Output PATH-based of all genes
	all_PATHsim_fn = file.path(my_dirs$spms_dir,"gsim_PATH_all.rds")
	if( !file.exists(all_PATHsim_fn) ){
		
		# Import Pathways
		cat(sprintf("%s: Import REACTOME ...\n",date()))
		react = get_PATH(my_dirs = my_dirs)
		dim(react); react[1:5,1:3]
		
		cat(sprintf("%s: Get PATH-based gene_sim ...\n",date()))
		tmp_out = get_full_geneSim_path(
			use_genes = names(tb1),path = react)
		gene_sim_PATH = 1 - tmp_out; rm(tmp_out)
		saveRDS(gene_sim_PATH,all_PATHsim_fn)
	}
	gene_sim_PATH = readRDS(all_PATHsim_fn)
	
	# Output ME-based of all genes
	all_MEsim_fn = file.path(my_dirs$spms_dir,"gsim_ME_all.rds")
	if( !file.exists(all_MEsim_fn) ){
		cat(sprintf("%s: Get ME-based gene_sim ...\n",date()))
		strat = "all"
		TYPE = "cond"
		tmp_gene_mut = 0
		while(TRUE){
			use_genes2 = names(tb1)[which(tb1 >= tmp_gene_mut)]
			cat(sprintf("minGM_thres = %s, nGenes = %s ...\n",
				tmp_gene_mut,length(use_genes2)))
			if( length(use_genes2) <= 300 ){ # need to set an upper limit
				break
			}
			tmp_gene_mut = tmp_gene_mut + 1
		}
		tmp_out = get_full_geneSim_ME(use_genes = use_genes2,
			mSPM = matSPM,strat = strat,TYPE = TYPE)
		gene_sim_ME = tmp_out$GS; rm(tmp_out)
		saveRDS(gene_sim_ME,all_MEsim_fn)
	}
	gene_sim_ME = readRDS(all_MEsim_fn)
	
	# If output file exists, we're done
	gsim_fn = file.path(my_dirs$spms_dir,
		sprintf("gsim_%s.rds",min_gene_mut))
	if( file.exists(gsim_fn) ){
		gsim = readRDS(gsim_fn)
		return(gsim)
	}
	
	# Decide on gene set
	use_genes = names(tb1)[which(tb1 >= min_gene_mut)]
	print(sprintf("Num genes = %s",length(use_genes)))
	
	if( length(use_genes) <= 3 ) return(NULL)
	
	# GO-based geneSim
	cat(sprintf("%s: Get GO-based geneSim ...\n",date()))
	int_GO_genes = intersect(use_genes,rownames(gene_sim_GO))
	gene_sim_GO = gene_sim_GO[int_GO_genes,int_GO_genes]
	
	# PATH-based geneSim
	cat(sprintf("%s: Get PATH-based geneSim ...\n",date()))
	int_PATH_genes = intersect(use_genes,rownames(gene_sim_PATH))
	gene_sim_PATH = gene_sim_PATH[int_PATH_genes,int_PATH_genes]
	# out = get_full_geneSim_path(use_genes = use_genes,
		# path = react)
	# gene_sim_PATH = 1 - out; rm(out)
	
	# ME-based geneSim
	cat(sprintf("%s: Get ME-based geneSim ...\n",date()))
	ME_genes = intersect(use_genes,rownames(gene_sim_ME))
	gene_sim_ME = gene_sim_ME[ME_genes,ME_genes]
	
	# Combine geneSims into array
	all_genes = unique(c(rownames(gene_sim_GO),
		rownames(gene_sim_PATH),
		rownames(gene_sim_ME)))
	num_genes = length(all_genes); num_genes

	fin = array(data = NA,dim = c(num_genes,num_genes,3),
		dimnames = list(all_genes,all_genes,c("GO","PATH","ME")))
	# str(fin)
	
	# GO
	genes_go = intersect(rownames(gene_sim_GO),all_genes)
	diag(fin[,,"GO"]) = 1
	fin[genes_go,genes_go,"GO"] = gene_sim_GO[genes_go,genes_go]
	rm(genes_go)
	any(is.na(fin[,,"GO"]))
	
	# PATH
	genes_path = intersect(rownames(gene_sim_PATH),all_genes)
	diag(fin[,,"PATH"]) = 1
	fin[genes_path,genes_path,"PATH"] = gene_sim_PATH[genes_path,genes_path]
	rm(genes_path)
	any(is.na(fin[,,"PATH"]))
	
	# ME
	genes_me = intersect(rownames(gene_sim_ME),all_genes)
	diag(fin[,,"ME"]) = 1
	fin[genes_me,genes_me,"ME"] = gene_sim_ME[genes_me,genes_me]
	rm(genes_me)
	any(is.na(fin[,,"ME"]))
	
	# Output
	saveRDS(fin,gsim_fn)
	
	return(fin)
	
}


# ----------
# Real Data: Optimal Transport Analysis
# ----------
comb_geneSim = function(WTS,geneSims){
	if(FALSE){
		WTS = rep(1,3)
		WTS = c(1,0,1)
		WTS = c(0,1,0)
		geneSims = gsim
	}
	
	# geneSims[1:5,1:5,]
	# any(geneSims[,,"PATH"] == 0)
	
	genes = rownames(geneSims); # genes
	ngenes = length(genes); ngenes
	fin = matrix(0,ngenes,ngenes)
	diag(fin) = 1
	fin = smart_names(fin,ROW = genes,COL = genes)
	
	for(ii in seq(ngenes)){
	for(jj in seq(ngenes)){
		if( ii < jj ){
			# ii = 1; jj = 5
			geneSims[ii,jj,]
			nWTS = WTS * (1 * !is.na(geneSims[ii,jj,]))
			if( sum(nWTS) == 0 ){
				nWTS = rep(0,3)
			} else {
				nWTS = nWTS / sum(nWTS)			
			}
			nWTS
			vec_gs = geneSims[ii,jj,]
			vec_gs[is.na(vec_gs)] = 0
			tmp_genesim = sum(vec_gs * nWTS)
			fin[ii,jj] = tmp_genesim
		}
	}}
	fin[lower.tri(fin)] = t(fin)[lower.tri(fin)]
	# fin[1:5,1:5]
	
	return(fin)
}

#' @title calc_full_DIST
#' @param my_dirs A list generated by \code{setdirs()}
#'	of directories specific to a cancer type
#' @param rmWGA Boolean set to TRUE to exclude WGA samples
#' @param min_gene_mut A non-negative integer to filter
#'	genes, based on mutation frequency, to include in distance 
#'	calculations
#' @param LAMBDA A positive numeric value or \code{Inf} to
#'	indicate unbalanced or balanced optimal transport distance
#'	calculation, respectively.
#' @param ncores A positive integer to specify the number of
#'	threads to decrease computational runtime.
#' @export
calc_full_DIST = function(my_dirs,rmWGA = TRUE,
	min_gene_mut,LAMBDA,ncores = 1){
	
	if(FALSE){
		rmWGA = TRUE
		min_gene_mut = 15
		# RULE = c("highLAM-lowTMB","highLAM-highTMB")[1]
		BAL = c(TRUE,FALSE)[1]
		LAMBDAs = c(1,1)
		
	}
	
	RULE = "highLAM-lowTMB"
	dataset = my_dirs$dataset
	
	# Set reference genes
	fav_genes = c("TP53","ATM","KRAS","BRAF","NRAS",
		"ERBB2","APC","SOX9","TTN","FRAS1")
	
	# Get geneSims
	gsim = run_full_geneSim(my_dirs = my_dirs,
		min_gene_mut = min_gene_mut,rmWGA = TRUE)
	str(gsim)
	sgenes = intersect(fav_genes,dimnames(gsim)[[1]]); sgenes
	round(gsim[sgenes,sgenes,],3)
	
	# Get SPM matrix
	mat_fn = file.path(my_dirs$spms_dir,sprintf("%s_mat.rds",dataset))
	matSPM = readRDS(mat_fn)
	# summary(colSums(matSPM))
	
	# Get weights
	wts = matrix(c(
		1,1,1, # equal
		1,0,0, # one approach
		0,1,0,
		0,0,1,
		1,1,0, # two approaches
		1,0,1,
		0,1,1),ncol = 3,byrow = TRUE)
	wts = t(apply(wts,1,function(xx) xx / sum(xx)))
	wts
	
	for(ww in seq(nrow(wts))){
		# ww = 4
		if( ww %in% c(1,5,6,7) ) next
		
		ww_str = paste(round(wts[ww,],3),collapse="-"); ww_str
		cat(sprintf("%s: Weights = (%s) ...\n",date(),ww_str))
		
		# Out filename
		LAMBDA_2 = LAMBDA
		BAL = ifelse(is.infinite(LAMBDA),TRUE,FALSE)
		if( BAL ) LAMBDA_2 = 1
		rds_fn = file.path(my_dirs$ithOT_dir,
			sprintf("OT_GM.%s_%s_BAL.%s_L1.%s_L2.%s_ww%s.rds",
			min_gene_mut,RULE,BAL,LAMBDA_2,LAMBDA_2,ww))
		if( file.exists(rds_fn) ) next
		
		# Define current geneSim
		tmp_gsim = comb_geneSim(WTS = wts[ww,],geneSims = gsim)
		# dim(tmp_gsim)
		# gsim[sgenes,sgenes,]
		# tmp_gsim[sgenes,sgenes]
		
		# Gene intersection
		int_genes = intersect(rownames(tmp_gsim),rownames(matSPM))
		length(int_genes)
		
		# Get gene COST
		COST = 1 - tmp_gsim[int_genes,int_genes]; rm(tmp_gsim)
		round(COST[sgenes,sgenes],5)
		
		# Define current sMUT
		sMUT = matSPM[int_genes,]
		sMUT = sMUT[,colSums(sMUT) > 0]
		dim(sMUT)
		
		if( any(is.na(COST)) ) stop("NA in COST")
		if( !all(diag(COST) == 0) ) stop("diag(COST) != 0")
		
		# Run OT
		if( !all(rownames(COST) == rownames(sMUT)) ) stop("genes not ordered")
		if( !all(colnames(COST) == rownames(sMUT)) ) stop("genes not ordered")
		EPS = 1e-3
		conv = 1e-3
		highLAM_lowMU = ifelse(RULE == "highLAM-lowTMB",TRUE,FALSE)
		subj_names = colnames(sMUT)#[1:40]
		ot_out = run_myOTs(ZZ = sMUT[,subj_names],COST = COST,
			EPS = EPS,LAMBDA1 = LAMBDA_2,LAMBDA2 = LAMBDA_2,
			balance = BAL,conv = conv,max_iter = 3e3,ncores = ncores,
			show = FALSE,show_iter = 50)
		
		# Save output
		saveRDS(ot_out,rds_fn)
		
	}
	
	return(NULL)
	
}


# ----------
# Real Data: Regression Analysis
# ----------
RDA_prepCLIN = function(my_dirs){
	
	dataset = my_dirs$dataset
	panc = down_PanCan(my_dirs = my_dirs)
	names(panc)
	DAT = panc$dat; # rm(panc)
	dim(DAT); DAT[1:3,]
	
	if( TRUE ){ # Get hg38 SPMs for batch variables
		uu = readRDS(file.path(my_dirs$spms_dir,"raw_spms.rds"))
		uu = unique(uu[,c("ID","Tumor_Sample_Barcode")])
		rownames(uu) = NULL
		uu$barcode_TSS = sapply(uu$Tumor_Sample_Barcode,
			function(xx) strsplit(xx,"-")[[1]][2],USE.NAMES = FALSE)
			table(uu$barcode_TSS)
		uu$barcode_plate = sapply(uu$Tumor_Sample_Barcode,
			function(xx) strsplit(xx,"-")[[1]][6],USE.NAMES = FALSE)
			table(uu$barcode_plate)
		uu$barcode_analyte = sapply(uu$Tumor_Sample_Barcode,
			function(xx) strsplit(strsplit(xx,"-")[[1]][5],"")[[1]][3],
			USE.NAMES = FALSE)
			table(uu$barcode_analyte)
		uu$barcode_SEQCENTER = sapply(uu$Tumor_Sample_Barcode,
			function(xx) strsplit(xx,"-")[[1]][7],USE.NAMES = FALSE)
			table(uu$barcode_SEQCENTER)
		# uu[1:3,]
		DAT = smart_merge(DAT,uu,all.x = TRUE)
		rm(uu)
	}
	if( TRUE ){ # Polish variables
		DAT = name_change(DAT,"ajcc_pathologic_tumor_stage","tumor_stage")
		DAT$tumor_stage[DAT$tumor_stage %in% paste0("Stage 0",c("","A","B","C"))] = "0"
		DAT$tumor_stage[DAT$tumor_stage %in% paste0("Stage I",c("","A","B","C"))] = "I"
		DAT$tumor_stage[DAT$tumor_stage %in% paste0("Stage II",c("","A","B","C"))] = "II"
		DAT$tumor_stage[DAT$tumor_stage %in% paste0("Stage III",c("","A","B","C"))] = "III"
		DAT$tumor_stage[DAT$tumor_stage %in% paste0("Stage IV",c("","A","B","C"))] = "IV"
		DAT$tumor_stage[DAT$tumor_stage %in% c("[Unknown]","[Discrepancy]")] = NA
		smart_table(DAT$tumor_stage)
		
		DAT = name_change(DAT,"histological_type","hist_type")
		DAT$hist_type[DAT$hist_type == "[Discrepancy]"] = NA
		
		DAT$hist_grade = DAT$histological_grade
		DAT$hist_grade[DAT$hist_grade == "[Discrepancy]"] = NA
		
		DAT$clin_stage = DAT$clinical_stage
		DAT$clin_stage = gsub("Stage ","",DAT$clin_stage)
		DAT$clin_stage[DAT$clin_stage %in% c("I","IA","IA1","IA2","IB",
			"IB1","IB2","IC")] = "I"
		DAT$clin_stage[DAT$clin_stage %in% c("II","IIA",
			"IIA1","IIA2","IIB","IIC")] = "II"
		DAT$clin_stage[DAT$clin_stage %in% c("III","IIIA",
			"IIIB","IIIC","IIIC1","IIIC2")] = "III"
		DAT$clin_stage[DAT$clin_stage %in% c("IV","IVA","IVB")] = "IV"
		
		DAT$BWA[DAT$BWA == "NA"] = NA
		smart_table(DAT$BWA)
		
		DAT$X20[DAT$X20 == "NA"] = NA
		DAT$X20 = as.numeric(DAT$X20)
		
		# Age
		DAT = name_change(DAT,"age_at_initial_pathologic_diagnosis","age")
		DAT$age_2 = DAT$age^2
		
		# Race
		DAT$race[DAT$race %in% c("AMERICAN INDIAN OR ALASKA NATIVE","ASIAN",
			"BLACK OR AFRICAN AMERICAN")] = "OTHER"
		DAT$race[grepl("NATIVE HAWA",DAT$race)] = "OTHER"
		DAT$race[is.na(DAT$race) | DAT$race %in% c("[Not Evaluated]","[Unknown]")] = NA
		smart_table(DAT$race)
		
		# Initial pathologic dx
		DAT$ipath_year = DAT$initial_pathologic_dx_year
		
	}
	
	if( dataset == "BLCA" ){
		DAT$tumor_stage[DAT$tumor_stage %in% c("I","II")] = "I/II"
		
		DAT$hist_grade = DAT$histological_grade
		DAT$hist_grade[DAT$hist_grade %in% c("[Unknown]","Low Grade")] = "Unknown/LowGrade"
		
		DAT$BWA[DAT$BWA %in% c("VN:0.5.7-6","VN:0.5.9-r16")] = "VN:0.5.REF"
		
		DAT$ipath_year[DAT$ipath_year %in% c("1999","2000","2001",
			"2002","2003","2004","2005","2006")] = "1999-2006"
		DAT$ipath_year[DAT$ipath_year %in% c("2007","2008","2009")] = "2007-2009"
		
		return(DAT)
		
	} else if( dataset == "BRCA" ){
		DAT$tumor_stage[DAT$tumor_stage %in% c("Stage X","IV")] = "IV/X"
		
		DAT$ipath_year[DAT$ipath_year %in% c("1988","1989",
			"1990","1991","1992","1993","1994","1995","1996")] = "1988-1996"
		DAT$ipath_year[DAT$ipath_year %in% c("1997","1998",
			"1999","2000")] = "1997-2000"
		DAT$ipath_year[DAT$ipath_year %in% c("2001","2002",
			"2003","2004","2005")] = "2001-2005"
		DAT$ipath_year[DAT$ipath_year %in% c("2006","2007",
			"2008","2009")] = "2006-2009"
		DAT$ipath_year[DAT$ipath_year %in% c("2010","2011",
			"2012","2013")] = "2010-2013"
		
		DAT$hist_type[!grepl("Infiltrating Ductal",DAT$hist_type)
			& !is.na(DAT$hist_type)] = "OtherCarc"
		
		DAT$SEQUENCER[grepl("Illumina HiSeq",DAT$SEQUENCER)] = "Illumina HiSeq"
		DAT$PAM50[DAT$PAM50 %in% c("Her2","Normal")] = "Her2,Normal"
		
		return(DAT)
		
	} else if( dataset == "CESC" ){
		DAT$ipath_year[DAT$ipath_year %in% c("1994","1995",
			"1996","1997","1998","1999","2000","2001","2002",
			"2003")] = "1994-2003"
		DAT$ipath_year[DAT$ipath_year %in% c("2004","2005",
			"2006","2007","2008","2009")] = "2004-2009"
		DAT$ipath_year[DAT$ipath_year %in% c("2010","2011",
			"2012","2013")] = "2010-2013"
		
		DAT$hist_grade[DAT$hist_grade %in% c("G1","G2",
			"GX")] = "G1/G2/GX"
		DAT$hist_grade[DAT$hist_grade %in% c("G3","G4")] = "G3/G4"
		
		DAT$hist_type[!grepl("Cervical Squa",DAT$hist_type)] = "NonCervSqua"
		
		DAT$clin_stage[DAT$clin_stage %in% c("II","III",
			"IV")] = "II/III/IV"
		
		DAT$BWA[DAT$BWA %in% c("VN:0.5.7-6","VN:0.5.9-r16")] = "nonVN:0.5.9"
		
		return(DAT)
		
	} else if( dataset == "COAD" ){
		DAT$hist_type[DAT$hist_type == "Colon Adenocarcinoma"] = "CA"
		# DAT$hist_type[DAT$hist_type == "Colon Mucinous Adenocarcinoma"] = "CMA"
		DAT$hist_type[grepl("Mucinous",DAT$hist_type)] = "MucAdeno"
		
		DAT$tumor_status2 = DAT$tumor_status
		DAT$tumor_status2[DAT$tumor_status2 == "[Discrepancy]"] = "WITH TUMOR"
		
		# Double check meaning of this variable: related to treatment? why CYT?
		DAT$ipath_year[DAT$ipath_year %in% c("1998","1999","2000","2001","2002")] = "1998-2002"
		DAT$ipath_year[DAT$ipath_year %in% c("2003","2004","2005","2006","2007","2008")] = "2003-2008"
		DAT$ipath_year[DAT$ipath_year %in% c("2009","2010")] = "2009-2010"
		DAT$ipath_year[DAT$ipath_year %in% c("2011","2012","2013")] = "2011-2013"
		
		return(DAT)
	
	} else if( dataset == "ESCA" ){
		DAT$tumor_stage[DAT$tumor_stage %in% c("I","II")] = "I/II"
		DAT$tumor_stage[DAT$tumor_stage %in% c("III","IV")
			| is.na(DAT$tumor_stage)] = "III/IV/NA"
		
		DAT$ipath_year[DAT$ipath_year %in% c("1998","1999",
			"2000","2001","2003")] = "1998-2003"
		DAT$ipath_year[DAT$ipath_year %in% c("2004","2005",
			"2006","2007","2008","2009")] = "2004-2009"
		DAT$ipath_year[DAT$ipath_year %in% c("2010","2011")] = "2010-2011"
		DAT$ipath_year[DAT$ipath_year %in% c("2012","2013")] = "2012-2013"
		
		DAT$hist_grade[DAT$hist_grade %in% c("G1","G2")] = "G1/G2"
		
		DAT$hist_type[grepl("Adeno",DAT$hist_type)] = "EsoAdeno"
		DAT$hist_type[grepl("Squam",DAT$hist_type)] = "EsoSquam"
		
		DAT$clin_stage[DAT$clin_stage %in% c("I","II")] = "I/II"
		
		return(DAT)
		
	} else if( dataset == "GBM" ){
		DAT$ipath_year[DAT$ipath_year %in% c("1989","1990",
			"1991","1992","1993","1994","1995",
			"1996","1997","1998","1999","2000","2001",
			"2002","2003","2004","2006","2007","2008")] = "1989-2008"
		DAT$ipath_year[DAT$ipath_year %in% c("2011","2012","2013")] = "2011-2013"
		
		DAT$BWA[DAT$BWA %in% c("VN:0.5.7-6","VN:0.5.9-r16")] = "VN:0.5.7-0.5.9"
		
		return(DAT)
	} else if( dataset == "HNSC" ){
		DAT$tumor_stage[DAT$tumor_stage %in% c("I","II")] = "I/II"
		
		DAT$ipath_year[DAT$ipath_year %in% c("1992","1993","1994",
			"1995","1996","1997","1998","1999","2000","2001","2002")] = "1992-2002"
		DAT$ipath_year[DAT$ipath_year %in% c("2003","2004","2005","2006")] = "2003-2006"
		DAT$ipath_year[DAT$ipath_year %in% c("2007","2008","2009")] = "2007-2009"
		DAT$ipath_year[DAT$ipath_year %in% c("2012","2013")] = "2012-2013"
		
		
		DAT$clin_stage[DAT$clin_stage %in% c("I","II")] = "I/II"
		DAT$clin_stage[DAT$clin_stage %in% c("IVA","IVB","IVC")] = "IV"
		
		DAT$hist_grade[DAT$hist_grade %in% c("G1","GX")] = "G1/GX"
		DAT$hist_grade[DAT$hist_grade %in% c("G3","G4")] = "G3/G4"
		
		return(DAT)
		
	} else if( dataset == "KIRC" ){
		DAT$tumor_stage[DAT$tumor_stage %in% c("I","II")] = "I/II"
		
		DAT$ipath_year[DAT$ipath_year %in% c("1998","1999","2000",
			"2001","2002")] = "1998-2002"
		DAT$ipath_year[DAT$ipath_year %in% c("2003","2004","2005")] = "2003-2005"
		DAT$ipath_year[DAT$ipath_year %in% c("2006","2007","2008",
			"2009","2010","2011","2012","2013")] = "2006-2013"
		
		DAT$hist_grade[DAT$hist_grade %in% c("G1","G2","GX")] = "G1/G2/GX"
		
		DAT$BWA[DAT$BWA %in% c("VN:0.5.9-r16","VN:0.5.9-tpx",
			"VN:0.6.2-r126")] = "VN:0.5.9-0.6.2"
		
		DAT$KIT[grepl("Nimblegen",DAT$KIT)] = "Nimblegen"
		
		return(DAT)
		
	} else if( dataset == "LGG" ){
		DAT$ipath_year[DAT$ipath_year %in% c("1992","1993",
			"1994","1995","1996","1997","1998")] = "1992-1998"
		DAT$ipath_year[DAT$ipath_year %in% c("1999","2000",
			"2001","2002","2003","2004","2005")] = "1999-2005"
		DAT$ipath_year[DAT$ipath_year %in% c("2006","2007",
			"2008","2009")] = "2006-2009"
		DAT$ipath_year[DAT$ipath_year %in% c("2010","2011",
			"2012","2013")] = "2010-2013"
		
		DAT$BWA[DAT$BWA %in% c("VN:0.5.7-6","VN:0.5.9-r16")] = "VN:0.5.7-0.5.9"
		
		return(DAT)
		
	} else if( dataset == "LIHC" ){
		DAT$tumor_stage[DAT$tumor_stage %in% c("III","IV")] = "III/IV"
		
		DAT$ipath_year[DAT$ipath_year %in% c("1995","1996",
			"1997","1998","1999","2000","2001","2002")] = "1995-2002"
		DAT$ipath_year[DAT$ipath_year %in% c("2003","2004",
			"2005","2006","2007")] = "2003-2007"
		DAT$ipath_year[DAT$ipath_year %in% c("2008","2009",
			"2010")] = "2008-2010"
		DAT$ipath_year[DAT$ipath_year %in% c("2011","2012",
			"2013")] = "2011-2013"
		
		DAT$hist_grade[DAT$hist_grade %in% c("G3","G4")] = "G3/G4"
		
		DAT$hist_type[grepl("Fibrola",DAT$hist_type)
			| grepl("Mixed",DAT$hist_type)] = "nonHepatocell"
		
		return(DAT)
		
	} else if( dataset == "LUAD" ){
		DAT$tumor_stage[DAT$tumor_stage %in% c("III","IV")] = "III/IV"
		
		DAT$ipath_year[DAT$ipath_year %in% c("1991","1992",
			"1993","1994","1995","1996","1997","1998","1999",
			"2000","2001")] = "1991-2001"
		DAT$ipath_year[DAT$ipath_year %in% c("2002","2003",
			"2004","2005")] = "2002-2005"
		DAT$ipath_year[DAT$ipath_year %in% c("2006","2007")] = "2006-2007"
		DAT$ipath_year[DAT$ipath_year %in% c("2008","2009")] = "2008-2009"
		DAT$ipath_year[DAT$ipath_year %in% c("2010","2011")] = "2010-2011"
		DAT$ipath_year[DAT$ipath_year %in% c("2012","2013")] = "2012-2013"
		
		return(DAT)
		
	} else if( dataset == "LUSC" ){
		DAT$tumor_stage[DAT$tumor_stage %in% c("III","IV")] = "III/IV"
		
		DAT$ipath_year[DAT$ipath_year %in% c("1992","1993",
			"1995","1996","1997","1998","1999","2000")] = "1992-2000"
		DAT$ipath_year[DAT$ipath_year %in% c("2001","2002")] = "2001-2002"
		DAT$ipath_year[DAT$ipath_year %in% c("2003","2004")] = "2003-2004"
		DAT$ipath_year[DAT$ipath_year %in% c("2005","2006")] = "2005-2006"
		DAT$ipath_year[DAT$ipath_year %in% c("2007","2008")] = "2007-2008"
		DAT$ipath_year[DAT$ipath_year %in% c("2009","2010")] = "2009-2010"
		DAT$ipath_year[DAT$ipath_year %in% c("2012","2013")] = "2012-2013"
		
		DAT$BWA[DAT$BWA %in% c("VN:0.5.7","VN:0.5.7-6")] = "VN:0.5.7-"
		
		return(DAT)
		
	} else if( dataset == "OV" ){
		DAT$ipath_year[DAT$ipath_year %in% c("1992","1993",
			"1994","1995","1996")] = "1992-1996"
		DAT$ipath_year[DAT$ipath_year %in% c("1997","1998")] = "1997-1998"
		DAT$ipath_year[DAT$ipath_year %in% c("2007","2008",
			"2009","2010","2011","2012","2013")] = "2007-2013"
		
		DAT$hist_grade[DAT$hist_grade %in% c("G1","G2","GB","GX")] = "G1/G2/GB/GX"
		DAT$hist_grade[DAT$hist_grade %in% c("G3","G4")] = "G3/G4"
		
		DAT$clin_stage[DAT$clin_stage %in% c("I","II")] = "I/II"
		
		return(DAT)
		
	} else if( dataset == "PAAD" ){
		DAT$tumor_stage[DAT$tumor_stage %in% c("III","IV")] = "III/IV"
		
		DAT$ipath_year[DAT$ipath_year %in% c("2001","2007",
			"2008")] = "2001-2008"
		DAT$ipath_year[DAT$ipath_year %in% c("2009","2010")] = "2009-2010"
		
		DAT$hist_grade[DAT$hist_grade %in% c("G3","G4","GX")] = "G3/G4/GX"
		
		DAT$hist_type[!grepl("Ductal Type",DAT$hist_type)
			& !is.na(DAT$hist_type)] = "Non-Ductal"
		
		return(DAT)
		
	} else if( dataset == "SARC" ){
		DAT$ipath_year[DAT$ipath_year %in% c("1994","1995",
			"1996","1997","1998","2000","2002","2003","2004",
			"2005")] = "1994-2005"
		DAT$ipath_year[DAT$ipath_year %in% c("2006","2007",
			"2008")] = "2006-2008"
		DAT$ipath_year[DAT$ipath_year %in% c("2009","2010")] = "2009-2010"
		DAT$ipath_year[DAT$ipath_year %in% c("2011","2012",
			"2013")] = "2011-2013"
		
		DAT$hist_type[grepl("Desmoid",DAT$hist_type)
			| grepl("Giant cell",DAT$hist_type)
			| grepl("Sheath Tumor",DAT$hist_type)
			| grepl("synovial; poor",DAT$hist_type)
			| grepl("Biphasic",DAT$hist_type)
			| grepl("Monophasic",DAT$hist_type)
			| grepl("Myxofib",DAT$hist_type)
			] = "DGSsBM"
		DAT$hist_type[grepl("Undifferentiated",DAT$hist_type)] = "Undiff"
		
		return(DAT)
		
	} else if( dataset == "SKCM" ){
		DAT$tumor_stage[DAT$tumor_stage %in% c("0","I")] = "I new"
		DAT$tumor_stage[DAT$tumor_stage %in% c("I/II NOS","II")] = "II new"
		DAT$tumor_stage[DAT$tumor_stage %in% c("III","IV")] = "III/IV"
		
		DAT$ipath_year[DAT$ipath_year %in% c("1978","1979",
			"1982","1984","1985","1986","1987","1989","1990",
			"1991","1992","1993","1994","1995","1996")] = "1978-1996"
		DAT$ipath_year[DAT$ipath_year %in% c("1997","1998",
			"1999","2000","2001")] = "1997-2001"
		DAT$ipath_year[DAT$ipath_year %in% c("2002","2003",
			"2004")] = "2002-2004"
		DAT$ipath_year[DAT$ipath_year %in% c("2005","2006")] = "2005-2006"
		DAT$ipath_year[DAT$ipath_year %in% c("2007","2008",
			"2009")] = "2007-2009"
		DAT$ipath_year[DAT$ipath_year %in% c("2010","2011",
			"2012","2013")] = "2010-2013"
		
		return(DAT)
		
	} else if( dataset == "STAD" ){
		DAT$tumor_stage[DAT$tumor_stage %in% c("I","II")] = "I/II"
		
		DAT$ipath_year[DAT$ipath_year %in% c("1996","2000",
			"2001","2002","2003","2004","2005","2006","2007")] = "1996-2007"
		DAT$ipath_year[DAT$ipath_year %in% c("2008","2009")] = "2008-2009"
		DAT$ipath_year[DAT$ipath_year %in% c("2012","2013")] = "2012-2013"
		
		DAT$hist_grade[DAT$hist_grade %in% c("G1","G2","GX")] = "G1/G2/GX"
		
		DAT$hist_type[grepl("Not Otherwise Specified",DAT$hist_type)] = "NotOtherSpec"
		DAT$hist_type[grepl("Diffuse Type",DAT$hist_type)
			| grepl("Papillary Type",DAT$hist_type)] = "Diff/Pap"
		DAT$hist_type[grepl("Signet Ring",DAT$hist_type)
			| grepl("Mucinous Type",DAT$hist_type)] = "Signet.Mucinous"
		DAT$hist_type[grepl("Intestinal",DAT$hist_type)] = "Intest,Tubular"
		
		DAT$BWA[DAT$BWA %in% c("VN:0.5.7","VN:0.5.7-6")] = "VN:0.5.7_type"
		DAT$BWA[DAT$BWA %in% c("VN:0.5.7_type","VN:0.5.9-r16")] = "non_VN:0.5.9-tpx"
		
		return(DAT)
		
	} else if( dataset == "UCEC" ){
		DAT$ipath_year[DAT$ipath_year %in% c("1995","1996",
			"1997","1999","2000","2001","2002","2003","2004",
			"2005","2006")] = "1995-2006"
		DAT$ipath_year[DAT$ipath_year %in% c("2007","2008",
			"2009")] = "2007-2009"
		DAT$ipath_year[DAT$ipath_year %in% c("2010","2011",
			"2012","2013")] = "2010-2013"
		
		DAT$hist_grade[DAT$hist_grade %in% c("G1","G2",
			"High Grade")] = "G1/G2/HighG"
		
		DAT$hist_type[grepl("Endometrioid endo",DAT$hist_type)] = "Endometrioid"
		DAT$hist_type[grepl("Mixed serous",DAT$hist_type)
			| grepl("Serous endo",DAT$hist_type)] = "Serous.Endo"
		
		DAT$clin_stage[DAT$clin_stage %in% c("I","II")] = "I/II"
		
		DAT$KIT[grepl("Nimblegen.SQEZ",DAT$KIT)] = "Nimblegen.SQEZ"
		
		return(DAT)
		
	}
	
	stop(sprintf("Write code for %s dataset",dataset))
	
}
RDA_calcTMB = function(DAT,mSPM){
	if(FALSE){
		DAT = dat; mSPM = matSPM
		
	}
	
	cat(sprintf("%s: Get TMB/log10_TMB variable ...\n",date()))
	int_subj = intersect(DAT$ID,colnames(mSPM))
	length(int_subj)
	DAT = DAT[which(DAT$ID %in% int_subj),]
	DAT = DAT[match(int_subj,DAT$ID),]
	mSPM = mSPM[,int_subj]
	DAT$TMB = colSums(mSPM)
	DAT$log10_TMB = log10(1 + DAT$TMB)
	DAT$log10_TMB_2 = DAT$log10_TMB^2
	DAT$TMB_bin = ifelse(DAT$log10_TMB > 2.5,1,0)
	DAT$int_TMB = DAT$log10_TMB * DAT$TMB_bin
	
	# Append top mutated genes
	tab = rowSums(mSPM)
	tab = tab[names(sort(tab,decreasing = TRUE)[1:10])]
	DAT = cbind(DAT,t(mSPM[names(tab),]))
	
	return(DAT)
}
RDA_getCYT = function(my_dirs,DAT){
	if(FALSE){
		dataset = "COAD"
		DAT = dat
	}
	
	cat(sprintf("%s: Get CYT variable ...\n",date()))
	cyt_fn = file.path(my_dirs$panc_dir,"mmc1 CYT.xlsx")
	cyt = readxl::read_excel(cyt_fn,sheet = 1)
	cyt = smart_df(cyt)
	cyt = name_change(cyt,"PatientID","ID")
	cyt = name_change(cyt,"Cytolytic.Activity","CYT")
	smart_table(cyt$HistoSubtype)
	smart_table(cyt$cancer)
	dim(cyt); cyt[1:3,]
	
	int_subj = intersect(cyt$ID,DAT$ID)
	length(int_subj)
	
	cyt = cyt[which(cyt$ID %in% int_subj),]
	smart_table(cyt$HistoSubtype)
	dim(cyt); cyt[1:3,]
	# dim(cyt2)
	
	DAT = smart_rmcols(DAT,"CYT")
	DAT = smart_merge(DAT,cyt[,c("ID","CYT")],all = TRUE)
	DAT$log10_CYT = log10(DAT$CYT)
	dim(DAT); DAT[1:3,]
	
	return(DAT)
}
RDA_makeCOVAR = function(DAT,MODEL){
	if(FALSE){
		DAT = dat
		MODEL = c("age","gender","tumor_stage","log10_TMB")
	}
	
	if( length(MODEL) == 1 && MODEL == "" ){
		XX = matrix(1,nrow = nrow(DAT))
		rownames(XX) = DAT$ID
		return(XX)
	}
	
	# Check all variables included
	if( !all(MODEL %in% names(DAT)) ){
		stop("Some variables missing from DAT")
	}
	
	# Create XX matrix
	MODEL_2 = c(MODEL,"IDnum")
	DAT$IDnum = seq(nrow(DAT))
	tmp_formula = as.formula(sprintf("~ %s",paste(MODEL_2,collapse=" + ")))
	COVARS = model.matrix(tmp_formula,data = DAT)
	XX = apply(COVARS,2,function(xx) as.numeric(xx))
	rownames(XX) = DAT$ID[XX[,"IDnum"]]
	XX = XX[,colnames(XX) != "IDnum"]
	XX[1:3,]
	any(is.na(XX))
	
	return(XX)
}
RDA_getOT_DIST = function(my_dirs,KWS = NULL,SUBJS){
	
	all_ot_fns = list.files(my_dirs$ithOT_dir,pattern = "rds")
	if( !is.null(KWS) ){
		for(KW in KWS){
			# KW = KWS[1]
			all_ot_fns = all_ot_fns[grepl(KW,all_ot_fns)]
			if( length(all_ot_fns) == 0 )
				stop("vector empty")
		}
	}
	
	# cat(sprintf("%s: Test all 7 gene_sim approaches!\n",date()))
	all_ot_fns = all_ot_fns[!grepl("ww1.rds",all_ot_fns)]
	all_ot_fns = all_ot_fns[!grepl("ww5.rds",all_ot_fns)]
	all_ot_fns = all_ot_fns[!grepl("ww6.rds",all_ot_fns)]
	all_ot_fns = all_ot_fns[!grepl("ww7.rds",all_ot_fns)]
	
	length(all_ot_fns); all_ot_fns # [1:5]
	if( length(all_ot_fns) == 0 ) stop("No OT files")
	
	# Get dimnames
	all_subjs = SUBJS
	num_subjs = length(all_subjs); num_subjs
	ARR = array(data = NA,
		dim = c(num_subjs,num_subjs,length(all_ot_fns)),
		dimnames = list(all_subjs,all_subjs,
			all_ot_fns))
	dim(ARR)
	
	for(fn in all_ot_fns){
		# fn = all_ot_fns[1]; fn
		tmp_rds = readRDS(file.path(my_dirs$ithOT_dir,fn))
		tmp_names = rownames(tmp_rds$DIST)
		setdiff(all_subjs,tmp_names)
		setdiff(tmp_names,all_subjs)
		tmp_names = intersect(tmp_names,all_subjs)
		length(tmp_names)
		ARR[tmp_names,tmp_names,fn] = tmp_rds$DIST[tmp_names,tmp_names]
		
		rm(tmp_rds)
	}
	# ARR[1:3,1:3,]
	
	return(ARR)
}
RDA_getEUC_DIST = function(mSPM,THRES = c(1,2,5,10,15,20,30,40,50,60)){
	if(FALSE){
		mSPM = matSPM
		THRES = tmp_GM
	}
	
	vnames = paste0("EUC.GM.",THRES)
	cat(sprintf("%s: Calculate euclidean distances ...\n",date()))
	# Calculate at different gene_freq thresholds (20,30,40)
	cnt = 1
	for(thres in THRES){
		# thres = THRES[1]; thres
		
		cat(sprintf("%s: GM = %s ...\n",date(),thres))
		mSPM_2 = mSPM
		
		# Subset genes mutated X times
		mSPM_2 = mSPM_2[rowSums(mSPM_2) >= thres,,drop = FALSE]
		
		# Remove subjects with no mutated genes
		mSPM_2 = mSPM_2[,colSums(mSPM_2) > 0,drop = FALSE]
		
		dim(mSPM); dim(mSPM_2)
		# tmp_dist2 = Rcpp_calc_EUC(ZZ = t(mSPM_2))
		tmp_dist = as.matrix(dist(t(mSPM_2),diag = TRUE,upper = TRUE))
		dim(tmp_dist); tmp_dist[1:3,1:3]
		
		if( cnt == 1 ){
			subjs = colnames(mSPM)
			nsubj = length(subjs)
			n_thres = length(THRES)
			
			EUC = array(data = NA,
				dim = c(nsubj,nsubj,n_thres),
				dimnames = list(subjs,subjs,vnames))
		}
		
		int_subj = intersect(colnames(tmp_dist),subjs)
		EUC[int_subj,int_subj,vnames[cnt]] = tmp_dist[int_subj,int_subj]
		cnt = cnt + 1
	}
	
	return(EUC)
	
}
new_ANA_COV = function(dataset,rm_subtype = FALSE,rm_CNA = FALSE){
	model_0 = c("age","BMI","BWA","clin_stage",
		"gender","height","hist_grade",
		"hist_type","ipath_year","KIT",
		"log10_TMB","ploidy","purity",
		"race","tumor_stage","weight")
	
	# Baseline covariates
	if( dataset == "BLCA" ){
		model = model_0
	} else if( dataset == "BRCA" ){
		extra_vars = c("PAM50")
		model = c(model_0,extra_vars)
		# "TP53"
	} else if( dataset == "CESC" ){
		model = model_0
		model = model[!(model %in% c("BMI","weight"))]
	} else if( dataset == "COAD" ){
		model = model_0
		model = model[!(model %in% c("weight","height","BMI","BWA"))]
	} else if( dataset == "ESCA" ){
		model = model_0
		model = model[!(model %in% c("weight","race"))]
	} else if( dataset == "GBM" ){
		model = model_0
	} else if( dataset == "HNSC" ){
		model = model_0
		# "TP53"
	} else if( dataset == "KIRC" ){
		extra_vars = c("barcode_SEQCENTER")
		model = c(model_0,extra_vars)
		model = model[!(model %in% c("BWA"))]
	} else if( dataset == "LGG" ){
		model = model_0
		# "IDH1","TP53","ATRX"
	} else if( dataset == "LIHC" ){
		model = model_0
		model = model[!(model %in% c("weight","height"))]
	} else if( dataset == "LUAD" ){
		model = model_0
		# "TP53"
	} else if( dataset == "LUSC" ){
		model = model_0
	} else if( dataset == "PAAD" ){
		model = model_0
		# "KRAS"
	} else if( dataset == "SARC" ){
		model = model_0
	} else if( dataset == "SKCM" ){
		model = model_0
	} else if( dataset == "STAD" ){
		model = model_0
		# "TP53"
	} else if( dataset == "UCEC" ){
		model = model_0
	} else {
		stop(sprintf("No covariates set for %s dataset",dataset))
		model = c("age","gender","tumor_stage","clin_stage",
			"hist_type","hist_grade","ipath_year","BWA",
			"SEQUENCER","KIT","CENTER","WGA_STATUS",
			"barcode_TSS","barcode_plate","barcode_analyte",
			"barcode_SEQCENTER","purity","ploidy",
			"log10_TMB","log10_TMB_2")
		
	}
	
	# Screen and remove tumor subtype info
	#		like tumor_grade, clin_stage, hist_type, hist_grade
	covars = list(model = model)
	if( rm_subtype ){
		for(OUT in names(covars)){
			tmp_vec = covars[[OUT]]
			tmp_vec = tmp_vec[!grepl("tumor_stage",tmp_vec)]
			tmp_vec = tmp_vec[!grepl("clin_stage",tmp_vec)]
			tmp_vec = tmp_vec[!grepl("hist_type",tmp_vec)]
			tmp_vec = tmp_vec[!grepl("hist_grade",tmp_vec)]
			tmp_vec = tmp_vec[!grepl("PAM50",tmp_vec)]
			covars[[OUT]] = tmp_vec
		}
	}
	if( rm_CNA ){
		for(OUT in names(covars)){
			tmp_vec = covars[[OUT]]
			tmp_vec = tmp_vec[!grepl("purity",tmp_vec)]
			tmp_vec = tmp_vec[!grepl("ploidy",tmp_vec)]
			covars[[OUT]] = tmp_vec
		}
	}
	
	return(covars)
}
get_fNULL_MOD = function(my_dirs,YY,keep_vars = NULL,
	rm_subtype = FALSE,rm_CNA = FALSE,add_genes = FALSE,
	strata = "all",verbose = TRUE){
	
	dataset = my_dirs$dataset
	dat_clin_fn = file.path(my_dirs$out_dir,"dat_clin.rds")
	if( !file.exists(dat_clin_fn) ){
		dat = RDA_prepCLIN(my_dirs = my_dirs)
		
		# Import matSPM
		mat_fn = file.path(my_dirs$spms_dir,
			sprintf("%s_mat.rds",dataset))
		matSPM = readRDS(mat_fn)
		tab = rowSums(matSPM)
		sort(tab,decreasing = TRUE)[1:10]
		
		# Calc TMB/MUT
		dat = RDA_calcTMB(DAT = dat,mSPM = matSPM); rm(matSPM)
		dat = RDA_getCYT(my_dirs = my_dirs,DAT = dat)
		rownames(dat) = dat$ID
		saveRDS(dat,dat_clin_fn)
	}
	DAT = readRDS(dat_clin_fn)
	dim(DAT)
	
	TYPE = "CONT"
	if( YY %in% c("OS","PFI","DSS","DFI") ){
		YY2 = sprintf("%s.time",YY)
		TYPE = "SURV"
		DAT = DAT[!is.na(DAT[,YY]),]
	} else if( YY %in% c("CYT") ){
		YY2 = "log10_CYT"
	}
	DAT = DAT[!is.na(DAT[,YY2]),]
	if( nrow(DAT) == 0 ) return(NULL)
	
	if( strata == "nHM" ){
		DAT = DAT[DAT$TMB_bin == 0,]
	}
	
	covars = new_ANA_COV(dataset = dataset,
		rm_subtype = rm_subtype,rm_CNA = rm_CNA); # covars
	start_model = covars$model; # start_model
	
	curr_model = start_model
	curr_model = unique(curr_model)
	curr_model = intersect(curr_model,names(DAT))
	cnt_uniq = apply(DAT[,curr_model],2,function(xx) length(unique(xx[!is.na(xx)])))
	cnt_uniq = cnt_uniq[cnt_uniq > 1]; # cnt_uniq
	curr_model = curr_model[curr_model %in% names(cnt_uniq)]
	curr_model = unique(c(curr_model,keep_vars))
	prop_nonNA = apply(DAT[,curr_model],2,function(xx) mean(!is.na(xx))); prop_nonNA
	prop_nonNA_thres = 0.7
	if( any(prop_nonNA < prop_nonNA_thres) ){
		rm_vars = names(prop_nonNA)[prop_nonNA < prop_nonNA_thres]; rm_vars
		curr_model = curr_model[!(curr_model %in% rm_vars)]
	}
	if( dataset %in% c("CESC","KIRC") ){
		curr_model = curr_model[curr_model != "KIT"]
	}
	
	if( add_genes ){
		mut_thres = 0.20
		tmp_genes = names(DAT)
		tmp_genes = tmp_genes[!(tmp_genes %in% c("CYT","log10_CYT"))]
		tmp_genes = tail(tmp_genes,n = 10)
		tmp_genes
		prop_mut = apply(DAT[,tmp_genes],2,function(xx) mean(xx)); prop_mut
		prop_mut = prop_mut[prop_mut >= mut_thres]; prop_mut
		if( length(prop_mut) > 0 ) curr_model = c(curr_model,names(prop_mut))
	}
	
	curr_model
	PVAL_thres = 1.5e-1; PVAL_thres
	
	while(TRUE){
		# curr_model = c(curr_model,"height")
		# curr_model = c("age","height","hist_grade","log10_TMB","ploidy","tumor_stage")
		if( !is.null(keep_vars) && all(curr_model %in% keep_vars) ) break
		if( length(curr_model) == 0 ) break
		
		if( TYPE == "CONT" ){
			glm_out = tryCatch(lm(as.formula(sprintf("%s ~ %s",
				YY2,paste(curr_model,collapse=" + "))),data = DAT),
				error = function(ee){NULL})
			if( is.null(glm_out) ) print("OLS is null")
			
			if( is.null(glm_out) ){
				DAT2 = DAT[,curr_model]
				for(vv in curr_model){
					# vv = curr_model[11]; vv
					nonNA_idx = !is.na(DAT2[,vv])
					DAT2 = DAT2[nonNA_idx,,drop = FALSE]
				}
				cnt_uniq = apply(DAT2,2,function(xx) length(unique(xx))); cnt_uniq
				if( any(cnt_uniq == 1) ){
					rm_var = names(cnt_uniq)[cnt_uniq == 1][1]; rm_var
					curr_model = curr_model[curr_model != rm_var]
					cat(sprintf("OUT = %s; Rm var = %s ...\n",YY,rm_var))
					next
				}
				stop("Check variables")
			}
			# summary(glm_out)
			
			cov_test = drop1(glm_out,.~.,test = "F")
			cov_test = as.matrix(cov_test)
			cov_test = cov_test[!is.na(cov_test[,ncol(cov_test)]),,drop = FALSE]
			# cov_test
		}
		if( TYPE == "SURV" ){
			glm_out = tryCatch(coxph(formula(sprintf("Surv(%s.time,%s) ~ %s",
				YY,YY,paste(curr_model,collapse = " + "))),data = DAT),
				error = function(ee){NULL},warning = function(ww){NULL})
			glm_out
			
			if( is.null(glm_out) ){
				DAT2 = DAT[,c(YY2,YY,curr_model)]
				for(vv in names(DAT2)){
					# vv = curr_model[11]; vv
					nonNA_idx = !is.na(DAT2[,vv])
					DAT2 = DAT2[nonNA_idx,,drop = FALSE]
				}
				cnt_uniq = apply(DAT2,2,function(xx) length(unique(xx))); cnt_uniq
				if( any(cnt_uniq == 1) ){
					rm_var = names(cnt_uniq)[cnt_uniq == 1][1]; rm_var
					curr_model = curr_model[curr_model != rm_var]
					cat(sprintf("OUT = %s; Rm var = %s ...\n",YY,rm_var))
					next
				} else {
					glm_out = suppressWarnings(coxph(formula(sprintf("Surv(%s.time,%s) ~ %s",
						YY,YY,paste(curr_model,collapse = " + "))),data = DAT))
					# Rm variable with largest pvalue
					tmp_mat = summary(glm_out)$coefficients; tmp_mat
					if( any(tmp_mat[,ncol(tmp_mat)] > 0.95) ){
						tmp_PVAL = tmp_mat[,ncol(tmp_mat)]; tmp_PVAL
						rm_var0 = names(tmp_PVAL)[which.max(tmp_PVAL)]; rm_var0
						tmp_vec = sapply(curr_model,function(xx) grepl(sprintf("^%s",xx),rm_var0)); tmp_vec
						rm_var = names(tmp_vec)[tmp_vec]; rm_var
						if( length(rm_var) > 1 ) stop("inspect this")
						curr_model = curr_model[curr_model != rm_var]
						cat(sprintf("OUT = %s; Rm var = %s ...\n",YY,rm_var))
						next
					}
					
				}
				stop("Check variables")
			}
			
			cov_test = tryCatch(car::Anova(glm_out,type = "III"),
				warning = function(ww){NULL}); cov_test
			if( is.null(cov_test) ){
				# drop variable with largest p-value
				cov_test = suppressWarnings(car::Anova(glm_out,type = "III"))
				rm_var = rownames(cov_test)[which.max(cov_test[,ncol(cov_test)])]
				curr_model = curr_model[curr_model != rm_var]
				next
			}
			
			cov_test = as.matrix(cov_test)
		}
		
		# Remove vars to keep
		cov_test2 = cov_test[!(rownames(cov_test) %in% keep_vars),,drop = FALSE]
		# rm(cov_test)
		# print(cov_test2)
		if( verbose ) cat(sprintf("%s; ",paste(curr_model,collapse = " + ")))
		if( nrow(cov_test2) == 0 ) break
		
		# Find var to remove
		if( any(cov_test2[,ncol(cov_test2)] > PVAL_thres) ){
			tmp_PVAL = cov_test2[,ncol(cov_test2)]
			names(tmp_PVAL) = rownames(cov_test2); tmp_PVAL
			# rm(cov_test2)
			rm_var = names(tmp_PVAL)[which.max(tmp_PVAL)]; rm_var
			# print(cov_test)
			cat(sprintf("OUT = %s; Rm var = %s ...\n",YY,rm_var))
			curr_model = curr_model[curr_model != rm_var]
		} else {
			# End if all p-values less than threshold
			break
		}
		
	}
	if( verbose ){
		cat("\n")
		print(cov_test)
		if( TYPE == "SURV" ){
			DAT2 = DAT[,c(YY2,YY,curr_model)]
			for(vv in names(DAT2)){
				# vv = curr_model[11]; vv
				nonNA_idx = !is.na(DAT2[,vv])
				DAT2 = DAT2[nonNA_idx,,drop = FALSE]
			}
			print(smart_table(DAT2[[YY]]))
		}
		# print(smart_table(!is.na(DAT[[YY2]])))
	}
	
	return(list(model = curr_model,DAT = DAT))
	
}
RDA_REG_final = function(my_dirs,DAT,YY,OT,EUC,
	MODEL,nPERM = 5e3){
	
	if(FALSE){
		DAT = mod_out$DAT; YY = YY_2; OT = ot;
		EUC = euc; MODEL = mod_out$model; nPERM = nPERM
		
	}
	
	all_LABS = c("EUC","OT")
	
	# Subset to GMs >= 2 and GMs <= 30
	nms = dimnames(OT)[[3]]
	nms = nms[!grepl("GM.60",nms)]
	nms = nms[!grepl("GM.50",nms)]
	nms = nms[!grepl("GM.40",nms)]
	OT = OT[,,nms,drop = FALSE]
	
	nms = dimnames(EUC)[[3]]
	nms = nms[!grepl("GM.60",nms)]
	nms = nms[!grepl("GM.50",nms)]
	nms = nms[!grepl("GM.40",nms)]
	EUC = EUC[,,nms,drop = FALSE]
	
	# Check missingness
	cat(sprintf("%s: Check missingness ...\n",date()))
	int_subjs = intersect(DAT$ID,dimnames(OT)[[1]])
	length(int_subjs)
	OT = OT[int_subjs,int_subjs,,drop = FALSE]
	DAT = DAT[which(DAT$ID %in% int_subjs),]
	for(nm in dimnames(OT)[[3]]){
		# nm = dimnames(OT)[[3]][1]; nm
		if( !any(is.na(OT)) ) break
		keep_subjs = dimnames(OT)[[1]]
		miss_subjs = keep_subjs[is.na(diag(OT[,,nm]))]
		keep_subjs = keep_subjs[!(keep_subjs %in% miss_subjs)]
		length(keep_subjs)
		OT = OT[keep_subjs,keep_subjs,,drop = FALSE]
	}
	
	int_subjs = intersect(DAT$ID,dimnames(OT)[[1]])
	int_subjs = intersect(int_subjs,dimnames(EUC)[[1]])
	length(int_subjs)
	EUC = EUC[int_subjs,int_subjs,,drop = FALSE]
	DAT = DAT[which(DAT$ID %in% int_subjs),]
	for(nm in dimnames(EUC)[[3]]){
		# nm = dimnames(EUC)[[3]][1]; nm
		if( !any(is.na(EUC)) ) break
		keep_subjs = dimnames(EUC)[[1]]
		miss_subjs = keep_subjs[is.na(diag(EUC[,,nm]))]
		keep_subjs = keep_subjs[!(keep_subjs %in% miss_subjs)]
		length(keep_subjs)
		EUC = EUC[keep_subjs,keep_subjs,,drop = FALSE]
	}
	
	# Create XX matrix
	TYPE = "CONT"
	DAT = DAT[!is.na(DAT[,YY]),]
	if( YY %in% c("OS","PFI","DSS","DFI") ){
		YY2 = sprintf("%s.time",YY)
		DAT = DAT[!is.na(DAT[,YY2]),]
		TYPE = "SURV"
	}
	XX = RDA_makeCOVAR(DAT = DAT,MODEL = MODEL)
	dim(XX)
	
	# Intersect/Align subjects
	int_subjs = intersect(rownames(XX),dimnames(OT)[[1]])
	int_subjs = intersect(int_subjs,dimnames(EUC)[[1]])
	length(int_subjs)
	OT 	= OT[int_subjs,int_subjs,,drop = FALSE]
	XX 	= XX[int_subjs,,drop = FALSE]
	DAT = DAT[which(DAT$ID %in% int_subjs),]
	DAT = DAT[match(int_subjs,DAT$ID),]
	EUC = EUC[int_subjs,int_subjs,,drop = FALSE]
	
	cnt = 1; tot = length(all_LABS)
	cat(sprintf("%s: Total Distances = %s ...\n",date(),tot))
	res = c()
	
	for(LAB in all_LABS){
		# LAB = all_LABS[2]
		
		# Get kernels
		KK = list()
		if( LAB == "EUC" ){
			for(nam in dimnames(EUC)[[3]]){
				tmp_dist = EUC[,,nam]
				KK[[nam]] = MiRKAT::D2K(D = tmp_dist)
				rm(tmp_dist)
			}
		} else if( LAB == "OT" ){
			for(nam in dimnames(OT)[[3]]){
				tmp_dist = OT[,,nam]
				KK[[nam]] = MiRKAT::D2K(D = tmp_dist)
				rm(tmp_dist)
			}
			
			tmp_LAB = names(KK)
			tmp_LAB = gsub("highLAM-lowTMB_","",tmp_LAB)
			tmp_LAB = gsub("BAL.FALSE_L1.0.5_L2.0.5","LAM.0.5",tmp_LAB)
			tmp_LAB = gsub("BAL.FALSE_L1.1_L2.1","LAM.1",tmp_LAB)
			tmp_LAB = gsub("BAL.FALSE_L1.5_L2.5","LAM.5",tmp_LAB)
			tmp_LAB = gsub("BAL.TRUE_L1.1_L2.1","LAM.Inf",tmp_LAB)
			tmp_LAB = gsub("_ww2.rds","_GO",tmp_LAB)
			tmp_LAB = gsub("_ww3.rds","_PATH",tmp_LAB)
			tmp_LAB = gsub("_ww4.rds","_ME",tmp_LAB)
			tmp_LAB
			names(KK) = tmp_LAB
		}
		
		# Count subjects
		Ns = nrow(DAT); Ns
		Ne = ifelse(YY %in% c("OS","PFI","DSS","DFI"),
			sum(DAT[,YY]),Ns); Ne
		if( Ne < 50 ){
			print(table(DAT[,YY]))
			VARS = c(MODEL,YY)
			if( YY %in% c("OS","PFI","DSS","DFI") ) VARS = c(VARS,sprintf("%s.time",YY))
			print(apply(DAT[,VARS],2,function(xx) mean(!is.na(xx))))
			stop("Check model variable missingness")
		}
		
		XX_2 = NULL
		if( ncol(XX) > 1 ){
			XX_2 = XX[,-1,drop = FALSE]
		}
		
		if( TYPE == "SURV" ){ # Run survival CoxPH regression
			cout = coxph(Surv(DAT[,YY2],DAT[,YY]) ~ .,
				data = smart_df(XX_2))
			RESI = as.numeric(cout$residuals)
			names(RESI) = int_subjs
			fit_perm = kernTEST(RESI = RESI,KK = KK,nPERMS = nPERM,
				iter1 = 2.5e2,iter2 = 5e3,verbose = TRUE)
			
		} else { # Run continuous regression
			
			# Fit null model and get residuals
			lm_out = lm(DAT[,YY] ~ .,data = smart_df(XX_2))
			RESI = as.numeric(lm_out$residuals)
			names(RESI) = int_subjs
			fit_perm = kernTEST(RESI = RESI,KK = KK,nPERMS = nPERM,
				iter1 = 2.5e2,iter2 = 5e3,verbose = TRUE)
			
		}
		
		if( LAB == "EUC" ){
			
			# Omnibus result
			res = rbind(res,smart_df(TYPE = TYPE,
				OUT = YY,MODEL = paste(MODEL,collapse = ","),
				LAB = sprintf("%s_omnibus",LAB),
				PVAL_perm = fit_perm$omni_PVAL,Ns = Ns,Ne = Ne))
			
			# Per kernel results
			res = rbind(res,smart_df(TYPE = TYPE,
				OUT = YY,MODEL = paste(MODEL,collapse = ","),
				LAB = names(KK),PVAL_perm = as.numeric(fit_perm$PVALs),
				Ns = Ns,Ne = Ne))
			
		} else if( LAB == "OT" ){
			
			# Omnibus result
			res = rbind(res,smart_df(TYPE = TYPE,
				OUT = YY,MODEL = paste(MODEL,collapse = ","),
				LAB = sprintf("%s_omnibus",LAB),
				PVAL_perm = fit_perm$omni_PVAL,Ns = Ns,Ne = Ne))
			
			# Per kernel results
			res = rbind(res,smart_df(TYPE = TYPE,
				OUT = YY,MODEL = paste(MODEL,collapse = ","),
				LAB = tmp_LAB,PVAL_perm = as.numeric(fit_perm$PVALs),
				Ns = Ns,Ne = Ne))
			
		}
		
		rm(fit_perm,KK)
	}
	
	return(res)
	
}

#' @title ANA_TEST_final
#' @param my_dirs A list generated by \code{setdirs()}
#'	of directories specific to a cancer type
#' @param nPERM A positive integer specifying
#'	the number of permutations to calculate
#' 	permutation-based individual kernel and omnibus
#'	p-values
#' @export
ANA_TEST_final = function(my_dirs,nPERM = 1e5){
	
	dataset = my_dirs$dataset
	setwd(my_dirs$ithOT_dir)
	
	# Get all completed OT thresholds
	OT_fns = list.files(my_dirs$ithOT_dir)
	OT_fns = OT_fns[!grepl("ww1.rds",OT_fns)];
	OT_fns = OT_fns[!grepl("ww5.rds",OT_fns)]
	OT_fns = OT_fns[!grepl("ww6.rds",OT_fns)]
	OT_fns = OT_fns[!grepl("ww7.rds",OT_fns)]
	# OT_fns
	GMs = sapply(OT_fns,function(xx){
		# xx = OT_fns[1]
		xx2 = strsplit(xx,"_")[[1]][2]
		xx2 = gsub("GM.","",xx2)
		xx2
	},USE.NAMES = FALSE) # we expect 3 gene sims, 4 lambdas
	tab = table(GMs); tab = tab[tab == 12]; tab
	GMs = sort(as.integer(names(tab))); GMs
	
	# Import SPMs, calculate EUC dist for various thresholds
	mat_fn = file.path(my_dirs$spms_dir,
		sprintf("%s_mat.rds",dataset))
	matSPM = readRDS(mat_fn)
	
	# Calculate euclidean distances
	euc = RDA_getEUC_DIST(mSPM = matSPM,THRES = GMs)
	rm(matSPM)
	
	# Get all OT
	ot = RDA_getOT_DIST(my_dirs = my_dirs,
		SUBJS = dimnames(euc)[[1]])
	GMs = sapply(dimnames(ot)[[3]],function(xx){
		xx2 = strsplit(xx,"_")[[1]][2]
		xx2 = gsub("GM.","",xx2)
		xx2
	},USE.NAMES = FALSE) # we expect 3 gene sims, 4 lambdas
	tab = table(GMs); tab
	
	# Regression constants
	all_strat = "all"
	if( dataset %in% c("COAD") ) all_strat = c("all","nHM")
	all_YY = c("CYT","OS","PFI","DSS")
	
	## Create dat_clin_KWS.rds image file
	dat_clin_fn = file.path(my_dirs$out_dir,"dat_clin.rds")
	if( file.exists(dat_clin_fn) ) unlink(dat_clin_fn)
	
	reg_fn = file.path(my_dirs$reg_dir,"reg.rds")
	if( file.exists(reg_fn) ) unlink(reg_fn)
	fin_reg = c()
	for(YY in all_YY){
		# YY = all_YY[1]; YY
		YY_2 = YY
		if( YY == "CYT" ){
			YY_2 = "log10_CYT"
		}
		if( YY == "OS" && dataset %in% c("CESC","SARC") ) next
		if( YY == "PFI" && dataset %in% c("CESC") ) next
		if( YY == "DSS" && dataset %in% c("BRCA","CESC","COAD","ESCA",
			"KIRC","LIHC","SARC","UCEC") ) next
	for(strat in all_strat){
		# strat = all_strat[1]
		cat(sprintf("%s: YY = %s; STRAT = %s ...\n",date(),YY,strat))
		rm_subtype = FALSE; rm_CNA = FALSE
		
		# Set null model
		mod_out = get_fNULL_MOD(my_dirs = my_dirs,YY = YY,
			keep_vars = "log10_TMB",rm_subtype = rm_subtype,
			rm_CNA = rm_CNA,add_genes = FALSE,strata = strat,
			verbose = TRUE)
		if( is.null(mod_out) ) next
		
		# Run kernel regression
		tmp_reg = RDA_REG_final(my_dirs = my_dirs,DAT = mod_out$DAT,
			YY = YY_2,OT = ot,EUC = euc,MODEL = mod_out$model,
			nPERM = nPERM)
		
		# Append
		tmp_df = smart_df(STRAT = strat,tmp_reg)
		print(tmp_df)
		fin_reg = rbind(fin_reg,tmp_df)
		rm(tmp_df)
		
	}}
	print(fin_reg[which(fin_reg$LAB %in% c("EUC","OT_omnibus")),])
	saveRDS(fin_reg,reg_fn)
	unlink(dat_clin_fn)
	
	return(NULL)
	
}

#' @title new_ANA_AGG
#' @param my_dirs A list generated by \code{setdirs()}
#'	of directories specific to a cancer type
#' @export
new_ANA_AGG = function(my_dirs){
	
	res_fn = file.path(my_dirs$fin_dir,"pan_kern.rds")
	
	# Gather all tumor types
	datasets = c("BLCA","BRCA","CESC","COAD",
							"ESCA","GBM","HNSC","KIRC",
							"LGG","LIHC","LUAD","LUSC",
							"PAAD","SARC","SKCM","STAD",
							"UCEC")
	
	if( !file.exists(res_fn) ){
		res = c()
		for(dataset in datasets){
			# dataset = datasets[1]; dataset
			
			cat(sprintf("%s ",dataset))
			my_dirs1 = setdirs(git_dir = my_dirs$git_dir,
				work_dir = my_dirs$work_dir,dataset = dataset,
				verbose = FALSE)
			
			tum_fns = list.files(my_dirs1$reg_dir)
			tum_fns = tum_fns[grepl(".rds$",tum_fns)]
			# tum_fns = tum_fns[grepl("^reg_",tum_fns)]
			if( length(tum_fns) == 0 ) next
			tum_res = c()
			
			for(fn in tum_fns){
				# fn = tum_fns[1]
				tmp_df = readRDS(file.path(my_dirs1$reg_dir,fn))
				tum_res = rbind(tum_res,tmp_df); rm(tmp_df)
			}
			
			res = rbind(res,smart_df(dataset = dataset,tum_res))
			rm(tum_res)
		}
		cat("\n")
		saveRDS(res,res_fn)
	}
	
	# Import
	res = readRDS(res_fn)
	res$YY = -log10(res$PVAL_perm)
	res$GM = NA
	idx = grep("EUC.GM.",res$LAB)
	res$GM[idx] = as.integer(gsub("EUC.GM.","",res$LAB[idx]))
	idx = grep("OT_GM.",res$LAB)
	res$GM[idx] = as.integer(sapply(res$LAB[idx],function(xx){
		# xx = res$LAB[idx][1]; xx
		gsub("GM.","",strsplit(xx,"_")[[1]][2])
	},USE.NAMES = FALSE))
	
	# Plotting
	dim(res); res[1:3,]
	
	# Get covariates selected per null model
	if( TRUE ){
		res2 = res[grepl("log10_TMB",res$MODEL),]
		res2 = res2[which(res2$GM <= 30 & res2$OUT != "DFI"),]
		ures = unique(res2[,c("dataset","STRAT","OUT","MODEL")]); # ures
		all_vars = sapply(ures$MODEL,function(xx)
			strsplit(xx,",")[[1]],USE.NAMES = FALSE)
		all_vars = sort(unique(unlist(all_vars)))
		all_vars = all_vars[all_vars != "log10_TMB"]
		all_vars
		for(vv in all_vars){
			ures[[vv]] = ifelse(grepl(vv,ures$MODEL),1,0)
		}
		dim(ures)
		ures[1:5,]
		ures$OUT2 = ures$OUT
		ures$OUT2[grepl("CYT",ures$OUT2)] = "CYT"
		MAT = as.matrix(ures[,all_vars])
		colnames(MAT)[grepl("SEQCENTER",colnames(MAT))] = "CENTER"
		MAT = MAT[,sort(colnames(MAT))]
		rownames(MAT) = sprintf("%s + %s",ures$dataset,ures$OUT2)
		idx = which(ures$dataset == "COAD")
		rownames(MAT)[idx] = sprintf("%s + %s:%s",ures$dataset,ures$OUT2,ures$STRAT)[idx]
		dim(MAT); MAT[1:4,]
		png(file.path(my_dirs$fin_dir,"null_models.png"),units = "px",
			height = 2000,width = 1600,res = 250,type = "cairo",pointsize = 20)
		bb = smart_heatmap(MAT = MAT,main = "Null Models",
			width = c(0.7,3,1),height = c(0.15,0.3,3,1),
			GRID = list(GRID = TRUE,lwd = 0.5),clustRC = rep(FALSE,2),
			nodePar_row = list(xx_shift = 0,xx = 0,adj = 0,lab.cex = 0.75),
			nodePar_col = list(xx_shift = -0.02,srt = 45,adj = 0,lab.cex = 0.65),
			make_key = FALSE,vec_cols = c("white","red")); rm(bb)
		dev.off()
		
	}
	
	# Get TMB distribution and num genes per cutoff
	if( TRUE ){
		
		# Num genes per cutoff
		cnts = c(2,5,10,15,20,30)
		mat = matrix(NA,length(datasets),length(cnts),
			dimnames = list(datasets,sprintf("GM.%s",cnts)))
		for(dataset in datasets){
			# dataset = datasets[1]
			my_dirs = setdirs(dataset = dataset,
				verbose = FALSE)
			tab = make_tb1(my_dirs = my_dirs,
				rmWGA = TRUE,show = FALSE)
			# tab[1:10]
			mat[dataset,] = sapply(cnts,function(xx) sum(tab >= xx))
		}
		mat = formatC(x = mat,format = "d",big.mark = ",")
		mat = smart_df(mat)
		names(mat) = gsub("GM.","$\\\\geq$",names(mat))
		mat[["Abbreviation"]] = rownames(mat)
		tmpd = name_TUMORTYPES()
		tmpd = tmpd[match(mat[["Abbreviation"]],tmpd$V1),]
		mat[["Tumor Type"]] = tmpd$V2
		rownames(mat) = NULL
		mat = mat[,c(ncol(mat),ncol(mat)-1,seq(ncol(mat)-2))]
		mat
		cap = "The number of genes at various"
		cap = sprintf("%s mutation frequency cutoffs",cap)
		cap = sprintf("%s for each tumor type.",cap)
		print_latex_table(DATA = mat,add_table = TRUE,
			my_align = "llrrrrrr",label = "gene_mut_freq",
			fontsize = "footnotesize",caption = cap)
		
		
		# TMB distribution
		dat = c()
		for(dataset in datasets){
			# dataset = datasets[1]
			my_dirs = setdirs(dataset = dataset,
				verbose = FALSE)
			mat_fn = file.path(my_dirs$spms_dir,
				sprintf("%s_mat.rds",my_dirs$dataset))
			matSPM = readRDS(mat_fn)
			dim(matSPM); matSPM[1:5,1:4]
			dat = rbind(dat,smart_df(dataset = dataset,
				TMB = colSums(matSPM)))
		}
		
		my_theme = theme(legend.position = c("none","bottom")[1],
			panel.background = element_blank(),
			panel.border = element_rect(color = "black",fill = NA,size = 1),
			# panel.spacing = unit(0.75,"lines"),
			panel.grid.major = element_line(colour = "grey50",
				size = 0.5,linetype = "dotted"),
			panel.spacing.x = unit(1.0,"lines"),
			# panel.spacing.y = unit(1.00,"lines"),
			# legend.text = element_text(size = 40),
			# axis.text.x = element_text(size = 35),
			# strip.text = element_text(size = 40),
			# strip.text.y = element_text(angle = 0,hjust = 0)
			text = element_text(size = 45)
			)
		
		TMB = NULL
		gg = ggplot(data = dat,aes(x = log10(TMB))) +
			geom_histogram(bins = 50) + facet_wrap(~ dataset) +
			xlab("log10(Tumor Mutation Burden)") + 
			ylab("Frequency") + my_theme
		
		png_fn = file.path(my_dirs$fin_dir,"log10_TMB.png")
		ggsave(png_fn,plot = gg,device = "png",width = 35,height = 25,
			units = "in")
		rm(gg)
		
		
	}
	
	# Aggregate scatter plot
	if( TRUE ){
		res2 = res[grepl("log10_TMB",res$MODEL),]
		table(res$GM)
		# res2 = res2[which(res2$GM <= 30),]
		# res2 = res2[which(res2$GM %in% sprintf(">= %s",c(5,10,15,20,30,40))),]
		# res2 = res2[which(!(res2$dataset == "COAD" & res2$STRAT == "all")),]
		res2 = res2[which(res2$LAB %in% c("EUC_omnibus","OT_omnibus")),]
		dim(res); dim(res2); res2[1:10,]
		
		dat = res2[which(grepl("EUC",res2$LAB)),]
		dat = name_change(dat,"YY","YY_EUC")
		dat = smart_rmcols(dat,c("PVAL_perm","GM","LAB"))
		dat = smart_merge(dat,smart_rmcols(res2[which(grepl("OT_",res2$LAB)),],
			c("PVAL_perm","GM","LAB")))
		dat = name_change(dat,"YY","YY_OT")
		dim(dat); dat[1:10,]
		
		my_theme = theme(legend.position = c("none","right","bottom")[2],
			text = element_text(size = 30),
			axis.text = element_text(size = 20),
			panel.background = element_blank(),
			# panel.spacing = unit(0.75,"lines"),
			panel.grid.major = element_line(colour = "grey50",
				size = 0.5,linetype = "dotted"))
		
		# table(dat$dataset)
		dat$dataset2 = dat$dataset
		idx = which(dat$dataset == "COAD" & dat$STRAT == "nHM")
		dat$dataset2[idx] = sprintf("%s.%s",dat$dataset,dat$STRAT)[idx]
		tmp_lev = sort(unique(dat$OUT)); # tmp_lev
		tmp_lev = tmp_lev[c(2,3,4,1)]; # tmp_lev
		dat$OUT = factor(dat$OUT,levels = tmp_lev,
			labels = c("Cytolytic Activity","Overall Survival",
				"Progression-Free Interval","Disease-Specific Survival"))
		
		thres = -log10(0.05)
		YY_EUC = YY_OT = OUT = dataset2 = NULL
		gg = ggplot(data = dat,mapping = aes(x = YY_EUC,y = YY_OT,group = dataset2)) +
			geom_point(alpha = 0.85,size = 5,aes(color = dataset2,shape = dataset2)) +
			geom_abline(intercept = 0,slope = 1) + # xlim(rr) + ylim(rr) +
			geom_segment(size = 2,linetype = 2,color = "red",aes(x = 0,y = thres,
				xend = thres,yend = thres)) +
			geom_segment(size = 2,linetype = 2,color = "red",aes(x = thres,y = thres,
				xend = thres,yend = 0)) +
			scale_color_manual(values = rep(smart_colors(9),2)) +
			scale_shape_manual(values = sort(rep(c(16,17),9))) +
			# facet_wrap(~ OUT,scales = "free") + 
			facet_wrap(~ OUT) + 
			xlab("-log10(Euclidean Omnibus)") + 
			ylab("-log10(Optimal Transport Omnibus)") + labs(color = "Dataset",shape = "Dataset") +
			guides(colour = guide_legend(override.aes = list(size = 5))) +
			my_theme
		png_fn = file.path(my_dirs$fin_dir,"overall_omnibus_scatter.png")
		ggsave(png_fn,plot = gg,device = "png",width = 13,height = 8,
			units = "in"); # rm(gg)
		
	}
	
	# KIRC/LGG/LUAD plot
	if( TRUE ){
		res[1:3,]
		res2 = res[which( ( (res$dataset == "KIRC" & res$OUT == "PFI")
			| (res$dataset == "LGG" & res$OUT %in% c("log10_CYT","OS"))
			| (res$dataset == "LUAD" & res$OUT == "DSS") )
			& !grepl("omnibus",res$LAB)
			& res$GM <= 30 & res$GM > 1),]
		dim(res2); res2[1:5,]
		
		my_theme = theme(legend.position = "bottom",
			text = element_text(size = 40),
			plot.title = element_text(hjust = 0.5),
			panel.grid.major = element_line(colour = "grey50",
				size = 0.5,linetype = "dotted"),
			axis.title.x = element_text(vjust = -0.5),
			axis.text.x = element_text(size = 32),
			panel.background = element_blank(),
			panel.spacing.x = unit(0.75,"lines"),
			strip.text.y = element_text(angle = 0,hjust = 0))
		
		res2 = res2[order(res2$GM),]
		res2$OUT = ifelse(res2$OUT == "log10_CYT","Cytolytic~Activity",
			ifelse(res2$OUT == "OS","Overall~Survival",
			ifelse(res2$OUT == "PFI","Progression-Free~Interval",
				"Disease-Specific~Survival")))
		uOUT = sort(unique(res2$OUT)); uOUT
		mOUT = c("Cytolytic~Activity","Overall~Survival",
			"Progression-Free~Interval","Disease-Specific~Survival")
		mOUT = mOUT[mOUT %in% uOUT]; mOUT
		res2$OUT = factor(res2$OUT,levels = mOUT)
		
		res2$LAM = sapply(res2$LAB,function(xx){
			if( grepl("EUC",xx) ){
				"EUC"
			} else {
				# xx = res2$LAB[2]; xx
				LAM = gsub("LAM.","",strsplit(xx,"_")[[1]][3])
				sprintf("OT (lambda = %s)",LAM)
			}
		},USE.NAMES = FALSE)
		res2$GS = sapply(res2$LAB,function(xx){
			if( grepl("EUC",xx) ){
				"EUC"
			} else if( grepl("GO",xx) ){
				"GO"
			} else if( grepl("PATH",xx) ){
				"PATH"
			} else if( grepl("ME",xx) ){
				"ME"
			} else {
				xx2 = strsplit(xx,"_")[[1]]
				xx2 = gsub(".rds","",xx2[length(xx2)])
				xx2 = gsub("ww","",xx2)
				ifelse(grepl("2",xx2),"GO",
					ifelse(grepl("3",xx2),"PATH","ME"))
			}
		},USE.NAMES = FALSE)
		
		max_yy = max(c(-log10(0.05),res2$YY)) + 0.5
		if( any(is.infinite(max_yy)) ){
			cat(sprintf("Skip %s b/c perm_pval = 0\n",dataset))
		}
		max_yy
		
		lev_LAM = sort(unique(res2$LAM)); lev_LAM
		res2$LAM = factor(res2$LAM,levels = lev_LAM[c(1,2,3,4,5)],
			labels = c("Euclidean",expression(paste("OT (",lambda," = 0.5)")),
				expression(paste("OT (",lambda," = 1.0)")),
				expression(paste("OT (",lambda," = 5.0)")),
				# expression(paste("OT (",lambda," = ",infinity,")"))
				"OT (Balanced)"))
		lev_GS = sort(unique(res2$GS)); lev_GS
		res2$GS = factor(res2$GS,levels = lev_GS[c(1,2,4,3)],
			labels = c("Euclidean","GO","PATH","ME"))
		res2[1:3,]
		
		res2$OUT2 = sprintf("%s + %s",res2$dataset,res2$OUT); table(res2$OUT2)
		GM = YY = GS = OUT2 = LAM = NULL
		gg = ggplot(data = res2,mapping = aes(x = factor(GM),y = YY,fill = GS)) +
			geom_bar(stat = "identity",position = position_dodge(),width = 0.8) +
			geom_hline(yintercept = -log10(0.05),size = 1.5,linetype = 2) +
			#geom_hline(data = tmp_df,size = 1.2,linetype = 2,
			#	mapping = aes(yintercept = YY,color = factor(GM_0))) +
			xlab("Gene Mutation Frequency Cutoff") + # ggtitle(dataset) +
			ylab("-log10(Permutation P-value)") + labs(fill = "Gene Similarity Method") + 
			facet_grid(OUT2 ~ LAM,labeller = label_parsed) +
			ylim(c(0,max_yy)) + guides(color = "none") + my_theme
		# gg
		png_fn = file.path(my_dirs$fin_dir,"top_tumor_regs.png")
		ggsave(png_fn,plot = gg,device = "png",width = 30,height = 20,units = "in")
		
		
	}
	
	if( TRUE ){ # Table of sample sizes
		ures = unique(res[,c("dataset","STRAT","OUT","Ns","Ne")])
		ures2 = unique(res[,c("dataset","STRAT")])
		rownames(ures2) = NULL
		ures2$CYT = ""; ures2$OS = ""; ures2$PFI = ""; ures2$DSS = ""
		for(ii in seq(nrow(ures2))){
			# ii = 2
			dd = ures2$dataset[ii]
			ss = ures2$STRAT[ii]
			idx1 = (ures$dataset == dd & ures$STRAT == ss)
			idx2 = (ures2$dataset == dd & ures2$STRAT == ss)
			ures[idx1,]
			ures2[idx2,]
			
			if( length(grep("CYT",ures$OUT[idx1])) > 0 ){
				idx3 = idx1 & grepl("CYT",ures$OUT)
				ures2$CYT[idx2] = sprintf("%s",
					ures$Ns[idx3])
			}
			
			for(tt in c("OS","PFI","DSS")){
				# tt = c("OS","PFI","DSS")[1]; tt
				if( length(grep(tt,ures$OUT[idx1])) > 0 ){
					idx3 = idx1 & grepl(tt,ures$OUT)
					ures2[[tt]][idx2] = sprintf("%s (%s)",
						ures$Ns[idx3],ures$Ne[idx3])
				}
			}
		}
		ures2 = name_change(ures2,"dataset","Tumor")
		ures2 = name_change(ures2,"STRAT","Strata")
		ures2$Strata = ifelse(ures2$Strata == "all","ALL","nHM")
		ures2$ABBV = apply(ures2[,c("Tumor","Strata")],1,function(xx){
			ifelse(xx[2] == "nHM",sprintf("%s-%s",xx[1],xx[2]),xx[1])
		})
		ures2
		
		tmp_df = name_TUMORTYPES()
		names(tmp_df) = c("Tumor","Name")
		ures2 = smart_merge(ures2,tmp_df)
		ures2
		
		print_latex_table(DATA = ures2[,c("Name","ABBV","CYT","OS","PFI","DSS")],
			add_table = TRUE,my_align = "llllll",fontsize = "tiny",repeat_VARS = "Name",
			caption = sprintf("Number of subjects and number of observed events"))
		
	}
	
	# Check on regression sample sizes
	sapply(unique(res$dataset),function(xx) sort(unique(res$Ne[res$dataset == xx])))
	
	if( TRUE ){ # Check top significant tumor types
		
		my_theme = theme(legend.position = "bottom",
			text = element_text(size = 40),
			plot.title = element_text(hjust = 0.5),
			panel.grid.major = element_line(colour = "grey50",
				size = 0.5,linetype = "dotted"),
			axis.title.x = element_text(vjust = -0.5),
			axis.text.x = element_text(size = 32),
			panel.background = element_blank(),
			panel.spacing.x = unit(0.75,"lines"),
			strip.text.y = element_text(angle = 0,hjust = 0))
		
		# In-depth euclidean vs optimal transport
		for(dataset in datasets){
			# dataset = c("BLCA","COAD","KIRC","LGG","LIHC","LUAD","PAAD")[2]
			cat(sprintf("%s: dataset = %s ...\n",date(),dataset))
			# if( dataset == "COAD" ) next
			res2 = res[which(res$dataset == dataset
				& res$OUT != "DFI" & res$GM <= 30 & res$GM > 1
				),]
			res2 = res2[order(res2$GM),]
			res2$OUT = ifelse(res2$OUT == "log10_CYT","Cytolytic~Activity",
				ifelse(res2$OUT == "OS","Overall~Survival",
				ifelse(res2$OUT == "PFI","Progression-Free~Interval",
					"Disease-Specific~Survival")))
			uOUT = sort(unique(res2$OUT)); uOUT
			mOUT = c("Cytolytic~Activity","Overall~Survival",
				"Progression-Free~Interval","Disease-Specific~Survival")
			mOUT = mOUT[mOUT %in% uOUT]
			res2$OUT = factor(res2$OUT,levels = mOUT)
			#res2$LAB = factor(res2$LAB,levels = c("EUC","OT_omnibus"),
			#	labels = c("Euclidean","Optimal Transport (Omnibus)"))
			# dim(res2); res2[1:5,]
			
			res2$LAM = sapply(res2$LAB,function(xx){
				if( grepl("EUC",xx) ){
					"EUC"
				} else {
					# xx = res2$LAB[2]; xx
					LAM = gsub("LAM.","",strsplit(xx,"_")[[1]][3])
					sprintf("OT (lambda = %s)",LAM)
				}
			},USE.NAMES = FALSE)
			res2$GS = sapply(res2$LAB,function(xx){
				if( grepl("EUC",xx) ){
					"EUC"
				} else if( grepl("GO",xx) ){
					"GO"
				} else if( grepl("PATH",xx) ){
					"PATH"
				} else if( grepl("ME",xx) ){
					"ME"
				} else {
					xx2 = strsplit(xx,"_")[[1]]
					xx2 = gsub(".rds","",xx2[length(xx2)])
					xx2 = gsub("ww","",xx2)
					ifelse(grepl("2",xx2),"GO",
						ifelse(grepl("3",xx2),"PATH","ME"))
				}
			},USE.NAMES = FALSE)
			
			max_yy = max(c(-log10(0.05),res2$YY)) + 0.5
			if( any(is.infinite(max_yy)) ){
				cat(sprintf("Skip %s b/c perm_pval = 0\n",dataset))
				next
			}
			
			lev_LAM = sort(unique(res2$LAM)); # lev_LAM
			res2$LAM = factor(res2$LAM,levels = lev_LAM[c(1,2,3,4,5)],
				labels = c("Euclidean",expression(paste("OT (",lambda," = 0.5)")),
					expression(paste("OT (",lambda," = 1.0)")),
					expression(paste("OT (",lambda," = 5.0)")),
					# expression(paste("OT (",lambda," = ",infinity,")"))
					"OT (Balanced)"))
			lev_GS = sort(unique(res2$GS)); # lev_GS
			res2$GS = factor(res2$GS,levels = lev_GS[c(1,2,4,3)],
				labels = c("Euclidean","GO","PATH","ME"))
			GM = YY = OUT = LAM = STRAT = NULL
			gg = ggplot(data = res2,mapping = aes(x = factor(GM),y = YY,fill = GS)) +
				geom_bar(stat = "identity",position = position_dodge(),width = 0.8) +
				geom_hline(yintercept = -log10(0.05),size = 1.5,linetype = 2) +
				#geom_hline(data = tmp_df,size = 1.2,linetype = 2,
				#	mapping = aes(yintercept = YY,color = factor(GM_0))) +
				xlab("Gene Mutation Frequency Cutoff") + 
				ggtitle(dataset) +
				ylab("-log10 Permutation P-value") + labs(fill = "Gene Similarity Method") + 
				ylim(c(0,max_yy)) + guides(color = "none") + my_theme
			# dataset = unique(res2$dataset)
			
			if( dataset == "COAD" ){
				gg = gg + facet_nested(OUT ~ LAM + STRAT,labeller = label_parsed)
				width = 45
			} else {
				# gg = gg + facet_grid(out2 ~ LAM,labeller = label_parsed)
				gg = gg + facet_grid(OUT ~ LAM,labeller = label_parsed)
				width = 30
			}
			# gg
			
			png_fn = file.path(my_dirs$fin_dir,sprintf("indepth_%s.png",dataset))
			ggsave(png_fn,plot = gg,device = "png",width = width,height = 20,units = "in")
			rm(gg)
			
		}
		
	}
	
	return(NULL)
	
}




#' @importFrom smartr smart_reqNames smart_table
#'	smart_merge name_change smart_df smart_rmcols
#'	smart_progress smart_names smart_mkdir smart_heatmap
#'	smart_colors print_latex_table
#' @importFrom ggplot2 aes element_blank element_line theme
#'	element_rect element_text facet_grid facet_wrap ylim unit
#'	geom_abline geom_bar geom_histogram geom_hline geom_point
#'	geom_segment geom_text ggplot ggsave ggtitle guide_legend
#'	guides label_parsed labs position_dodge scale_color_manual
#'	xlab ylab scale_shape_manual
#' @importFrom ggh4x facet_nested
#' @importFrom survival coxph
#' @importFrom stats dhyper fisher.test glm pnorm
#'	quantile as.formula dist drop1 formula lm
#'	model.matrix
#' @importFrom car Anova
#' @importFrom grDevices dev.off png
#' @importFrom utils str head tail
#' @importFrom readxl read_excel
#' @importFrom GOSemSim godata mgeneSim
#' @importFrom ROKET run_myOTs kernTEST
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

