

#' enumerate all TCGA tumor sites (anatomic site codes)
#' @examples
#' head(tcga_sites())
#' @return character vector
#' @export
tcga_sites = function() c("acc", "blca", "brca", "cesc", 
 "chol", "coad", "dlbc", "esca", "fppp", "gbm", "hnsc",
 "kich", "kirc", "kirp", "laml", "lgg", "lihc", "luad", 
 "lusc", "meso", "ov", "paad", "pcpg", "prad", "read", 
 "sarc", "skcm", "stad", "tgct", "thca", "thym", 
 "ucec", "ucs", "uvm")




#' map assay names used in GDC to shorter (improvised) names 
#' @examples
#' gdmap = gdc_assay_map()
#' summary(nchar(gdmap))
#' head(gdmap)
#' @return character vector
#' @export
gdc_assay_map = function() c(biospecimen_clin="clin__bio__nationwidechildrens_org__Level_1__biospecimen__clin",
clin="clin__bio__nationwidechildrens_org__Level_1__clinical__clin",
cna="cna__illuminahiseq_dnaseqc__hms_harvard_edu__Level_3__segmentation__seg",
meth27="methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data",
meth450="methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data",
mir_ga_gene="mirnaseq__illuminaga_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data",
mir_hiseq_gene="mirnaseq__illuminaga_mirnaseq__bcgsc_ca__Level_3__miR_isoform_expression__data",
mir_hiseq_gene="mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data",
mir_hiseq_iso="mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_isoform_expression__data",
rppa="protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data",
rnaseq_exon_expr="rnaseq__illuminahiseq_rnaseq__unc_edu__Level_3__exon_expression__data",
rnaseq_gene="rnaseq__illuminahiseq_rnaseq__unc_edu__Level_3__gene_expression__data",
rnaseq_spljunc="rnaseq__illuminahiseq_rnaseq__unc_edu__Level_3__splice_junction_expression__data",
rsem_genes="rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data",
rsem_genes_norm="rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data",
rsem_isoforms="rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms__data",
rsem_isoforms_norm="rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms_normalized__data",
rnaseq_exon_quant="rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__exon_quantification__data",
rnaseq_junc_quant="rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__junction_quantification__data",
snp_seg_hg18="snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg18__seg",
snp_seq_hg19="snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg19__seg",
snp_seg_hg18_minus_germline="snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg18__seg",
snp_seg_hg19_minus_germline="snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg",
tx_agilent="transcriptome__agilentg4502a_07_3__unc_edu__Level_3__unc_lowess_normalization_gene_level__data")

# private functions
assay_cands = gdc_assay_map() 
dcc_prefix = function() "gs://firecloud-tcga-open-access/tcga/dcc/"
assay_template = function(site, assay=assay_cands["rsem_genes_norm"]) paste0(dcc_prefix(), site, "/", assay)
assay_src = function(site, assay=assay_cands["rsem_genes_norm"]) {
    system(paste0("gsutil ls ", assay_template(site=site, assay=assay)), intern=TRUE)
}
allsites = tcga_sites()

#' retrieve GDC text files from AnVIL buckets
#' @param site character(1) anatomic site of tumor
#' @param brief_assay character(1) element of names(gdc_assay_map())
#' @param assay character(1) full name of assay in GDC nomenclature
#' @param suffix character(1) used to annotate folder name where text files are localized
#' @param localize_dry logical(1) set value of parameter `dry` for AnVIL::localize call; default FALSE implies that data will be retrieved
#' @param choice numeric(1) function fails if AnVIL has multiple versions
#' of assay quantifications in place, and will report in order the names
#' of versions present; set `choice` to the position in the
#' failure report of the version to be used
#' @return invisibly, the name of the folder where text files
#' were localized
#' @note When successful, this function populates a new folder
#' with text files drawn from GCP buckets managed in the AnVIL
#' project.  A tempdir() could be used, this may be default in future,
#' so this function may be redesigned in the near future.
#' @examples
#' if (interactive()) { # will egress
#' curd = getwd()
#' td = tempdir()
#' setwd(td)
#' nn = localize_tumor_site(site="acc", brief_assay="rppa",
#'  localize_dry = TRUE)
#' dir(nn)
#' unlink(nn, recursive=TRUE)
#' setwd(curd)
#' }
#' @export
localize_tumor_site = function(site, brief_assay = "rsem_genes_norm", 
                               assay = assay_cands[brief_assay], suffix="_terra_txt", localize_dry=FALSE, choice=NULL) {
    dir2use = paste0(site, "_", brief_assay, "_", suffix)
    if (dir.exists(dir2use)) stop("target folder already exists, please rename")
    dir.create(dir2use)
    src = assay_src(site=site, assay=assay)
    if (length(src)!=1) {
        message("multiple assay versions present")
        print(src)
        if (is.null(choice)) stop("multiple assay folders found and choice is NULL, please set choice")
        src = src[choice]
        }
    localize(src, dir2use, dry=localize_dry)
    invisible(dir2use)
}

#' Create an annotated matrix of quantifications for a specified TCGA tumor site and assay
#' @importFrom utils read.delim
#' @import AnVIL
#' @inheritParams localize_tumor_site
#' @return a matrix with rownames determined by column 1 of
#' the first discovered text file (skipping first line)
#' and colnames determined by the identifier present in the
#' first line of each text file
#' @note If `localize_dry` is TRUE, NULL is returned.
#' @examples
#' if (interactive()) { # will egress
#'  curd = getwd()
#'  td = tempdir()
#'  setwd(td)
#'  tst = build_tcga_mat("acc", "rppa", suffix="accrppatst")
#'  tst[1:4,1:4]
#'  fis = dir()
#'  nf = grep("accrppatst", fis, value=TRUE)
#'  if (length(nf)==1) unlink(nf, recursive=TRUE)
#'  setwd(curd)
#' }
#' @export
build_tcga_mat = function(site, brief_assay = "rsem_genes_norm", 
             assay = assay_cands[brief_assay], 
             localize_dry=FALSE,
             suffix="_terra_txt", choice=NULL) {
 newdir = localize_tumor_site(site=site, brief_assay=brief_assay, assay=assay, localize_dry=localize_dry, suffix=suffix, choice=choice)
 if (localize_dry) {
   message("localize_dry is TRUE, returning NULL")
   return(NULL)
   }
 allfi = dir(newdir,full.names=TRUE)
 one_subj = read.delim(allfi[1], skip=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)
 #head(one_subj)
 dat = matrix(NA, nrow=nrow(one_subj), ncol=length(allfi))
 rownames(dat) = one_subj[,1] # consider colClasses for next
 for (i in 1:length(allfi)) dat[,i] = read.delim(allfi[i], skip=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)[,2]
 grabid = function(x) strsplit(readLines(x, n=1), "\t")[[1]][2]
 #grabid(allfi[1])
 cn = vapply(allfi, grabid, character(1))
 colnames(dat) = cn
 dat
}
