#' @title Convert transcript IDs between different databases
#' @description Preprocessing, mapping/converting, searching SNP for data
#' @name convert_transcriptID
#'
#' @param dat a dataframe including Transcript_version, Gene, exon, Nucleotide, AA_changes
#' @param db `EnsDb` object.
#' @return a new dataset with converting information
#' @export
#' @examples
#' library(EnsDb.Hsapiens.v86)
#' db=EnsDb.Hsapiens.v86
#' dat<-read.csv(system.file("extdata",
#'                           "variant_list_test.csv",
#'                           package = "CvtransctiptID"),
#'               stringsAsFactors = FALSE, encoding = "UTF-8", row.names = NULL, sep = ",")
#' new_dat<-convert_transcriptID(dat, db)
#'
convert_transcriptID <- function(dat, db){

  #preprocess data
  str<-as.array(dat$Nucleotide_changes)
  str_split<-as.integer(unlist(stringr::str_extract_all(str, "[0-9]+")))
  dat$start_loc<-str_split
  dat$ref<-substr(dat$Nucleotide_changes, 1, 1)
  dat$alt<-stringr::str_sub(dat$Nucleotide_changes, -1, -1)
  dat<-dat[order(dat$Transcript_version, dat$start_loc),]

  #Mapping protein coordinates to transcript coordinates
  ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  dat_ensmbl_id<-biomaRt::getBM(attributes=c("refseq_mrna", "ensembl_transcript_id", "hgnc_symbol"), filters = "refseq_mrna",
                       values = dat$Transcript_version, mart= ensembl, uniqueRows = FALSE)
  dat_ensmbl_id<-dat_ensmbl_id[order(dat_ensmbl_id$refseq_mrna),]
  newvalues <- factor(dat_ensmbl_id$ensembl_transcript_id)
  dat$em_id <- newvalues[ match(dat$Transcript_version, dat_ensmbl_id$refseq_mrna) ]

  #search SNP
  loc<-dat$start_loc
  tx_name<-dat$em_id
  cds <- IRanges::IRanges(start = loc, width = 1, name = tx_name)
  cds_tx <- ensembldb::cdsToTranscript(cds, db)
  cds_tx_gn<-ensembldb::transcriptToGenome(cds_tx, db)
  cds_tx_gn_df<-as.data.frame(cds_tx_gn)

  #subset dat with only detectable transcripts
  dat2<-dat[dat$em_id %in% cds_tx_gn_df$group_name == T, ]

  #return final result
  results<-cbind(dat2, cds_tx_gn_df)
  results
}
