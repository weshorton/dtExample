#' tmpClassify
#' @description
#' Group data_dt into 3 levels each for up/down expression:
#' low sig: 0.1 > p >= 0.05
#' medium sig: 0.05 > p >= 0.01
#' high sig: 0.01 > p
#' @param data_dt proteomics data.table with differential expression results
#' @param col1_v column in data_dt indicating the fold change results (logFC for TMT; Z-score for Silac)
#' @param col2_v column in data_dt indicating the significance value (FDR for TMT; BH_adjusted for Silac)
#' @param lfc_v value for cut-off for fold change. Will be included in plot as vertical dotted line
#' @param pval_v vector of length 2! first value indicates the threshold for 'high' significance, second for 'low'
#' @param newName_v name of new column
#' @import data.table
#' @details
#' Add a new column in data_dt summarizing the direction and level of each entry.
#' @return data_dt with one new column
#' @export

tmpClassify <- function(data_dt, col1_v = "logFC", col2_v = "FDR", lfc_v = 0.5, pval_v = c(0.01, 0.05), newName_v = "diffExp") {

  
  print(truelength(data_dt))
  data_dt[get(col1_v) > lfc_v & get(col2_v) >= pval_v[2] & get(col2_v) < 0.1, (newName_v) := "upLow"]
  data_dt[get(col1_v) < -lfc_v & get(col2_v) >= pval_v[2] & get(col2_v) < 0.1, (newName_v) := "downLow"]
  
  data_dt[get(col1_v) > lfc_v & get(col2_v) >= pval_v[1] & get(col2_v) < pval_v[2], (newName_v) := "upMed"]
  data_dt[get(col1_v) < -lfc_v & get(col2_v) >= pval_v[1] & get(col2_v) < pval_v[2], (newName_v) := "downMed"]
  
  data_dt[get(col1_v) > lfc_v & get(col2_v) >= pval_v[1] & get(col2_v) < pval_v[2], (newName_v) := "upHigh"]
  data_dt[get(col1_v) < -lfc_v & get(col2_v) >= pval_v[1] & get(col2_v) < pval_v[2], (newName_v) := "downHigh"]
  
  data_dt[is.na(get(newName_v)), (newName_v) := "NO"]
  return(data_dt)
}
