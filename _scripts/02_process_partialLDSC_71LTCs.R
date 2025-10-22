#############################################
####                                     ####
####      GEMINI, process results        ####
####      partial rg - BMI               ####
####                                     ####
####                                     ####
####      Ninon Mounier                  ####
####      02/08/2023                     ####
####                                     ####
#############################################


# assumed you are: 
#  1. running in the repo directory
#  2. have downloaded the "GEMINI GWAS" to directory GWASs_GEMINI/ 
#     get from Zenodo https://doi.org/10.5281/zenodo.14284046
#  3. downloaded the BMI summary statstics to a provided directory
#  4. run the partialLDSC analysis


library(tidyverse)
library(openxlsx)
library(MetBrewer)
library(plotly)
library(igraph)


#### xlsx tables ####

# reads in _data/[name]_difference.tsv
# writes nice xlxs: file _results/_tables/[name].xlsx

make_table <- function(name){
  
  all_res = read_tsv(paste0("_data/", name, "_difference.tsv"))
  
  # restrict our analyses to the non-obesity pairs (but keep the ones with non-significant rg)
  # add rg p-value and fdr
  all_res %>%
    mutate(rg.P = 2*pnorm(-abs(rg/rg.SE)),
           rg.FDR = p.adjust(rg.P, method = "fdr"), .before = partial_rg) -> all_res
  # table(all_res$rg.FDR<0.05) 
  # FALSE  TRUE 
  #   953  1532 
  
  all_res %>%
    mutate(diff.FDR = p.adjust(diff.P, method = "fdr")) -> all_res
  # nrow(all_res) 
  # 2485
  
  
  table(all_res$diff.FDR<0.05)
  # FALSE  TRUE (200 blocks, 71 conditions + BMI)
  #   1123  1362 
 
  
  # add partial rg p-value and fdr (for all pairs)
  all_res %>% 
    mutate(partial_rg.P = 2*pnorm(-abs(partial_rg/partial_rg.SE)), 
           partial_rg.FDR = p.adjust(partial_rg.P, method = "fdr"), .before=rg_cov) -> res_rg 
  
  # pairs with non-significant rg but significant difference
  res_rg %>% filter(rg.FDR > 0.05, diff.FDR<0.05) %>% nrow # 341
  # pairs with significant rg but significant difference
  res_rg %>% filter(rg.FDR < 0.05, diff.FDR<0.05) %>% nrow # 1021
  
  # pairs with non-significant rg, significant difference, and significant partial rg
  res_rg %>% filter(rg.FDR > 0.05, diff.FDR<0.05, partial_rg.FDR<0.05) %>% nrow # 33
  
  # pairs with significant rg, significant difference, and stronger partial rg
  res_rg %>% filter(rg.FDR < 0.05, diff.FDR < 0.05, rg*diff.T<0) %>% nrow # 120
  
  
  # pairs with significant rg, significant difference, and non-significant partial rg
  res_rg %>% filter(rg.FDR < 0.05, diff.FDR < 0.05, partial_rg.FDR > 0.05) %>% nrow # 161
  
  
  
  
  #### xlsx file with results ####
  
  
  ## create workbook
  wb <- createWorkbook() 
  
  ## ceate worksheet
  addWorksheet(wb, "all results", gridLines = TRUE)
  
  
  ## add results
  writeDataTable(wb, 1, res_rg %>% arrange(diff.P), tableStyle = "TableStyleLight8") 
  
  
  # column widths
  setColWidths(wb, 1, cols = 1:2, widths = 25)
  setColWidths(wb, 1, cols = 1:14, widths = 13)
  setColWidths(wb, 1, cols = 11, hidden = TRUE)
  
  ## number of digit for decimal numbers
  dec2 <- createStyle(numFmt = "0.00")
  dec4 <- createStyle(numFmt = "0.0000")
  dec5 <- createStyle(numFmt = "0.00000")
  sci <- createStyle(numFmt = "SCIENTIFIC")
  
  
  addStyle(wb, 1, style = dec4, cols = c(3,4), rows = 2:(nrow(res_rg)+1), gridExpand = TRUE)
  addStyle(wb, 1, style = sci, cols = c(5,6), rows = 2:(nrow(res_rg)+1), gridExpand = TRUE)
  addStyle(wb, 1, style = dec4, cols = c(7,8), rows = 2:(nrow(res_rg)+1), gridExpand = TRUE)
  addStyle(wb, 1, style = sci, cols = c(9,10), rows = 2:(nrow(res_rg)+1), gridExpand = TRUE)
  addStyle(wb, 1, style = dec5, cols = c(11), rows = 2:(nrow(res_rg)+1), gridExpand = TRUE)
  addStyle(wb, 1, style = dec2, cols = c(12), rows = 2:(nrow(res_rg)+1), gridExpand = TRUE)
  addStyle(wb, 1, style = sci, cols = c(13,14), rows = 2:(nrow(res_rg)+1), gridExpand = TRUE)
  
  
  # color rules
  n_signif = sum(res_rg %>% pull(diff.FDR)<0.05)
  # diff.FDR < 0.05 -> column 14 green
  signif <- createStyle(bgFill = "#a5d4ac")
  
  conditionalFormatting(wb, 1, cols = 14, rows = 2:(nrow(res_rg)+1), rule = "$N2<0.05",
                        style = signif)
  
  # diff.FDR < 0.05 & rg.FDR> 0.05 & partial_rg.FDR < 0.05  -> column 7:10 grey
  fully_explained <- createStyle(bgFill = "#dedede")
  
  conditionalFormatting(wb, 1, cols = 7:10, rows = 2:(nrow(res_rg)+1), rule = "$J2>0.05",
                        style = fully_explained)
  
  # diff.FDR < 0.05 & partial_rg stronger  -> column 12 bold
  stronger_partial <- createStyle(bgFill = "#ece0ff", textDecoration = "Bold")
  
  conditionalFormatting(wb, 1, cols = 12, rows = 2:(nrow(res_rg)+1), rule = "$C2*$L2<0",
                        style = stronger_partial)
  
  # diff.FDR < 0.05 & rg not significant  -> column 3:6 grey
  new_signal <- createStyle(bgFill = "#dedede")
  conditionalFormatting(wb, 1, cols = 3:6, rows = 2:(nrow(res_rg)+1), rule = "$F2>0.05",
                        style = new_signal)
  
  
  
  
  ## borders
  bord_right <- createStyle(border = "right")
  bottom = nrow(res_rg)+1
  addStyle(wb, 1, style = bord_right, cols = 2, rows = 1:bottom, gridExpand = FALSE, stack=TRUE)
  addStyle(wb, 1, style = bord_right, cols = 6, rows = 1:bottom, gridExpand = FALSE, stack=TRUE)
  addStyle(wb, 1, style = bord_right, cols = 10, rows = 1:bottom, gridExpand = FALSE, stack=TRUE)
  addStyle(wb, 1, style = bord_right, cols = 14, rows = 1:bottom, gridExpand = FALSE, stack=TRUE)
  
  bord_bottom <- createStyle(border = "bottom")
  addStyle(wb, 1, style = bord_bottom, cols = 1:ncol(res_rg), rows = bottom, gridExpand = FALSE, stack=TRUE)
  
  
  freezePane(wb, 1,
             firstRow = TRUE,
             firstActiveCol = 3)
  
  
  
  # add a README sheet
  readme = c(paste0("This table contains the results for ",length(unique(c(res_rg$condition.1, res_rg$condition.2))), " conditions (", nrow(res_rg), " pairs)."),
             "Raw and partial correlations (adjusted for BMI genetics) are compared for all pairs.",
             " ",
             "A:B : pair of conditions",
             "C:F : raw genetic correlations (if greyed out, fdr < 0.05)",
             "G:J : partial genetic correlations (if greyed out, fdr < 0.05)",
             "L:N : testing for difference between the two:",
             paste0("if green in column N, significant difference (fdr < 0.05 - ", n_signif, " pairs)"),
             "if purple in column L, partial correlation stronger than raw correlation"
  )
  
  addWorksheet(wb, "README")
  writeData(wb, 2, readme) 
  
  
  
  
  ## save workbook
  saveWorkbook(wb, paste0("_results/", name, ".xlsx"), overwrite = TRUE)
}

make_table("71LTCs_070823_200b")  # 71 conditions - 200 blocks (MAIN ONE USED FOR EVERYTHING ELSE)

#make_table("71LTCs_070823_500b")  # 71 conditions - 500 blocks
#make_table("71LTCs_070823_1000b") # 71 conditions - 1000 blocks
#make_table("71LTCs_070823_1500b") # 71 conditions - 1500 blocks
#make_table("71LTCs_070823_2000b") # 71 conditions - 2000 blocks
#make_table("71LTCs_070823_2500b") # 71 conditions - 2500 blocks
#make_table("71LTCs_070823_3000b") # 71 conditions - 3000 blocks
#make_table("71LTCs_070823_3500b") # 71 conditions - 3500 blocks
#make_table("71LTCs_070823_4000b") # 71 conditions - 4000 blocks
#make_table("71LTCs_070823_4500b") # 71 conditions - 4500 blocks

make_table("71LTCs_070823_200b_2015") # 71 conditions - 200 blocks (2015 BMI data)

make_table("71LTCs_070823_200b_obesity") # 71 conditions - 200 blocks (GEMINI obesity data)


#### interactive scatter plots ####

## plot (using plotly)
# reads in file _results/_tables/[name]_BMI.xlsx
# rg vs partial_rg
# color depending on which ones are significant
# shape if significant difference upper triangle is stronger, lower triangle if weaker
# saves figure as  _results/_figures/[name]_BMI.html
# saves figure as  _results/_figures/[name]_BMI_signif.html

# pairs = all (plot all pairs)
# pairs = signif (plot significant pairs)
# pairs = both (both plots)
make_scatterplot <- function(name, pairs = "both"){
  
  if(!pairs %in% c("all", "signif", "both")) stop("pairs should be \"all\", \"signif\" or \"both\".")
  
  res_file = paste0("_results/", name, ".xlsx")
  if(!file.exists(res_file)) stop("results file should first be created using make_table().")
  res_rg = read.xlsx(res_file)
  
  res_rg %>% 
    # color
    mutate(significance = case_when(
      rg.FDR > 0.05 & partial_rg.FDR > 0.05 ~ "no correlation significant",
      rg.FDR < 0.05 & partial_rg.FDR > 0.05 ~ "only raw correlation significant",
      rg.FDR > 0.05 & partial_rg.FDR < 0.05 ~ "only partial correlation significant",
      rg.FDR < 0.05 & partial_rg.FDR < 0.05 ~ "both significant"),
      # shape
      difference = case_when(
        diff.FDR > 0.05 ~ "no significant difference",
        diff.FDR < 0.05 &  rg*diff.T<0 ~ "partial correlation stronger",
        diff.FDR < 0.05 &  rg*diff.T>0 ~ "partial correlation weaker"),
      pair = paste0(condition.1, "~", condition.2)) -> res_plot
  
  
  theme_set(theme_bw())
  my_colors = met.brewer("Egypt", 4)
  
  if(pairs %in% c("all", "both")){
    res_plot %>%
      ggplot(aes(x=rg, y=partial_rg, color=significance, shape=difference, text=pair)) + 
      geom_abline(slope=1, intercept = 0) +
      geom_point() +
      scale_shape_manual(values=c(1, 2, 6))+ 
      scale_color_manual(values=my_colors) +
      labs(x="raw genetic correlation", y="partial genetic correlation\n(adjusted for BMI)",
           color="", shape="")-> p
    
    
    ggplotly(p) -> myplot
    # to remove () from the legend
    for (i in 1:length(myplot$x$data)){
      if (!is.null(myplot$x$data[[i]]$name)){
        myplot$x$data[[i]]$name = stringr::str_remove(myplot$x$data[[i]]$name, "\\(")
        myplot$x$data[[i]]$name = stringr::str_remove(myplot$x$data[[i]]$name, "\\)")
        myplot$x$data[[i]]$name = stringr::str_replace(myplot$x$data[[i]]$name, ",", ",\\\n")
        
        
      }
    }
    
    htmlwidgets::saveWidget(
      widget = myplot, #the plotly object
      file = paste0("_results/", name, ".html"), #the path & file name
      selfcontained = TRUE #creates a single html file
    )
  }


  # same figure, only pairs with significant difference
  if(pairs %in% c("signif", "both")){
    res_plot %>%
      filter(diff.FDR<0.05) %>%
      ggplot(aes(x=rg, y=partial_rg, color=significance, shape=difference, text=pair)) + 
      geom_abline(slope=1, intercept = 0) +
      geom_point() +
      scale_shape_manual(values=c(2, 6))+ 
      scale_color_manual(values=my_colors) +
      labs(x="raw genetic correlation", y="partial genetic correlation\n(adjusted for BMI)",
           color="", shape="")-> p_signif
    
    
    ggplotly(p_signif) -> myplot
    for (i in 1:length(myplot$x$data)){
      if (!is.null(myplot$x$data[[i]]$name)){
        myplot$x$data[[i]]$name = stringr::str_remove(myplot$x$data[[i]]$name, "\\(")
        myplot$x$data[[i]]$name = stringr::str_remove(myplot$x$data[[i]]$name, "\\)")
        myplot$x$data[[i]]$name = stringr::str_replace(myplot$x$data[[i]]$name, ",", ",\\\n")
        
        
      }
    }
    
    htmlwidgets::saveWidget(
      widget = myplot, #the plotly object
      file = paste0("_results/", name, "_signif.html"), #the path & file name
      selfcontained = TRUE #creates a single html file
    )
  }
}

# only do it for 200 blocks
make_scatterplot("71LTCs_070823_200b")
make_scatterplot("71LTCs_070823_200b_2015")

