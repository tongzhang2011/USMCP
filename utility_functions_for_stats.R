
library(dplyr)
# region1 <- "Cluster_1"
# region2 <- "Cluster_2"
# vplot(d1, "Cluster_1", "Cluster_2")


vplot <- function(d1,region1, region2) {
  ROI1 <- region1
  ROI2 <- region2
  
  my_title <- paste(ROI2, ROI1, sep = " / ")
  
  df1 <-  dplyr::select(d1, ends_with(ROI1))
  df2 <-  dplyr::select(d1, ends_with(ROI2))
  
  df <- cbind(df1, df2)
  
  df$p_value <- t.test_v2_df(df1, df2)
  df$FC <- FC(df1, df2)
  

  
  valid_p <- as.numeric(table(is.na(df$p_value))[1])
  my_title <- paste(my_title,  
                    paste0("(n = ", valid_p, ")"))
  
  p_cut <- 0.01
  FC_cut <- 1
  
  #  p_cut <- 0.01
  # FC_cut <- 0.5
  
  
  
  d_up <- filter(df, p_value < p_cut, FC > FC_cut)
  d_dn <- filter(df, p_value < p_cut, FC < -FC_cut)
  
  Increase <- dim(d_up)[1]
  Decrease <- dim(d_dn)[1]
  
  point_size <- 1
  
  # create a plot
  with(df,plot(FC, -log(p_value,10),pch = 20,
               xlab = "",
               ylab = "",
               cex.axis = 1, font.axis = 1,
               cex.lab = 0.8, font.lab = 1,
               cex = point_size))
  
  xmin <- par("usr")[1]
  xmax <- par("usr")[2]
  ymin <- par("usr")[3]
  ymax <- par("usr")[4]
  
  xrange <- xmax - xmin
  yrange <- ymax - ymin
  xlab_pos_1 <- xmin + 0.2*xrange
  xlab_pos_2 <- xmin + 0.8*xrange
  ylab_pos <- ymin + 0.9*yrange
  
  
  with(subset(df,p_value > p_cut),points(FC, -log(p_value,10),pch = 20, col = "gray",
                                         cex = point_size))
  with(subset(df,p_value < p_cut  & FC > FC_cut ),points(FC, -log(p_value,10),pch = 20, col = "red",
                                                         cex = point_size))
  with(subset(df,p_value < p_cut  & FC < -FC_cut ),points(FC, -log(p_value,10),pch = 20, col = "blue",
                                                          cex = point_size))
  abline(h = -log(p_cut,10),lty =2 ,col = "black")
  abline(v = FC_cut, lty =2 ,col = "red")
  abline(v = -FC_cut, lty =2 ,col = "blue")
  box(lwd = 2)
  text(xlab_pos_1, ylab_pos, paste("Decrease: ", Decrease, sep = ""),
       cex = 1.2, col = "blue")
  text(xlab_pos_2, ylab_pos, paste("Increase: ", Increase, sep = ""),
       cex = 1.2, col = "red")
  title(main = my_title)
  
}




p_value_plot <- function(d1,region1, region2) {
  ROI1 <- region1
  ROI2 <- region2
  
  my_title <- paste(ROI2, ROI1, sep = " / ")
  
  df1 <-  dplyr::select(d1, ends_with(ROI1))
  df2 <-  dplyr::select(d1, ends_with(ROI2))
  
  df <- cbind(df1, df2)
  
  df$p_value <- t.test_v2_df(df1, df2)
  # df$FC <- FC(df1, df2)
  
  
  valid_p <- as.numeric(table(is.na(df$p_value))[1])
  my_title <- paste(my_title,  
                    paste0("(n = ", valid_p, ")"))
  
  hist(df$p_value, 
       breaks = seq(0, 1, length.out = 100),
       col = c("blue", rep("gray", 99)),
       xlab = "",
       ylab = "",
       main = my_title)
  
}






sig_list <- function(d1, region1, region2) {
  ROI1 <- region1
  ROI2 <- region2
  
  my_title <- paste(ROI2, ROI1, sep = "_vs_")
  
  df1 <-  dplyr::select(d1, ends_with(ROI1))
  df2 <-  dplyr::select(d1, ends_with(ROI2))
  
  df <- cbind(df1, df2)
  
  df$p_value <- t.test_v2_df(df1, df2)
  df$FC <- FC(df1, df2)
  
  
  valid_p <- as.numeric(table(is.na(df$p_value))[1])
  my_title <- paste(my_title,  
                    paste0("(n = ", valid_p, ")"))
  
  p_cut <- 0.01
  FC_cut <- 1
  
  df_sig <- filter(df,p_value < p_cut)
  df_sig <- filter(df_sig,abs(FC) > FC_cut)
  
  write.csv(df_sig, paste(my_title, "sig_proteins_FC_1_Mapping.csv", sep = "_"))
  
}

####################################################################################calcluate p and FC
# step 1: define t.test and FC for this case
###########################################
#a t.test function that requires at least 2 valid numbers from both groups
at_least_number <- 3

t.test_at_least_2_v2 <- function(num1, num2) {
  if( sum(!is.na(num1)) >= at_least_number &  sum(!is.na(num2)) >= at_least_number ) {
    my_p.value <- t.test(num1,num2,na.rm=TRUE,var.equal=TRUE)$p.value
  } else {  my_p.value <- NA
  }
  return( my_p.value)
}

#apply this funtion to two dataframe.
#df1 contains the first group;
#df2 contains the second group;
t.test_v2_df <- function(df1,df2) {
  n <- dim(df1)[1]
  my_p_value <- vector(length = n)
  
  for (i in 1:n) {
    my_p_value[i] <- t.test_at_least_2_v2(as.numeric(df1[i,]), as.numeric(df2[i,]))
  }
  return(my_p_value)
}
#end of the function
##########################################

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#similarly, FC also requires at least two 
#valid values from each comparing group
fc_at_least_2  <- function(num1, num2) {
  if( sum(!is.na(num1)) >= at_least_number &  sum(!is.na(num2)) >= at_least_number ) {
    fc <-  mean(num2, na.rm=TRUE) - mean(num1, na.rm=TRUE)
  } else { fc <- NA}
  return( fc)
}

#define a function use two dataframes as input
#i.e: do the calculation at the datafrmae level
FC <- function(df1, df2) {
  n <- dim(df1)[1]
  my_fc <- vector(length = n)
  
  for (i in 1:n) {
    my_fc[i] <- fc_at_least_2(as.numeric(df1[i,]), as.numeric(df2[i,]))
  }
  return(my_fc)
}
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

vplot_one_vs_all_other <- function(d1,region2) {

  ROI2 <- region2
  
  my_title <- paste(ROI2, "Others", sep = " / ")
  
  df1 <-  dplyr::select(d1, !ends_with(ROI2))  # df1: all others 
  df2 <-  dplyr::select(d1, ends_with(ROI2))
  
  df <- cbind(df1, df2)
  
  df$p_value <- t.test_v2_df(df1, df2)
  df$FC <- FC(df1, df2)
  
  
  
  valid_p <- as.numeric(table(is.na(df$p_value))[1])
  my_title <- paste(my_title,  
                    paste0("(n = ", valid_p, ")"))
  
  p_cut <- 0.01
  FC_cut <- 0.5
  
  #  p_cut <- 0.01
  # FC_cut <- 0.5
  
  
  
  d_up <- filter(df, p_value < p_cut, FC > FC_cut)
  d_dn <- filter(df, p_value < p_cut, FC < -FC_cut)
  
  Increase <- dim(d_up)[1]
  Decrease <- dim(d_dn)[1]
  
  point_size <- 1
  
  # create a plot
  with(df,plot(FC, -log(p_value,10),pch = 20,
               xlab = "",
               ylab = "",
               cex.axis = 1, font.axis = 1,
               cex.lab = 0.8, font.lab = 1,
               cex = point_size))
  
  xmin <- par("usr")[1]
  xmax <- par("usr")[2]
  ymin <- par("usr")[3]
  ymax <- par("usr")[4]
  
  xrange <- xmax - xmin
  yrange <- ymax - ymin
  xlab_pos_1 <- xmin + 0.2*xrange
  xlab_pos_2 <- xmin + 0.8*xrange
  ylab_pos <- ymin + 0.9*yrange
  
  
  with(subset(df,p_value > p_cut),points(FC, -log(p_value,10),pch = 20, col = "gray",
                                         cex = point_size))
  with(subset(df,p_value < p_cut  & FC > FC_cut ),points(FC, -log(p_value,10),pch = 20, col = "red",
                                                         cex = point_size))
  with(subset(df,p_value < p_cut  & FC < -FC_cut ),points(FC, -log(p_value,10),pch = 20, col = "blue",
                                                          cex = point_size))
  abline(h = -log(p_cut,10),lty =2 ,col = "black")
  abline(v = FC_cut, lty =2 ,col = "red")
  abline(v = -FC_cut, lty =2 ,col = "blue")
  box(lwd = 2)
  text(xlab_pos_1, ylab_pos, paste("Decrease: ", Decrease, sep = ""),
       cex = 1.2, col = "blue")
  text(xlab_pos_2, ylab_pos, paste("Increase: ", Increase, sep = ""),
       cex = 1.2, col = "red")
  title(main = my_title)
  
}

p_value_plot_one_vs_all_others <- function(d1, region2) {
  ROI2 <- region2
  
  my_title <- paste(ROI2, "Others", sep = " / ")
  
  df1 <-  dplyr::select(d1, !ends_with(ROI2))  # df1: all others 
  df2 <-  dplyr::select(d1, ends_with(ROI2))
  
  df <- cbind(df1, df2)
  
  df$p_value <- t.test_v2_df(df1, df2)

  # df$FC <- FC(df1, df2)
  
  
  valid_p <- as.numeric(table(is.na(df$p_value))[1])
  my_title <- paste(my_title,  
                    paste0("(n = ", valid_p, ")"))
  
  hist(df$p_value, 
       breaks = seq(0, 1, length.out = 100),
       col = c("blue", rep("gray", 99)),
       xlab = "",
       ylab = "",
       main = my_title)
  
}




sig_list_one_vs_all_others <- function(d1, region2) {
  ROI2 <- region2
  
  my_title <- paste(ROI2, "Others", sep = "_vs_")
  
  df1 <-  dplyr::select(d1, !ends_with(ROI2))  # df1: all others 
  df2 <-  dplyr::select(d1, ends_with(ROI2))
  
  df <- cbind(df1, df2)
  
  df$p_value <- t.test_v2_df(df1, df2)
  df$FC <- FC(df1, df2)
  
  
  valid_p <- as.numeric(table(is.na(df$p_value))[1])

  
  p_cut <- 0.01
  FC_cut <- 0.5
  
  df_sig <- filter(df,p_value < p_cut)
  df_sig <- filter(df_sig,abs(FC) > FC_cut)
  
  write.csv(df_sig, paste(my_title, "sig_proteins_FC_1.csv", sep = "_"))
  
}




sig_list_one_vs_all_others_only_up_proteins <- function(d1, region2) {
  ROI2 <- region2
  
  my_title <- paste(ROI2, "Others", sep = "_vs_")
  
  df1 <-  dplyr::select(d1, !ends_with(ROI2))  # df1: all others 
  df2 <-  dplyr::select(d1, ends_with(ROI2))
  
  df <- cbind(df1, df2)
  
  df$p_value <- t.test_v2_df(df1, df2)
  df$FC <- FC(df1, df2)
  
  
  valid_p <- as.numeric(table(is.na(df$p_value))[1])
  
  
  p_cut <- 0.01
  FC_cut <- 1
  
  df_sig <- filter(df,p_value < p_cut)
 # df_sig <- filter(df_sig,abs(FC) > FC_cut)
  
  # only up-regulated proteins
  df_sig <- filter(df_sig,FC > FC_cut)
  
  write.csv(df_sig, paste(my_title, "sig_proteins_FC_1_up_only.csv", sep = "_"))
  
}












