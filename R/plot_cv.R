plot_cv <- function(res, cols=NULL, add=TRUE){
  library("Hmisc")
  xlab = expression(Log(lambda))
  ylab = "deviance"#res$name
  if(is.null(cols)){
    cols = 1:NCOL(res$cvm)
  }
  
  if( is.null(dim(res$cvm))){
    print("cols")
    print(cols)
    Hmisc::errbar(log(res$lambda), res$cvm, res$cvup, res$cvlo, type = "b",add = FALSE, 
                  xlab = xlab, ylab = ylab, col = cols[1], errbar.col = cols[1])
  } else {
    nmeth <- ncol(res$cvm)
    for(i in 1:nmeth ){
      add_to_plot <- (add & i > 1)
      if(add){
        ylim <- c(min(res$cvlo), max(res$cvup))
        xlim <- c(min(log(res$lambdas)), max(log(res$lambdas)))
        Hmisc::errbar(log(res$lambdas[,i]), res$cvm[,i], res$cvup[,i], res$cvlo[,i], type = "b", 
                      xlab = xlab, ylab = ylab, col = cols[i], errbar.col = cols[i], add = add_to_plot)
        legend('bottomright', col = cols, legend = res$alphas, lwd = 1)
      } else {
        Hmisc::errbar(log(res$lambdas[,i]), res$cvm[,i], res$cvup[,i], res$cvlo[,i], type = "b", xlab = xlab, 
                      ylab = ylab, col = cols[i], errbar.col = cols[i])
      }
    }
  }
}