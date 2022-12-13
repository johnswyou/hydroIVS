#' Perform Input Variable Selection
#'
#' This function performs input variable selection for a given input-output
#' dataset.
#'
#' @param y target vector \[N x 1\]
#' @param x input matrix \[N x D\]
#' @param ivsm input variable selection (IVS) method \[string:'none','boruta','rrf','ea_cmi_htc','ea_cmi_tol','knn_cmi_tol','knn_cmi_bi_tol',pmis_bic','pcis_bic'\]
#' @param ivs_param parameters for input variable selection method
#' \itemize{
#' \item rrf
#' \itemize{
#' \item ivs_param\[1\] = regularization parameter \[0 < scalar <= 1\]
#' \item ivs_param\[2\] = ntrees \[integer > 1)\] usually between 128 and 1000
#' }
#' \item ea_cmi
#' \itemize{
#' \item ivs_param\[1\] =  Hampel test criterion threshold \[scalar > 0\]
#' }
#' \item ea_cmi_tol
#' \itemize{
#' \item ivs_param\[1\] = CMI/MI threshold \[0 < scalar <= 1\]
#' \item ivs_param\[2\] = no. of nearest neighbours \[1 < integer < sample size-1\]
#' }
#' \item knn_cmi_tol
#' \itemize{
#' \item ivs_param\[1\] = CMI/MI threshold \[0 < scalar <= 1\]
#' \item ivs_param\[2\] = no. of nearest neighbours \[1 < integer < sample size-1\]
#' }
#' }
#'
#' @return list that contains indices and names of selected input variables
#'
#' @references
#'
#' Quilty, J., J. Adamowski, B. Khalil, and M. Rathinasamy (2016), Bootstrap rank-
#' ordered conditional mutual information (broCMI): A nonlinear input variable
#' selection method for water resources modeling, Water Resour. Res., 52,
#' doi:10.1002/2015WR016959.
#'
#' Quilty, J., and J. Adamowski (2018), Addressing the incorrect usage of wavelet-
#' based hydrological and water resources forecasting models for real-world
#' applications with best practices and a new forecasting framework, J. Hydrol.,
#' doi:10.1016/j.jhydrol.2018.05.003.
#'
#' @author John Quilty
#' @export

#  Date Created:
#
#  Sep. 10, 2018
#
#  Date(s) Modified:
#
#  Sep. 26, 2018 - Included 'ntrees' in 'rrf' as a specificable parameter
#  Oct. 29, 2018 - Included the 'score' (importance) of each selected feature

ivsIOData <- function(y,x,ivsm,ivs_param){

  names_inputs = colnames(x)
  inputs = as.matrix(x)
  n_inputs = ncol(inputs) # number of inputs
  tmp_ind = seq(1,n_inputs) # temporary sequence

  switch(ivsm,

         none={ # NO ivs

           sel_inputs = tmp_ind
           scores = rep(0, length(sel_inputs)) # same for all variables

         },

         boruta={

           ivs = Boruta::Boruta(x=inputs,y=y) # Boruta IVS
           sel_inputs = tmp_ind[Boruta::attStats(ivs)$decision=="Confirmed"] # selected inputs
           scores = colMeans(ivs[["ImpHistory"]][,sel_inputs]) # mean importance across different Boruta runs

         },

         rrf={

           ivs = RRF::RRF(x=inputs,y=y,
                          ntree=ivs_param[1], #500, # usually set between 128 and 1000, 500 is typical
                          mtry = max(floor(n_inputs/3), 1), # usually set to max(floor(ncol(x)/3), 1)
                          flagReg=1, # set to 0 for standard RF (we want RRF for IVS though...)
                          coefReg=ivs_param[2]) # set (0,1], smaller values force sparsity
           sel_inputs = ivs$feaSet # selected inputs
           scores = (ivs$importance / max(ivs$importance))[sel_inputs] # normalized importance measure

         },

         ea_cmi_htc={

           ivs = EAcmi_framework_htc(x=inputs,y=y,ivs_param,silent=TRUE)
           sel_inputs = ivs$Input # selected inputs
           scores = ivs$CMI # CMI during forward selection

         },

         ea_cmi_tol={

           ivs = EAcmi_framework_tol(x=inputs,y=y,ivs_param,silent=TRUE)
           sel_inputs = ivs$Input # selected inputs
           scores = ivs$CMI # CMI during forward selection

         },

         knn_cmi_tol={

           ivs = KNNcmi_framework_tol(x=inputs,y=y,
                                      ivs_param[1],k=ivs_param[2],
                                      silent=TRUE)
           sel_inputs = ivs$Input # selected inputs
           scores = ivs$CMI # CMI during forward selection

         },

         knn_cmi_bi_tol={

           ivs = KNNcmi_bi_framework_tol(x=inputs,y=y,
                                         thresh=ivs_param[1],
                                         k=ivs_param[2],silent=TRUE)
           sel_inputs = ivs$Input # selected inputs
           scores = ivs$CMI # CMI during forward selection

         },

         pmis_bic={

           ivs = PMI_grnn_framework_bic(x=inputs,y=y,
                                        ybin=FALSE, silent=TRUE)
           sel_inputs = ivs$Input # selected inputs
           scores = ivs$PMI # PMI during forward selection

         },

         pcis_bic={

           ivs = PCIS_framework_bic(x=inputs,y=y,
                                    silent=TRUE)
           sel_inputs = ivs$Input # selected inputs
           scores = ivs$PC # Partial Correlation during forward selection

         }

  )

  names_sel_inputs = names_inputs[sel_inputs] # names of selected inputs

  return(list(
    sel_inputs = sel_inputs,
    names_sel_inputs = names_sel_inputs,
    scores = scores
  ))

}
