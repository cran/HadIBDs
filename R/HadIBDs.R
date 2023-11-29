#'Incomplete Block Designs using Hadamard Matrix (HadIBDs)
#'
#' @param v is expressed as product of (4t_i-1), where t_i = 2^x,(i=1,2,...) and (x = 0,1,2...)
#'@importFrom utils combn
#' @return This function generates an IBD based on modified Hadamard matrices or their Kronecker product along with the Parameters, Information matrix, Average variance factor and Canonical efficiency factor of the generated design.
#' @references 1) R.C. Bose, K.R. Nair (1939). Partially balanced incomplete block designs, Sankhya 4, 337-372. https://www.jstor.org/stable/40383923.
#'
#' 2) M.N. VARTAK (1955). On an application of Kronecker product of matrices to statistical
#' designs,The Annals of Mathematical Statistics 26, 420-438.

#' @export
#'
#' @examples
#' library(HadIBDs)
#' Hadamard_to_IBDs(9)
Hadamard_to_IBDs<-function(v){
  if(v<7){
    v_less_7<-TRUE
    v<-7
    message(paste("As the entered value of  v is not the product of (4t_i-1) form, (i=1,2,...) and t_i = 2^x, (x = 0,1,2...); design will not exist. The design for the nearest possible parameter is for v =",v))
    cat("\n")
  }else{
    v_less_7<-FALSE
  }
  factors_we_need<-function(v){
    v=prod(v)
    temp=prod(v)
    store=c()
    i=2
    while(i<=temp){
      if(temp%%i==0){
        store=c(store,i)
        temp=v/prod(store)
        i=1
      }
      i=i+1
    }
    ########
    store=sort(store)
    return(store)
  }
  #########
  steps<-NULL
  for(i in 1:50){
    steps<-c(steps,kronecker(c(-1,1),i))
  }
  set<-c(v,(v+steps))
  #if_3_5<-c(3,5)
  #################
  ################final factors we got
  for(value in set){
    ###case 1 : All are (4t-1) form
    #value=33

    check_35<-c(3,5)
    if(all((factors_we_need(value)+1)%%4==0) || all(check_35 %in% factors_we_need(value))){
      for(k in (factors_we_need(value)+1)){
        # Factor wise factors_to_be_needed_checking
        not_needed<-0
        if(length(setdiff(factors_we_need(k),2))!=0){  #only 3 and 5 containing vectors excluding
          #not_needed<-0
          is_there_three<-length(which(factors_we_need(value)==3))
          is_there_five<-length(which(factors_we_need(value)==5))
          if(is_there_five>=1 && is_there_three>=1) { # 15 adding portion
            #is_there_three>=is_there_five || is_there_three<=is_there_five
            add_15<-rep(15,min(is_there_three,is_there_five))
            rest<-prod(factors_we_need(value))/(prod(add_15))
            final_factors<-c(rest,add_15)
            if(all((final_factors+1)%%4==0)){
              not_needed<-1
              break
            }
          }
        }


        if(length(setdiff(unique(factors_we_need(prod(factors_we_need(value)+1))),2))==0){
          final_factors<-factors_we_need(value)
          not_needed=1
          break
        }

        if(length(setdiff(unique(factors_we_need(prod(factors_we_need(value))+1)),2))==0){
          final_factors<-factors_we_need(value)
          not_needed=1
          break
        }
      }
      if(not_needed==1){
        break
      }

    }
  }
  ####################### Direct or indirect
  if(length(setdiff(unique(factors_we_need(prod(final_factors)+1)),2))==0){
    final_factors<-prod(final_factors)
  }

  ##################################

  ########################%%%%%%%%%%%#######################$$$$$$$$$$$$$$$
  ###################### Error message
  if(prod(final_factors)!=v && v_less_7==FALSE){
    message(paste("As the entered value of  v is not the product of (4t_i-1) form, (i=1,2,...) and t_i = 2^x, (x = 0,1,2...); design will not exist. The design for the nearest possible parameter is for v =",prod(final_factors)))
    cat("\n")
  }
  ##################
  ############Hadamard we need of order power of 2
  Hadamard_generate<-function(order){
    times<-c()
    i=2
    while(i<=order ){
      if(order%%2==0){
        times<-c(times,2)
        order=order/2
        i=2
      }
    }
    times<-length(times)
    H2<-matrix(c(1,1,1,-1),nrow=2,byrow=T)
    #############
    kronecker_product<-function(mat1,mat2){
      return(kronecker(mat1,mat2))
    }
    ##################
    final_layer<-rep(2,times)
    HadamardMatrix<-1
    for(k in 1:length(final_layer)){
      HadamardMatrix<-kronecker_product(HadamardMatrix,H2)
    }
    return(HadamardMatrix)
  }
  ###############################################################
  hadamard_of_all_final_factors<-list()
  for(l in (final_factors)){
    hadamard_of_all_final_factors<-append(hadamard_of_all_final_factors,list(Hadamard_generate(l+1)[-1,-1]))
  }
  ####################Resultant design
  kronecker_product<-function(mat1,mat2){
    return(kronecker(mat1,mat2))
  }
  final_hadamard_we_need<-Reduce(kronecker_product,hadamard_of_all_final_factors)
  #############################
  ###############################create hadamard now
  products<-final_hadamard_we_need

  ################make incidence matrix
  products[products==-1]<-0
  incidence_matrix=products
  ############
  ####################incidence to design
  design_should_be<-NULL
  #design_should_be<-t(apply(incidence_matrix,1,function(x) as.which(x==1)))
  for(i in 1:nrow(incidence_matrix)){
    design_should_be<-rbind(design_should_be,t(which(incidence_matrix[i,]==1)))
    #print(c(sum(incidence_matrix[i,]),i))
  }
  ######################
  v=max(design_should_be)
  ################## replication matrix
  rep=length(which(design_should_be==design_should_be[1,1]))
  rep_matrix=diag(rep,nrow=nrow(design_should_be))
  ####################block size matrix
  k_matrix=diag(ncol(design_should_be),nrow=nrow(design_should_be),ncol=nrow(design_should_be))
  ################### c matrix
  c_matrix=rep_matrix-(incidence_matrix)%*%solve(k_matrix)%*%t(incidence_matrix)
  #####################################
  eig=eigen(c_matrix)$values
  eig = eig[eig>(10^(-10))]
  eig=round(eig,3)
  ############variance and avg variance
  ###########combn function
  p_matrix<-matrix(0,nrow=choose(v,2),ncol=v)
  comb_mat<-t(combn(v,2))
  p_matrix[c(comb_mat[,1],comb_mat[,2])]<-c(1,-1)
  ########## Variance covariance part
  variances<-(p_matrix)%*%MASS::ginv(c_matrix)%*%t(p_matrix)
  #######variance part
  var<-diag(variances)
  ########## avg variance
  Avg_var<-mean(var)
  ########CEF
  #CEF=rep_matrix[1,1]/2*Avg_var
  cef=(1/mean(1/eig))*(1/rep_matrix[1,1])
  ##################message
  if(length(unique(eig))==1){
    list_print=list('BIB_Design' = design_should_be,'Number of Treatments (v)'=v,'Number of Blocks (b)'=nrow(design_should_be),'Number of Replications(r)'=rep_matrix[1,1],'Block Size (k)'=ncol(design_should_be), 'C_Matrix' = round(c_matrix,4), 'Variance_Factor'=round(Avg_var,4), 'Cannonical_Efficiency_Factor'=round(cef,4))
    return(list_print)
  }else{
    list_print=list("PBIB_Design" = design_should_be,"Number of Treatments (v)"=v,"Number of Blocks (b)"=nrow(design_should_be),"Number of Replications(r)"=rep_matrix[1,1],"Block Size (k)"=ncol(design_should_be),'C_Matrix' = round(c_matrix,4), 'Variance_Factor'=round(Avg_var,4), 'Cannonical_Efficiency_Factor'=round(cef,4))
    return(list_print)
  }
}
