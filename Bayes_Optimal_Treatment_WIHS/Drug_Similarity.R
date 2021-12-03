############################################# Similarity Matrix #################################################


eta <- 0.5
Classes <- c("NRTI", "NNRTI", "PI", "INSTI", "EI", "PE")
num_classes <- length(Classes)


Regimen_Reconstruct <- function(Z, num_classes){
  
  # reconstruct the regimen 
  # divide drugs into different classes and sort them alphabetically
  # args: Z: regimen used by subject i at visit j 
  #       num_classes: number of classes
  # returns: regimen: regimen with drugs divided into different classes
  
  # drug class index
  class_index <- sapply(abbr_group[Z], function(x) which(Classes == x))
  
  # divide drugs into different classes
  regimen <- sapply(1:num_classes, function(x) Z[which(class_index == x)])
  return_list <- list(regimen=regimen, class_index=class_index)
  return(return_list)
}



Third_Layer_Similarity <- function(Zx, Zy, num_classes, eta){
  
  # similarity in the third layer 
  # in this layer, their child nodes are leaves
  # we only need to count the number of same drugs
  # args: Zx, Zy: two regimens (after reconstruction), 
  #       num_classes: number of classes, eta: tuning parameter
  # returns: rho_third: similarities for each nodes
  #          simi_third: total similarities in third layer     
  
  # similarities for each nodes
  rho_third <- sapply(1:num_classes, function(x) sum(!is.na(match(Zx[[x]], Zy[[x]]))))
  
  # similarity in third layer
  simi_third <- eta*sum(rho_third)
  return_list <- list(rho_third=rho_third, simi_third=simi_third)
  return(return_list)
}



Second_Layer_Similarity <- function(Zx, Zy, num_classes, rho_third, eta){
  
  # similarity in the second layer 
  # if the same parts of two trees contain the same number of elements
  # the corresponding nodes of two trees on the second layer share the same subset tree structure
  # args: Zx, Zy: two regimens (after reconstruction), num_classes: number of classes, 
  #       rho_third: similarties for each nodes in third layer, eta: tuning parameter
  # returns: rho_second: similarities for each nodes
  #          simi_second: total similarities in second layer   
  
  # similarities for each nodes
  rho_second <- sapply(1:num_classes, function(x) 
    if (length(Zx[[x]])>0 & length(Zx[[x]])==length(Zy[[x]])){ eta*(1+eta)^(rho_third[x]) })
  rho_second[sapply(rho_second, is.null)] <- 0
  rho_second <- unlist(rho_second)
  
  # similarity in second layer
  simi_second <- sum(rho_second)
  return_list <- list(rho_second=rho_second, simi_second=simi_second)
  return(return_list)
}



First_Layer_Similarity <- function(Zx, Zy, Zx_class_index, Zy_class_index, num_classes, rho_second, eta){
  
  # similarity in the first layer
  # for the root, only need to check whether the drugs contained 
  # in these two regimens belong to the same classes
  # args: Zx, Zy: two regimens (after reconstruction), 
  #       Zx_class_index, Zy_class_index: class index for two regimens (after reconstruction), 
  #       num_classes: number of classes, 
  #       rho_second: similarties for each nodes in second layer, eta: tuning parameter
  # returns: simi_first: similarities in first layer   
  
  # indicator for usage of each drug classes
  Zx_ind <- 1:num_classes %in% Zx_class_index
  Zy_ind <- 1:num_classes %in% Zy_class_index
  
  # only if Zx and Zy contain the same drug classes
  if (all(Zx_ind == Zy_ind)){
    simi_first <- eta*prod(1+rho_second)
  }else{ simi_first <- 0 }
  
  return(simi_first)
}



Drug_Similarity <- function(Z, Zk){
  
  # drug similarity between regimen Z and kernel regimens Zk
  # args: Z: ART regimen, 
  #       Zk: a number of D representative regimens in the kernel regression
  # returns: drug_simi      
  
  D <- dim(Zk)[1] # number of representative regimens
  
  if (length(Z) == 0){
    drug_simi <- rep(0, D)
  }
  else{
    drug_simi <- rep(NA, D)
    # reconstruct regimen Zx
    Zx_list <- Regimen_Reconstruct(Z, num_classes) 
    Zx <- Zx_list$regimen; Zx_class_index <- Zx_list$class_index
    
    for (d in 1:D){
      similarity <- 0
      # reconstruct regimen Zy
      Zy_list <- Regimen_Reconstruct(as.vector(na.omit(Zk[d,])), num_classes)
      Zy <- Zy_list$regimen; Zy_class_index <- Zy_list$class_index
      
      # calculate similarities in three layers
      # third layer similarities 
      thrid_layer <- Third_Layer_Similarity(Zx, Zy, num_classes, eta)
      rho_third <- thrid_layer$rho_third
      similarity <- similarity + thrid_layer$simi_third
      
      # second layer similarities 
      second_layer <- Second_Layer_Similarity(Zx, Zy, num_classes, rho_third, eta)
      rho_second <- second_layer$rho_second
      similarity <- similarity + second_layer$simi_second
      
      # first layer similarities 
      similarity <- similarity + First_Layer_Similarity(Zx, Zy, Zx_class_index, Zy_class_index, num_classes, rho_second, eta)
      
      # drug similarities 
      drug_simi[d] <- similarity
    } 
  }
  
  return(drug_simi)
}