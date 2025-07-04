#' Transform fragments file
#'
#' This function transforms fragments file associating each fragment to a single-cell to a fragments file associating fragments to metacells
#'
#'
#' @param input_file path to original fragments file
#' @param output_name name of the supercell fragments file (whithout extension)
#' @param output_path path where the new fragments file will be saved
#' @param membership membership vector 
#' @param n_skip number of lines to skip in fragments file, (by default 0)
#'
#' @details
#' 
#'
#' @return 
#'
#'@examples
#'\dontrun{
#' 
#'}
#' @export

transform_fragment_file = function(input_file, output_name, output_path, membership, n_skip = 0){
  
  source_python(paste0(system.file(package = "SuperCellMultiome"),"/python_functions.py"))
  
  matching_names=data.frame(single_cell_names=names(membership),
                            super_cell_names=sapply(membership,function(j) paste0("SC_",j)))
  rownames(matching_names)=matching_names$single_cell_names
  write.table(matching_names,file=paste0(output_path,output_name,"_matching_names.txt"),row.names = F, col.names = T,quote = F,sep="\t")
  
  update_fragments(input_file,
                   output_name,
                   paste0(output_path,output_name,"_matching_names.txt"),
                   output_path,n_skip)
}

#' Transform fragments file parallel
#'
#' This function transforms fragments file associating each fragment to a single-cell to a fragments file associating fragments to metacells
#'
#'
#' @param input_file path to original fragments file
#' @param output_name name of the supercell fragments file (whithout extension)
#' @param output_path path where the new fragments file will be saved
#' @param membership membership vector 
#' @param tmp_path directory path for tempory files 
#' @param nb_cl number of cores to use. (default \code{nb_cl = parallel::detectCores()-2}) 
#' @param nb_row_split row number for split files
#' @param prefixMC prefix of metacell name in membership 
#' @param returnOutputFileName wheter to return the new fragment file name 
#'
#' @details
#' 
#'
#' @return 
#'
#'@examples
#'\dontrun{
#' 
#'}
#' @importFrom foreach %dopar%
#' @export


transform_fragment_file_parallel = function(input_file, 
                                            membership, 
                                            output_name = NULL, 
                                            output_path = NULL, 
                                            tmp_path = "./tmp/", 
                                            nb_cl = NULL,
                                            nb_row_split = "10000000",
                                            prefixMC = "SC_",
                                            bgzip_path = NULL,
                                            tabix_path = NULL,
                                            returnOutputFileName = T){
  
  if (is.null(bgzip_path)) {
    bgzip_command = "bgzip"
  } else {
    bgzip_command = bgzip_path
  }
  
  if (is.null(tabix_path)) {
    tabix_command = "tabix"
  } else {
    tabix_command = tabix_path
  }
  
  if (is.null(output_name)) {
    output_name <- paste0("MC_",fs::path_file(input_file))
  }
  
  if (is.null(output_path)) {
    output_path <- fs::path_dir(input_file)
  }
  
  output_path <- fs::path_abs(output_path)
  
  full_output_name <- paste0(output_path,"/",output_name)
  full_output_name_tsv <- fs::path_ext_remove(full_output_name) 
  
  dir.create(output_path,showWarnings = F,recursive = T)
  dir.create(tmp_path,showWarnings = F,recursive = T)
  
  if (file.exists(full_output_name)) {
    print(paste0(full_output_name, " aggregated framgent files already exists, it will be overwritten"))
  }
  
  if(is.null(nb_cl)){
    nb_cl = parallel::detectCores()-2
  } 
  cl <- parallel::makeCluster(nb_cl)
  doParallel::registerDoParallel(cl)
  
  init_path = getwd()
  setwd(tmp_path)
  
  matching_names=data.frame(single_cell_names=names(membership),
                            super_cell_names=sapply(membership,function(j) paste0(prefixMC,j)))
  
  rownames(matching_names)=matching_names$single_cell_names
  write.table(matching_names,file=paste0(output_name,"_matching_names.txt"),row.names = F, col.names = T,quote = F,sep="\t")
  
  # Fragment file split    
  print("Start fragments file split")
  system(command = paste0("gunzip -c ", input_file," | split -l ",nb_row_split," - frags_subset"))
  list_fragments = list.files(path = ".",pattern = "frags_subset")
  print("Start parallel update")
  foreach::foreach(i = 1:length(list_fragments),.packages = c("data.table")) %dopar% {
    setwd(tmp_path)
    if(i == 1){
      library(data.table)
      tmp = data.table::fread(paste0(list_fragments[i]))
    }else{
      tmp = data.table::fread(paste0(list_fragments[i]))
    } 
    tmp$V4 = matching_names[tmp$V4,"super_cell_names"] 
    
    data.table::fwrite(tmp,paste0("up_",list_fragments[i]),quote = F,row.names = F, col.names = F,sep="\t")
  }
  parallel::stopCluster(cl)
  #setwd(init_path)
  print("Start concat")
  system(command = paste0("cat ","up_* > ",full_output_name_tsv))
  print("End concat, start bgzip")
  # system(command = paste0(bgzip_command, " -h"))
  # print(paste0("bgzip command ",bgzip_command))
  system(command = paste0(bgzip_command, " -f ",full_output_name_tsv))
  print("End bgzip, start tabix")
  system(command = paste0(tabix_command, " -f -p bed ",full_output_name))
  print("End tabix")
  system(command = paste0("rm ","*frags_subset*"))
  setwd(init_path)
  
  if(returnOutputFileName) {
    return(full_output_name)
  }
  
}




#' AggregateFragmentFile parallel
#'
#' This function transforms fragments file associating each fragment to a single-cell to a fragments file associating fragments to metacells
#'
#'
#' @param input_file path to original fragments file
#' @param output_name name of the supercell fragments file (whithout extension)
#' @param output_path path where the new fragments file will be saved
#' @param membership membership vector 
#' @param tmp_path directory path for tempory files 
#' @param nb_cl number of cores to use. (default \code{nb_cl = parallel::detectCores()-2}) 
#' @param nb_row_split row number for split files
#' @param prefixMC prefix of metacell name in membership 
#' @param returnOutputFileName wheter to return the new fragment file name 
#'
#' @details
#' 
#'
#' @return 
#'
#'@examples
#'\dontrun{
#' 
#'}
#' @importFrom foreach %dopar%
#' @export


AggregateFragmentFile = function(input_file, 
                                 membership, 
                                 output_name = NULL, 
                                 output_path = NULL, 
                                 tmp_path = "./tmp/", 
                                 nb_cl = NULL,
                                 nb_row_split = "10000000",
                                 #prefixMC = "SC_",
                                 bgzip_path = NULL,
                                 tabix_path = NULL,
                                 returnOutputFileName = T){
  
  if (is.null(bgzip_path)) {
    bgzip_command = "bgzip"
  } else {
    bgzip_command = bgzip_path
  }
  
  if (is.null(tabix_path)) {
    tabix_command = "tabix"
  } else {
    tabix_command = tabix_path
  }
  
  if (is.null(output_name)) {
    output_name <- paste0("MC_",fs::path_file(input_file))
  }
  
  if (is.null(output_path)) {
    output_path <- fs::path_dir(input_file)
  }
  
  output_path <- fs::path_abs(output_path)
  
  full_output_name <- paste0(output_path,"/",output_name)
  full_output_name_tsv <- fs::path_ext_remove(full_output_name) 
  
  dir.create(output_path,showWarnings = F,recursive = T)
  dir.create(tmp_path,showWarnings = F,recursive = T)
  
  if (file.exists(full_output_name)) {
    print(paste0(full_output_name, " aggregated framgent files already exists, it will be overwritten"))
  }
  
  if(is.null(nb_cl)){
    nb_cl = parallel::detectCores()-2
  } 
  cl <- parallel::makeCluster(nb_cl)
  doParallel::registerDoParallel(cl)
  
  # init_path = getwd()
  # setwd(tmp_path)
  
  matching_names=data.frame(single_cell_names=names(membership),
                            super_cell_names=membership)
  
  rownames(matching_names)=matching_names$single_cell_names
  write.table(matching_names,file=paste0(tmp_path,"/",
                                                output_name,"_matching_names.txt"),
                                         row.names = F, col.names = T,quote = F,sep="\t")
  
  # Fragment file split    
  print("Start fragments file split")
  system(command = paste0("gunzip -c ", input_file," | split -l ",nb_row_split," - ",tmp_path,"/frags_subset"))
  list_fragments = list.files(path = tmp_path,pattern = "frags_subset",full.names = T)
  print("Start parallel update")
  foreach::foreach(i = 1:length(list_fragments),.packages = c("data.table")) %dopar% {
    # setwd(tmp_path)
    if(i == 1){
      library(data.table)
      tmp = data.table::fread(paste0(list_fragments[i]))
    }else{
      tmp = data.table::fread(paste0(list_fragments[i]))
    } 
    tmp <- tmp[tmp$V4 %in% rownames(matching_names),]
    tmp$V4 = matching_names[tmp$V4,"super_cell_names"]
    
    data.table::fwrite(tmp,paste0(list_fragments[i],"_up"),quote = F,row.names = F, col.names = F,sep="\t")
  }
  parallel::stopCluster(cl)
  #setwd(init_path)
  print("Start concat")
  system(command = paste0("cat ",tmp_path,"/*_up > ",full_output_name_tsv))
  print("End concat, start bgzip")
  # system(command = paste0(bgzip_command, " -h"))
  # print(paste0("bgzip command ",bgzip_command))
  system(command = paste0(bgzip_command, " -f ",full_output_name_tsv))
  print("End bgzip, start tabix")
  system(command = paste0(tabix_command, " -f -p bed ",full_output_name))
  print("End tabix")
  system(command = paste0("rm ",tmp_path,"/*frags_subset*"))
  # setwd(init_path)
  
  if(returnOutputFileName) {
    return(full_output_name)
  }
  
}
