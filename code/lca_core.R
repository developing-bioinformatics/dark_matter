function(results, parallel=FALSE, nclus=4) {
  require(dplyr)
  require(multidplyr)
  
  sub_lca = function(x) {
    nombre = vector()
    taxnames = c('superkingdom',
                 'phylum',
                 'order',
                 'family',
                 'genus',
                 'species')
    shortnames = apply(x[, taxnames], 2, unique)
    countshnames = sapply(shortnames, length)
    numcount = countshnames == 1
    lastuni = tail(names(shortnames[numcount == T]), n = 1)
    nombre = as.data.frame(x[1, which(colnames(x) == lastuni)])
    newtax <- as.list(ifelse(countshnames == 1, shortnames, NA))
    
    #if(length(nombre) == 0){ 
    if(all(is.na(newtax))) {
      ret = x %>%
        mutate(last_common = NA) %>%
        mutate(level_lca = NA) %>%
        mutate(superkingdom = NA) %>%
        mutate(phylum = NA) %>%
        mutate(class = NA) %>%
        mutate(order = NA) %>%
        mutate(family = NA) %>%
        mutate(genus = NA) %>%
        mutate(species = NA)
    } else {
      ret = x %>%
        mutate(last_common = as.character(nombre[[1]])) %>%
        mutate(level_lca = lastuni) %>%
        mutate(superkingdom = newtax$superkingdom) %>%
        mutate(phylum = newtax$phylum) %>%
        mutate(class = newtax$class) %>%
        mutate(order = newtax$order) %>%
        mutate(family = newtax$family) %>%
        mutate(genus = newtax$genus) %>%
        mutate(species = newtax$species)
    }
    return(ret)
  }
  if(parallel == TRUE){
    cluster <- new_cluster(nclus)
    cluster_library(cluster, 'dplyr')
    cluster_copy(cluster, 'sub_lca')
    #split groups across multiple CPU cores
    prepare = results %>%
      group_by(QueryID) %>%
      partition(cluster) %>%  #split groups into cluster units
      do({sub_lca(.)}) %>%
      collect()
  } else {
    prepare = results %>%
      group_by(QueryID) %>%
      do({sub_lca(.)})
  }
  return(prepare)
}