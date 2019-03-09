get_10Xgenemapping <- function(object, assay = "RNA", gtfpath) {
  # Map the gene names according to the gtf used for alignment, and
  # return a dataframe of the genes and their corresponding FBgn.
  # The Flybase gene number might be out-of-date though.
  ## Extract the names of all detected genes of an assay in a Seuratv3 object
  genes <- get_10Xgenenames(object, assay = assay)

  # Load .gtf file used for cellranger count
  gene_mapping <- rtracklayer::import(gtfpath)
  gene_mapping <- as.data.frame(gene_mapping)
  gene_mapping <- dplyr::filter(gene_mapping, type == "gene",
                                gene_name %in% genes)
  gene_mapping <- dplyr::select(gene_mapping, gene_id, gene_name)
  gene_mapping$gene_id <- update_fbid(gene_mapping$gene_id)
  result <- gene_mapping$gene_name
  names(result) <- gene_mapping$gene_id
  return(result)
}

update_fbid <- function (id) {
  # Convert out-dated FBid to current version
  ## Load lookup table
  if (!file.exists(fbid_obj)) {
    generate_fbid_version_table(fbid_table)
  }
  version_dict <- readRDS(fbid_obj)
  ## Find FBid that need conversion
  index_update <- id %in% names(version_dict)
  ## Convert with lookup table while keeping the order
  result <- sapply(seq(length(id)), function (x) {
    if (index_update[x]) {
      return(version_dict[id[x]])
    }
    return(id[x])
  })
  result <- unname(result)
  return(result)
}

generate_fbid_version_table <- function(path) {
  ## Read the FBgn annotation table from FlyBase
  id_table <- data.table::fread(input = path, header = TRUE,
                                sep = "\t", stringsAsFactors = FALSE,
                                data.table = TRUE)
  print(contable_ver)
  query <- id_table$`secondary_FBgn#(s)`
  names(query) <- id_table$`primary_FBgn#`
  lookup <- lapply(names(query), function(x) {
    old_id <- strsplit(query[[x]], split = ",")[[1]]
    if (length(old_id) > 0) {
      result <- rep(x, length(old_id))
      names(result) <- old_id
      return(result)
    }
  })
  lookup <- lookup[sapply(lookup, function(x) !is.null(x))]
  lookup <- unlist(lookup)
  saveRDS(object = lookup, file = fbid_obj)
  print("The version table is sucessfully generated and saved.")
  return(0)
}

generate_flybase_sym <- function (path) {
  # Read the tsv file of Flybase Synonym
  # http://flybase.org/cgi-bin/get_static_page.pl?file=bulkdata7.html&title=Current%20Release
  # and generate a named vector, with which you can use alias/name to query FBid
  # Usage: named_vector <- generate_flybase_sym([the tsv file])
  # named_vector[*alias/name*] gives its corresponding flybase id
  fbsym <- data.table::fread(file = path, sep = "\t", header = T,
                             stringsAsFactors = F, data.table = F, skip = 4, quote = "")
  # Keep only Drosophila genes
  colnames(fbsym)[1] <- "primary_fbid"
  fbsym <- filter(fbsym, organism_abbreviation == "Dmel")
  fbsym <- fbsym[grepl("^FBgn", fbsym$primary_fbid), ]

  # Generate a named vector with *aliases as names* and
  # *fbid as values*
  alias_dict_list <- apply(fbsym, 1, function(x) {
    alias <- strsplit(x["symbol_synonym(s)"], ",")[[1]]
    if (length(alias) > 0) {
      dict <- rep(x["primary_fbid"], length(alias))
      names(dict) <- normalize_genename(alias)
      return(dict)
    }
  })
  alias_dict <- unlist(alias_dict_list)

  # Since unlist appends names of the list items to names of vector items
  # For example, if ChAT is #10000 in the file, before unlist() the name:value would be ChAT:FBid
  # After unlist it would be 10000.ChAT:FBid
  # I'll remove the list item names (the "10000." part) for they disrupt mapping
  names(alias_dict) <- gsub("^[0-9]*\\.", "", names(alias_dict))
  symbol_dict <- fbsym$primary_fbid

  # Generate a named vector with *symbols as names* and
  # *fbid as values*
  names(symbol_dict) <- normalize_genename(fbsym$current_symbol)
  dict_obj <- list("symbol" = symbol_dict,
                   "alias" = alias_dict)
  return(dict_obj)
}

categorize_gene <- function(name, flybase_sym) {
  # Separate gene names to 2 categories
  # 1. Matching current symbol
  # 2. Matching aliases

  ## Remove mitochondrial genes
  name <- name[!grepl("^mt:", name)]
  name <- normalize_genename(name)
  ## Here I suppose that genes with names matching official symbols do not need conversion
  in_symbol <- name[name %in% names(flybase_sym[["symbol"]])]
  names(in_symbol) <- flybase_sym[["symbol"]][in_symbol]
  not_in_symbol <- name[!name %in% names(flybase_sym[["symbol"]])]
  alias <- not_in_symbol[not_in_symbol %in% names(flybase_sym[["alias"]])]
  names(alias) <- flybase_sym[["alias"]][alias]

  ## Create warnings if there's multiple mapping for aliases
  for (item in not_in_symbol) {
    matching_num <- sum(names(flybase_sym[["alias"]]) == item)
    if (matching_num > 1) {
      warning(paste0("'", denormalize_genename(item), "' is matched with multiple aliases. ",
                     "Please check the FB id generated here to determine if the mapping is correct."))
    }
  }

  ## Remind the user about non-convertible gene names
  unmapped <- setdiff(not_in_symbol, alias)
  if (length(unmapped) > 0) {
    unmapped <- paste(unmapped, collapse = ", ")
    warning(paste0("Please note the following genes are not found in the database and won't be processed: ",
                   denormalize_genename(unmapped), "."))
  }
  genelist <- list("symbol" = in_symbol, "alias" = alias)
  return(genelist)
}

convert_fbid <- function(name, flybase_sym) {
  # Convert a list of fly gene to flybase id
  ## Get a list of dataframes with a conversion table with Flybase IDs and symbol/alias
  genelist <- categorize_gene(name, flybase_sym = flybase_sym)
  fbid <- unique(c(names(genelist[["symbol"]]), names(genelist[["alias"]])))
  return(fbid)
}

convert_10Xgenename <- function(name, flybase_sym, dict10X) {
  # Convert a list of fly gene names to make it compatible
  # with the 10X gene names in the Seurat object

  ## Generate a named list from the Seurat object
  ## work as a dictionary for conversion

  ## Convert input names to Fbid
  input_fbid <- convert_fbid(name, flybase_sym)
  ## Throw away genes that are not present in 10X
  input_fbid <- input_fbid[input_fbid %in% names(dict_10X)]
  name_10Xformat <- dict_10X[input_fbid]
  name_10Xformat <- normalize_genename(name_10Xformat)
  return(name_10Xformat)
}