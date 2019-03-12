#' Make Gene Names Legal Variable Names
#'
#' It is common that gene names are saved as names of rows or columns of an
#' expression matrix. In R, this could create some problems because some gene
#' names are not legal variable names and thus will be converted when they are
#' saved as column names (which R assumes that they should be legal variable
#' names). This automatic confusion might interfere with name matching.
#' \code{normalize_genename()} converts gene names to make them legal variable
#' names, thereby make sure the gene names compared are consistently after
#' conversion. \code{normalize_genename()} will add a prefix (gn_) to gene
#' names to make gene names starting with a number legal, and then calls
#' \code{\link[base]{make.names}()} to convert illegal punctuations.
#' @param gene a character vector containing gene names
#'
#' @return a character vector containing prefixed and converted gene names
#' @export
#'
#' @seealso \code{\link{denormalize_genename}()}
#'
#' @examples
#' normalize_genename("128up")
normalize_genename <- function(gene) {
  # Make gene names legal colnames which can be safely reversed by removing the "gn_" prefix
  genelist <- paste0("gn_", gene)
  genelist <- make.names(genelist, unique = FALSE)
  return(genelist)
}


#' Remove Prefix Added When Making Gene Names Legal Variable Names
#'
#' \code{denormalize_genename()} removes the prefix
#' \code{\link{normalize_genename}()} adds, and make the gene names more
#' readable for humans.
#' @inheritParams normalize_genename
#'
#' @return a character vector containing converted gene names with prefix
#' removed
#' @export
#' @examples
#' gene <- normalize_genename("128up")
#' # [1] "gn_128up"
#' denormalize_genename(gene)
#' # [1] "128up"
denormalize_genename <- function(gene) {
  name <- substr(start = 4, stop = nchar(gene), gene)
  return(name)
}


#' Make an Reference for Gene Name Conversion
#'
#' \code{make_database()} loads gene name / ID conversion tables from online
#' databases (e.g., FlyBase) and a user-provided GTF file that is used to
#' annotate the features of the NGS dataset of interest (e.g., the GTF file
#' STAR or \code{cellranger count} uses), and makes conversion tables from the
#' loaded data.
#'
#' @inheritParams prepare_database
#' @return a list containing multiple coversion tables and meta data
#' @export
#'
#' @examples
#' dummypath <- system.file("extdata", "dummy.gtf",
#' package = "genofeatutil")
#' dmeldb <- make_database(species = "dmel",
#'                            gtf.path = dummypath)
make_database <- function(species = "dmel", gtf.path, version = NULL) {
  db <- prepare_database(species = species,
                         gtf.path = gtf.path,
                         version = version)

  # Generate conversion vectors for FBgn update
  db <- generate_fbid_version_table(db)

  # Generate conversion vectors from the synonym table
  db <- generate_flybase_sym(db)
  db <- generate_gene_mapping(db)

  return(db)
}


#' Prepare an Unprocessed Database for Gene Name Conversion
#'
#' \code{prepare_database()} loads gene name / ID conversion tables from online
#' databases (e.g., FlyBase) and a user-provided GTF file that is used to
#' annotate the features of the NGS dataset of interest (e.g., the GTF file
#' STAR or \code{cellranger count} uses) into memory for later processing by
#' \code{\link{make_database}()}.
#'
#' @param species a character string describing the species of the genes.
#' (Options: "\code{dmel}")
#' @param gtf.path a character string containing the path to the GTF file
#' containing the input gene names (e.g., the GTF files used by
#' \code{cellranger count})
#' @param version a character specifying the version desired
#' (e.g., "FB2019_01".)
#'
#' @return a list containing data frames and character strings of meta data
prepare_database <- function(species = "dmel", gtf.path, version = NULL) {
  if (!species %in% c("dmel", "test")) {
    stop("prepare_database does not support ", species, " now.")
  }
  # Load GTF file (contains gene id and names)
  id_mapping <- rtracklayer::import(gtf.path)
  id_mapping <- as.data.frame(id_mapping)
  id_mapping <- id_mapping[id_mapping$type == "gene", ]

  # Load data from species-specific gene databases
  if (species == "dmel") {
    result <- fetch_flybase(version = version)
    version <- result[["version"]]
  }

  # Load data from local if in test environment
  if (species == "test") {
    fbgnpath <- system.file("extdata", "dummy_fbgn.tsv",
                            package = "genofeatutil")
    synopath <- system.file("extdata", "dummy_syno.tsv",
                            package = "genofeatutil")
    result <- fetch_flybase(paths = c(fbgnpath,
                                      synopath))
  }

  # Keep gene id and gene name from the GTF file
  gtf <- id_mapping[ , c("gene_id", "gene_name")]
  result[["gtf"]] <- gtf
  result[["metadata"]] <- c("species" = species,
                            "date" = date(),
                            "GTF_path" = gtf.path,
                            "FlyBase_ver" = version)

  return(result)
}


#' Generate Conversion Tables for Symbols and Aliases to FlyBase ID
#'
#' @param db a list generated by \code{\link{prepare_database}()}
#'
#' @return a list of 2 (a named vector that translate symbols to FB ids and
#' another that translates aliases to FB ids)
generate_flybase_sym <- function(db) {
  # Generate a a named vector for conversion based on Flybase Synonym
  # with which you can use alias/name to query FBid
  # Usage: named_vector <- generate_flybase_sym([the tsv file])
  # named_vector[*alias/name*] gives its corresponding flybase id
  fbsym <- db[["syno"]]

  # Keep only Drosophila genes
  colnames(fbsym)[1] <- "primary_fbid"
  fbsym <- fbsym[fbsym$organism_abbreviation == "Dmel", ]
  fbsym <- fbsym[grepl("^FBgn", fbsym$primary_fbid), ]

  # Generate a named vector with *aliases as names* and
  # *fbid as values*
  alias_dict_list <- list()
  for (row in seq(nrow(fbsym))) {
    alias <- strsplit(fbsym[row, "symbol_synonym(s)"], ",")[[1]]
    if (length(alias) > 0) {
      fbid <- fbsym[row, "primary_fbid"]
      dict <- rep(fbid, length(alias))
      names(dict) <- normalize_genename(alias)
      alias_dict_list[[fbid]] <- dict
    }
  }

  alias_dict <- unlist(alias_dict_list)

  # Since unlist appends names of the list items to names of vector items
  # For example, if ChAT is #10000 in the file, before unlist() the
  # name:value would be ChAT:FBid
  # After unlist it would be 10000.ChAT:FBid
  # I'll remove the list item names (the "10000." part) for they disrupt
  # mapping
  names(alias_dict) <- gsub("FBgn[0-9]*\\.", "", names(alias_dict))

  # Generate a named vector with *symbols as names* and
  # *fbid as values*
  symbol_dict <- fbsym$primary_fbid
  names(symbol_dict) <- normalize_genename(fbsym$current_symbol)

  # Add conversion vectors into db and remove raw data
  db[["symbol_dict"]] <- symbol_dict
  db[["alias_dict"]] <- alias_dict
  db <- db[names(db) != c("syno")]
  return(db)
}


#' Generate a Conversion Vector for FlyBase IDs from Older Versions
#'
#' \code{generate_fbid_version_table()} takes the FBgn ID annotation table
#' from FlyBase, and prepare a conversion vector.
#' \code{generate_fbid_version_table()} takes a list generated by
#' \code{\link{prepare_database}()} and updates it by adding the conversion
#' vector and removing the raw FBgn ID annotation table.
#' @param db a list generated by \code{prepare_database()}
#'
#' @return a list without FBgn ID annotation table but with a derived
#' conversion vector
generate_fbid_version_table <- function(db) {
  # Read the FBgn annotation table from FlyBase
  id_table <- db[["fbid"]]

  # Report the version used
  version <- db[["metadata"]]["FlyBase_ver"]
  message("Using FlyBase version: ", version, " to update FBgn.")

  # Generate a conversion vector
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

  # Remove raw data and return db with the conversion vector
  db[["id_dict"]] <- lookup
  db <- db[names(db) != "fbid"]
  return(db)
}


#' Generate a Mapping Vector to Translate Gene IDs
#'
#' \code{generate_gene_mapping()} takes gene names from a character vector, an
#' expression matrix / data frame, or a Seurat object, and generate a vector
#' that can be used to translate gene IDs to the gene names provided to
#' \code{generate_gene_mapping()}. The purpose of this translation is to make
#' sure different datasets are consistent in the gene names used during
#' analysis, and to thereby ensure comparsion between datasets or querying are
#' accurate.
#' @param db a list generated by \code{\link{prepare_database}()}
#' @param ... other arguments that are passed to \code{\link{get_gene_names}()}
#' (see below)
#' @inheritParams get_gene_names
#'
#' @return a named character vector, in which the values are gene names and
#' the names are gene ids
generate_gene_mapping <- function(db) {
  # Check if db is legit
  if (!"gtf" %in% names(db)) {
    stop(paste("The database list seems to be wrong. Please make sure that",
               "you generated it by prepare_database() before using gene",
               "name conversion functions."))
  }
  id_mapping <- db[["gtf"]]

  # Update FBgn
  if (db[["metadata"]]["species"] %in% c("dmel", "test")) {
    id_mapping$gene_id <- update_fbgn(id_mapping$gene_id, db = db)
  }

  result <- id_mapping$gene_name
  names(result) <- id_mapping$gene_id
  db[["to_name_dict"]] <- result
  db <- db[names(db) != "gtf"]
  return(db)
}


#' Update FlyBase Gene Numbers from an Older Version
#'
#' \code{update_fbgn()} takes a character vector of FlyBase gene numbers and a
#' list from \code{\link{make_database}()}, and convert the vector according
#' to the list.
#' @param id a character vector containing FlyBase gene numbers
#' @param db a list from \code{make_database()}
#'
#' @return a character vector containing updated FlyBase gene numbers
#' @export
#'
#' @examples
#' dummypath <- system.file("extdata", "dummy.gtf",
#' package = "genofeatutil")
#' dmeldb <- make_database(species = "dmel",
#'                            gtf.path = dummypath)
#' update_fbgn("FBgn0032045", db = dmeldb)
update_fbgn <- function (id, db) {
  # Convert out-dated FBid to current version
  ## Load lookup table
  if ("id_dict" %in% names(db)) {
    version_dict <- db[["id_dict"]]
  } else {
    stop(paste("The database list seems to be wrong. Please make sure that",
               "you generated it by prepare_database() before using gene",
               "name conversion functions."))
  }

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


#' Convert Gene Names to FlyBase Gene Numbers
#'
#' \code{convert_gene_to_fbgn()} takes a vector of gene names and query the
#' synonym table from FlyBase to convert them to a vector of FlyBase gene
#' numbers. Due to potential multiple mapping (one alias corresponding to many
#' FlyBase gene numbers, like Cha can be ChAT or Charaff), it's important to
#' note that when there's warning message, the vector
#' /code{convert_gene_to_fbgn()} returns will not retain the same order as the
#' input vector of gene names.
#'
#' @param genes a character vector containing gene names
#' @param db a list from \code{make_database()}
#'
#' @return a character vector containing FlyBase gene numbers
#' @export
#'
#' @examples
#' # Prepare the database (using dummy local data)
#' dummypath <- system.file("extdata", "dummy.gtf",
#' package = "genofeatutil")
#' testdb <- make_database(species = "test",
#' gtf.path = dummypath)
#'
#' # Convert gene names
#' convert_gene_to_fbgn("CG31610", db = testdb)
convert_gene_to_fbgn <- function(genes, db) {
  # Separate gene names to 2 categories and convert to FBgn
  # 1. Matching current symbol
  # 2. Matching aliases

  ## Reading from the db
  if (!"symbol_dict" %in% names(db) | !"alias_dict" %in% names(db)) {
    stop(paste("The database list seems to be wrong. Please make sure that",
               "you generated it by prepare_database() before using gene",
               "name conversion functions."))
  }
  symbol_dict <- db[["symbol_dict"]]
  alias_dict <- db[["alias_dict"]]


  ## Remove mitochondrial genes
  if (length(genes[!grepl("^mt:", genes)]) > 0) {
    unordered <- TRUE
  }
  genes <- genes[!grepl("^mt:", genes)]
  genes <- normalize_genename(genes)

  ## Here I suppose that genes with names matching official symbols do
  ## not need conversion
  index_symbol <- genes %in% names(symbol_dict)
  genes[index_symbol] <- symbol_dict[genes[index_symbol]]

  # Find the names that are not official symbols but are documented aliases
  index_alias <-
    !(genes %in% names(symbol_dict)) & (genes %in% names(alias_dict))
  genes[index_alias] <- alias_dict[genes[index_alias]]

  # Keep only the FBgn of official symbols and the aliases and drop uncoverted
  # ones
  unmapped <- genes[!(index_symbol | index_alias)]
  genes <- genes[index_symbol | index_alias]

  ## Create warnings if there's multiple mapping for aliases
  for (item in genes[index_alias]) {
    matching_num <- sum(names(alias_dict) == item)
    if (matching_num > 1) {
      unordered <- TRUE
      warning(paste0("'", denormalize_genename(item),
                     "' is matched with multiple aliases. All the FBgn ",
                     "that correspond to it will be append to the end ",
                     "of the result.")
                     )

      # Create a vector of the rest of the alias:FBgn pair a non-unique
      # alias mapped to, and append that to the alias vector.
      multimap_id <- alias_dict[names(alias_dict) == item][-1]
      alias <- append(x = genes,
                      values =  multimap_id)
    }
  }

  ## Remind the user about non-convertible gene names
  if (length(unmapped) > 0) {
    unordered <- TRUE
    unmapped <- paste(unmapped, collapse = ", ")
    warning(paste("Please note the following genes are not found in",
                  "the database and won't be processed:",
                   denormalize_genename(unmapped), "."))
  }
  if (unordered) {
    message(paste("Please note that the conversion result of",
                  "convert_gene_to_fbgn might not retain the same",
                  "order of the vector of gene names that are converted."))
  }
  return(genes)
}


#' Convert Gene Names or FlyBase Gene Numbers to Gene Names Specified in the
#' Database
#'
#' \code{convert_to_genename()} takes a character vector of gene names and
#' FlyBase gene numbers and uses the database list
#' \code{\link{make_database}()} generated to convert them to gene names that
#' are consistent with the version of the provided GTF file.
#'
#' @inheritParams distinguish_fbgn
#' @param db a list from \code{make_database()}
#' @param normalize a logical expression. If set as TRUE, it will call
#' \code{\link{normalize_genename}()} to make consistent the converted gene
#' names.
#'
#' @export
#' @examples
#' # Make databases from dummy data
#' dummypath <- system.file("extdata", "dummy.gtf",
#' package = "genofeatutil")
#' testdb <- make_database(species = "test",
#' gtf.path = dummypath)
#'
#' convert_to_genename("FBgn0086917", testdb)
convert_to_genename <- function(x, db, normalize = TRUE) {
  if (!"to_name_dict" %in% names(db)) {
    stop(paste("The database list seems to be wrong. Please make sure that",
               "you generated it by prepare_database() before using gene",
               "name conversion functions."))
  }
  genename_index <- !distinguish_fbgn(x)
  if (length(which(genename_index)) > 0) {
    # Convert gene names to FBgn
    converted_names <- convert_gene_to_fbgn(x[genename_index], db)
    if (length(converted_names) == length(x[genename_index])) {
      x[genename_index] <- converted_names
    } else {
      x <- c(x[!genename_index], converted_names)
      warning(paste("Because there were multiple mappings of the aliases,",
                    "please note that the length and  order of the output is",
                    "different from the input."))
    }
  }

  # Examine duplication
  dup_all <- Reduce(`|`, duplicated(x))
  if (dup_all) {
    x <- unique(x)
    warning(paste("There are duplications of genes in your input, and",
                  "duplicated items are removed. As a result, the order and",
                  "length of output won't be the same as the input."))
  }

  # Convert to gene name
  dict <- db[["to_name_dict"]]
  mapped_index <- x %in% names(dict)
  map_all <- Reduce(`&`, mapped_index)
  if(map_all) {
    result <- dict[x]
  } else {
    mapped <- x[mapped_index]
    result <- dict[mapped]
    warning(paste("The following FlyBase gene numbers are not found in the",
                  "GTF file you loaded:",
                  paste(x[!mapped_index], collapse = ", ")))
  }

  if (normalize) {
    result <- normalize_genename(result)
  }
  return(result)
}

#' Distinguish FlyBase Gene Numbers from Gene Names
#'
#' @param x a character vector containing genes and FlyBase gene numbers
#'
#' @return a logical vector
distinguish_fbgn <- function(x) {
  if (!is.character(x)) {
    stop(paste("distinguish_fbgn() takes only a character vector of FlyBase",
               "gene numbers or gene names"))
  }
  FBgn_index <- grepl("^FBgn", x)
  return(FBgn_index)
}