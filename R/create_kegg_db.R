##' create KEGG.db package
##'
##' 
##' @title create_kegg_db
##' @param species one of KEGG supported species, e.g. hsa for human
##' @param author a string for the Author field in the generated package DESCRIPTION
##' @param maintainer a string for the Maintainer field in the generated package DESCRIPTION
##' @return KEGG.db package generated in working directory
##' @export
##' @author Guangchuang Yu and Ziru Chen
create_kegg_db <- function(species, author = NULL, maintainer = NULL) {
    packagedir <- tempfile() # tempdir() maynot empty

    ## skeleton
    prepare_pkg_skeleton(packagedir, author = author, maintainer = maintainer)

    ## sqlite
    sqlite_path <- paste(packagedir, "inst", "extdata", sep=.Platform$file.sep)
    prepare_kegg_db(species, sqlite_path)

    ## build pkg
    pkgbuild::build(packagedir, dest_path = ".")
}

##' Create KEGG.db package for microbiota
##'
##' This function queries KEGG organism list and builds KEGG.db for selected microbiota groups.
##'
##' @title create_kegg_db_microbiota
##' @param output_dir output directory to write the generated KEGG.db package
##' @param include microbiota groups to include
##' @param author a string for the Author field in the generated package DESCRIPTION
##' @param maintainer a string for the Maintainer field in the generated package DESCRIPTION
##' @return paths of the generated package tarballs
##' @export
create_kegg_db_microbiota <- function(
  output_dir = ".",
  include = c("bacteria", "archaea", "fungi", "viruses"),
  author = NULL,
  maintainer = NULL
) {
  include <- unique(as.character(include))
  valid_include <- c("bacteria", "archaea", "fungi", "viruses")
  include <- intersect(include, valid_include)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  org_url <- "https://rest.kegg.jp/list/organism"
  org <- read_kegg_table(org_url)
  if (is.null(org) || nrow(org) == 0) stop("Failed to download KEGG organism list.")

  if (ncol(org) < 4) {
    org[[4]] <- ""
  }
  colnames(org)[1:4] <- c("kegg_taxid", "species", "name", "lineage")

  keep <- rep(FALSE, nrow(org))
  if ("bacteria" %in% include) keep <- keep | grepl("^Prokaryotes;Bacteria", org$lineage)
  if ("archaea" %in% include) keep <- keep | grepl("^Prokaryotes;Archaea", org$lineage)
  if ("fungi" %in% include) keep <- keep | grepl("^Eukaryotes;Fungi", org$lineage)
  if ("viruses" %in% include) keep <- keep | grepl("^Viruses", org$lineage)

  species <- unique(org$species[keep])
  species <- species[nzchar(species)]

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(output_dir)

  before <- list.files(".", pattern = "^KEGG\\.db_.*\\.tar\\.gz$", full.names = TRUE)
  create_kegg_db(species, author = author, maintainer = maintainer)
  after <- list.files(".", pattern = "^KEGG\\.db_.*\\.tar\\.gz$", full.names = TRUE)

  setdiff(after, before)
}


prepare_pkg_skeleton <- function(packagedir, author = NULL, maintainer = NULL) {
    .fcp <- function(..., todir, file) {
        file.copy(from = system.file("KEGG.db", ..., file, package = "createKEGGdb"),
                  to = paste(todir, file, sep = .Platform$file.sep))
    }

    if(!dir.exists(packagedir)) {
        dir.create(packagedir)
    }

    ## to store sqlite
    sqlite_path <- paste(packagedir, "inst", "extdata", sep=.Platform$file.sep)
    if(!dir.exists(sqlite_path)){
        dir.create(sqlite_path,recursive = TRUE)
    }

    R_src <- paste(packagedir, "R", sep=.Platform$file.sep)
    if(!dir.exists(R_src)){
        dir.create(R_src,recursive = TRUE)
    }

    .fcp("R", todir = R_src, file = "zzz.R")
    .fcp(todir = packagedir, file = "DESCRIPTION")
    .fcp(todir = packagedir, file = "LICENSE")
    .fcp(todir = packagedir, file = "NAMESPACE")

    if (!is.null(author) || !is.null(maintainer)) {
        desc_path <- file.path(packagedir, "DESCRIPTION")
        desc_lines <- readLines(desc_path, warn = FALSE)
        if (!is.null(author)) {
            desc_lines <- sub("^Author:.*$", paste0("Author: ", author), desc_lines)
        }
        if (!is.null(maintainer)) {
            desc_lines <- sub("^Maintainer:.*$", paste0("Maintainer: ", maintainer), desc_lines)
        }
        writeLines(desc_lines, con = desc_path)
    }
}


##' @importFrom magrittr %<>%
read_kegg_table <- function(url) {
  tryCatch(
    utils::read.delim(
      url,
      sep = "\t",
      header = FALSE,
      quote = "",
      stringsAsFactors = FALSE,
      fill = TRUE,
      comment.char = ""
    ),
    error = function(e) NULL
  )
}

kegg_list <- function(db, organism = NULL) {
  url <- if (is.null(organism)) {
    paste0("https://rest.kegg.jp/list/", db)
  } else {
    paste0("https://rest.kegg.jp/list/", db, "/", organism)
  }

  df <- read_kegg_table(url)
  if (is.null(df) || nrow(df) == 0 || ncol(df) < 2) return(NULL)
  df[, 1:2, drop = FALSE]
}

kegg_link <- function(organism, db) {
  url <- paste0("https://rest.kegg.jp/link/", organism, "/", db)
  df <- read_kegg_table(url)
  if (is.null(df) || nrow(df) == 0 || ncol(df) < 2) return(NULL)
  df[, 1:2, drop = FALSE]
}

get_path2name <- function(species){
  safe_kegg_list <- function(db, organism = NULL) {
    if (is.null(organism)) return(kegg_list(db))
    kegg_list(db, organism)
  }

  if (length(species) == 1) {
    keggpathid2name.df <- safe_kegg_list("pathway", species)
  } else {
    keggpathid2name.list <- vector("list", length(species))
    names(keggpathid2name.list) <- species
    for (i in species) {
      keggpathid2name.list[[i]] <- safe_kegg_list("pathway", i)
    }
    keggpathid2name.list <- keggpathid2name.list[!vapply(keggpathid2name.list, is.null, logical(1))]
    if (length(keggpathid2name.list) == 0) {
      keggpathid2name.df <- NULL
    } else {
      keggpathid2name.df <- do.call(rbind, keggpathid2name.list)
      rownames(keggpathid2name.df) <- NULL
    }
  }

  if (is.null(keggpathid2name.df) || nrow(keggpathid2name.df) == 0) {
    data.frame(path_id = character(0), path_name = character(0), stringsAsFactors = FALSE)
  } else {
    keggpathid2name.df[,2] <- sub("\\s-\\s[a-zA-Z ]+\\(\\w+\\)$", "", keggpathid2name.df[,2])
    colnames(keggpathid2name.df) <- c("path_id","path_name")
    keggpathid2name.df
  }
}


##' @importFrom magrittr %<>%
download.organisms.KEGG <- function(organism) {
  keggpathid2extid.df <- kegg_link(organism, "pathway")

  if (is.null(keggpathid2extid.df) || nrow(keggpathid2extid.df) == 0) {
    write(paste(Sys.time(), "Pathway data of", organism, "is null."), stderr())
    return(NULL)
  }

  message(paste0(Sys.time(), " Getting KEGG data of ", organism, "."))
  keggpathid2extid.df[,1] %<>% gsub("[^:]+:", "", .)
  keggpathid2extid.df[,2] %<>% gsub("[^:]+:", "", .)
  colnames(keggpathid2extid.df) <- c("pathway_id","gene_or_orf_id")
  message(paste(Sys.time(), "KEGG data of", organism, "has been downloaded."))
  keggpathid2extid.df
}


get_organisms_list <- function(db){
  organisms <- kegg_list(db)
  if (is.null(organisms) || nrow(organisms) == 0) return(character(0))
  as.character(organisms[,2])
}


##' @importFrom RSQLite dbDriver
##' @importFrom RSQLite dbConnect
##' @importFrom RSQLite dbWriteTable
##' @importFrom RSQLite dbDisconnect
prepare_kegg_db <- function(organisms, sqlite_path) {
  dbfile <- file.path(sqlite_path, "KEGG.sqlite")
  unlink(dbfile)
  ###################################################
  ### create database
  ###################################################
  drv <- dbDriver("SQLite")
  db <- dbConnect(drv, dbname=dbfile)

  if (length(organisms) >= 1 && any(organisms == "all")) {
    organisms <- get_organisms_list("organism")
  }

  use_map_pathway_names <- length(organisms) > 50

  dbWriteTable(
    conn = db,
    name = "pathway2gene",
    value = data.frame(pathway_id = character(0), gene_or_orf_id = character(0), stringsAsFactors = FALSE),
    row.names = FALSE
  )

  pathway2gene_count <- 0L
  pathway2name_count <- 0L

  map_names <- NULL
  if (use_map_pathway_names) {
    dbWriteTable(
      conn = db,
      name = "pathway2name",
      value = data.frame(path_id = character(0), path_name = character(0), stringsAsFactors = FALSE),
      row.names = FALSE
    )

    kegg_map_pathway <- kegg_list("pathway")
    if (!is.null(kegg_map_pathway) && nrow(kegg_map_pathway) > 0) {
      kegg_map_pathway[, 2] <- sub("\\s-\\s[a-zA-Z ]+\\(\\w+\\)$", "", kegg_map_pathway[, 2])
      map_ids <- gsub("[^:]+:", "", kegg_map_pathway[, 1])
      map_digits <- sub("^[^0-9]+", "", map_ids)
      map_names <- as.character(kegg_map_pathway[, 2])
      names(map_names) <- map_digits
    }
  } else {
    KEGGPATHID2NAME <- get_path2name(organisms)
    dbWriteTable(conn = db, "pathway2name", KEGGPATHID2NAME, row.names = FALSE)
    pathway2name_count <- nrow(KEGGPATHID2NAME)
  }

  for (organism in organisms) {
    KEGGPATHID2EXTID <- download.organisms.KEGG(organism)
    if (!is.null(KEGGPATHID2EXTID)) {
      dbWriteTable(conn = db, "pathway2gene", KEGGPATHID2EXTID, row.names = FALSE, append = TRUE)
      pathway2gene_count <- pathway2gene_count + nrow(KEGGPATHID2EXTID)
      message(paste(Sys.time(), "KEGG data of", organism, "has been added to the sqlite database."))

      if (use_map_pathway_names && !is.null(map_names)) {
        path_ids <- unique(KEGGPATHID2EXTID[, 1])
        path_digits <- sub("^[^0-9]+", "", path_ids)
        KEGGPATHID2NAME_org <- data.frame(
          path_id = path_ids,
          path_name = unname(map_names[path_digits]),
          stringsAsFactors = FALSE
        )
        dbWriteTable(conn = db, "pathway2name", KEGGPATHID2NAME_org, row.names = FALSE, append = TRUE)
        pathway2name_count <- pathway2name_count + nrow(KEGGPATHID2NAME_org)
      }
    }
  }
  
  ###################################################
  ### append the metadata
  ###################################################
  metadata <- rbind(c("PATHNAMESOURCENAME", "KEGG PATHWAY"),
                    c("PATHNAMESOURCEURL", "ftp://ftp.genome.jp/pub/kegg/pathway"),
                    c("PATHNAMESOURCEDATE", format(Sys.Date(), "%Y%m%d")),
                    c("KEGGSOURCENAME", "KEGG GENOME"),
                    c("KEGGSOURCEURL", "ftp://ftp.genome.jp/pub/kegg/genomes"),
                    c("KEGGSOURCEDATE", format(Sys.Date(), "%Y%m%d")),
                    c("GOEXTSOURCEDATE", "2015-Sepec2go27"),
                    c("GOEXTSOURCENAME", "Gene Ontology External Link"),
                    c("GOEXTSOURCEURL", "http://www.geneontology.org/external2go"),
                    c("Db type", "KEGGDB"),
                    c("DBSCHEMA", "KEGG_DB"),
                    c("DBSCHEMAVERSION", "2.1"))
  
  metadata <- as.data.frame(metadata)
  colnames(metadata) <- c("name", "value") #makeAnnDbPkg规定的
  dbWriteTable(conn = db, "metadata", metadata, row.names=FALSE)
  
  map.counts <- rbind(c("pathway2name", pathway2name_count),
                      c("pathway2gene", pathway2gene_count))
  map.counts <- as.data.frame(map.counts)
  colnames(map.counts) <- c("map_name","count")
  dbWriteTable(conn = db, "map_counts", map.counts, row.names=FALSE)

  map.metadata <- rbind(c("ENZYMEID2GO","Gene Ontology External Link","http://www.geneontology.org/external2go","2015-Sepec2go27"),
                        c("GO2ENZYMEID","Gene Ontology External Link","http://www.geneontology.org/external2go","2015-Sepec2go27"),
                        c("EXTID2PATHID","KEGG GENOME","ftp://ftp.genome.jp/pub/kegg/genomes","2011-Mar15"),
                        c("PATHID2EXTID","KEGG GENOME","ftp://ftp.genome.jp/pub/kegg/genomes","2011-Mar15"),
                        c("PATHNAME2ID","KEGG PATHWAY","ftp://ftp.genome.jp/pub/kegg/pathway",format(Sys.Date(),"%Y%m%d")),
                        c("PATHID2NAME","KEGG PATHWAY","ftp://ftp.genome.jp/pub/kegg/pathway",format(Sys.Date(),"%Y%m%d")))
  map.metadata <- as.data.frame(map.metadata)
  colnames(map.metadata) <- c("map_name","source_name","source_url","source_date")
  dbWriteTable(conn = db, "map_metadata", map.metadata, row.names=FALSE)
  
  dbDisconnect(db)
  invisible(dbfile)
}


utils::globalVariables(".")

