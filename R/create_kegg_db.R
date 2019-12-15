##' create KEGG.db package
##'
##' 
##' @title create_kegg_db
##' @param species one of KEGG supported species, e.g. hsa for human
##' @return KEGG.db package generated in working directory
##' @export
##' @author Guangchuang Yu and Ziru Chen
create_kegg_db <- function(species) {
    packagedir <- tempfile() # tempdir() maynot empty
    dir.create(packagedir)

    sqlite_path <- paste(packagedir, "inst", "extdata", sep=.Platform$file.sep)
    if(!dir.exists(sqlite_path)){
        dir.create(sqlite_path,recursive = TRUE)
    }

    R_src <- paste(packagedir, "R", sep=.Platform$file.sep)
    if(!dir.exists(R_src)){
        dir.create(R_src,recursive = TRUE)
    }

    ## file.copy(from = system.file("KEGG.db", "R", "zzz.R", package = "createKEGGdb"),
    ##           to = paste(R_src, "zzz.R", sep = .Platform$file.sep))

    fcp("R", todir = R_src, file = "zzz.R")
    fcp(todir = packagedir, file = "DESCRIPTION")
    fcp(todir = packagedir, file = "LICENSE")
    fcp(todir = packagedir, file = "NAMESPACE")

    prepare_kegg_db(species, sqlite_path)

    pkgbuild::build(packagedir, dest_path = ".")
}


fcp <- function(..., todir, file) {
    file.copy(from = system.file("KEGG.db", ..., file, package = "createKEGGdb"),
              to = paste(todir, file, sep = .Platform$file.sep))
}


##' @importFrom magrittr %<>%
get_path2name <- function(){
  keggpathid2name.df <- clusterProfiler:::kegg_list("pathway")
  keggpathid2name.df[,1] %<>% gsub("path:map", "", .)
  colnames(keggpathid2name.df) <- c("path_id","path_name")
  return(keggpathid2name.df)
}


##' @importFrom magrittr %<>%
download.organisms.KEGG <- function(organism) {
  keggpathid2extid.df <- clusterProfiler:::kegg_link(organism, "pathway")
  if (is.null(keggpathid2extid.df)){
    write(paste(Sys.time(),"Pathway data of",organism,"is null."), stderr())
  }else{
    message(paste0(Sys.time()," Getting KEGG data of ",organism,"."))
    keggpathid2extid.df[,1] %<>% gsub("[^:]+:", "", .)
    keggpathid2extid.df[,2] %<>% gsub("[^:]+:", "", .)
    colnames(keggpathid2extid.df) <- c("pathway_id","gene_or_orf_id")
    message(paste(Sys.time(),"KEGG data of",organism,"has been downloaded."))
  }
  return(keggpathid2extid.df)
}


get_organisms_list <- function(db){
  organisms <- clusterProfiler:::kegg_list(db)
  organisms_list <- as.character(organisms[,2])
  return(organisms_list)
}


##' @importFrom clusterProfiler download_KEGG
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

  KEGGPATHID2NAME <- get_path2name()
  ###################################################
  ### put the pathway2name data into the tables 
  ###################################################
  dbWriteTable(conn = db, "pathway2name", KEGGPATHID2NAME, row.names=FALSE)
  
  if (length(organisms) == 1){
    if(organisms == "all"){
      organisms <- get_organisms_list("organism")
    }
  }
  for(organism in organisms){
    KEGGPATHID2EXTID <- download.organisms.KEGG(organism)
    if(!is.null(KEGGPATHID2EXTID)){
      ###################################################
      ### put the pathway2gene data into the tables 
      ###################################################
      # 数据是直接添加进去，不会自动去重
      dbWriteTable(conn = db, "pathway2gene", KEGGPATHID2EXTID, row.names=FALSE,append = TRUE)
      message(paste(Sys.time(),"KEGG data of",organism,"has been added to the sqlite database."))
      
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
  
  map.counts <- rbind(c("pathway2name", nrow(KEGGPATHID2NAME)),
                      c("pathway2gene", nrow(KEGGPATHID2EXTID)))
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
