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

    ## skeleton
    prepare_pkg_skeleton(packagedir)

    ## sqlite
    sqlite_path <- paste(packagedir, "inst", "extdata", sep=.Platform$file.sep)
    prepare_kegg_db(species, sqlite_path)

    ## build pkg
    pkgbuild::build(packagedir, dest_path = ".")
}


prepare_pkg_skeleton <- function(packagedir) {
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
}


##' @importFrom clusterProfiler download_KEGG
##' @importFrom RSQLite dbDriver
##' @importFrom RSQLite dbConnect
##' @importFrom RSQLite dbWriteTable
##' @importFrom RSQLite dbDisconnect
prepare_kegg_db <- function(species, sqlite_path) {
    dbfile <- file.path(sqlite_path, "KEGG.sqlite")
    unlink(dbfile)

###################################################
### download data 
###################################################
    ## KEGG species list
    ## https://www.genome.jp/kegg/catalog/org_list.html
    ## 
    ## http://rest.kegg.jp/list/organism

    kegg <- clusterProfiler::download_KEGG(species)
    KEGGPATHID2NAME <- kegg$KEGGPATHID2NAME
    colnames(KEGGPATHID2NAME) <- c("path_id", "path_name")
    KEGGPATHID2NAME$path_id <- sub(species, "", KEGGPATHID2NAME$path_id)

    KEGGPATHID2EXTID <- kegg$KEGGPATHID2EXTID
    colnames(KEGGPATHID2EXTID) <- c("pathway_id","gene_or_orf_id")

###################################################
### create database
###################################################
    ## Create the database file

    drv <- dbDriver("SQLite")
    db <- dbConnect(drv, dbname=dbfile)
    ## dbDisconnect(db)

###################################################
### put the data into the tables 
###################################################
    dbWriteTable(conn = db, "pathway2name", KEGGPATHID2NAME, row.names=FALSE)
    dbWriteTable(conn = db, "pathway2gene", KEGGPATHID2EXTID, row.names=FALSE)
    ## dbListTables(db)
    ## test <- dbReadTable(conn = db,"pathway2name")
    ## head(test)
    ## test <- dbReadTable(conn = db,"pathway2gene")
    ## head(test)

###################################################
### append the metadata
###################################################
#dbSendQuery(conn = db,"drop table if exists metadata")
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
    
    ## dbListTables(db)
    ## dbListFields(conn = db, "metadata")
    ## dbReadTable(conn = db,"metadata")
    
    
    map.metadata <- rbind(c("ENZYMEID2GO","Gene Ontology External Link","http://www.geneontology.org/external2go","2015-Sepec2go27"),
                          c("GO2ENZYMEID","Gene Ontology External Link","http://www.geneontology.org/external2go","2015-Sepec2go27"),
                          c("EXTID2PATHID","KEGG GENOME","ftp://ftp.genome.jp/pub/kegg/genomes","2011-Mar15"),
                          c("PATHID2EXTID","KEGG GENOME","ftp://ftp.genome.jp/pub/kegg/genomes","2011-Mar15"),
                          c("PATHNAME2ID","KEGG PATHWAY","ftp://ftp.genome.jp/pub/kegg/pathway",format(Sys.Date(),"%Y%m%d")),
                          c("PATHID2NAME","KEGG PATHWAY","ftp://ftp.genome.jp/pub/kegg/pathway",format(Sys.Date(),"%Y%m%d")))
    map.metadata <- as.data.frame(map.metadata)
    colnames(map.metadata) <- c("map_name","source_name","source_url","source_date")
    dbWriteTable(conn = db, "map_metadata", map.metadata, row.names=FALSE)
    
    
    ## dbListTables(db)
    ## dbListFields(conn = db, "map_metadata")
    ## dbReadTable(conn = db,"map_metadata")
    dbDisconnect(db)
    invisible(dbfile)
}
