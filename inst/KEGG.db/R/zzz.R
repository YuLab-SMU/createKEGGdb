datacache <- new.env(hash=TRUE, parent=emptyenv())

KEGG <- function() showQCData("KEGG", datacache)
KEGG_dbconn <- function() dbconn(datacache)
KEGG_dbfile <- function() dbfile(datacache)
KEGG_dbschema <- function(file="", show.indices=FALSE) dbschema(datacache, file=file, show.indices=show.indices)
KEGG_dbInfo <- function() dbInfo(datacache)

.onLoad <- function(libname, pkgname)
{
	    ## Connect to the SQLite DB
	    dbfile <- system.file("extdata", "KEGG.sqlite", package=pkgname, lib.loc=libname)
    assign("dbfile", dbfile, envir=datacache)
        dbconn <- dbFileConnect(dbfile)
        assign("dbconn", dbconn, envir=datacache)
	    ## Create the AnnObj instances
	    ann_objs <- createAnnObjs.SchemaChoice("KEGG_DB", "KEGG", "KEGG", dbconn, datacache)
	    mergeToNamespaceAndExport(ann_objs, pkgname)
	        packageStartupMessage(AnnotationDbi:::annoStartupMessages("KEGG.db"))
}

.onUnload <- function(libpath)
{
	    dbFileDisconnect(KEGG_dbconn())
}
