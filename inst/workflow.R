library("DiagrammeR")


grViz("digraph createKEGGdb {
rankdir = TD
node [shape = box, style=filled]
layout = dot
compound =true
#color = crimson

download [label='Query KEGG pathway for specific species']
database [label='Pack KEGG data into a sqlite file']
pkg0 [label='KEGG.db package skeleton']
pkg [label='Build KEGG.db package']

download -> database -> pkg
pkg0 -> pkg

}") -> x

yyplot::gv2file(x, file = 'diagram.png' )
