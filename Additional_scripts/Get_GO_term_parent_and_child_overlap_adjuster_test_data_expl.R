##### tester

library("VennDiagram")
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

### deal with the problem that topographically close go terms from different enrichement analyses will not overlap in a venn diagram
###### example data using response to caffeine

### hierarchy:  
#"GO:0010243" response to organonitrogen compound ## grand parent
# GO:0043279 response to alkaloid ## parent
# GO:0014074 response to purine-containing compound ## parent
# GO:0031000 response to caffeine ### term
# GO:0071313 cellular response to caffeine ## child

### problem set (first term in list is related to all)
### one unrealted term overlaps in all set. 
### Want 2 terms in the middle
set1 = c("GO:0010243","GO:0043279","GO:0002181", "GO:0008039", "GO:0001522","GO:0042689") ## has 2 terms from caffeine
set2 = c("GO:0031000","GO:0002181", "GO:0046692", "GO:0007394","GO:0046845")
set3 = c("GO:0071313","GO:0002181", "GO:0090090", "GO:0032543","GO:0008045")

### plot 
venn.plot.1 <- venn.diagram(
list(set1 = set1, set2 = set2, set3 = set3), filename = NULL,
                            cat.col = c( "black",   "red",     "green3"),
                            fill=c("black",   "red",     "green3"), margin = 0.2)

grid.arrange(gTree(children=venn.plot.1),ncol = 1) ### only one term in the middle

#### run python Get_GO_term_parent_and_child_overlap_adjuster.py with test data:
## python Get_GO_term_parent_and_child_overlap_adjuster.py  -g /Users/dparker/Documents/Gen_BioInf/obo_GO_file/go.obo -T

## read in output
setwd("/Users/dparker/Desktop")
set1z = as.data.frame(read.table("GO_overlap_adj_set1.txt", header = T))
set2z = as.data.frame(read.table("GO_overlap_adj_set2.txt", header = T))
set3z = as.data.frame(read.table("GO_overlap_adj_set3.txt", header = T))



venn.plot.1z <- venn.diagram(
list(set1z = set1z$set1, set2a = set2z$set2, set3a = set3z$set3), filename = NULL,
                            cat.col = c( "black",   "red",     "green3"),
                            fill=c("black",   "red",     "green3"), margin = 0.2)

grid.arrange(gTree(children=venn.plot.1),gTree(children=venn.plot.1z),ncol = 2) ### to in the middle - Go term numbers are the same





