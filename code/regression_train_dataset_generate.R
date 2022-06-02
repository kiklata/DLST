#regression train dataset generate

marker = c('PTPRC','CD45','CD3D','CD3E','CD4','CD8A','CD8B','CD19','CD79A',
'MS4A1','IGHG1','MZB1','SDC1','CD68','CD163','CD14','FGFBP2',
'FCG3RA','CX3CR1','KLRB1','NCR1','TPSAB1','TPSB2','C1QA','C1QB',
'S100A9','S100A8','MMP19','LAMP3','IDO1','IDO2','CCL19','CCR7',
'CD1E','CD1C','CLEC9A', 'CD141')

library(dplyr)
load("A_norm_exp.Rdata")
A_gene = colnames(A_norm_exp)
load("B_norm_exp.Rdata")
B_gene = colnames(B_norm_exp)
load("C_norm_exp.Rdata")
C_gene = colnames(C_norm_exp)
load("D_norm_exp.Rdata")
D_gene = colnames(D_norm_exp)
load("E_norm_exp.Rdata")
E_gene = colnames(E_norm_exp)
load("F_norm_exp.Rdata")
F_gene = colnames(F_norm_exp)
load("G_norm_exp.Rdata")
G_gene = colnames(G_norm_exp)
load("H_norm_exp.Rdata")
H_gene = colnames(H_norm_exp)

intersect_Gene = Reduce(intersect,list(A_gene,B_gene,C_gene,D_gene,E_gene,F_gene,G_gene,H_gene))

seq_marker = marker[marker %in% intersect_Gene]

A_norm_marker = select(A_norm_exp,seq_marker) 
B_norm_marker = select(B_norm_exp,seq_marker) 
C_norm_marker = select(C_norm_exp,seq_marker) 
D_norm_marker = select(D_norm_exp,seq_marker) 
E_norm_marker = select(E_norm_exp,seq_marker) 
F_norm_marker = select(F_norm_exp,seq_marker) 
G_norm_marker = select(G_norm_exp,seq_marker) 
H_norm_marker = select(H_norm_exp,seq_marker) 

intersect_Gene = Reduce(intersect,list(A_gene,B_gene,C_gene,D_gene,E_gene,F_gene,G_gene,H_gene))

immuno_marker_exp = rbind(A_norm_marker,B_norm_marker,C_norm_marker,D_norm_marker,
      E_norm_marker,F_norm_marker,G_norm_marker,H_norm_marker)
immuno_marker_exp$ID = rownames(immuno_marker_exp)
save(immuno_marker_exp,file = 'immuno_marker.Rdata')

