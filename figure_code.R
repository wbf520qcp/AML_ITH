########################################################################
#This script contains code for analysis and figures for heterogenity paper
########################################################################
library(Seurat);library(ComplexHeatmap);library(symphony);library(pheatmap)
library(dplyr);library(tidydr);library(fgsea);library(enrichR)
library(rstatix)  #for wilcox test
library(cowplot)  #for theme_cowplot
library(Rmagic)   #for magic to impute missing values
library(ggcorrplot);library("ggvenn")
library("Hmisc") #calculate p valye of correlation matrix
library(edgeR)   #for cpm conversion
library(ggpubr)  #for ggboxplot
library(smoother) #for gaussian smoothing in ordered heatmap
library(viridis) #for crispr visualization
library(RColorBrewer)  #for brewer.pal color
library(ggrastr)  #for geom_point_rast
library(circlize)  #for colorRamp2
library(riverplot);library(tidyverse)
library(cluster);library(factoextra);library(EnhancedVolcano);library(data.table)
library(UpSetR)  #for upset plot
library(hypeR)  # for geneset enrichment (hypergenomic)
library(patchwork)  #for plot_layout in pathway visualization
library(msigdbr)
library(readxl)
library(gplots)
library(SCEVAN)
library(survivalAnalysis)  #for multivariate cox model
library(survminer);library(survival)
library(destiny)  #for diffusion map
library(GSVA)
library(gtsummary)  #for clinical table
library("easier") #for EaSleR to predict immune response
library(SCP)
library(viper)
library(data.table)
library(factoextra) #fviz_nbclust:Dertermining and Visualizing the Optimal Number of Clusters 
library(NMF)
library(glmnet) #for LASSO regression
library(pROC) #these 5 packages are used for figure 5 ML to select genes for signatures
library(caret)
library(e1071)
library(glmnet)
library(lars)
library("faux")
library("DataExplorer")
library("caret")
library("randomForest");library("VennDiagram")
library(corrplot)

#------------------------------------------------------------------------------
############################################################################
#Figure 1: study population; Project symphony reference, lymphoid, all AML, 
#          deg between cell types, PCA plot, SCENIC
############################################################################
#------------------------------------------------------------------------------
#############################################
#Study design is plotted in biorender
############################################
#oncoprint for study population
metadata=read_excel("Research/AML/info_files/patients_info/Copy of Metadata_Hemepath 4.14.2022.xlsx",sheet = "Clean_metadata")
metadata1=read_excel("Research/AML/info_files/patients_info/Copy of Metadata_Hemepath 4.14.2022.xlsx",sheet = "Clinical")
metadata=as.data.frame(metadata[metadata$Bofei_AML_cellAnalysis=="Include in Analysis",])
metadata=metadata[startsWith(metadata$SampleID,"PT"),]
plot.data=metadata[,c(1,25,29:32)]
plot.data=merge(metadata1[,c("SampleID","AML Status at BL")],plot.data,by="SampleID")
rownames(plot.data)=plot.data$SampleID
plot.data[plot.data$`AML Status at BL`!="Denovo",2]="Secondary"
plot.data=plot.data[,2:ncol(plot.data)];plot.data=t(plot.data)
rownames(plot.data)=c("Status","Karyotype","Del5/5q","Del7/7q","Trisomy8","Del17/17p")#,"EZH2","TP53","FLT3","NPM1","RAS","SF3B1","IDH","DTA"

######process the MDL_Pos_Mutations column to plot oncoprint; mutations in separate columns are not correct
mut=as.data.frame(str_split_fixed(metadata$MDL_Pos_Mutations, ',',Inf))
mut=as.data.frame(add_column(mut, Patient = metadata$SampleID, .before = 1))
mut[18,3]="TP53"
mut_new=mut[,1:2]
for (i in colnames(mut)[3:ncol(mut)]) {
  temp=mut[,c("Patient",i)]
  names(temp)=names(mut_new)
  mut_new=rbind(mut_new,temp)
}
mut_new=reshape2::dcast(mut_new,Patient~V1)
rownames(mut_new)=mut_new$Patient
mut_new=mut_new[c(20,1:19),-c(1,2)]
plot.data=rbind(as.data.frame(plot.data),as.data.frame(t(mut_new)))




lgd1 = Legend(labels = c("Mutation", "No Mutation", "Not Tested"), legend_gp = gpar(fill = c("brown4", "lightgray", "snow4")), title = "Mutation")
lgd2 = Legend(labels=c("Diploid", "Complex", "Not Complex"), legend_gp = gpar(fill=c("rosybrown","palevioletred", "lightpink")), title="Karyotype")
lgd3 = Legend(labels=c("Del5/5q", "Del7/7q", "Trisomy8","Del17/17p"), legend_gp=gpar(fill=c("peachpuff1", "steelblue3",  "indianred2","tan2")), title= "Cytogenetic abnormalities")
lgd4= Legend(labels=c("Secondary", "Denovo"), legend_gp = gpar(fill=c("wheat3", "cornflowerblue")), title="AML Status")
pd = packLegend(lgd4,lgd2, lgd3, lgd1,  direction = "horizontal", gap = unit("2", "cm"))

png("Research/manuscript/cell_of_origin/Figure1/study_cohort.png", res=200, width=10, height=7.5,units="in")
Heatmap(as.matrix(plot.data), col = c("1" = "brown4", "0"="lightgray", "NP" = "snow4", "Diploid"="rosybrown",
                        "Complex"="palevioletred", "Not Complex"="lightpink","del7"="steelblue3","del7q"="steelblue3",
                           "del5q"="peachpuff1", "Trisomy8"="indianred2", "del17_17p"="tan2", "No"="lightgray",
                           "Secondary"= "wheat3", "Denovo"="cornflowerblue"), 
        column_title = "Patients",column_title_gp = gpar(fontsize = 18), row_title=" ", 
        row_split=c(rep("A", 6), rep("B",27)),border=FALSE, width=unit(7, "in"),
        height=unit(4,"in"), rect_gp = gpar(col = "white", lwd = 2), show_heatmap_legend = FALSE)
draw(pd, just=c("center","bottom"), y=unit(0.8, "cm"))
dev.off()

rm(plot.data,metadata,lgd1,lgd2,lgd3,lgd4,i,temp,mut_new,mut,metadata.metadata1)


#------------------------------------------------------------------------------
#AML cell fraction by cytogenetic group
#------------------------------------------------------------------------------
allcell=readRDS("/Users/bwang8/Research/AML/AML_subsets/allcell_new_noharmony_update.RDS")
allcell$cyto3=allcell$cyto4
allcell$cyto3[allcell$cyto3 %in% c("Diploid-mono","Diploid-Nonmono")]='Diploid'

meta=allcell@meta.data[,c("orig.ident","cyto3","class2")]
fraction=meta %>% group_by(orig.ident,class2)  %>% summarise(count=n()) 
fraction=fraction %>% group_by(orig.ident) %>% mutate(Percent = round(count/sum(count),4))
fraction=as.data.frame(fraction); fraction=fraction[fraction$class2=="AML",]
fraction=merge(fraction,distinct(meta[,c('orig.ident',"cyto3")],orig.ident,.keep_all = TRUE),by='orig.ident')
fraction$cyto3=factor(fraction$cyto3,levels = c("Diploid","Del5.5q","Del7.7q","Double del"))

p=ggplot(fraction,aes(x=cyto3,y=Percent)) +
  geom_boxplot(aes(fill=cyto3),width=0.5,color='gray40',outlier.shape=NA,size=0.5)+
  geom_jitter(fill='gray80',color='black',shape=21,width = 0.2,size=1,stroke = 0.5)+
  ylab("Fraction of AML cells") + theme_classic() + 
  scale_fill_manual(values=c("indianred2", "darkseagreen3", "gold2","deepskyblue"))+
  theme(legend.position="none",
        axis.title.y = element_text(size = 12),axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12,hjust = 0.5),axis.text.y = element_text(size = 12))

kruskal.test(fraction$Percent~fraction$cyto3)
p + annotate(geom="text",x=2.5,y=0.9,label="kruskal-Wallis, p = 0.5416",cex=4)
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_Supp/supp_AMLfraction_cyto3_box.pdf",height=12,width=10, units="cm")
rm(p,fraction,meta)


############################################
#oncoprint of patient cohort
metadata <- read_excel("Research/AML/info_files/patients_info/Copy of Metadata_Hemepath 4.14.2022.xlsx",sheet = "Clean_metadata")
clinical=as.data.frame(read_excel("Research/AML/info_files/patients_info/Copy of Metadata_Hemepath 4.14.2022.xlsx",sheet = "Clinical"))
metadata=as.data.frame(metadata[metadata$Bofei_AML_cellAnalysis=="Include in Analysis",])
metadata=metadata[metadata$Disease=="AML",]
plot.data <- metadata[!(metadata$SampleID %like% "^N"),c(1,25,29,30,32,36:41,43,44,46:48,50)]
plot.data$AML_status=clinical$`AML Status at BL`
plot.data$Gender=clinical$Gender
plot.data$AML_status=ifelse(plot.data$AML_status=="Denovo",'Denovo','Secondary')
rownames(plot.data)=plot.data$SampleID
plot.data=plot.data[,2:ncol(plot.data)];plot.data=as.data.frame(t(plot.data))
rownames(plot.data)=c("Karyotype","Del5/5q","Del7/7q","Del17/17p","EZH2","TP53","FLT3","NPM1","RAS","SF3B1","IDH1","IDH2","DNMT3A","TET2","ASXL1","ELN2017","Status","Gender")

mut=as.data.frame(str_split_fixed(metadata$MDL_Pos_Mutations, ',',Inf))
mut=as.data.frame(add_column(mut, Patient = metadata$SampleID, .before = 1))
mut[18,3]="TP53"
mut_new=mut[,1:2]
for (i in colnames(mut)[3:ncol(mut)]) {
  temp=mut[,c("Patient",i)]
  names(temp)=names(mut_new)
  mut_new=rbind(mut_new,temp)
}
mut_new=reshape2::dcast(mut_new,Patient~V1)
rownames(mut_new)=mut_new$Patient
mut_new=mut_new[c(20,1:19),-c(1,2,23)]##subset(mut_new, select=-c(None,Patient))
plot.data=rbind(plot.data[c("Gender","Karyotype","Del5/5q","Del7/7q","Del17/17p"),],as.data.frame(t(mut_new)),plot.data[c("ELN2017","Status"),])




lgd1 = Legend(labels = c("Mutation", "No Mutation"), legend_gp = gpar(fill = c("purple", "lightgray")), title = "Mutation")
lgd2 = Legend(labels=c("Diploid", "Complex", "Not Complex"), legend_gp = gpar(fill=c("indianred2", "gold2", "lightgreen")), title="Karyotype")
lgd3 = Legend(labels=c("Alternation", "No alternation"), legend_gp=gpar(fill=c("deepskyblue","lightgray")), title= "Chromosome Alternations")
lgd4= Legend(labels=c("Favorable", "Intermediate", "Adverse"), legend_gp = gpar(fill=c("darkseagreen3","khaki", "paleturquoise")), title="ELN 2017")
lgd5= Legend(labels=c("Secondary", "Denovo"), legend_gp = gpar(fill=c("orange", "seashell")), title="AML Status")
lgd6 = Legend(labels = c("Male", "Female"), legend_gp = gpar(fill = c("deepskyblue", "coral")), title = "Gender")
pd = packLegend(lgd2, lgd3, lgd1,lgd4,lgd5,lgd6, direction = "horizontal", gap = unit("1", "cm"))


pdf("Research/manuscript/intratumor_heterogeneity/ITH_Supp/oncoprint_cohort1.pdf", width=8, height=12)
Heatmap(as.matrix(plot.data), col = c("Diploid"="indianred2","Complex"="gold2", "Not Complex"="lightgreen",
             "del7"="steelblue3","del7q"="steelblue3","del5q"="steelblue3", "del17_17p"="steelblue3","No"="lightgray",
                           "1" = "purple", "0"="lightgray", 
             "Male"="deepskyblue","Female"="coral",
                           "Favorable"="darkseagreen3", "Intermediate"="khaki", "Adverse"="paleturquoise", 
                           "Secondary"= "orange", "Denovo"="seashell1"), 
        row_title=" ", row_split=c(rep("A", 5), rep("B",26), rep("C",2)),border=FALSE, width=unit(5, "in"),
        height=unit(6,"in"), rect_gp = gpar(col = "black", lwd = 1), show_heatmap_legend = FALSE)
draw(pd, just=c("center","bottom"), y=unit(4, "cm"))
dev.off()
rm(plot.data,lgd1,lgd2,lgd3,lgd4,lgd5,lgd6,mat,pd,clinical,i,mut,mut_new,temp,metadata)








############################################
#diffusion map
sce <- as.SingleCellExperiment(amlnew)
aml.dm=DiffusionMap(sce,k=5,verbose=T,n_pcs=20,distance = "euclidean")
dpt.aml <- DPT(aml.dm)
amlnew@meta.data$dpt_pca20 <- rank(dpt.aml$dpt)
df <- data.frame(DC1 = eigenvectors(aml.dm)[, 1], DC2 = eigenvectors(aml.dm)[, 2], 
                 DC3 = eigenvectors(aml.dm)[, 3], DC4 = eigenvectors(aml.dm)[, 4],
                 cell_type = amlnew$class2,cytogenetic=amlnew$cyto3,coty2=amlnew$cyto2)
ggplot(df, aes(x = DC1, y = DC2, colour = cytogenetic)) + geom_point(size=0.25)  + 
  xlab("Diffusion component 1") + ylab("Diffusion component 2") +theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=6))) + 
  theme(legend.text=element_text(size=12), legend.title=element_text(size=16),
        axis.title=element_text(size=14,face="bold"),axis.text=element_text(size=14)) +
  scale_colour_manual(values = c("indianred2", "darkseagreen3", "gold2","deepskyblue"))
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure1/diffusion_cyto.pdf",height=10,width=15, units="cm")


#------------------------------------------------------------------------------
# UMAP of AML cells by cytogenetics
cyto3.colors <- c("Diploid"="indianred2","Del5.5q"="chartreuse1","Del7.7q"="gold2","Double deletion"="deepskyblue")

DimPlot(amlnew,group.by="cyto3",cols=cyto3.colors,raster=F)+
  NoAxes()+ggtitle("")+
  theme(legend.title=element_text(size = 0),
        legend.text=element_text(size = 8),
        legend.key.height = unit(0.25,"cm"),
        legend.key.width = unit(0.25,"cm"))+NoAxes()
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure1/umap_cyto.pdf",height=15,width=18, units="cm")




#------------------------------------------------------------------------------
#DEG of cytogenetic groups 
#------------------------------------------------------------------------------
#heatmap
gene.heat=c()
for (i in c('diploid',"del5","del7","double")){
  deg=read.csv(paste0("Research/manuscript/intratumor_heterogeneity/ITH_figure2/DEG_cyto/",i,".csv"))
  deg=deg[!(grepl(pattern="^(MT|AJ|AL|AC|LIN|AP)", deg$X)),]
  deg=deg[(deg$p_val_adj<0.05 & deg$pct.1>0.1),]
  deg=deg[order(deg$avg_log2FC,decreasing = T),]
  gene.heat=c(gene.heat,deg$X[1:10])
}

mat <- amlnew@assays$RNA@data[gene.heat,]
group.list=amlnew@meta.data
group.list=data.frame(cyto=group.list$cyto3,id=rownames(group.list),patient=group.list$orig.ident)
group.list=group.list[order(match(group.list$cyto,c("Diploid","Del5.5q","Del7.7q","Double deletion"))),]
mat=as.matrix(mat[,group.list$id])
mat=t(scale(t(mat)))
column_ha=HeatmapAnnotation(df=data.frame(
  cyto=factor(group.list[,1],levels = c("Diploid","Del5.5q","Del7.7q","Double deletion"))),
  annotation_name_side = "left",
  annotation_name_gp= gpar(fontsize = 10),show_legend = c(TRUE),
  col = list(cyto=c("Diploid"="indianred2", "Del5.5q"="darkseagreen3", "Del7.7q"="gold2","Double deletion"="deepskyblue")))

pdf("Research/manuscript/intratumor_heterogeneity/ITH_figure1/cyto3_deg_heatmap.pdf",width=10, height=10)
p=Heatmap(mat,name="Expression",col=colorRamp2(c(-2, 0, 2), c("#2166AC","white","#B2182B")),
          cluster_columns=FALSE, cluster_rows = FALSE, column_names_rot = 45, top_annotation = column_ha,
          column_title=NULL,show_heatmap_legend = T,row_names_side = "left", show_column_names = FALSE,
          row_split = c(rep("1",10),rep("2",10),rep("3",10),rep("4",10)),row_gap = unit(1, "mm"),
          row_names_gp = gpar(fontsize = 9),column_names_gp = gpar(fontsize = 12))
draw(p, annotation_legend_side = 'right',heatmap_legend_side="right")
dev.off()
rm(i,p,mat,group.list,deg,gene.heat,column_ha)




#############################################################################
#run cytosig on NIH web to predict cytokines in each CG group 
#############################################################################
#prepare input file
human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
BM_labels <- getBM(attributes = c("hgnc_symbol","ensembl_gene_id","entrezgene_id", 'chromosome_name', "start_position","end_position","transcript_start","transcript_end"),filters = "hgnc_symbol",values = VariableFeatures(amlnew),mart = human)
BM_labels$length=BM_labels$transcript_end-BM_labels$transcript_start
BM_labels=BM_labels[order(BM_labels$length,decreasing = T),]
BM_labels=BM_labels %>% distinct(hgnc_symbol,.keep_all = TRUE)
gene.length=BM_labels[order(match(BM_labels$hgnc_symbol,VariableFeatures(amlnew))),"length"]

expri=amlnew@assays$RNA@counts
expri=expri[VariableFeatures(amlnew),]#expri=expri[BM_labels$hgnc_symbol,]

rpk <- apply(expri, 2, function(x) x/(gene.length/1000))
tpm <- apply(rpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)
tpm=log2(tpm/10 +1)
gm=groupMeans(tpm,amlnew$cyto4); gm=t(scale(t(gm)))
write.csv(gm,file="Research/manuscript/cell_of_origin/paper/Figure3/cytosig/CytoSig_aml5Cyto.csv")
rm(rpk,tpm,gm,expri,gene.length,human,BM_labels)
####-----------------------------------------------------------------------------
#cytosig is run on NIH website; use the result to plot heatmap
cytokine=read_excel("Research/manuscript/intratumor_heterogeneity/ITH_figure1/cytosig/heatmap_logTPM.xlsx")
cytokine=as.data.frame(cytokine)
rownames(cytokine)=cytokine$ID;cytokine=cytokine[,2:ncol(cytokine)]
cytokine=as.data.frame(t(cytokine))

p=Heatmap(scale(cytokine[c("Diploid","Del5.5q","Del7.7q","Double deletion"),]),
          col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
          cluster_columns=TRUE, cluster_rows = FALSE, column_names_rot = 90, column_title=NULL, 
          show_heatmap_legend = T,row_names_side = "left", show_column_names = TRUE,
          heatmap_legend_param = list(title = ""), width = ncol(cytokine)*unit(5, "mm"), 
          height = nrow(cytokine)*unit(5, "mm"),row_names_gp = gpar(fontsize = 15),column_names_gp = gpar(fontsize = 11))
pdf("Research/manuscript/intratumor_heterogeneity/ITH_figure1/cytosig_AML_cyto3_logTPM.pdf",width=11.5, height=2.5)
draw(p, annotation_legend_side = 'top',heatmap_legend_side="right")
dev.off()








#############################################
#############################################
del7_diploid=read.csv("Research/manuscript/intratumor_heterogeneity/Figure1/DEG/del7_vs_diploid.csv")
del7_diploid=del7_diploid[!(grepl(pattern="^(MT|AJ|AL|AC|LIN)",del7_diploid$X)),]
EnhancedVolcano(del7_diploid,lab = del7_diploid$X,x = 'avg_log2FC',y = 'p_val_adj',xlim = c(-3,3),
                title="", subtitle = "",caption="",axisLabSize=10,pCutoff=10e-15,FCcutoff=0.5,
                pointSize = 0.8,labSize = 3,legendLabSize = 10,legendIconSize = 4)
#drawConnectors = TRUE,widthConnectors = 0.5,max.overlaps=Inf)
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure1/deg_pathway/deg_del7_diploid_volcano.pdf",width=16,height =14,unit="cm")

###run gsea on each cyto group
deg=read.csv("/Users/bwang8/Research/manuscript/intratumor_heterogeneity/ITH_figure2/deg_pathway/diploid.csv")
deg$gene=deg$X
fgseaRes=run_fgsea(group1="diploid",group2 = "",deg,gmt_file  = "/Users/bwang8/Research/AML/info_files/AUCell_geneset/h.all.v7.5.1.symbols.gmt.txt",out.path = "Research/manuscript/intratumor_heterogeneity/ITH_figure2/deg_pathway")

#### combine NES in each group and plot heatmap
gsea.diploid=read.csv("Research/manuscript/intratumor_heterogeneity/ITH_figure1/deg_pathway/fgsea_diploid_plot/gseaHALLMARK_diploid.csv")
gsea.diploid=gsea.diploid[,c("pathway","NES")]
gsea.del5=gsea.del5[,c("pathway","NES")]
gsea.del7=gsea.del7[,c("pathway","NES")]
gsea.double=gsea.double[,c("pathway","NES")]
##contruct the matrix of NES for 4 cytogenetic groups
hallmark=merge(gsea.diploid,gsea.del5,all=TRUE,by="pathway")
hallmark=merge(hallmark,gsea.del7,by="pathway",all=TRUE)
hallmark=merge(hallmark,gsea.double,by="pathway",all=TRUE)
colnames(hallmark)=c('pathway','Diploid',"Del5/5q","Del7/7q","Double deletion")
rownames(hallmark)=hallmark$pathway;hallmark=hallmark[,2:5]
rownames(hallmark)=unlist(lapply(rownames(hallmark),function(x){substr(x,10,nchar(x))}))
hallmark=hallmark[c("HYPOXIA","ALLOGRAFT_REJECTION","ANGIOGENESIS","APICAL_JUNCTION",
                    "APOPTOSIS","UV_RESPONSE_DN","UV_RESPONSE_UP","COMPLEMENT","COAGULATION",
                    "G2M_CHECKPOINT","MYC_TARGETS_V1","MYC_TARGETS_V2","E2F_TARGETS", "MITOTIC_SPINDLE",
                    "INFLAMMATORY_RESPONSE","INTERFERON_ALPHA_RESPONSE","INTERFERON_GAMMA_RESPONSE"            
                    ,"ALLOGRAFT_REJECTION","TNFA_SIGNALING_VIA_NFKB","ESTROGEN_RESPONSE_EARLY", 
                    "KRAS_SIGNALING_UP","MTORC1_SIGNALING","PROTEIN_SECRETION","EPITHELIAL_MESENCHYMAL_TRANSITION"),]
# hallmark.all=hallmark.all[c("APICAL_JUNCTION",
#                              "ANGIOGENESIS","EPITHELIAL_MESENCHYMAL_TRANSITION",
#                              "UV_RESPONSE_DN","UV_RESPONSE_UP",
#                              "ALLOGRAFT_REJECTION","COMPLEMENT", "INFLAMMATORY_RESPONSE","INTERFERON_ALPHA_RESPONSE",                               "INTERFERON_GAMMA_RESPONSE","COAGULATION",
#                              "OXIDATIVE_PHOSPHORYLATION",
#                              "APOPTOSIS","HYPOXIA","PROTEIN_SECRETION",  
#                              "E2F_TARGETS","G2M_CHECKPOINT","MYC_TARGETS_V1","MYC_TARGETS_V2","P53_PATHWAY", 
#                              "MITOTIC_SPINDLE", 
#                              "ESTROGEN_RESPONSE_EARLY","ESTROGEN_RESPONSE_LATE","MTORC1_SIGNALING",
#                              "PI3K_AKT_MTOR_SIGNALING","TGF_BETA_SIGNALING","TNFA_SIGNALING_VIA_NFKB"),]

pdf("Research/manuscript/intratumor_heterogeneity/ITH_figure1/deg_pathway/gsea_hallmark_4group_heatmap.pdf",width=6,height =8)
Heatmap(as.matrix(hallmark),name="NES",col=colorRamp2(c(-2, 0, 2), c("#2166AC","white","#B2182B")),
        cluster_columns=FALSE, cluster_rows = FALSE, #column_names_rot = 90, #width = ncol(gm)*unit(20, "mm"),
        column_title=NULL,show_heatmap_legend = T,row_names_side = "left", show_column_names = TRUE,
        column_split = c(rep("",1),rep(" ",1),rep("  ",1),rep("   ",1)),
        column_gap = unit(1, "mm"),heatmap_legend_param = list(title = ""),
        width = ncol(hallmark)*unit(10, "mm"), height = nrow(hallmark)*unit(5, "mm"),
        row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 11))
dev.off()
rm(hallmark,gsea.del5,gsea.del7,gsea.diploid,gsea.double,deg)


######validate the proliferation pathway enrichment in bulk data
bulk.hallmark=read.csv("Research/manuscript/cell_of_origin/Figure1_new/scoring/TCGA_Abbas_BEAT_ssGSEA_tpm_Bofei.csv")
bulk.meta=read.csv("Research/manuscript/cell_of_origin/paper/Figure1/bulk_3sets_meta_Bofei.csv")
colnames(bulk.hallmark)[1]="Sample"
df=merge(bulk.meta[,c('Sample',"CG_group")],bulk.hallmark[,c("Sample","HALLMARK_G2M_CHECKPOINT","HALLMARK_E2F_TARGETS","HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2")],by="Sample")
df$CG_group[df$CG_group %in% c("Diploid-mono", "Diploid-Nonmono")]="Diploid"
df$CG_group=factor(df$CG_group,levels = c('Diploid',"Del5.5q","Del7.7q","Double del","Other"))

kruskal.test(df$HALLMARK_MYC_TARGETS_V1~df$CG_group)
p=ggplot(df,aes(x=CG_group,y=HALLMARK_MYC_TARGETS_V2)) +
  geom_boxplot(aes(fill=CG_group),width=0.5,color='gray40',outlier.shape=NA,size=0.2)+
  geom_jitter(fill='gray80',color='black',shape=21,width = 0.2,size=1,stroke = 0.2)+
  ylab("HALLMARK MYC TARGETS V2") + theme_classic() + #ylim(0,0.2) +
  scale_fill_manual(values=c("indianred2", "chartreuse1", "gold2","deepskyblue","mediumpurple1"))+
  theme(legend.position="none",
        axis.title.y = element_text(size = 12),axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12,angle = 30,hjust = 1),axis.text.y = element_text(size = 12))

p + annotate(geom="text",x=3,y=0.55,label="kruskal-Wallis, p = 0.012",cex=5)
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure1/bulk_E2F_cyto3_box.pdf",width=12,height =12,unit="cm")




#------------------------------------------------------------------------------
#test all drugs by cytogenetic groups,select significant ones, heatmap
#------------------------------------------------------------------------------
drug.auc=read.csv("Research/AML/info_files/public_AMLdata/Beat_AML/BEAT_drug_AUC.csv")
colnames(drug.auc)[1]="Sample"
drug.auc=merge(bulk.meta[,c("Sample","CG_group")],drug.auc,by="Sample")
drug.auc=drug.auc[drug.auc$CG_group != "Other",]
drug.auc$CG_group[drug.auc$CG_group %in% c("Diploid-mono", "Diploid-Nonmono")]="Diploid"
drug.auc$CG_group1=drug.auc$CG_group
drug.auc$CG_group1[drug.auc$CG_group1 != 'Diploid']="Nondiploid"

drug.index=c()
for (i in 3:(ncol(drug.auc)-1)){
  test=kruskal.test(drug.auc[,i]~drug.auc$CG_group)
  if(test$p.value < 0.05){
    print(test$p.value)
    drug.index=c(drug.index,colnames(drug.auc)[i])}
}
drug.index1=c()
for (i in 3:(ncol(drug.auc)-1)){
  test=wilcox.test(drug.auc[drug.auc$CG_group1=="Diploid",i],drug.auc[drug.auc$CG_group1=="Nondiploid",i])
  if(test$p.value < 0.05){drug.index1=c(drug.index1,colnames(drug.auc)[i])}
}


gm=groupMeans(t(drug.auc[,drug.index]),drug.auc$CG_group)
gm=as.data.frame(na.omit(gm))
gm=gm[,c("Diploid","Del5.5q","Del7.7q","Double del")]
pdf("Research/manuscript/intratumor_heterogeneity/ITH_figure1/beat_drug_cyto3.pdf",width=10,height=6)
Heatmap(scale(t(gm)),col=colorRamp2(c(-2, 0, 2), c("#2166AC","white","#B2182B")),#switch color
        #colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
        cluster_columns=TRUE, cluster_rows = FALSE, column_names_rot = 90, column_title=NULL, 
        show_heatmap_legend = T,row_names_side = "left", show_column_names = TRUE,
        heatmap_legend_param = list(title = ""), width = ncol(gm)*unit(20, "mm"), 
        height = nrow(gm)*unit(4, "mm"),row_names_gp = gpar(fontsize = 15),column_names_gp = gpar(fontsize = 14))
dev.off()
rm(gm,drug.index,drug.index1,drug.auc,test,i)

##for each drug, test one CG group with all others; manually label * on the heatmap above
for (i in drug.index[1:4]){
  temp=drug.auc[,c("CG_group",i)]
  for (j in unique(drug.auc$CG_group)){
    temp$CG_group1=temp$CG_group
    temp$CG_group1[temp$CG_group1!=j]='Other'
    test=pairwise.wilcox.test(temp[[i]], temp$CG_group1)
    if (test$p.value < 0.05){
      print(paste0(i," ",j," ",test$p.value))
    }
  }
}



#------------------------------------------------------------------------------
#correlation of drug sensitivity with hallmark score and cell type abundance
#------------------------------------------------------------------------------
##code to process sample name in bulk.hallmark is in Figure1 code of inflammation paper
bulk.hallmark=read.csv("Research/manuscript/cell_of_origin/Figure1_new/scoring/TCGA_Abbas_BEAT_ssGSEA_tpm_Bofei.csv")
bulk.hallmark=bulk.hallmark[,1:51];colnames(bulk.hallmark)[1]="Sample"
bulk.hallmark=bulk.hallmark[(bulk.hallmark$Sample %in% drug.auc$Sample),]
colnames(bulk.hallmark)[2:51]=unlist(lapply(colnames(bulk.hallmark)[2:51],function(x){substr(x,10,nchar(x))}))
cor.mat=data.frame(matrix(ncol=122,nrow=50))
rownames(cor.mat)=colnames(bulk.hallmark)[2:51];colnames(cor.mat)=colnames(drug.auc)[3:124]
#rownames(cor.mat)=unlist(lapply(rownames(cor.mat),function(x){substr(x,10,nchar(x))}))
p.mat=cor.mat
for (i in colnames(bulk.hallmark)[2:51]){
  for (j in colnames(drug.auc)[3:124]){
    cor.mat[i,j]=cor(bulk.hallmark[,i],drug.auc[,j],use = "complete.obs")
    p.mat[i,j]=cor.test(bulk.hallmark[,i],drug.auc[,j])$p.value
  }
}

#ggcorrplot(t(cor.mat),method="circle",tl.srt=90)
pdf("Research/manuscript/intratumor_heterogeneity/ITH_figure1/beat_drug_hallmark_corr.pdf",width=15,height =12)
corrplot(as.matrix(cor.mat),p.mat = as.matrix(p.mat),method = "circle",insig = "pch",pch.cex = 0.5,col = colorRampPalette(c("blue","white","red"))(200),tl.offset=0,tl.cex = 0.6)
dev.off()


##abundance of cell types (deconvolution) with drug response
drug.auc=merge(bulk.deconv,drug.auc,by="Sample")
corr=rcorr(as.matrix(drug.auc[,c(2:6,10:131)]))
correlation=as.data.frame(corr$r[1:5,6:127]);pvalue=corr$P[1:5,6:127]
col.index=c()
for (i in 1:ncol(correlation)){
  if(any(pvalue[,i]<0.05)){col.index=c(col.index,i)}
}
correlation=correlation[,col.index];pvalue=pvalue[,col.index]
col.index=c()
for (i in 1:ncol(correlation)){
  correlation[pvalue[,i]>0.05,i]=NA
  if(sum(!is.na(correlation[abs(correlation[,i])>0.3,i]))>0){col.index=c(col.index,i)}
}
correlation=correlation[,col.index];pvalue=pvalue[,col.index]

pdf("Research/manuscript/intratumor_heterogeneity/ITH_figure2/celltype_drug_resistance_corr.pdf",width=10,height=6)
pheatmap(correlation, cluster_rows=FALSE, cluster_cols=FALSE,cellwidth = 15, cellheight = 20,na_col="white",
         scale="none",show_colnames=T,show_rownames=T, fontsize_row = 15, fontsize_col = 12)
dev.off()
rm(corr,correlation,drug.auc,pvalue)









############################################################################
#Figure 2: Symphony projection, celltype distribution
#          deg between cell types, PCA plot, SCENIC
############################################################################
# Symphony projection:code in heterogenity_paper.R;plotBasic
p = ref.in %>%dplyr::sample_frac(1L) %>% # permute rows randomly
  ggplot(aes(x =UMAP1,y = UMAP2),label= TRUE) + xlab("UMAP1")+ylab("UMAP2")+
  geom_point_rast(aes(col = get("BioClassification")), size = 0.4, stroke = 0.4, shape = 16)
p = p + theme_bw() +labs(color = "Cell type")+theme(legend.position="right") + theme_classic()+
  theme(legend.text = element_text(size=12), legend.title=element_text(size=12),
        axis.title = element_text(size=12), axis.text = element_text(size=12)) +
  guides(colour = guide_legend(override.aes = list(size = 4))) + guides(alpha = 'none')
p = p + scale_color_manual(values = type.colors1)
p=p+ annotate(geom="text",x=label.coor$UMAP1.median,y=label.coor$UMAP2.median,label=label.coor$BioClassification,cex=4) 
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure2/Symphony/reference.pdf",height=12,width=16, units="cm")


normal=readRDS("/Users/bwang8/Research/AML/normal/normal456_annotated.RDS")
normal.lympho=subset(normal,subset=class2 %in% c("CLP","lympho"))
##plot_Basic functions is in heterogeneity_paper.R
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure2/Symphony/lymphoid_3healthy.pdf",height=12,width=16.8, units="cm")
### project cells on AML Andy's map, put in supplementary
### code is in Symphony_mapping_Andy.R
###Reference map and projections of amlnew are saved at the following path
#/Users/bwang8/Research/manuscript/intratumor_heterogeneity/ITH_figure2/symphony_AndyRef






#------------------------------------------------------------------------------
# cell fraction in AML
#------------------------------------------------------------------------------
meta=amlnew@meta.data[,c("orig.ident","cyto3","class2","Andy.prediction")]
show_col(ColAssign(letters[1:20]))

fraction=meta %>% group_by(cyto3,class2)  %>% summarise(count=n()) 
fraction=fraction %>% group_by(cyto3) %>% mutate(Percent = round(count/sum(count),4))
fraction=as.data.frame(fraction)
fraction$orig.ident=factor(fraction$orig.ident,levels = c("PT12A","PT19A",'PT20A',"PT23A","PT25A","PT30A","PT32A","PT9A", "PT13A","PT17A","PT22A","PT26A","PT14A","PT15A","PT21A","PT27A","PT28A","PT10A","PT16A","PT29A"))
#fraction$cyto3=factor(fraction$cyto3,levels = c("Diploid","Del5.5q","Del7.7q","Double deletion"))
fraction$class2=factor(fraction$class2,levels = c("HSC-like","CMP.LMPP-like","CLP-like","GMP-like","lympho-like", "Mono-like", "Baso-like","Erythroid-like","DC-like"))
#fraction$Andy.prediction=factor(fraction$Andy.prediction,levels = c("LSPC-Quiescent","LSPC-Primed","LSPC-Cycle","GMP-like","ProMono-like", "Mono-like","cDC-like"))
ggplot(fraction, aes(fill=class2, y=Percent, x=cyto3)) + theme_classic() + xlab("") +
  geom_bar(position="stack",stat="identity") + labs(fill='Cell type') +
  scale_x_discrete(labels=c("Diploid (n=7)","Del5.5q (n=5)","Del7.7q (n=5)","Double del (n=3)")) +
  scale_fill_manual(values = brewer.pal(n = 10, name = 'Set3'))+
  theme(plot.title = element_text(hjust = 0.5),axis.title.y = element_text(size = 16),
        axis.text.y =element_text(size=12), axis.text.x =element_text(angle = 90,vjust = 0.5,size=12))
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure2/celltype_cyto3_bar.pdf",height=12,width=10, units="cm")



####same stacked barplot for Andy's prediction
ggplot(fraction, aes(fill=Andy.prediction, y=Percent, x=cyto3)) + theme_classic() + xlab("") +
  geom_bar(position="stack",stat="identity") + labs(fill='Cell type') +
  scale_x_discrete(labels=c("Diploid (n=7)","Del5.5q (n=5)","Del7.7q (n=5)","Double del (n=3)")) +
  scale_fill_manual(values=c("#8DD3C7","#1F78B4","#FFFFB3","#FB8072","#BEBADA","#FDB462","#D9D9D9"))+
  theme(plot.title = element_text(hjust = 0.5),axis.title.y = element_text(size = 16),
        axis.text.y =element_text(size=12), axis.text.x =element_text(angle = 90,vjust = 0.5,size=12))
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure2/LSPC_cyto3_bar.pdf",height=12,width=10, units="cm")


#####pie chart of all AML cells
#####Used in the supplementary
type.colors= c("HSC-like" = "#8DD3C7","CMP.LMPP-like" = "#FFFFB3","CLP-like" = "#BEBADA","GMP-like" = "#FB8072", "lympho-like" = "#80B1D3","Mono-like" = "#FDB462", "Baso-like" = "#B3DE69", "Erythroid-like" = "#FCCDE5", "DC-like" = "#D9D9D9")
meta=amlnew@meta.data
fraction=meta %>% group_by(class2)  %>% summarise(count=n()) %>% mutate(Percent = round(count/sum(count),4))
fraction=as.data.frame(fraction)
fraction$class2=factor(fraction$class2,levels = names(type.colors))
df <- fraction %>% mutate(csum = rev(cumsum(rev(count))), 
                          pos = count/2 + lead(csum, 1),pos = if_else(is.na(pos), count/2, pos))

df[2,5]=19800;df[4,5]=100;df[5,5]=1800;df[1,5]=7200;df[9,5]=9800
df[7,5]=48000.0;df[3,5]=30000;df[6,5]=15500;df[8,5]=10500

ggplot(fraction, aes(x="", y=count, fill=class2)) +geom_bar(stat="identity", width=1, color="white") + 
  coord_polar("y", start=0,direction = -1) +scale_fill_manual(values = type.colors) + theme_void() + 
  theme(legend.position="right",legend.text=element_text(size=12), legend.title=element_blank()) +   
  geom_label_repel(data = df,aes(y = pos, label = paste0(Percent *100, "%")),size = 3, 
                   nudge_x = c(1,0.7,0,1,0,0,0,0.8,0), show.legend=FALSE) 
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure2/class2_prop_pie.pdf",width=15,height =10,unit="cm")
rm(meta,fraction,df,type.colors)


#------------------------------------------------------------------------------
# dotplot of normal cell type markers in each AML cell type
# demonstrate the phentypical annotation of Symphony
#------------------------------------------------------------------------------
#group_dot_plot function in R script dot_plot.R
feat_stem=c("KIT","CD44","CD34")
feat_gmp=c("ELANE","AZU1","MPO")
feat_lympho=c("CD3E","CD19","CD79A","NCAM1")
feat_mono=c("CD14","ITGAX","S100A12")
feat_ery=c("HBD","HBB","GATA2");feat_dc=c("CD1C","PLD4","CST3")

allcell_dot <- group_dot_plot(amlnew, group = "class2",
           features=list(feat_stem,feat_gmp,feat_lympho,feat_mono,feat_ery,feat_dc))
leg <- get_legend(dot_plot_leg(amlnew, feat_stem))
pdf("Research/manuscript/intratumor_heterogeneity/ITH_figure2/class2_normal_marker.pdf", height=8, width=5)
ggarrange(allcell_dot)
dev.off()

dot_plot_leg(amlnew, rev(feat_stem), group = "class2")
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure2/class2_normal_marker1.pdf",width=10,height =15,unit="cm")
rm(feat_stem,feat_gmp,feat_lympho,feat_ery,feat_mono,feat_dc,allcell_dot,leg)






#------------------------------------------------------------------------------
# UMAP of AML cells highlight by cell type and cytogenetics
#------------------------------------------------------------------------------
celltype_2.colors <- c(
  "HSC-like"="#8DD3C7","CMP.LMPP-like"="#FFFFB3","CLP-like"="#BEBADA","GMP-like"="#FB8072",
  "lympho-like"="#80B1D3", "Mono-like"="#FDB462", "Baso-like"="#B3DE69" ,"Erythroid-like"="#FCCDE5",
  "DC-like"="#D9D9D9")
amlnew$class2 <- factor(amlnew$class2,levels=names(celltype_2.colors))

amlnew <- RunUMAP(amlnew,reduction = "harmony",seed.use = 123,dims=1:20,umap.method='uwot',min.dist=0.01,spread=1, a= 0.8, b=5,verbose=FALSE)
DimPlot(amlnew,group.by="class2",cols=celltype_2.colors,raster=F)+
  NoAxes()+ggtitle("")+
  theme(legend.title=element_text(size = 0),
        legend.text=element_text(size = 8),
        legend.key.height = unit(0.25,"cm"),
        legend.key.width = unit(0.25,"cm"))+NoAxes()
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure2/AMLcell_UMAP.pdf",height=15,width=18, units="cm")


###color by each cytogenetics
ref.in=data.frame(amlnew@reductions$umap@cell.embeddings)
ref.in$cyto=amlnew$cyto3;colnames(ref.in)[1:2]=c("UMAP1", "UMAP2")
query.in = ref.in[ref.in$cyto=="Double deletion",]#c("indianred2", "chartreuse1", "gold2","deepskyblue")
p4=ggplot()+geom_point(data=ref.in, aes(x=UMAP1, y=UMAP2),color = 'gray',size=0.3) +
  geom_point(data=query.in %>%dplyr::sample_frac(1L),aes(x = UMAP1, y = UMAP2),color = "deepskyblue",size = 0.3, stroke = 0.8, shape = 16) + theme_classic() + theme(legend.position="none") +
  theme(legend.text = element_text(size=10), legend.title=element_text(size=10),
        axis.title = element_text(size=10), axis.text = element_text(size=10)) 

grid.arrange(p1, p2, p3, p4, nrow = 2, ncol=2)
g <- arrangeGrob(p1, p2, p3, p4, nrow=2, ncol=2)
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_Supp/supp_UMAP_cyto3_box.pdf",height=15,width=15, units="cm",g)

rm(g,ref.in,query.in,p1,p2,p3,p4)




#------------------------------------------------------------------------------
# DEG between cell types in class2 and class3, provide table top 30 genes in each 
# as signatures in supplementary
#------------------------------------------------------------------------------
deg.class3=read.csv("Research/manuscript/intratumor_heterogeneity/ITH_figure2/deg_celltype/deg.class3.csv")
deg.class3=deg.class3[deg.class3$p_val_adj<0.5,]
deg.class3=deg.class3[!(grepl(pattern="^(MT|AJ|AL|AC|LIN|AP)", deg.class3$gene)),]
deg.class3=deg.class3 %>% arrange(desc(avg_log2FC)) %>% group_by(cluster) %>% slice(1:30)
#Reduce(rbind,by(deg.class3,deg.class3["cluster"],head,n = 30))  #same as above, different methods
deg.class3$timevar=c(rep(1:30,9))
deg.class3=reshape(as.data.frame(deg.class3[,6:8]), idvar = "cluster", timevar = "timevar", direction = "wide")
write.csv(deg.class3,file = "Research/manuscript/intratumor_heterogeneity/ITH_Supp/class3_signature.csv",row.names = FALSE)






#------------------------------------------------------------------------------
#cluster membership heatmap with Andy's LSPC states
#------------------------------------------------------------------------------
a <- table(amlnew$class2,amlnew$Andy.prediction)
b <- t(a/rowSums(a))
b=b[c("LSPC-Quiescent","LSPC-Primed","LSPC-Cycle","GMP-like","ProMono-like", "Mono-like","cDC-like"),c("HSC-like","CMP.LMPP-like","CLP-like","GMP-like","lympho-like", "Mono-like", "Baso-like","Erythroid-like","DC-like")]
p <- pheatmap::pheatmap(mat = b, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), border_color = "black", cluster_cols = FALSE, cluster_rows = FALSE,cellheight=20, cellwidth = 20)
pdf("Research/manuscript/intratumor_heterogeneity/ITH_figure2/Andy_cyto3_membership.pdf",width=12, height=12)
Heatmap(t(b),col=colorRamp2(c(0,1), c("peachpuff","#B2182B")),
        cluster_columns=FALSE, cluster_rows = FALSE, column_names_rot = 90, column_title=NULL, 
        show_heatmap_legend = T,row_names_side = "left", show_column_names = TRUE,
        heatmap_legend_param = list(title = ""), width = ncol(b)*unit(8, "mm"), 
        height = nrow(b)*unit(12, "mm"),row_names_gp = gpar(fontsize = 15),column_names_gp = gpar(fontsize = 14))
dev.off()








#------------------------------------------------------------------------------
#test HSC-like are phenotypical functional
#------------------------------------------------------------------------------
df=amlnew@meta.data[,c("class2","Eppert_LSC_higher_than_HSC","Eppert_LSC.R_Up","NG.AML.LSC","LSC17","KMT2A.r_LSC")]
gm=groupMeans(t(df[,2:6]),df$class2)
gm=gm[,c("HSC-like","CMP.LMPP-like","CLP-like","GMP-like","lympho-like", "Mono-like", "Baso-like","Erythroid-like","DC-like")]
pdf("/Users/bwang8/Research/manuscript/intratumor_heterogeneity/ITH_figure2/LSC_score_class2_new.pdf",width=8,height =8)
Heatmap(t(scale(t(gm))),
        col=colorRamp2(c(-2, 0, 2), c("#2166AC","white","#B2182B")),
        #colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
        cluster_columns=FALSE, cluster_rows = FALSE, column_names_rot = 90, column_title=NULL, 
        show_heatmap_legend = T,row_names_side = "left", show_column_names = TRUE,
        heatmap_legend_param = list(title = ""), width = ncol(gm)*unit(8, "mm"), 
        height = nrow(gm)*unit(8, "mm"),row_names_gp = gpar(fontsize = 15),column_names_gp = gpar(fontsize = 14))
dev.off()






#------------------------------------------------------------------------------
#SCENIC: transcription factors for cell type are lineage specific
# code to process  scenicoptions and select top10 TFs in each cell type is in scenic_new.R 
#------------------------------------------------------------------------------
setwd("/Users/bwang8/Research/manuscript/cell_of_origin/SCENIC_cells0.1frac_class2")
scenicOptions <- readRDS("./int/scenicOptions_final.Rds")
aml.sub=readRDS("/Users/bwang8/Research/manuscript/cell_of_origin/SCENIC_cells0.1frac_class2/amlsub_0.1frac_class2.RDS")
cellInfo=aml.sub@meta.data
pdf("/Users/bwang8/Research/manuscript/intratumor_heterogeneity/ITH_figure2/class2_scenic_heatmap1.pdf",width=12,height =6)
Heatmap(t(regulonActivity_byCellType_Scaled[unique(tf),]),
        row_names_side = "left",column_names_rot = 90,
        name="Regulon activity",cluster_columns = F,cluster_rows = F,
        col=colorRamp2(c(-2, 0, 2), c("#2166AC","white","#B2182B")),
        #column_split = c(rep("1",5),rep("2",5),rep("3",5),rep("4",5),rep("5",5),rep("6",5),rep("7",5),rep("8",5),rep         ("9",5)),column_gap = unit(1, "mm"),
        height = unit(40, "mm"),width = unit(200, "mm"),row_names_gp = gpar(fontsize = 15))

dev.off()

rm(cellInfo,scenicOptions,aml.sub,regulonActivity_byCellType,regulonActivity_byCellType_Scaled,AUC.Mat,i,j,tf,test.x,sel.TFs,sel.TFs2,regulonAUC)







############################################################################
#Figure 3: ITH by cell types: diffusion map, scenic,pathways
#         
############################################################################
#------------------------------------------------------------------------------
#take the HSC/GMP of each cytogenetic group and run DEG
#Also run erythroid as a control.It is more committed and evenly distributed, expect to show less variety
#------------------------------------------------------------------------------
hsc=subset(amlnew,subset=class2 %in% c("HSC-like"))
hsc=analyze_seurat_harmony(hsc)
deg.hsc=FindAllMarkers(hsc,logfc.threshold=0, min.pct=0)
write.csv(deg.hsc,file="Research/manuscript/intratumor_heterogeneity/ITH_figure3/HSC_analysis/deg_hsc_all.csv",row.names = F)
##GMP
gmp=subset(amlnew,subset=class2 %in% c("GMP-like"))
gmp=analyze_seurat_harmony(gmp)
deg.gmp=FindAllMarkers(gmp,logfc.threshold=0, min.pct=0)
write.csv(deg.gmp,file="Research/manuscript/intratumor_heterogeneity/ITH_figure3/GMP_analysis/deg_gmp_all.csv",row.names = F)
deg.sub=deg.eryth[deg.eryth$cluster=='Diploid',]
fgseaRes=run_fgsea(group1="Diploid",group2 = "",deg.sub,gmt_file  = "Research/AML/info_files/AUCell_geneset/h.all.v7.5.1.symbols.gmt.txt",out.path = "Research/manuscript/intratumor_heterogeneity/ITH_figure3/Erythroid_analysis")



#### combine NES in each group and plot heatmap; use p.adj < 0.2 instead of 0.25 to reduce the total pathways retained
gsea.diploid=read.csv("Research/manuscript/intratumor_heterogeneity/ITH_figure3/HSC_analysis/gsea_diploid.csv")
gsea.diploid=gsea.diploid[gsea.diploid$padj < 0.2,]
gsea.diploid=gsea.diploid[,c("pathway","NES")]
gsea.del5=gsea.del5[,c("pathway","NES")]
gsea.del7=gsea.del7[,c("pathway","NES")]
gsea.double=gsea.double[,c("pathway","NES")]
##contruct the matrix of NES for 4 cytogenetic groups
hallmark=merge(gsea.diploid,gsea.del5,all=TRUE,by="pathway")
hallmark=merge(hallmark,gsea.del7,by="pathway",all=TRUE)
hallmark=merge(hallmark,gsea.double,by="pathway",all=TRUE)
colnames(hallmark)=c('pathway','Diploid.HSC',"Del5/5q.HSC","Del7/7q.HSC","Double_deletion.HSC")


##do the same for GMP cells and create hallmark1
#hallmark.all=read.csv("Research/manuscript/intratumor_heterogeneity/ITH_figure3/pheno_hallmark.csv")
hallmark.all=merge(hallmark.hsc,hallmark.lmpp,by="pathway",all=TRUE)
rownames(hallmark.all)=hallmark.all$pathway;hallmark.all=hallmark.all[,2:9]
rownames(hallmark.all)=unlist(lapply(rownames(hallmark.all),function(x){substr(x,10,nchar(x))}))

hallmark.all=hallmark.all[c("APICAL_JUNCTION","APICAL_SURFACE",
                            "ADIPOGENESIS","ANGIOGENESIS","EPITHELIAL_MESENCHYMAL_TRANSITION","MYOGENESIS",  
                            "PANCREAS_BETA_CELLS",
                            "DNA_REPAIR", "UV_RESPONSE_DN","UV_RESPONSE_UP",
                            "ALLOGRAFT_REJECTION","COAGULATION","COMPLEMENT", "INFLAMMATORY_RESPONSE",   
                            "INTERFERON_ALPHA_RESPONSE","INTERFERON_GAMMA_RESPONSE","IL6_JAK_STAT3_SIGNALING",
                            "FATTY_ACID_METABOLISM","CHOLESTEROL_HOMEOSTASIS","GLYCOLYSIS","HEME_METABOLISM",         
                            "OXIDATIVE_PHOSPHORYLATION","XENOBIOTIC_METABOLISM",
                            "APOPTOSIS","HYPOXIA","PROTEIN_SECRETION","UNFOLDED_PROTEIN_RESPONSE",  
                            "REACTIVE_OXYGEN_SPECIES_PATHWAY",
                            "E2F_TARGETS","G2M_CHECKPOINT","MYC_TARGETS_V1","MYC_TARGETS_V2","P53_PATHWAY", 
                            "MITOTIC_SPINDLE", 
                            "ESTROGEN_RESPONSE_EARLY","ESTROGEN_RESPONSE_LATE","ANDROGEN_RESPONSE","NOTCH_SIGNALING", 
                            "IL2_STAT5_SIGNALING","KRAS_SIGNALING_UP","MTORC1_SIGNALING","KRAS_SIGNALING_DN",
                            "PI3K_AKT_MTOR_SIGNALING","TGF_BETA_SIGNALING","TNFA_SIGNALING_VIA_NFKB"),]

pdf("Research/manuscript/intratumor_heterogeneity/ITH_figure3/gsea_hallmark_4groupHSC.LMPP.GMP.Eryth_heatmap.pdf",width=15,height =12)
Heatmap(as.matrix(hallmark.all),name="NES",col=colorRamp2(c(-2, 0, 2), c("#2166AC","white","#B2182B")),
        cluster_columns=FALSE, cluster_rows = FALSE, #column_names_rot = 90, #width = ncol(gm)*unit(20, "mm"),
        column_title=NULL,show_heatmap_legend = T,row_names_side = "left", show_column_names = TRUE,
        column_split = c(rep("1",4),rep("2",4),rep("3",4),rep("4",4)),
        row_split = c(rep("1",2),rep("2",5),rep("3",3),rep("4",7),rep("5",6),rep("6",5),rep("7",6),rep("8",11)),
        column_gap = unit(3, "mm"),row_gap = unit(3, "mm"),heatmap_legend_param = list(title = ""),
        width = ncol(hallmark.all)*unit(8, "mm"), height = nrow(hallmark.all)*unit(5, "mm"),
        row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 11))
dev.off()
rm(hallmark,gsea.del5,gsea.del7,gsea.diploid,gsea.double,deg)



#------------------------------------------------------------------------------
#scatter plot show 3 Dimensions of pathways, replace the heatmap above
#------------------------------------------------------------------------------
df=read_excel("Research/manuscript/intratumor_heterogeneity/ITH_figure3/pheno_function_hallmark.xlsx",sheet = "pheno_function")
df=as.data.frame(df)
df$NES_pheno=as.numeric(df$NES_pheno);df$NES_function=as.numeric(df$NES_function)
df=df[!(is.na(df$NES_pheno) | is.na(df$NES_function)),]
df$label=NA
df$label[abs(df$NES_pheno)>2 | abs(df$NES_function)>2]=df$pathway[abs(df$NES_pheno)>2 | abs(df$NES_function)>2]
df$label=unlist(lapply(df$label,function(x){substr(x,10,nchar(x))}))
df$Cyto=factor(df$Cyto,levels = c("Diploid","Del5","Del7","Double_del"))
ggplot(data=df, aes(x = NES_pheno, y = NES_function, color=Cyto,label=label)) + #
  geom_point(stroke = 0, alpha = 1,shape = 16,size=3) + theme_classic() +
  geom_text_repel(color="black",size=2,segment.size=0.1,point.padding=NA,max.overlaps=40) +
  #scale_color_gradient2(midpoint=0, low="steelblue4", mid="#f0f0f0",high="#a50026",space ="Lab") +
  scale_color_manual(name="CG",values = c("indianred2", "chartreuse1", "gold2","deepskyblue")) +
  geom_hline(yintercept=0, linetype="dashed", color = "indianred")+
  geom_vline(xintercept=0, linetype="dashed", color = "indianred")+
  xlab("NES score of phenotypical axis") + ylab("NES score of functional axis")+
  labs(color = "Cytogenetic groups")
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure3/hallmark_NES_function+pheno_scatter.pdf",height=12,width=15, units="cm")





#------------------------------------------------------------------------------
#Use 50 hallmark to cluster cytogenetic groups, may apply PCA/NMF (not included)
#------------------------------------------------------------------------------
df=amlnew@meta.data[,c("class2",names(pathways))]
df=df[df$class2=="CMP.LMPP-like",]
df=df[,2:51]
pca = prcomp(df, center = TRUE,rank. = 50) #, scale = TRUE
summary(pca)
transform = as.data.frame(pca$x[,1:2])
transform$cellID=rownames(transform)
transform=merge(transform,amlnew@meta.data[,c("cellID","cyto3")],by="cellID")
transform$cyto3=factor(transform$cyto3,levels = c('Diploid',"Del5.5q","Del7.7q","Double deletion"))
ggplot(transform, aes(PC1, PC2, color=cyto3)) + 
  geom_point(aes(color=cyto3),size=0.4) + #stat_ellipse() + 
  xlab("PC1 (27.71%)") +ylab("PC2 (15.73%)") +theme_classic()+
  scale_color_manual(name="CG",values = c("indianred2", "chartreuse1", "gold2","deepskyblue")) +
  theme(axis.text = element_text(size = 12),axis.title = element_text(size=12),legend.text = element_text(size = 12))
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure3/GMP_analysis/pca_hallmark_aucell.pdf",height=10,width=13, units="cm")


#PCA on 50 features didn't clearly seperate cytogenetic groups, apply LASSO to remove
df$cellID=rownames(df)
df=merge(df,amlnew@meta.data[,c("cellID","cyto3")],by="cellID")
rownames(df)=df$cellID; df=df[,-1]
y=as.matrix(df[,51])
y=factor(y,levels = c("Diploid","Del5.5q","Del7.7q","Double deletion"))
y=unclass(y)
x=as.matrix(df[,-51])
cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "multinomial")
plot(cv.lasso)
fit=glmnet(x,y,family="multinomial",alpha=1,lambda=cv.lasso$lambda.min)
path.keep=names(pathways)
for (i in 1: length(coef(fit))){
  path.keep=intersect(path.keep,row.names(coef(fit)[[i]])[which(coef(fit)[[i]]!=0)][-1])
}
df=df[,1:50]
df=df[,path.keep]
pca = prcomp(df, center = TRUE,rank. = ncol(df)) #, scale = TRUE
summary(pca)
transform = as.data.frame(pca$x[,1:2])
transform$cellID=rownames(transform)
transform=merge(transform,amlnew@meta.data[,c("cellID","cyto3")],by="cellID")
transform$cyto3=factor(transform$cyto3,levels = c('Diploid',"Del5.5q","Del7.7q","Double deletion"))
ggplot(transform, aes(PC1, PC2, color=cyto3)) + 
  geom_point(aes(color=cyto3),size=0.4) + #stat_ellipse() + 
  xlab("PC1 (28.42%)") +ylab("PC2 (17.65%)") +theme_classic()+
  scale_color_manual(name="CG",values = c("indianred2", "chartreuse1", "gold2","deepskyblue")) +
  theme(axis.text = element_text(size = 12),axis.title = element_text(size=12),legend.text = element_text(size = 12))
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure3/LMPP_analysis/pca_hallmark_aucell_LASSO48.pdf",height=10,width=13, units="cm")
rm(df,cv.lasso,fit,pca,x,y,i,path.keep,transform)


###NMF: doesn't work well
nmf.out=nmf(df,50)
coef=nmf.out@fit@H
coef=as.data.frame(t(coef))
coef$cellID=rownames(coef)
coef=merge(coef,amlnew@meta.data[,c("cellID","cyto3")],by="cellID")

ggplot(coef, aes(V1, V2, color=cyto3)) + theme_classic()+
  geom_point(aes(color=cyto3),size=0.4) + #stat_ellipse() + 
  xlab("PC1 (28.42%)") +ylab("PC2 (17.65%)") +
  scale_color_manual(name="CG",values = c("indianred2", "chartreuse1", "gold2","deepskyblue")) +
  theme(axis.text = element_text(size = 12),axis.title = element_text(size=12),legend.text = element_text(size = 12))



#------------------------------------------------------------------------------
#Calculate shannon entropy to represent ITH in bulk AML 
#------------------------------------------------------------------------------
##example of measuring ITH in TCGA. Other datasets are calculated in the same way. 
raw.count=as.data.frame(GetAssayData(amlnew,slot = "counts",assay="RNA"))
human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
BM_labels <- getBM(attributes = c("hgnc_symbol","ensembl_gene_id","entrezgene_id", 'chromosome_name', "start_position","end_position","transcript_start","transcript_end"),filters = "hgnc_symbol",values = rownames(raw.count),mart = human)
BM_labels$length=BM_labels$transcript_end-BM_labels$transcript_start
BM_labels=BM_labels[order(BM_labels$length,decreasing = T),]
BM_labels=BM_labels %>% distinct(hgnc_symbol,.keep_all = TRUE)
raw.count=raw.count[BM_labels$hgnc_symbol,]
gene.length=BM_labels[order(match(BM_labels$hgnc_symbol,rownames(raw.count))),"length"]

rpk <- apply(raw.count, 2, function(x) x/(gene.length/1000))
#normalize by the sample size using rpk values
tpm <- apply(rpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)
p_i=tpm/apply(tpm,2,sum)
ith=apply(p_i,2,function(x){-sum(x[x!=0] * log2(x[x!=0]))})
ith=data.frame(cellID=names(ith),Shanno=ith)
###
ith=read.csv("Research/manuscript/intratumor_heterogeneity/ITH_figure3/sc_ith.csv")
ith=merge(ith,amlnew@meta.data[,c("cellID",'cyto3',"LSC17")],by="cellID")
p=ggplot(ith,aes(x=cyto3,y=Shanno)) + 
  geom_boxplot(aes(fill=cyto3),width=0.5,color='gray40',outlier.shape=NA,size=0.2)+
  #geom_jitter(fill='gray80',color='black',shape=21,width =0.2,size=0.25,stroke = 0.2)+
  ylab("Shanno Entropy") + theme_classic() + ylim(6,9) +
  scale_fill_manual(values =c("indianred2", "darkseagreen3", "gold2","deepskyblue"))+
  theme(legend.position="none",
        axis.title.y = element_text(size = 10),axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10))
kruskal.test(ith$Shanno~ith$cyto3)
p+annotate(geom="text",x=2.5,y=9,label="kruskal-Wallis, p < 2.2e-16",cex=4)
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure3/sc_ith_box.pdf",width=10,height =10,unit="cm")

p=ggplot(ith,aes(x=cyto3,y=Shanno,fill=cyto3)) + theme_classic()+
  geom_violin() + ylab("Shanno Entropy")+ ylim(6,9) +
  scale_fill_manual(values =c("indianred2", "darkseagreen3", "gold2","deepskyblue"))+
  theme(axis.text.x = element_text(size = 7),axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),axis.title.y = element_text(size = 10),legend.position="none")+ 
  stat_summary(fun=mean, geom="point", shape=18,size=2, color="red")
p+annotate(geom="text",x=2.5,y=9,label="kruskal-Wallis, p < 2.2e-16",cex=4)
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure3/sc_ith_violin.pdf",width=10,height =10,unit="cm")


###correlate ITH with LSC17
ggscatter(ith,add="reg.line",conf.int = FALSE,
          cor.coef = TRUE,cor.coef.size = 3,cor.method = "pearson",ylab = "Shanno Entropy",
          x= "LSC17", y= "Shanno",fill = "cyto3",shape = 21, size = 0.2,
          palette = c("indianred2", "darkseagreen3", "gold2","deepskyblue")) + ylim(6,9) +
  theme(legend.position = "right", axis.title = element_text(size = 10),axis.text = element_text(size = 10))





#------------------------------------------------------------------------------
#SCENIC regulon activity across cytogenetics in each AML cell type 
# Merge together show heatmap like pathways
#------------------------------------------------------------------------------
###process the regulon activity, select top TFs in each CG group in each cell type, merge the regulon activity df
write.xlsx(regulon.eryth, file="Research/manuscript/intratumor_heterogeneity/ITH_figure3/regulon.xlsx", sheetName="Eryth", append = TRUE,row.names=TRUE)

tf.all=merge(regulon.hsc,regulon.lmpp,by="TF",all=TRUE)
tf.all=merge(tf.all,regulon.gmp,by="TF",all=TRUE)
tf.all=merge(tf.all,regulon.eryth,by="TF",all=TRUE)
rownames(tf.all)=tf.all$TF;tf.all=tf.all[,2:17]
colnames(tf.all)=unlist(lapply(colnames(tf.all),function(x){substr(x,1,nchar(x)-2)}))
tf.all=tf.all[unique(c(regulon.hsc$TF,regulon.lmpp$TF,regulon.gmp$TF,regulon.eryth$TF)),]

pdf("Research/manuscript/intratumor_heterogeneity/ITH_figure3/scenic_4group_4cell_heatmap.pdf",width=16,height =10)
Heatmap(as.matrix(t(tf.all)),name="Regulon activity",col=colorRamp2(c(-2, 0, 2), c("#2166AC","white","#B2182B")),
        cluster_columns=FALSE, cluster_rows = FALSE, #column_names_rot = 90, #width = ncol(gm)*unit(20, "mm"),
        column_title=NULL,show_heatmap_legend = T,row_names_side = "left", show_column_names = TRUE,
        #column_split = c(rep("1",4),rep("2",4),rep("3",4),rep("4",4)),
        row_split = c(rep("1",4),rep("2",4),rep("3",4),rep("4",4)),
        column_gap = unit(3, "mm"),row_gap = unit(3, "mm"),heatmap_legend_param = list(title = ""),
        height = ncol(tf.all)*unit(6, "mm"), width = nrow(tf.all)*unit(4, "mm"),
        row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 11))
dev.off()

rownames(regulon.eryth)=regulon.eryth[,1];regulon.eryth=regulon.eryth[,2:5]
pdf("/Users/bwang8/Research/manuscript/intratumor_heterogeneity/ITH_figure3/Erythroid_analysis/scenic_top5.pdf",width=6,height =10)
Heatmap(as.matrix(regulon.eryth),name="Regulon activity",col=colorRamp2(c(-2, 0, 2), c("#2166AC","white","#B2182B")),
        cluster_columns=FALSE, cluster_rows = FALSE, #column_names_rot = 90, #width = ncol(gm)*unit(20, "mm"),
        column_title=NULL,show_heatmap_legend = T,row_names_side = "left", show_column_names = TRUE,
        #column_split = c(rep("1",4),rep("2",4),rep("3",4),rep("4",4)),
        #row_split = c(rep("1",4),rep("2",4),rep("3",4),rep("4",4)),
        column_gap = unit(3, "mm"),row_gap = unit(3, "mm"),heatmap_legend_param = list(title = ""),
        width = ncol(regulon.eryth)*unit(6, "mm"), height = nrow(regulon.eryth)*unit(4, "mm"),
        row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 11))
dev.off()



############################################################################
#Figure 5: LSPC signature and phenotypic signature
#          generate a mixed signature combining both axis
#          discuss HOPX that is contained in all signatures
############################################################################
signature=load_genesets("/Users/bwang8/Research/AML/info_files/AUCell_geneset/Dicklab_Genesets.gmt.txt")
deg.class2=read.csv("Research/manuscript/intratumor_heterogeneity/ITH_figure1/deg_celltype/deg.class2.csv")
hsc=deg.class2[deg.class2$cluster=="HSC-like",]
hsc=hsc[order(hsc$avg_log2FC,decreasing = T),]
hsc=hsc[!(grepl(pattern="^(MT|AJ|AL|AC|LIN|AP)", hsc$gene)),]
lspc.hsc.union=union(signature$`LSPC-Quiescent`,union(signature$`LSPC-Primed-Top100`,signature$`LSPC-Cycle-Top100`))
lspc.hsc.union=union(lspc.hsc.union,hsc[1:100,'gene'])
lspc.hsc.overlap=intersect(signature$`LSPC-Quiescent`,intersect(signature$`LSPC-Primed-Top100`,hsc[1:100,'gene']))


####feature selection methods in machine learning to reduce the size of signatures while keep its prognostic value
dat=bulk.tpm[lspc.hsc.union,]
dat=dat[rownames(dat) %in% rownames(bulk.tpm),]
dat=as.data.frame(t(dat));dat$Sample=rownames(dat)
dat=merge(dat,bulk.meta[,c("Sample","OS_days")])
rownames(dat)=dat$Sample
dat= dplyr::select(dat,-c("Sample"))
dat=dat[!is.na(dat$OS_days),]
#dat=dat[,c(row.names(coef(fit))[which(coef(fit)!=0)][-1],"OS_days")]
folds <- createFolds(y=dat[,262],k=10) #create folds
fc=as.numeric()
model_pre=as.numeric()
#intest=createDataPartition(y=dat[,1],p=0.25,list=F)
###ten fold cross validation
for (i in 1:10){
  fold_test=dat[folds[[i]],]
  fold_train=dat[-folds[[i]],]
  model=glm(OS_days~.,data=fold_train,family=binomial(link='logit'))
  pre=predict(model,newdata=fold_test,type='response')
  fc=append(fc,fold_test[,48])
  model_pre=append(model_pre,as.numeric(pre))
}

###LASSO cox
dat=bulk.tpm[lspc.hsc.union,]
dat=dat[rownames(dat) %in% rownames(bulk.tpm),]
dat=as.data.frame(t(dat));dat$Sample=rownames(dat)
dat=merge(dat,bulk.meta[,c("Sample","OS_days","Dead_Alive")])
rownames(dat)=dat$Sample
dat= dplyr::select(dat,-c("Sample"))
dat=dat[!is.na(dat$OS_days),]
dat$Dead_Alive=ifelse(dat$Dead_Alive=='Alive',0,1)
dat=dat[dat$OS_days>0,]
y=as.matrix(dat[,c(262,263)]);colnames(y)=c("time","status")
x=as.matrix(dat[,-c(262,263)])
#fit2=coxph(as.formula(paste0("Surv(OS_days, Dead_Alive==1)", paste0(colnames(x), collapse=" + "), sep=" ~ ")), dat)
#fit2=coxph(Surv(OS_days, Dead_Alive==1) ~ ., dat)
cv.lasso2 <- cv.glmnet(x, y, alpha = 1, family = "cox")
fit2=glmnet(x,y,family="cox",alpha=1,lambda=cv.lasso2$lambda.min)
coef(fit2)
row.names(coef(fit2))[which(coef(fit2)!=0)][-1]

##LASSO logistic
dat=bulk.tpm[lspc.hsc.union,]
dat=dat[rownames(dat) %in% rownames(bulk.tpm),]
dat=as.data.frame(t(dat));dat$Sample=rownames(dat)
dat=merge(dat,bulk.meta[,c("Sample","OS_days")])
rownames(dat)=dat$Sample
dat= dplyr::select(dat,-c("Sample"))
dat=dat[!is.na(dat$OS_days),]
dat$OS_days=ifelse(dat$OS_days>median(dat$OS_days),0,1)
y=as.matrix(dat[,262])
x=as.matrix(dat[,-262])
cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial")
#plot(cv.lasso)
fit=glmnet(x,y,family="binomial",alpha=1,lambda=cv.lasso$lambda.min)
coef(fit)
row.names(coef(fit))[which(coef(fit)!=0)][-1]


######Random forest to rank logistic signature
control <- rfeControl(functions = rfFuncs, # random forest
                      method = "repeatedcv", # repeated cv
                      repeats = 5, # number of repeats
                      number = 10) # number of folds
result_rfe1 <- rfe(x = x[,row.names(coef(fit.beat))[which(coef(fit.beat)!=0)][-1]], 
                   y = as.factor(y), sizes = c(1:30),rfeControl = control, metric="Accuracy")

result_rfe1
predictors(result_rfe1)
row.names(varImp(result_rfe1))[1:30]


####------------------------------------------------------------------------------
#CCS42 level in bulk and single cell data
amlnew$Andy.prediction=factor(amlnew$Andy.prediction,levels = c("LSPC-Quiescent","LSPC-Primed","LSPC-Cycle","GMP-like","ProMono-like", "Mono-like","cDC-like"))
VlnPlot(amlnew,features = "CCR42",group.by = "class2",pt.size = 0) +
  #scale_fill_manual(values=c("#8DD3C7","#1F78B4","#FFFFB3","#FB8072","#BEBADA","#FDB462","#D9D9D9"))+
  scale_fill_manual(values =brewer.pal(n = 10, name = 'Set3'))+ ylab("CCS42") +
  stat_summary(fun=median, geom="point", shape=18,size=2, color="red") +
  theme(axis.text.x=element_text(size=10,angle = 90,vjust = 0.5),axis.title.y = element_text(size = 12),
        axis.text.y=element_text(size = 12),axis.title.x = element_blank(),legend.position="none",
        legend.text=element_text(size = 12),legend.title=element_text(size = 12))
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure5/CCS42_class2.pdf",width=10,height =10,unit="cm")
  

df=bulk.meta[,c("LSPC_HSC_LASSO_logit","CG_group","FAB_Simple")]
df=df[df$FAB_Simple %in% c("M0","M1","M2","M4","M5"),]
kruskal.test(bulk.meta$LSPC_HSC_LASSO_logit~bulk.meta$FAB_Simple)
ggplot(df,aes(x=FAB_Simple,y=LSPC_HSC_LASSO_logit)) +
  geom_boxplot(aes(fill=FAB_Simple),width=0.5,color='gray40',outlier.shape=NA,size=0.1)+
  geom_jitter(fill='gray80',color='black',shape=21,width =0.2,size=1,stroke = 0.2)+
  ylab("CCS42") + theme_classic() + #ylim(0,0.2) +
  annotate(geom="text",x=3,y=1.3,label="kruskal-Wallis, p < 2.2e-16",cex=6) +
  scale_fill_manual(values=c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd"))+
  theme(legend.position="none",axis.title.y = element_text(size = 12),axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12))
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure5/CCS42_bulk_fab.pdf",height=10,width=10, units="cm")






####------------------------------------------------------------------------------
#score all different signatures. LSPC signatures is insignificant for survival, put as supp
bulk.tpm=read.csv("Research/AML/info_files/public_AMLdata/Combined_tpm_beat_tcga_abbas_Bofei.csv")#process colname
rownames(bulk.tpm)=bulk.tpm$X;bulk.tpm=bulk.tpm[,2:971]
gsva_score=gsva(as.matrix(bulk.tpm), list(LSPC.quiescent=signature$`LSPC-Quiescent`,LSPC.primed=signature$`LSPC-Primed-Top100`,LSPC.cycle=signature$`LSPC-Cycle-Top100`,HSC.like100=hsc[1:100,'gene'],lspc.hsc.union=lspc.hsc.union,quiescent.hsc.overlap=intersect(signature$`LSPC-Quiescent`,hsc[1:100,'gene'],Eppert_LSC=signature$Eppert_LSC,Fares_HSC=signature$Fares_HSC_common),primed.hsc.overlap=intersect(signature$`LSPC-Primed-Top100`,hsc[1:100,'gene'])),method="ssgsea",abs.ranking=F)
gsva_score=as.data.frame(t(gsva_score))
gsva_score$Sample=rownames(gsva_score)
bulk.meta=merge(bulk.meta,gsva_score,by="Sample")


survival=bulk.meta[,c("Sample","OS_days","Dead_Alive","LSPC.quiescent", "LSPC.primed","LSPC.cycle", "HSC.like100",  "lspc.hsc.union","quiescent.hsc.overlap","primed.hsc.overlap","LSPC_HSC_LASSO_logit","LSPC_HSC_LASSO_cox","LSPC_HSC_LASSO_RF")]
survival[(!is.na(survival$OS_days) & survival$OS_days> 5*365 ),"Dead_Alive"]="Alive"
survival[(!is.na(survival$OS_days) & survival$OS_days> 5*365 ),"OS_days"]= 5*365
#survival=survival[survival$Cohort=="BEAT2",]
fit=survfit(Surv(OS_days,Dead_Alive == "Dead") ~ LSPC_HSC_LASSO_logit>quantile(LSPC_HSC_LASSO_logit,0.5), survival)
p=ggsurvplot(fit,data=survival,pval=T,palette = c("blue","red"),xlab="Survival Time (Years)",
             xscale=365,xlim=c(0,5*365),break.x.by = 1*365,surv.median.line = "hv",
             legend.title="lspc.hsc.lasso_logistic_RF", risk.table = F,#title="BEAT2, >0.5",
             legend.labs = c("Below Median", "Above Median"),pval.coord = c(1000,0.8),pval.size =10) 
pdf("Research/manuscript/intratumor_heterogeneity/ITH_figure5/LSPC_HSC.LASSO_logit_survival.pdf",width=12, height=10)
p$plot+theme(axis.text.x=element_text(size=28),axis.title.y = element_text(size = 27),
             axis.text.y=element_text(size = 28),axis.title.x = element_text(size = 27),
             legend.text=element_text(size = 25),legend.title=element_text(size = 26)) 
dev.off()
rm(deg.class2,hsc,signature,lspc.hsc.union,survival,fit,p)

#------------------------------------------------------------------------------
##correlation with LSC, HSC signatures etc. may put a cor matrix in supp as validation
pdf("Research/manuscript/intratumor_heterogeneity/ITH_Supp/CCS42_cor_bulk.pdf", width = 10, height = 10)
corrplot(cor(bulk.meta[,c("LSPC_HSC_LASSO_logit","LSC17","Eppert_LSC", "Fares_HSC","Eppert_HSC")]),
         method = 'circle',type="upper",
         col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),tl.col = "black") +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 10),legend.text = element_text(size = 10)) 
#amlnew@meta.data[,c("CCR42","LSC17","Eppert_LSC.R_Up","Eppert_HSC.R_Up", "Fares_HSC")]
dev.off()

#------------------------------------------------------------------------------
####venn plot to bring up HOPX and CD69
VENN.LIST=list(LSPC.Quiescent=signature$`LSPC-Quiescent`,LSPC.Primed=signature$`LSPC-Primed-Top100`,HSC.like=hsc, CCS42=lsc_hsc$LSC_HSC_LASSO_logistic)
venn.plot <- venn.diagram(VENN.LIST, NULL,color="white",
                          fill=c("RoyalBlue", "Salmon","Green","gold2"),#alpha=c(0.5,0.5,0.5),
                          cat.cex=1.5,# main="Venn diagram of three feature selection methods",
                          resolution=800,cex=2.55,main.cex = 2, main.fontface = 2)
pdf("Research/manuscript/intratumor_heterogeneity/ITH_Supp/LSPC_HSC_CCS_Venn.pdf",width=12,height =10)
grid.draw(venn.plot)
dev.off()

rm(fit,venn.plot,p,VENN.LIST,survival)



#HOPX and CD69 survival
survival=merge(survival,data.frame(CD69=as.numeric(bulk.tpm["CD69",]),HOPX=as.numeric(bulk.tpm["HOPX",]),Sample=as.character(colnames(bulk.tpm))),by="Sample")
survival$logCD69=log2(survival$CD69 + 1);survival$logHOPX=log2(survival$HOPX + 1)
fit=survfit(Surv(OS_days,Dead_Alive == "Dead") ~ HOPX>quantile(HOPX,0.5), survival)
p=ggsurvplot(fit,data=survival,pval=T,palette = c("blue","red"),xlab="Survival Time (Years)",
             xscale=365,xlim=c(0,5*365),break.x.by = 1*365,surv.median.line = "hv",
             legend.title="HOPX", risk.table = F,#title="BEAT2, >0.5",
             legend.labs = c("Below Median", "Above Median"),pval.coord = c(1000,0.8),pval.size =10) 
pdf("Research/manuscript/intratumor_heterogeneity/ITH_Supp/HOPX_survival.pdf",width=12, height=10)
p$plot+theme(axis.text.x=element_text(size=28),axis.title.y = element_text(size = 27),
             axis.text.y=element_text(size = 28),axis.title.x = element_text(size = 27),
             legend.text=element_text(size = 25),legend.title=element_text(size = 26)) 
dev.off()


####HOPX expression compare to healthy in bulk and single cell
bulk.temp=read.csv("Research/AML/info_files/public_AMLdata/MetaData_All_v4.csv")
bulk.temp=bulk.temp[bulk.temp$Healthy_Disease %in% c("Healthy, Individual CD34+","Healthy, pooled CD34+"),]
health_cd34=data.frame(Sample=bulk.temp$Sample,FAB_Simple="Healthy")
df=rbind(bulk.meta[,c('Sample',"FAB_Simple")],health_cd34)
df$group=ifelse(df$FAB_Simple=="Healthy","Healthy",'AML')
HOPX=as.data.frame(t(bulk.tpm["HOPX",]));HOPX$Sample=rownames(HOPX)
df=merge(df,HOPX,by="Sample");df$HOPX=log(df$HOPX+1)
df=df[df$FAB_Simple %in% c("M0","M1","M2","Healthy"),]

df %>% group_by("group") %>% pairwise_wilcox_test(HOPX ~ group, p.adjust.method = "BH") %>% add_significance("p.adj") 
ggplot(df,aes(x=group,y=HOPX)) +
  geom_boxplot(aes(fill=group),width=0.5,color='gray40',outlier.shape=NA,size=0.1)+
  geom_jitter(fill='gray80',color='black',shape=21,width =0.2,size=1,stroke = 0.2)+
  ylab("HOPX") + theme_classic() + #ylim(0,0.2) +
  annotate(geom="text",x=1.5,y=5.3,label="Wilcoxon, p =0.312",cex=4) +
  scale_fill_manual(values=c("brown1","Green"))+
  theme(legend.position="none",axis.title.y = element_text(size = 12),axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12))
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure5/HOPX_M012_health.pdf",height=10,width=8, units="cm")

rm(bulk.temp,df,health_cd34,HOPX)

###single cell
normal=readRDS("/Users/bwang8/Research/AML/normal/normal456_annotated.RDS")
aml.normal=merge(amlnew,y=normal,add.cell.ids=NULL,merge.data=T)
aml.normal=analyze_seurat_harmony(aml.normal)
aml.normal$group=ifelse(aml.normal$cellID %in% amlnew$cellID,"AML","Healthy")
VlnPlot(aml.normal,features = "HOPX",group.by = "group",pt.size = 0) +
  scale_fill_manual(values =c("brown1","Green"))+ ylab("HOPX") +
  #stat_summary(fun=median, geom="point", shape=18,size=2, color="red") +
  theme(axis.text.x=element_text(size=12,angle = 0,hjust = 0.5),axis.title.y = element_text(size = 12),
        axis.text.y=element_text(size = 12),axis.title.x = element_blank(),legend.position="none",
        legend.text=element_text(size = 12),legend.title=element_text(size = 12))
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure5/HOPX_AML_health456.pdf",width=8,height =10,unit="cm")





#------------------------------------------------------------------------------
####drug resistance of CCR42 and HOPX, CD69
drug.auc=read.csv("Research/AML/info_files/public_AMLdata/Beat_AML/BEAT_drug_AUC.csv")
colnames(drug.auc)[1]="Sample"
df=as.data.frame(t(bulk.tpm[c("HOPX","CD69"),]));df$Sample=rownames(df)
df=merge(bulk.meta[,c("Sample","LSPC_HSC_LASSO_logit")],df,by="Sample")
df=df[(df$Sample %in% drug.auc$Sample),];drug.auc=drug.auc[drug.auc$Sample %in% df$Sample,]
df$HOPX=log(df$HOPX+1);df$CD69=log(df$CD69+1)

cor.mat=data.frame(matrix(ncol=122,nrow=3))
rownames(cor.mat)=colnames(df)[2:4];colnames(cor.mat)=colnames(drug.auc)[2:123]
p.mat=cor.mat;index=c()
for (i in colnames(df)[2:4]){
  for (j in colnames(drug.auc)[2:123]){
    cor.mat[i,j]=cor(df[,i],-drug.auc[,j],use = "complete.obs")
    p.mat[i,j]=cor.test(df[,i],drug.auc[,j])$p.value
    if (any(p.mat[,j]<0.05,na.rm = TRUE)){index=c(index,j)}
  }
}

index=unique(index)
p.mat=p.mat[,index];cor.mat=cor.mat[,index]
pdf("Research/manuscript/intratumor_heterogeneity/ITH_figure5/CCS42_HOPX_drug_cor.pdf",width=12,height =3)
corrplot(as.matrix(cor.mat),p.mat = as.matrix(p.mat),method = "circle",insig = "pch",pch.cex = 0.5,col = colorRampPalette(c("blue","white","red"))(200),tl.offset=0,tl.cex = 0.6)
dev.off()

rm(df,p.mat,i,j,index,cor.mat)




#------------------------------------------------------------------------------
############################################
#Figure4 del7/7q analysis: This is actually figure 3, following the phenotypical and functional axes. 
#Enes works on it. Figure 4 is spatial by Chris 
############################################
del7=readRDS("Research/manuscript/intratumor_heterogeneity/ITH_figure4/del7.RDS")
########UMAP 
del7.colors <- c(`0` = '#1f77b4', `1` = '#ff7f0e')
del7$RNA_snn_res.0.1 <- factor(del7$RNA_snn_res.0.1,levels=names(del7.colors))
DimPlot(del7,group.by="RNA_snn_res.0.1",cols=del7.colors,raster=F)+
  NoAxes()+ggtitle("")+
  theme(legend.title=element_text(size = 0),
        legend.text=element_text(size = 8),
        legend.key.height = unit(0.25,"cm"),
        legend.key.width = unit(0.25,"cm"))+NoAxes()
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure4/del7_umap.pdf",height=15,width=16, units="cm")


########volcano plot of DEGs
deg.del7=read.csv("Research/manuscript/intratumor_heterogeneity/ITH_figure4/del7_FindAllMarkers.csv")
deg.del7=deg.del7[!(grepl(pattern="^(MT|AJ|AL|AC|LIN|AP)", deg.del7$gene)),]
deg.del7=distinct(deg.del7,gene, .keep_all=TRUE)

EnhancedVolcano(deg.del7,lab = deg.del7$gene,x = 'avg_log2FC',y = 'p_val_adj',xlim = c(-2,2),
                title="", subtitle = "",caption="",axisLabSize=12,pCutoff=10e-20,FCcutoff=0.5,
                titleLabSize = 10,legendLabSize = 10,legendIconSize = 3,
                selectLab = c('CD74','DUSP2',"CD52","IFITM3","HLA-DRB1","HLA-DMB","HLA-DQA1","AVP",'HBB','PRTN3','S100B', "CDK13","CDK6","ID3", "ARL5B","SKAP2","JAK1","MBNL1",'CD69',"ZNF460","SOX4"),
                pointSize = c(ifelse(deg.del7$avg_log2FC>1, 1, 0.5)),labSize = 2)
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure4/del7_deg_volcano.pdf",width=12,height =12,unit="cm")



####LSC signature in 2 clusters
df=del7@meta.data[,c("RNA_snn_res.0.1","LSC17","chr7","chr7p","chr7q")]
p=ggplot(df,aes(x=RNA_snn_res.0.1,y=chr7,fill=RNA_snn_res.0.1)) + theme_classic()+
  geom_violin() + scale_fill_manual(values = c('#1f77b4', '#ff7f0e')) + 
  ylab("Chr7 gene score")+ #ylim(9,12) +
  theme(axis.text.x = element_text(size = 7),axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),axis.title.y = element_text(size = 10),legend.position="none")+ 
  stat_summary(fun=mean, geom="point", shape=18,size=4, color="red")

wilcox.test(df[df$RNA_snn_res.0.1==0,"chr7"],df[df$RNA_snn_res.0.1==1,"chr7"])
p + annotate(geom="text",x=1.5,y=0.08,label="Wilcoxon, p < 2.2e-16",cex=5)
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure4/del7_chr7Score_violin.pdf",width=8,height =8,unit="cm")

###patient distribution across 2 clusters
meta.by.cyto = del7@meta.data %>% group_by(RNA_snn_res.0.1,orig.ident)  %>% summarise(count=n())
meta.by.cyto=meta.by.cyto %>% group_by(RNA_snn_res.0.1) %>% mutate(Percent = round(count/sum(count),4))
ggplot(meta.by.cyto, aes(fill=orig.ident, y=Percent, x=RNA_snn_res.0.1)) + theme_bw() + 
  xlab("") +labs(fill='Patient') + theme_classic()+
  geom_bar(position="stack",stat="identity")+ scale_fill_manual(values =brewer.pal(n = 10, name = 'Set3'))+
  #scale_x_discrete(labels=c("Diploid (n=7)","Del5.5q (n=5)","Del7.7q (n=5)","Double positive (n=3)")) +
  theme(axis.text.x=element_text(size=12), #text = element_text(face = "bold"),
        axis.title.y = element_text(size = 12),axis.text.y = element_text(size = 12),
        legend.text=element_text(size=10),legend.title=element_blank()) 
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure4/del7_pt_dist.pdf",width=8,height =8,unit="cm")



###cellular states in 2 clusters
df=del7@meta.data[,c("RNA_snn_res.0.1","orig.ident","class2","Andy.prediction")]
df$class2=factor(df$class2,levels = c("HSC-like","CMP.LMPP-like","CLP-like","GMP-like","lympho-like", "Mono-like", "Baso-like","Erythroid-like","DC-like"))
df$Andy.prediction=factor(df$Andy.prediction,levels = c("LSPC-Quiescent","LSPC-Primed","LSPC-Cycle","GMP-like","ProMono-like", "Mono-like","cDC-like"))
meta.by.cyto=df %>% group_by(RNA_snn_res.0.1,class2)  %>% summarise(count=n())
temp=as.data.frame(reshape(as.data.frame(meta.by.cyto[,c(1,2,3)]), timevar = "RNA_snn_res.0.1", idvar = "class2", direction = "wide"))
meta.by.cyto=meta.by.cyto %>% group_by(RNA_snn_res.0.1) %>% mutate(Percent = round(count/sum(count),4))

p=ggplot(meta.by.cyto, aes(fill=class2, y=Percent, x=RNA_snn_res.0.1)) + 
  theme_bw() + xlab("") +labs(fill='Cell type') +   theme_classic()+
  geom_bar(position="stack",stat="identity")+ scale_fill_manual(values =brewer.pal(n = 10, name = 'Set3'))+
  #scale_fill_manual(values=c("#8DD3C7","#1F78B4","#FFFFB3","#FB8072","#BEBADA","#FDB462","#D9D9D9"))+
  #scale_x_discrete(labels=c("Diploid (n=7)","Del5.5q (n=5)","Del7.7q (n=5)","Double positive (n=3)")) +
  theme(axis.text.x=element_text(size=12),axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text=element_text(size=10),legend.title=element_blank()) 

rownames(temp)=temp$class2;temp=temp[,2:3]
chisq.test(temp)
p + annotate(geom="text",x=1.5,y=1.05,label="Chi-square, p < 0.01",cex=4)
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure4/del7_class2_bar.pdf",width=10,height =10,unit="cm")



#------------------------------------------------------------------------------
#SCP package for DEGs and pathways
#------------------------------------------------------------------------------
#del7 <- RunDEtest(del7, group_by = 'RNA_snn_res.0.1', fc.threshold = 1, only.pos = F)
#del7_DEG <- del7@tools$DEtest_RNA_snn_res.0.1$AllMarkers_wilcox
#del7_DEG <- del7_DEG[with(del7_DEG, p_val_adj < 0.05),]
#VolcanoPlot(del7, group_by = 'RNA_snn_res.0.1', DE_threshold = 'p_val_adj < 0.05', nlabel = 10)
del7_FindAllMarkers <- read.csv("Research/manuscript/intratumor_heterogeneity/ITH_figure4/del7_FindAllMarkers.csv")

del7_FindAllMarkers$group1 <- del7_FindAllMarkers$cluster
del7_FindAllMarkers$group2 <- "others"
del7@tools$DEtest_RNA_snn_res.0.1$AllMarkers_wilcox <- del7_FindAllMarkers

del7 <- RunEnrichment(
  srt = del7, group_by = "RNA_snn_res.0.1", db = "Hallmark", species = "Homo_sapiens")
EnrichmentPlot(
  srt = del7, group_by = "RNA_snn_res.0.1", db = "Hallmark", plot_type = "bar",palette = "Set3",topTerm = 10)
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure4/del7_ORA.pdf", height = 14, width = 20)

del7 <- RunGSEA(
  srt = del7, group_by = "RNA_snn_res.0.1", db = "Hallmark", species = "Homo_sapiens",
  DE_threshold = "p_val_adj < 0.05")
GSEAPlot(srt = del7, group_by = "RNA_snn_res.0.1", db = "Hallmark")
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure4/del7_GSEA.pdf", height = 14, width = 20)


#------------------------------------------------------------------------------
##cell communications between 2 clusters
#------------------------------------------------------------------------------
data.input <-del7@assays$RNA@data
meta <- data.frame(ID=del7$cellID,labels=del7$RNA_snn_res.0.1)
meta$labels=ifelse(meta$labels==0,"More","Less")
# Create a CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels")
groupSize <- as.numeric(table(cellchat@idents))
# DB
CellChatDB <- CellChatDB.human 
# showDatabaseCategory(CellChatDB)
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","Cell-Cell contact","ECM-Receptor")) # use Secreted Signaling
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
# slow ...
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#############################
# Part II: Inference of cell-cell communication network
# slow ...
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

#############################
# Part III" Visualization
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions") #cellchat@net$weight
pathways.show.all <- cellchat@netP$pathways
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
pdf("Research/manuscript/intratumor_heterogeneity/ITH_figure4/cell_interaction_chord1.pdf",height=10,width=10)
#netVisual_chord_gene(cellchat, sources.use = "More", targets.use = "Less", lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cellchat, sources.use = "Less", targets.use = "More", lab.cex = 0.5,legend.pos.y = 30)
dev.off()

netVisual_bubble(cellchat, sources.use = c("More","Less"), targets.use = c("Less","More"), remove.isolate = FALSE)





#------------------------------------------------------------------------------
##deconvolve bulk with 2 clusters and check survival
#------------------------------------------------------------------------------
#create single cell reference matrix by subseting cells 
meta=del7@meta.data;meta$ID=rownames(meta)
meta=meta %>% group_by(RNA_snn_res.0.1) %>% sample_frac(size=0.2)
counts=del7@assays$RNA@data
counts=counts[,colnames(counts) %in% meta$ID];counts=counts[,meta$ID]
colnames(counts)=meta$RNA_snn_res.0.1
genes=counts[!(grepl(pattern="^(MT|AJ|AL|AC|LIN|AP)", rownames(counts))),]
fwrite(as.data.frame(genes),file="Research/manuscript/intratumor_heterogeneity/ITH_figure4/del7_20percent_cibersort_reference.txt",sep="\t",row.names=T,quote=T)
#subset for scenic
del7.sub=subset(del7,cellID %in% meta$ID)

results_abbas=read.csv("/Users/bwang8/Downloads/Abbas_del7_res0.1.csv")
results_abbas$Mixture=substr(results_abbas$Mixture,2,10)
results_abbas$Mixture=gsub("\\.","\\-",results_abbas$Mixture)
colnames(results_abbas)[1]="Sample"
for (i in 1:nrow(results_abbas)){
  results_abbas$group[i]=names(which.max(results_abbas[i,2:3]))
}
survival=merge(results_abbas[,c(1,7)],bulk.meta,by="Sample")
survival[(!is.na(survival$OS_days) & survival$OS_days> 5*365 ),"Dead_Alive"]="Alive"
survival[(!is.na(survival$OS_days) & survival$OS_days> 5*365 ),"OS_days"]= 5*365

fit=survfit(Surv(OS_days,Dead_Alive == "Dead") ~ group, survival)
p=survminer::ggsurvplot(fit,data=survival,pval = T,xscale=365,xlim=c(0,9*365),break.x.by = 3*365, xlab = "Year")
ggsave("")
survminer::pairwise_survdiff(Surv(OS.time,OS == 1) ~ group,data=phei,p.adjust.method = "BH")

res.cox <- coxph(Surv(OS.time,OS) ~ group, data = phei)

p=survminer::ggsurvplot(fit,data=allsurvival,pval = F,xscale=365,xlim=c(0,9*365),break.x.by = 3*365, xlab = "Year",
                        legend.labs = paste("Cluster", 1:4),risk.table = F,palette = c("indianred2", "chartreuse1", "gold2","deepskyblue"))
png("",res=200,width=12, height=10,units="in")
p$plot+annotate("text", x = 7.5*365, y = 0.8, label = "logrank p = 0.0084", cex=10, col="black", vjust=0, 
                hjust = 1.1, fontface=22) +
  theme(axis.text.x=element_text(size=32),axis.title.y = element_text(size = 27),
        axis.text.y=element_text(size = 32),axis.title.x = element_text(size = 27),
        legend.text=element_text(size = 23),legend.title=element_text(size = 23)) 
dev.off()




#------------------------------------------------------------------------------
##CNV score in 2 clusters; show as barplot in supp
#------------------------------------------------------------------------------
#calculate the inferCNV score for each cell and add to metadata of seurat objects
groupSum <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups) == nrow(mat))
  gm <- lapply(unique(groups), function(x) {
    if (sparse) {
      Matrix::colSums(mat[which(groups == x),, drop = F], na.rm = na.rm) /nrow(mat[which(groups == x),, drop = F])
    }
    else {
      colSums(mat[which(groups == x),, drop = F], na.rm = na.rm)/nrow(mat[which(groups == x),, drop = F])
    }
  }) %>% Reduce("rbind", .)
  rownames(gm) <- unique(groups)
  return(gm)
}
cnvscore.all=data.frame(matrix(ncol = 0, nrow = 22))
for (i in list.files("/Users/bwang8/Research/AML/cleaned_PTs/infercnvScore_20PTs_Mono")){
  infer_obj=readRDS(paste0("/Users/bwang8/Research/AML/cleaned_PTs/infercnvScore_20PTs_Mono/",i))
  gene.order=infer_obj@gene_order;gene.order$gene=rownames(gene.order)
  expr=infer_obj@expr.data-1
  cnvscore=groupSum(expr,gene.order$chr)
  cnvscore.all=cbind(cnvscore.all,cnvscore)
}
rm(gene.order,expr,cnvscore,infer_obj,i)

###Add infercnv score to metadata, plot frequency (proportion of cells) of gain/del by chr; then seperate by diploid/nondiploid and compute effect size (Hedge’s g)
cnvscore.all=as.data.frame(t(cnvscore.all[,colnames(amlnew)]))
colnames(cnvscore.all)=paste0(colnames(cnvscore.all),"_cnvScore")
amlnew=AddMetaData(amlnew, as.data.frame(cnvscore.all))


##proportion of gain/del on each chromosome
#Cells with an average absolute CNV level that was above the 0.1 quantile of the entire population were considered as potentially malignant
cnv.freq <- function(cnvscore,name,cut.gain,cut.del){
  #gain=unlist(lapply(colnames(cnvscore),function(x){sum(cnvscore[,x]>quantile(cnvscore.all[cnvscore.all[,x]>0,x],probs=0.1))/nrow(cnvscore)}))
  gain=unlist(lapply(colnames(cnvscore),function(x){sum(cnvscore[,x]>cut.gain)/nrow(cnvscore)}))
  del=unlist(lapply(colnames(cnvscore),function(x){sum(cnvscore[,x]<cut.del)/nrow(cnvscore)}))
  #del=unlist(lapply(colnames(cnvscore),function(x){sum(cnvscore[,x]<quantile(cnvscore.all[cnvscore.all[,x]<0,x],probs=0.9))/nrow(cnvscore)}))
  ##split by cyto2 and plot for diploid and nondiploid seperately
  #del=unlist(lapply(colnames(cnvscore.all[df$cyto2=="Non-Diploid",]),function(x){sum(cnvscore.all[df$cyto2=="Non-Diploid",x]<0)/nrow(cnvscore.all[df$cyto2=="Non-Diploid",])}))
  merged=rbind(data.frame(chr=paste0("chr",1:22),cyto=rep("gain",22),proportion=gain),
               data.frame(chr=paste0("chr",1:22),cyto=rep("del",22),proportion=-del))
  merged$chr=factor(merged$chr,levels = paste0("chr",1:22))
  merged$cyto=factor(merged$cyto,levels = c("gain","del"))
  #y.up=ceiling(quantile(merged$proportion[merged$proportion>0],probs=0.8) * 10) / 10
  #yup1=ceiling((max(merged$proportion[merged$proportion>0])-0.1) * 10) / 10
  #y.down=floor(quantile(merged$proportion[merged$proportion<0],probs=0.2) * 10) / 10
  #ydown1=floor((min(merged$proportion[merged$proportion<0])+0.1) * 10) / 10
  p=ggplot(merged,aes(x=chr,y=proportion,fill=cyto))+geom_bar(stat="identity") + theme_linedraw() +#+coord_flip()
    theme(legend.text = element_text(size=26), legend.title=element_blank(),panel.grid = element_blank(),
          axis.title.x = element_blank(), axis.text.x = element_text(size=25,angle = 90,vjust=0.5),
          legend.position="top",axis.title.y = element_text(size=28), 
          axis.text.y = element_text(size=28),axis.ticks.x = element_blank()) + 
    scale_y_continuous(name="Proportion cells",labels=c(0.8,0.4,"0",0.4,0.8),breaks=c(0.8,0.4,0,-0.4,-0.8),limits=c(-0.9,0.9)) +
    scale_fill_manual(values=c("gain"="#B2182B","del"="#2166AC"),labels=c('Gain', 'Deletion')) +
    guides(colour = guide_legend(override.aes = list(size = 4))) + guides(alpha = 'none')
  pdf(paste0("Research/manuscript/intratumor_heterogeneity/ITH_figure4/",name,"CNV_bar.pdf"),width=10,height=10)
  print(p)
  dev.off()
}
cnvscore=del7@meta.data[,694:715]
cnvscore=cnvscore[rownames(cnvscore) %in% del7@meta.data[del7$RNA_snn_res.0.1==1,"cellID"],]
cnv.freq(cnvscore,"del7_cluster1_",quantile(cnvscore[cnvscore[,8]>0,8],probs=0.7),quantile(cnvscore[cnvscore[,7]<0,7],probs=0.9))  #cutoff of cnv for cells 25,-25





############################################
########SCENIC 
del7.sub=readRDS("/Users/bwang8/Research/manuscript/intratumor_heterogeneity/ITH_figure4/SCENIC/del7_res0.1_20percent.RDS")
scenicOptions <- readRDS("./int/scenicOptions.Rds")
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
AUC.Mat <- getAUC(regulonAUC)
regulonActivity_byCellType <- sapply(split(rownames(del7.sub@meta.data), del7.sub@meta.data$RNA_snn_res.0.1),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
rownames(regulonActivity_byCellType_Scaled)=unlist(lapply(rownames(regulonActivity_byCellType_Scaled),function(x){str_split(x," ")[[1]][1]}))
rownames(regulonActivity_byCellType_Scaled)=unlist(lapply(rownames(regulonActivity_byCellType_Scaled),function(x){str_split(x,"_")[[1]][1]}))
rownames(AUC.Mat)=unlist(lapply(rownames(AUC.Mat),function(x){str_split(x," ")[[1]][1]}))
rownames(AUC.Mat)=unlist(lapply(rownames(AUC.Mat),function(x){str_split(x,"_")[[1]][1]}))
print("prepare heatmap")
x <- ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity",cluster_columns = F)
sel.TFs <- data.frame(TF=rownames(AUC.Mat),p=NA,row.names = rownames(AUC.Mat))
for (i in sel.TFs$TF) {
  test.x <- wilcox.test(AUC.Mat[i,del7.sub@meta.data$RNA_snn_res.0.1==0],AUC.Mat[i,del7.sub@meta.data$RNA_snn_res.0.1!=0])
  sel.TFs[i,"p"] <- test.x$p.value
}
sel.TFs <- sel.TFs[order(sel.TFs$p),]
sel.TFs$tumor_hi <- NA
sel.TFs$tumor_hi[sel.TFs$TF %in% rownames(regulonActivity_byCellType_Scaled)[regulonActivity_byCellType_Scaled[,"0"]>0]] <- "cluster0"
sel.TFs2 <- sel.TFs[!is.na(sel.TFs$tumor_hi),]
pdf("Research/manuscript/intratumor_heterogeneity/ITH_figure4/scenic_heatmap_top10.pdf",width=10, height=10)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[c(sel.TFs2$TF[1:10],sel.TFs[is.na(sel.TFs$tumor_hi),"TF"][1:10]),],name="Regulon activity",cluster_columns = F,cluster_rows = F)
dev.off()
ha = ComplexHeatmap::HeatmapAnnotation(Field = anno_simple(del7.sub$RNA_snn_res.0.1),annotation_name_side = "left")
Mat <- AUC.Mat[sel.TFs$TF[1:10],];Mat <- t(scale(t(Mat)))

ComplexHeatmap::Heatmap(Mat,cluster_columns = T,top_annotation = ha,show_column_names = F,clustering_method_columns="ward.D")







#------------------------------------------------------------------------------
############################################
#Hierarchical model; clonal architecture (SCEVAN)
############################################
cell.high=rownames(pt10a@meta.data[pt10a$chr17.7674950.A.G=="Mutant",])#chr17.7675238.T.C
cell.high=rownames(pt16a@meta.data[pt16a$chr17.7675148.G.T=="Mutant",])#chr1.114713909.G.T
DimPlot(pt16a,label = T,cells.highlight = cell.high)

table(pt10a$class1,pt10a$chr17.7674950.A.G)
table(pt10a$subclone,pt10a$chr17.7674950.A.G)


#####SCEVAN on AML cells only
pt16_scevan=readRDS("Research/manuscript/intratumor_heterogeneity/Figure3/SCEVAN/pt16A_AML_SCEVAN.RDS")
df=pt16_scevan@meta.data[,c("cellID","class2","subclone","class","Andy.prediction")]
temp=as.data.frame.matrix(table(df$class2,df$subclone))
colnames(temp)=c("subclone1","subclone2","subclone3","subclone4","subclone5")
#temp=as.data.frame(apply(temp,2,function(x){x/sum(x)}))
temp$celltype=rownames(temp)
type.colors= c("HSC-like" = "#8DD3C7","CMP.LMPP-like" = "#FFFFB3","CLP-like" = "#BEBADA","GMP-like" = "#FB8072", "lympho-like" = "#80B1D3","Mono-like" = "#FDB462", "Baso-like" = "#B3DE69", "Erythroid-like" = "#FCCDE5", "DC-like" = "#D9D9D9")
type.colors1=c("LSPC-Quiescent"="#8DD3C7","LSPC-Primed"="#FFFFB3","LSPC-Cycle"="#BEBADA","GMP-like"="#FB8072","ProMono-like"="#80B1D3","Mono-like"="#FDB462","cDC-like"="#B3DE69")
temp$celltype=factor(temp$celltype,levels = temp$celltype)

chisq.test(temp[,1:6])$p.value

composition_pie <- function(scevan_mtx,group.by,out.path,factor.level,type.colors){
  df=scevan_mtx@meta.data[,c("cellID","subclone","class",group.by)]
  temp=as.data.frame.matrix(table(df[,group.by],df$subclone))
  colnames(temp)=paste0("subclone",1:length(unique(df[!is.na(df$subclone),"subclone"])))
  temp$group=rownames(temp)
  temp$group=factor(temp$group,levels = factor.level)
  plot.panel=par(mfrow=c(2,ceiling(length(unique(df[!is.na(df$subclone),"subclone"]))/2)),mar=c(1,1,1,1))
  p=list()
  for (i in 1:length(unique(df[!is.na(df$subclone),"subclone"]))){
    p[[i]]=print(ggplot(temp, aes(x="", y=!!sym(paste0("subclone",i)), fill=group)) +
      geom_bar(stat="identity", width=1, color="white") + 
      coord_polar("y", start=0,direction = -1) +scale_fill_manual(values = type.colors) + theme_void() + 
      theme(legend.position="right",legend.text=element_text(size=12), legend.title=element_blank(),
            plot.title = element_text(hjust = 0.5)) + labs(title = paste0("subclone",i)))
  }
  pdf(paste0(out.path,".pdf"),width=15, height=15)
  do.call(grid.arrange,p)
  dev.off()
}
composition_pie(pt25_scevan,"Andy.prediction","Research/manuscript/intratumor_heterogeneity/Figure3/SCEVAN/AML/PT25A_AML/pt25A_sub_LSC_pie.",c("LSPC-Quiescent","LSPC-Primed","LSPC-Cycle","GMP-like","ProMono-like","Mono-like","cDC-like"),type.colors1)
composition_pie(pt25_scevan,"class2","Research/manuscript/intratumor_heterogeneity/Figure3/SCEVAN/AML/PT25A_AML/pt25A_sub_class2_pie.",c("HSC-like","CMP.LMPP-like","CLP-like","GMP-like","lympho-like", "Mono-like", "Baso-like","Erythroid-like","DC-like"),type.colors)
composition_pie(pt25_scevan,"Phase","Research/manuscript/intratumor_heterogeneity/Figure3/SCEVAN/AML/PT25A_AML/pt25A_sub_phase_pie.",c("G1","S","G2M"),c("G1"="blue","S"="green","G2M"="red"))

####mutations by subclone
mut_heat<- function(object,name){
  df=object@meta.data[,c(520:555)]
  df=df[,apply(df,2,function(x){!("Not performed" %in% x)})]
  df=df[apply(df,1,function(x){sum(x=="Mutant")>0}),]
  df$cellID=rownames(df)
  df=merge(df,object@meta.data[,c("cellID","class2","subclone","Andy.prediction")],by="cellID")
  df=df[order(df$subclone),]
  df=df[!is.na(df$subclone),]
  ht=Heatmap(t(df[,c(2:(ncol(df)-3))]),name="Mutation",col=c("Wild type"="gray","Mutant"="firebrick1"),
          column_split = df$subclone,gap = unit(1, "mm"), 
          column_title=NULL,show_heatmap_legend = TRUE,row_names_side = "left",
          top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "white"))),
          row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 9),show_column_names = FALSE)
  pdf(paste0("Research/manuscript/intratumor_heterogeneity/Figure3/SCEVAN/",name,"_AML/mutation_subclone.pdf"),width = 12,height=2)
  draw(ht)
  dev.off()
}
mut_heat(pt13_scevan,"PT13A")




####run SCEVAN on local computer, not work
pt9a=subset(amlnew,subset=orig.ident %in% "PT9A")
pt9a=analyze_seurat(pt9a)
saveRDS(pt9a,file="/Users/bwang8/Downloads/pt9A_AML.RDS")
##del7 AML cells
del7_scevan=readRDS("Research/manuscript/intratumor_heterogeneity/Figure3/SCEVAN/del7_sub_SCEVAN.RDS")





#------------------------------------------------------------------------------
#Calculate shannon entropy to represent ITH in bulk AML 
#------------------------------------------------------------------------------
bulk.deconv=read.csv("Research/manuscript/intratumor_heterogeneity/ITH_figure1/combined3_20percent_class3_result.csv")
bulk.meta=read.csv("Research/manuscript/cell_of_origin/paper/Figure1/bulk_3sets_meta_Bofei.csv")

##example of measuring ITH in TCGA. Other datasets are calculated in the same way. 
tcga.tpm=read.csv("/Users/bwang8/Library/CloudStorage/OneDrive-InsideMDAnderson/LabMembers/Resources/BulkRNA/AML_TPM/TCGA_AML_TPM.csv")
rownames(tcga.tpm)=tcga.tpm$Genes;tcga.tpm=tcga.tpm[,2:ncol(tcga.tpm)]
colnames(tcga.tpm)=gsub("\\.","\\-",colnames(tcga.tpm))
colnames(tcga.tpm)=unlist(lapply(colnames(tcga.tpm),function(x){substr(x,1,12)}))
tcga.tpm=tcga.tpm[rowSums(is.na(tcga.tpm)) != ncol(tcga.tpm), ]
p_i=tcga.tpm/apply(tcga.tpm,2,sum)
tcga.ith=apply(p_i,2,function(x){-sum(x[x!=0] * log2(x[x!=0]))})
tcga.ith=as.data.frame(t(tcga.ith));tcga.ith$Sample=rownames(tcga.ith)
colnames(tcga.ith)[1]="ITH_score"

###
bulk.ith=read.csv("Research/manuscript/intratumor_heterogeneity/ITH_figure1/combined3_ITH_score.csv")

temp=merge(beat2.ith,bulk.deconv,by="Sample")
temp$cell_type=apply(temp[,3:7],1,function(x){names(which.max(x))})
temp=merge(temp,bulk.meta[,c("Sample","Cohort","CG_group1","Cytogenetic_risk_incomplete","HALLMARK.IFNG.RESPONSE","LSC17","OS_days","Dead_Alive")],by="Sample")
temp=merge(temp,bulk.meta[,c(1,64:102)],by="Sample")
temp$cell_type=factor(temp$cell_type,levels = c("Primitive.like","GMP.like","Lympho.like","Committed.like","Erythroid.like"))
temp$CG_group1[temp$CG_group1 %in% c("Diploid-mono","Diploid-ONS","Diploid-Nonmono")]="Diploid"
temp$CG_group1=factor(temp$CG_group1,levels = c("Diploid","Del5.5q","Del7.7q","Double del","Other"))
#temp=temp[,colSums(is.na(temp[,15:53])) != nrow(temp[,15:53])]
temp$mutation=apply(temp[,16:ncol(temp)],1,function(x){sum(x=="Positive",na.rm=TRUE)})
temp$mutation[temp$mutation>4]="5 and more"
temp$Cytogenetic_risk_incomplete[temp$Cytogenetic_risk_incomplete %in% c("","Unknown")]="Unknown"

ggplot(temp,aes(x=mutation,y=ITH_score)) + 
  geom_boxplot(aes(fill=mutation),width=0.5,color='gray40',outlier.shape=NA,size=0.2)+
  geom_jitter(fill='gray80',color='black',shape=21,width =0.2,size=1.25,stroke = 0.2)+
  ylab("ITH") + theme_classic() + ylim(8,12) +
  #scale_fill_manual(values =c("indianred2", "darkseagreen3", "gold2","deepskyblue","darkslategray"))+
  theme(legend.position="none",
        axis.title.y = element_text(size = 12),axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12))
p=ggplot(temp,aes(x=cell_type,y=LSC17,fill=cell_type)) + theme_classic()+
  geom_violin() + scale_fill_manual(values = brewer.pal(n = 10, name = 'Set3')[c(1,4,5,6,8)]) + 
  ylab("ITH")+ ylim(9,12) +
  theme(axis.text.x = element_text(size = 7),axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),axis.title.y = element_text(size = 10),legend.position="none")+ 
  stat_summary(fun=mean, geom="point", shape=18,size=2, color="red")

wilcox_stat=temp %>% group_by("cell_type") %>% pairwise_wilcox_test(ITH_score~cell_type, p.adjust.method = "BH") %>% add_significance("p.adj") %>% add_xy_position(x="cell_type", fun="max", step.increase = 0.1)
colnames(wilcox_stat)[1]="CG_group1"
p + stat_pvalue_manual(wilcox_stat, label="p.adj.signif", tip.length = 0.01,hide.ns = T,bracket.size = 0.5, size = 4.2)

p+annotate(geom="text",x=3,y=11.7,label=paste0("kruskal-Wallis, p = ", signif(kruskal.test(temp$ITH_score~temp$cell_type)$p.value,3)),cex=4)
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure2/ITH_score_celltype_beat_violin.pdf",width=10,height =10,unit="cm")

temp$Cluster=factor(temp$Cluster,levels = c("Primitive","GMP","Intermediate","Mature"))


###correlate ITH with LSC17
ggscatter(temp,add="reg.line",conf.int = FALSE,
          cor.coef = TRUE,cor.coef.size = 3,cor.method = "pearson",ylab = "ITH",
          x= "LSC17", y= "ITH_score",fill = "cell_type",shape = 21, size = 2,
          palette = brewer.pal(n = 10, name = 'Set3')[c(1,4,5,6,8)]) + ylim(8,13) +
  theme(legend.position = "right", axis.title = element_text(size = 10),axis.text = element_text(size = 10))
ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure3/ITH_score_LSC17_corr_celltype.pdf",width=15,height =12,unit="cm")


######survival by cell type
survival=temp[,c("Sample","OS_days","Dead_Alive","cell_type","ITH_score")]
survival[(!is.na(survival$OS_days) & survival$OS_days> 5*365 ),"Dead_Alive"]="Alive"
survival[(!is.na(survival$OS_days) & survival$OS_days> 5*365 ),"OS_days"]= 5*365
#survival=survival[survival$Cohort=="BEAT2",]
fit=survfit(Surv(OS_days,Dead_Alive == "Dead") ~ ITH_score>quantile(ITH_score,0.5), survival)
p=ggsurvplot(fit,data=survival,pval=T,surv.median.line = "hv",
             palette = brewer.pal(n = 10, name = 'Set3')[c(1,4,5,6,8)],
             xlab="Survival Time (Years)",xscale=365,xlim=c(0,5*365),break.x.by = 1*365,
             legend.title="Parsimonious IFNG", risk.table = FALSE,#title="BEAT2, >0.5",
             legend.labs = c("Primitive.like","GMP.like","Lympho.like","Committed.like","Erythroid.like"),
             pval.coord = c(1000,0.8),pval.size = 10) 
pdf("Research/manuscript/intratumor_heterogeneity/ITH_figure3/survival_celltype_allbulk.pdf",width=12, height=10)
p$plot+theme(axis.text.x=element_text(size=20),axis.title.y = element_text(size = 20),
             axis.text.y=element_text(size = 20),axis.title.x = element_text(size = 20),
             legend.text=element_text(size = 15),legend.title=element_text(size = 16)) 
dev.off()


################
#correlation between celltype percentage and LSC17 score
meta=amlnew@meta.data
meta.by.cyto=meta %>% group_by(orig.ident,class2)  %>% summarise(count=n())
meta.by.cyto=meta.by.cyto %>% group_by(orig.ident) %>% mutate(Percent = round(count/sum(count),4))
meta.by.cyto$class2=factor(meta.by.cyto$class2,levels=c("HSC-like","CMP.LMPP-like","CLP-like","GMP-like","lympho-like", "Mono-like", "Baso-like","Erythroid-like","DC-like"))

composition=reshape(as.data.frame(meta.by.cyto[,c(1,2,4)]), timevar = "class2", idvar = "orig.ident", direction = "wide")
colnames(composition)=substr(colnames(composition),9,30);composition[is.na(composition)]=0
colnames(composition)[1]="Patient"; composition=composition[,c('Patient',"HSC-like","CMP.LMPP-like","CLP-like","GMP-like","lympho-like", "Mono-like", "Baso-like","Erythroid-like","DC-like")]
lsc17=groupMeans(t(amlnew$LSC17),groups = amlnew$orig.ident,sparse = FALSE)
lsc17=as.data.frame(t(lsc17));lsc17$Patient=rownames(lsc17)
colnames(lsc17)[1]="LSC17"
composition=merge(composition,lsc17,by="Patient")
rownames(composition)=composition$orig.ident;composition=composition[,c(2:ncol(composition))]

correlation <- rcorr(as.matrix(composition),type="pearson")

pdf("/Users/bwang8/Research/manuscript/intratumor_heterogeneity/ITH_figure1/celltype_composition_corr1.pdf",width=10,height =10)
heatmap.2(correlation$r,col=colorRampPalette(c("blue",'white',"red"))(20), key=TRUE,density.info="none",trace="none",
          keysize = 1,margins=c(10,10))
dev.off()
rm(lsc17,composition,meta,meta.by.cyto,correlation)







#------------------------------------------------------------------------------
#drug targets proliferation (not work well)
#------------------------------------------------------------------------------
###drugs target proliferation pathwaysin Malani
# Palbociclib, Abemaciclib,Ribociclib,AT7519, Dinaciclib,Alvocidib
# Serdemetan: MDM2 inhibitor in tp53 mutant
# Cytarabine:interferes with DNA synthesis
# PKI-402, Temsirolimus, Sirolimus, Everolimus: indirect, target signaling (mTOR) then regulate cell cycle
# Duvelisib, Omipalisib:PI3K
# Ruxolitinib,Baricitinib:JAK

###drugs target proliferation pathwaysin BeatAML
#Palbociclib, AT7519,Flavopiridol
#"BEZ235","GDC-0941": PI3K_mTOR
#"Ruxolitinib","TG101348":JAK
drug.auc=read.csv("Research/AML/info_files/public_AMLdata/Beat_AML/BEAT_drug_AUC.csv")
colnames(drug.auc)[1]="Sample"
df=merge(drug.auc[,c("Sample","Palbociclib","AT7519","Flavopiridol","BEZ235","GDC.0941","Ruxolitinib..INCB018424.","TG101348")],bulk.meta[,c("Sample","CG_group")],by="Sample")
df$CG_group[df$CG_group %in% c("Del5.5q","Double del")]="Del5/5q"
df$CG_group[df$CG_group != "Del5/5q"]="Other"
temp=reshape(df[,2:9], direction="long", varying = list(1:7), sep="")
temp=melt(temp, id = c("CG_group","time")) 
pdf("Research/manuscript/intratumor_heterogeneity/ITH_figure1/beat_cycle_drug_cyto.pdf",width = 15,height=10)
bwplot(Palbociclib ~CG_group | paste0("Car ", time), data = temp)
dev.off()



ggplot(df,aes(x=CG_group,y=TG101348))+#scale_y_continuous(trans=squish_trans(0,6),breaks = seq(6,15))+
  geom_boxplot(aes(fill=CG_group),width=0.5,color='gray40',outlier.shape=NA,size=0.2)+
  geom_jitter(fill='gray80',color='black',shape=21,width =0.2,size=1.25,stroke = 0.2)+
  ylab("Palbociclib") + theme_classic() + #ylim(0,0.2) +
  #scale_fill_manual(values=c("darkseagreen3","indianred2", "darkslategray", "gold2","deepskyblue","mediumpurple1"))+
  theme(legend.position="none",
        axis.title.y = element_text(size = 12),axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12,angle = 90),axis.text.y = element_text(size = 12))

df %>% group_by("CG_group") %>% pairwise_wilcox_test(TG101348 ~ CG_group, p.adjust.method = "BH") %>%
  add_significance("p.adj") 
#ggsave("Research/manuscript/intratumor_heterogeneity/ITH_figure2/beat_cycle_drug_cyto.pdf",height=15,width=15, units="cm")








#------------------------------------------------------------------------------
#correlate genetic mutations with cytogenetic groups and cell compositions
#------------------------------------------------------------------------------
temp=bulk.meta[,c(1,64:102,548)]
temp$CG_group1[temp$CG_group1 %in% c("Diploid-mono", "Diploid-Nonmono","Diploid-ONS")]="Diploid"
temp[is.na(temp)]=""
temp$mut_count=apply(temp[,2:40],1,function(x){sum(x=="Positive")})

temp=merge(temp,bulk.deconv[,1:6],by="Sample")
rownames(temp)=temp$Sample
colnames(temp)[2:40]=unlist(lapply(colnames(temp)[2:40], function(x){substr(x,1,nchar(x)-4)}))
temp=temp[order(match(temp[,"CG_group1"],c("Diploid","Del5.5q","Del7.7q","Double del","Other"))),]
##chi-square test of each mutation
chi.res=c()
for (i in 2:40){
  df=as.data.frame(temp[,c(41,i)])
  df=df%>% group_by(CG_group1,!!sym(colnames(temp)[i]))  %>% summarise(count=n()) %>% as.data.frame()
  df=as.data.frame(reshape(df, timevar = "CG_group1", idvar = colnames(temp)[i], direction = "wide"))
  df[is.na(df)]=0
  chi.res=c(chi.res,chisq.test(df[,2:6])$p.value)
}

chi.res[chi.res<=0.01]="**"
chi.res[chi.res<0.05 & chi.res>0.01]="*"
chi.res[chi.res>=0.05]="ns"

column_ha=HeatmapAnnotation(bar=anno_barplot(temp$mut_count,add_numbers = F,axis=F,border = F,gp=gpar(fill="darksalmon"),annotation_name_gp= gpar(fontsize = 0)))
pdf("Research/manuscript/intratumor_heterogeneity/ITH_figure1/bulk_cyto_mutation_heatmap.pdf",width = 10,height=9)
p=Heatmap(t(temp[,2:40]),name="Mutation",
          col=c("slategray2","steelblue","red"),#rect_gp=gpar("black"),
          column_split = c(rep(1,294),rep(2,36),rep(3,34),rep(4,14),rep(5,294)),gap = unit(1, "mm"), 
          top_annotation = column_ha, #right_annotation = row_ha,
          column_title=NULL,show_heatmap_legend = T,row_names_side = "left",
          row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 9),show_column_names = FALSE) +
  rowAnnotation(anno=anno_text(chi.res),annotation_name_gp= gpar(fontsize = 4))
draw(p, annotation_legend_side = 'top',heatmap_legend_side="right")
dev.off()

###can also merge with beat2 drug sensitivity samples 
#temp=subset(temp,select=-c(ASXL2,NF1,WT1,SMC1A,SMC3,SPI1,PHF6))
#temp=temp[order(match(temp[,"CG_group"],c("Diploid","Del5.5q","Del7.7q","Double del")),-temp$mut_count),] 

######oncoprint of selected genes and cell type composition
pdf("Research/manuscript/intratumor_heterogeneity/ITH_figure1/bulk_selected_mutation_heatmap.pdf",width = 15,height=2)
Heatmap(t(temp[,c("NPM1","TP53","DNMT3A","FLT3","FLT3ITD","IDH1","IDH2","KRAS","RUNX1","TET2")]),name="Mutation",
        col=c("slategray2","steelblue","red"),#rect_gp=gpar("black"),
        column_split = c(rep(1,294),rep(2,36),rep(3,34),rep(4,14),rep(5,294)),gap = unit(1, "mm"), 
        column_title=NULL,show_heatmap_legend = FALSE,row_names_side = "left",
        row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 9),show_column_names = FALSE) 
dev.off()

df3=as.data.frame(pivot_longer(temp[,c(1,41,43:47)], cols = names(temp[,c(1,41,43:47)])[3:7], names_to = "cell_type",values_to = "Percentage"))
df3$cell_type=factor(df3$cell_type,levels = c("Primitive.like","GMP.like","Lympho.like","Committed.like", "Erythroid.like"))
pdf("Research/manuscript/intratumor_heterogeneity/ITH_figure1/bulk_cell_composition.pdf",width = 15,height=12)
ggplot(df3, aes_string(fill="cell_type",y="Percentage",x="Sample"))+theme_classic()+
  geom_bar(position="stack",stat="identity",show.legend = FALSE) + xlab("")+ylab("Percentage") + 
  scale_fill_manual(values =brewer.pal(n = 10, name = 'Set3')[c(1,4,5,6,8)]) +
  #scale_x_discrete(labels=c("Diploid (n=7)","Del5.5q (n=5)","Del7.7q (n=5)","Double del (n=3)")) +
  theme(axis.text.x=element_blank(), #labs(fill='Cell type') +
        axis.title.y = element_blank(),axis.text.y = element_text(size = 22),
        legend.text=element_text(size=22),legend.title=element_blank()) 
dev.off()

















############################################
#Heatmap of Hallmark pathways across cytogenetic groups
pathway_score=read_excel("Research/manuscript/cell_of_origin/Figure1_new/Inflammation_cyto/avg_pathScore_byCyto3.xlsx",sheet="hallmark")
pathway_score=as.data.frame(pathway_score);rownames(pathway_score)=pathway_score$...1
p=pheatmap(t(pathway_score[,c("Diploid","Del5.5q","Del7.7q","Double del")]), cluster_rows=F, show_rownames=T,cluster_cols=F,show_colnames = T,scale="none",angle_col=90, fontsize_row=20,fontsize_col = 12,cellheight=20, cellwidth = 15)
png("/Users/bwang8/Research/manuscript/cell_of_origin/Figure1_new/fig1d_hallmark_cyto3.png",res=200,width=14, height=8,units="in")
print(p)
dev.off()


#Heatmap of Hallmark pathways across patients;hallmark is the AUCell score of hallmark pathways
df=data.frame(ID=rownames(amlnew@meta.data),patient=amlnew$orig.ident)
df=merge(df,hallmark[,c(2,38:87)],by.x = "ID",by.y = 'X')
colnames(df)[3:52]=gsub("AUCell_","",colnames(df)[3:52])
gm=groupMeans(t(df[,3:52]),df$patient)
gm=gm[,c(20,1:19)]
p=pheatmap(gm, cluster_rows=F, show_rownames=T,cluster_cols=F,show_colnames = T,scale="row",angle_col=90, fontsize_row=18,fontsize_col = 12,cellheight=15, cellwidth = 15)
png("/Users/bwang8/Research/manuscript/cell_of_origin/Figure1_new/hallmark_patient.png",res=200,width=12, height=12,units="in")
print(p)
dev.off()


############################################
#Density plot of representative pathways in each CG group
hallmark=read.csv("Research/manuscript/cell_of_origin/Figure1_new/scoring/AML_prog_HALLMARK_AUCell.csv")
df=data.frame(ID=rownames(amlnew@meta.data),cyto=amlnew$cyto3)
df=merge(df,hallmark[,c(2,38:87)],by.x = "ID",by.y = 'X')
colnames(df)[3:52]=gsub("AUCell_","",colnames(df)[3:52])
for (path in c("HALLMARK_NOTCH_SIGNALING","HALLMARK_MYC_TARGETS_V2","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_PI3K_AKT_MTOR_SIGNALING")){
  p=ggplot(df, aes_string(x=path,fill="cyto")) + geom_density(alpha=0.4)+ theme_classic() + ylab("Density")+
    scale_fill_manual(values=c("indianred2", "chartreuse1", "gold2","deepskyblue"))+ xlab(path)+
    theme(axis.text.x = element_text(size = 30),legend.position = 'top',axis.title.x=element_text(size=24),
          axis.text.y = element_text(size = 30, angle = 0),axis.title.y=element_text(size=28),legend.text = element_text(size = 28),legend.title = element_text(size = 32))
  png(paste0("Research/manuscript/cell_of_origin/Figure1_new/fig1e_",path,"cyto3_density.png"),res=200,width=12, height=10,units="in")
  print(p)
  dev.off()
  
  p=ggplot(df,aes_string(x="cyto", y=path,fill="cyto"))+geom_violin(violinwidth=1) +
    geom_boxplot(width=0.3, color="grey", alpha=0.2) + theme_classic() +
    scale_fill_manual(values=c("indianred2", "chartreuse1", "gold2","deepskyblue")) +
    scale_x_discrete(labels=c('Diploid', 'Del5.5q', 'Del7.7q', "Double del"))+
    theme(legend.position = "none", axis.title = element_text(size = 25),
          axis.text.x = element_text(size = 26,angle = 0,hjust = 0.5),axis.title.x = element_blank(),
          axis.text.y = element_text(size = 32),axis.ticks = element_blank()) 
  test=kruskal.test(df[[path]]~df$cyto)
  png(paste0("Research/manuscript/cell_of_origin/Figure1_new/fig1e_",path,"cyto3_boxViolin.png"),res=200,width=8, height=10,units="in")
  print(p+annotate(geom="text",x=2.5,y=max(df[[path]])+0.01,label=paste0("kruskal-Wallis, p=",signif(test$p.value, digits=3)),cex=8))
  dev.off()
}
rm(hallmark,df,path,p,test)

#############################################################################
#SCP style Pathway figure
deg.del7=deg.del7[!(grepl(pattern="^(MT|AJ|AL|AC|LIN|AP)", deg.del7$X)),]
common.gene=c(deg.diploid$X[1:10],deg.del5$X[c(1:8,10,13)],deg.del7$X[1:10],deg.double$X[c(1:2,4:11)])
mat <- amlnew@assays$RNA@data[common.gene,]
group.list=amlnew@meta.data
group.list=data.frame(cyto1=group.list[["cyto3"]],class2=group.list[["class2"]],id=rownames(group.list),patients=group.list$orig.ident)
group.list$cyto1=as.character(group.list$cyto1)
group.list$cyto1[group.list$cyto1=="Double deletion"]="Double del"
group.list=group.list[order(match(group.list$cyto1,c("Diploid","Del5.5q","Del7.7q","Double del"))),]
group.list$patients=factor(group.list$patients,levels = paste0('PT',c(9,10,12:17,19:23,25:30,32),"A"))
rownames(group.list)=group.list$id
mat=as.matrix(mat[,rownames(group.list)])
mat=t(scale(t(mat)))

column_ha=HeatmapAnnotation(df=data.frame(
  Cytogenetics=factor(group.list[,1],levels = c("Diploid","Del5.5q","Del7.7q","Double del")),
  Patients=factor(group.list[,4])), annotation_name_side = "left",#gp = gpar(col=c("white")),
  annotation_name_gp= gpar(fontsize = 10),show_legend = c(TRUE,TRUE),
  col = list(Cytogenetics=c("Double del"="deepskyblue","Del7.7q"="gold2","Del5.5q"="darkseagreen3","Diploid"="indianred2")))
png("Research/manuscript/cell_of_origin/Figure1_new/fig1f_deg_cyto_heatmap.png",res=200,width=10,height=10,units="in")
p=Heatmap(mat,name="Expression",col=colorRamp2(seq(-2, 2, length = 100), palette_scp(palette = "RdBu")),
          top_annotation = column_ha, cluster_columns=FALSE, cluster_rows = FALSE, 
          column_title=NULL,show_heatmap_legend = T,row_names_side = "left", show_column_names = FALSE,
          row_split = c(rep("",10),rep(" ",10),rep("  ",10),rep("   ",10)),row_gap = unit(1, "mm"),
          row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 10))
draw(p, annotation_legend_side = 'left',heatmap_legend_side="left")
dev.off()

rm(p,column_ha,mat,common.gene,group.list)







##############################################
#validation of representative pathways in each CG
############################################
#DEG, GSEA of each CG group
deg.del7=FindMarkers(amlnew,ident.1 = "Del7.7q",min.pct = 0,logfc.threshold = 0)
run_fgsea <- function(deg,gene.path,out.path,name){
  ranks <- setNames(deg$avg_log2FC, rownames(deg))
  geneset=load_genesets(gene.path)
  fgseaRes <- fgsea(geneset, ranks, minSize=1, maxSize=500)
  fgseaRes=as.data.frame(fgseaRes)
  fgseaRes=fgseaRes[fgseaRes$padj<0.25,]
  fgseaRes=fgseaRes[order(fgseaRes$NES,decreasing = T),]
  
  for (pathway in fgseaRes$pathway){
    p=plotEnrichment(geneset[[pathway]],ranks) + labs(title=pathway)+
      theme(axis.title = element_text(size = 25),
            axis.text.x = element_text(size = 26,angle = 0,hjust = 0.5),axis.title.x = element_blank(),
            axis.text.y = element_text(size = 32),axis.ticks = element_blank()) 
    png(paste0(out.path,name,"_other_",pathway,".png"),res=200,width=13, height=9,units="in")
    print(p)
    dev.off()
  }
  
  for (i in 1:nrow(fgseaRes)){
    fgseaRes$leadingEdge[i]=paste(fgseaRes$leadingEdge[[i]],collapse = ",")
  }
  fgseaRes$leadingEdge=sapply(fgseaRes$leadingEdge,as.character)
  if (dim(fgseaRes)[1]>0){
    write.csv(as.matrix(fgseaRes),file=paste0(out.path,"/gsea_",name,"_other",".csv"),row.names = FALSE)
  }
}
run_fgsea(deg.diploid,"/Users/bwang8/Research/AML/info_files/AUCell_geneset/c2.cp.reactome.v2022.1.Hs.symbols.gmt.txt","Research/manuscript/cell_of_origin/Figure1_new/GSEA/diploid/Reactome/","diploid")




############################################
#ssgsea score bulk data, select common genes and convert to tpm. May lose some genes but those are not important

gsvaScore.bulk=data.frame(matrix(ncol = 0, nrow = 970))
for (i in c("h.all.v7.5.1.symbols.gmt.txt","c2.cp.reactome.v2022.1.Hs.symbols.gmt.txt","c5.go.v7.5.1.symbols.gmt.txt")){
  geneset=load_genesets(paste0("/Users/bwang8/Research/AML/info_files/AUCell_geneset/",i))
  gsva_score=gsva(as.matrix(bulk.tpm), geneset,method="ssgsea",abs.ranking=F)
  gsva_score=as.data.frame(t(gsva_score))
  gsvaScore.bulk=cbind(gsvaScore.bulk,gsva_score)
}
rm(i,geneset,gsva_score)

############################################
#bulk pathway by CG group




bulk.hallmark=gsvaScore.bulk %>% dplyr::select(contains("HALLMARK"))
bulk.hallmark$Sample=rownames(bulk.hallmark)

bulk.hallmark=merge(bulk.meta[,c('Sample',"CG_group",'Sex')],bulk.hallmark,by='Sample')

gm=groupMeans(t(bulk.hallmark[,4:53]),bulk.hallmark$CG_group)
gm=gm[,c("Diploid","Del5.5q",'Del7.7q',"Double del","Other")];gm=as.data.frame(scale(t(gm)))
p=pheatmap(gm, cluster_rows=F, show_rownames=T,cluster_cols=F,show_colnames = T,scale="none",angle_col=90, fontsize_row=20,fontsize_col = 12,cellheight=20, cellwidth = 15)
png("Research/manuscript/cell_of_origin/Figure1_suppl/bulk_ssGSEA_cyto.png",res=200,width=14, height=7,units="in")
print(p)
dev.off()

bulk.hallmark$CG_group=factor(bulk.hallmark$CG_group,levels = c('Diploid', 'Del5.5q', 'Del7.7q', "Double del",'Other'))
p=ggplot(bulk.hallmark,aes_string(x="CG_group", y="HALLMARK_INTERFERON_GAMMA_RESPONSE",fill="CG_group"))+
  geom_violin(violinwidth=1) + geom_boxplot(width=0.3, color="grey", alpha=0.2) + theme_classic() +
  scale_fill_manual(values=c("indianred2", "chartreuse1", "gold2","deepskyblue","turquoise4")) +
  theme(legend.position = "none", axis.title = element_text(size = 25),
        axis.text.x = element_text(size = 26,angle = 0,hjust = 0.5),axis.title.x = element_blank(),
        axis.text.y = element_text(size = 32),axis.ticks = element_blank()) 
kruskal.test(bulk.hallmark[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]~bulk.hallmark$CG_group)
png("Research/manuscript/cell_of_origin/Figure1_suppl/bulk_ssGSEA_cyto_boxplot.png",res=200,width=9, height=10,units="in")
p+annotate(geom="text",x=3,y=0.58,label="kruskal-Wallis, p = 0.16",cex=8)
dev.off()
rm(bulk.hallmark,gm,p)


############################################
#expression of genes in bulk cohort; may exclude "Other"
expri=bulk.tpm[bulk.tpm$X %in% c("IFNG","CD74","BST2","IRF2","IRF1","IFITM1","IFNGR1","IFITM3"),]
expri=as.data.frame(t(expri));colnames(expri)=expri[1,]
expri=expri[2:971,]; expri$Sample=rownames(expri)
expri[,1:8] = as.data.frame(sapply(expri[1:8], as.numeric))
expri=merge(expri,bulk.meta,by="Sample")
expri$CG_group=factor(expri$CG_group,levels=c("Diploid","Del5.5q","Del7.7q","Double del","Other"))
expri=expri[order(match(expri$CG_group,c("Diploid","Del5.5q","Del7.7q","Double del","Other"))),]
#expri[expri$CG_group=="Diploid" & expri$NPM1_mut=="Positive","CG_group"]="Diplod(NPM1+)" #darkorchid"
p=ggplot(expri,aes_string(x="CG_group", y="BST2",fill="CG_group"))+
  geom_violin() + geom_boxplot(width=0.3, color="grey", alpha=0.2) + theme_classic() +
  scale_fill_manual(values=c("indianred2", "chartreuse1", "gold2","deepskyblue","turquoise4")) +
  theme(legend.position = "none", axis.title = element_text(size = 25),
        axis.text.x = element_text(size = 26,angle = 0,hjust = 0.5),axis.title.x = element_blank(),
        axis.text.y = element_text(size = 32),axis.ticks = element_blank()) 
kruskal.test(expri[["BST2"]]~expri$CG_group)
png("Research/manuscript/cell_of_origin/Figure1_suppl/bulk_3set_BST2_cyto.png",res=200,width=9, height=10,units="in")
p+annotate(geom="text",x=3,y=2000,label="kruskal-Wallis, p = 0.078",cex=8)
dev.off()

column_ha=HeatmapAnnotation(df=data.frame(NPM1=factor(expri$NPM1_mut,levels=c("Positive","Negative","NA")),
                                          FLT3.ITD=factor(expri$FLT3ITD_mut,levels=c("Positive","Negative","NA")),
                                          RUNX1=factor(expri$RUNX1_mut,levels=c("Positive","Negative","NA")),
                                          KIT=factor(expri$KIT_mut,levels=c("Positive","Negative","NA")),
                                          DNMT3A=factor(expri$DNMT3A_mut,levels=c("Positive","Negative","NA")),
                                          FLT3=factor(expri$FLT3_mut,levels=c("Positive","Negative","NA")),
                      CG=factor(expri$CG_group,levels = c("Diploid","Del5.5q","Del7.7q","Double del","Other"))),
                            annotation_name_side = "left",
                     annotation_name_gp= gpar(fontsize = 12),show_legend = c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE),
              annotation_legend_param = list(NPM1=list(title="Mutation",at = c("Positive", "Negative","NA"),labels = c("True", "False","NA"))),
                            col = list(NPM1 = c("Positive"="tomato","Negative"="steelblue2","NA"="gray"),
                                      FLT3.ITD = c("Positive"="tomato","Negative"="steelblue2","NA"="gray"),
                                      RUNX1 = c("Positive"="tomato","Negative"="steelblue2","NA"="gray"),
                                      DNMT3A = c("Positive"="tomato","Negative"="steelblue2","NA"="gray"),
                                      KIT = c("Positive"="tomato","Negative"="steelblue2","NA"="gray"),
                                      FLT3 = c("Positive"="tomato","Negative"="steelblue2","NA"="gray"),
                                      CG = c("Double del"="deepskyblue","Del7.7q"="gold2","Del5.5q"="chartreuse1","Diploid"="indianred2","Other"="turquoise4")))

rm(expri,p)












############################################
#bulk immune response
#CYT:cytolytic activity, Rooney et al, Cell, 2015; TLS:tertiary lymphoid structures, Cabrita et al., Nature, 2020;
#Chemokines: Chemokine signature, Messina et al., Nat. Sci. Rep., 2012
#Davoli_IS: Davoli immune signature, Davoli et al., Science 2017; 
#IFNy: IFNy  signature, Ayers_exp: IS Expanded immune signature, Tcell_inflamed: T-cell inflamed signature Ayers et al., JCI, 2017
#Roh_IS: Roh immune score, Roh et al., Sci. Transl. Med., 2017
#RIR: Repressed immune resistance, Jerby-Arnon et al., Cell, 2018
hallmarks_immuneResponse=c("CYT","Roh_IS","chemokines","Davoli_IS","IFNy","Ayers_expIS","Tcell_inflamed","RIR","TLS")
immune_response_scores=compute_scores_immune_response(RNA_tpm = bulk.tpm,selected_scores = hallmarks_immuneResponse)
cell_fractions <- compute_cell_fractions(RNA_tpm=bulk.tpm)
tf_activity <- compute_TF_activity(RNA_tpm=bulk.tpm)
pathway_activities <- compute_pathways_activity(RNA.counts=counts, remove_sig_genes_immune_response = TRUE)

tf_activity$Sample=rownames(tf_activity)
tf_activity=merge(bulk.meta[,c("Sample","CG_group")],tf_activity, by="Sample")
gm=groupMeans(t(tf_activity[,2:119]),)
gm=gm[,c('Diploid',"Del5.5q","Del7.7q","Double del")]
#code for heatmap is the same as SCENIC; heatmap for immune_response_scores is the same as tf_activity
#stacked barplot for cell fractions
df=cell_fractions[cell_fractions$CG_group != "Other",]
df=reshape(df, direction = "long", varying = list(names(df)[3:13]),v.names = "Percentage",
           idvar = c("Sample","CG_group"),timevar = "celltype",times = colnames(df)[3:13])

ggplot(df, aes_string(fill="celltype",y="Percentage",x="Sample")) + 
  theme_classic()+xlab("")+ylab("Percentage") + geom_bar(position="stack",stat="identity") + 
  scale_fill_manual(values=c(brewer.pal(n = 10, name = 'Set3'),"red"))+ theme_bw()+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.15))+ ylab('')+
  theme(axis.text.x=element_blank(), axis.title.y = element_text(size = 28), axis.ticks = element_blank(),
        axis.text.y = element_blank(),legend.position = "top",
        panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + 
  annotate(ymin = 1.005, ymax = 1.15,xmin = c(-Inf,294,330,364),xmax = c(294,330,364,Inf),geom = "rect",
           fill = c("indianred2","darkseagreen3", "gold2","deepskyblue")) #, 

rm(hallmarks_immuneResponse,immune_response_scores,tf_activity,cell_fractions,pathways_activity)


#deconvolution of bulk with LM22
treg.percent=(bulk.lm22 %>% group_by(CG_group) %>% summarise(count = sum(T.cells.regulatory..Tregs. >0.01)))$count / (bulk.lm22 %>% group_by(CG_group)  %>% summarise(count=n()))$count
treg.percent=data.frame(cg=c("Diploid","Del5.5q","Del7.7q","Double del"),treg=treg.percent)
treg.percent$cg=factor(treg.percent$cg,levels=c("Diploid","Del5.5q","Del7.7q","Double del"))
png("Research/manuscript/cell_of_origin/Figure1_suppl/bulk_3set_LM22_Treg.png",res=200,width=9, height=10,units="in")
ggplot(data=treg.percent, aes(x=cg, y=treg,fill=cg)) + geom_bar(stat="identity") + theme_classic() +
  scale_fill_manual(values=c("indianred2", "darkseagreen3", "gold2","deepskyblue")) +
  ylab("Percent of samples with Treg > 1%") +
  theme(axis.text = element_text(size=26),axis.title.y = element_text(size=28),legend.position = "none",  
        axis.title.x = element_blank(),legend.text = element_text(size=22))
dev.off()




############################################
#Nanostring STM/JCI paper
plot_vln <- function(gene){
  p=ggplot(nano.expri,aes_string(x="CG_group", y=gene,fill="CG_group"))+
    geom_violin(violinwidth=1) + geom_boxplot(width=0.1, color="grey", alpha=0.2) + theme_classic() +
    scale_fill_manual(values=c("indianred2", "chartreuse1", "gold2","deepskyblue","turquoise4")) +
    theme(legend.position = "none", axis.title = element_text(size = 25),
          axis.text.x = element_text(size = 26,angle = 0,hjust = 0.5),axis.title.x = element_blank(),
          axis.text.y = element_text(size = 32),axis.ticks = element_blank()) 
  return(p)
}
vln_plots <- map(c("IFNG","CD74","BST2","IRF2","IRF1","IFITM1","IFNGR1"), plot_vln)#"IFITM3",
pdf("Research/manuscript/cell_of_origin/Figure1_suppl/nano_gene_box.pdf",width = 10,height=10)
print(ggarrange(plotlist=vln_plots, ncol=1,nrow=1, legend="none"))
dev.off()

#plot heatmap of genes with clinical information as row anotation
nano.expri=as.data.frame(nano.expri);rownames(nano.expri)=nano.expri$Patient_Id
expri=nano.expri[,c("Patient_Id","IFNG","CD74","BST2","IRF2","IRF1","IFITM1","IFNGR1")]
expri=expri[expri$Patient_Id %in% nano.meta$Patient_Id,]
nano.meta=nano.meta[order(match(nano.meta$CG_group,c("Diploid","Del5.5q","Del7.7q","Double del","Other"))),]
expri=expri[order(match(expri$Patient_Id,nano.meta$Patient_Id)),]

column_ha=HeatmapAnnotation(df=data.frame(NPM1=factor(nano.meta$NPM1,levels=c(1,0)),
                                          FLT3.ITD=factor(nano.meta$`FLT3-ITD`,levels=c(1,0)),
                                          JAK2=factor(nano.meta$JAK2,levels=c(1,0)),
                                          KIT=factor(nano.meta$KIT,levels=c(1,0)),
                                          PML=factor(nano.meta$PML,levels=c(1,0)),
                                          FLT3.TKD=factor(nano.meta$`FLT3-TKD`,levels=c(1,0)),
              CG=factor(nano.meta$CG_group,levels = c("Diploid","Del5.5q","Del7.7q","Double del","Other"))),
                            annotation_name_side = "left",gp = gpar(col=c("white")),
                  annotation_name_gp= gpar(fontsize = 12),show_legend = c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE),
            annotation_legend_param = list(NPM1=list(title="Mutation",at = c("1", "0"),labels = c("True", "False"))),
                            col = list(NPM1 = c("1"="tomato","0"="steelblue2"),
                      FLT3.ITD = c("1"="tomato","0"="steelblue2"),PML = c("1"="tomato","0"="steelblue2"),
                        JAK2 = c("1"="tomato","0"="steelblue2"),KIT = c("1"="tomato","0"="steelblue2"),
                      FLT3.TKD = c("1"="tomato","0"="steelblue2"),
            CG = c("Double del"="deepskyblue","Del7.7q"="gold2","Del5.5q"="chartreuse1","Diploid"="indianred2","Other"="turquoise4")))

p=Heatmap(t(scale(expri[,2:8])),name="Expression",col=colorRamp2(seq(-4,4,0.1), colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(81)),
          top_annotation = column_ha, cluster_columns=FALSE, cluster_rows = FALSE, 
          column_title=NULL,show_heatmap_legend = T,row_names_side = "left", show_column_names = FALSE,
          row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 10),width = ncol(expri[,2:8])*unit(25, "mm"), height = nrow(expri[,2:8])*unit(0.5, "mm"))
png("Research/manuscript/cell_of_origin/Figure1_suppl/nano_gene_onco.png",res=200,width=9.8, height=5.5,units="in")
draw(p, annotation_legend_side = 'right',heatmap_legend_side="right")
dev.off()


#scatter plot between CD74 and BST2
png("Research/manuscript/cell_of_origin/Figure1_suppl/nano_CD74_BST2_cor.png",res=200,width=10, height=10,units="in")
ggscatter(nano.expri, x="BST2", y="CD74",add="reg.line",conf.int = TRUE,cor.coef = TRUE, cor.coef.size = 8,
          cor.method = "pearson",xlab = "BST2", ylab = "CD74",size=2) +
  theme(legend.position = "none", axis.title = element_text(size = 30),
        axis.text = element_text(size = 28),axis.ticks = element_blank())
dev.off()

rm(nano.cyto,gene,vln_plots,plot_vln,p,column_ha,nano.cyto)

############################################
#validate in Other scRNA AML data
vangalen2019=readRDS("Research/AML/info_files/VanGalen_Cell2019/RDS_and_MD/D0/vangalen_d0_aml.rds")
vangalen2019$cyto=ifelse(vangalen2019$haa_orig.ident=="AML328-D0","Del7.7q","Other")
vangalen2019$cyto1=ifelse(vangalen2019$haa_orig.ident %in% c("AML210A-D0","AML419A-D0","AML916-D0","AML329-D0","AML556-D0","AML475-D0","AML921A-D0"),"Diploid","Other")
#vangalen2019=subset(vangalen2019,subset = haa_orig.ident %in% c("AML3210A-D0","AML419A-D0","AML916-D0","AML921-D0","AML314--D0","AML475-D0","AML997-D0","AML328-D0",'AML329-D0',"AML556-D0")) #only del7 and diploid
p=ggplot(vangalen2019@meta.data,aes_string(x="cyto", y="HALLMARK_INTERFERON_GAMMA_RESPONSE",fill="cyto"))+
  geom_violin(violinwidth=1) + geom_boxplot(width=0.15, color="grey", alpha=0.2) + theme_classic() +
  scale_fill_manual(values=c("gold2","turquoise4")) +
  theme(legend.position = "none", axis.title = element_text(size = 26),
        axis.text.x = element_text(size = 30,angle = 0,hjust = 0.5),axis.title.x = element_blank(),
        axis.text.y = element_text(size = 32),axis.ticks = element_blank()) 
vangalen2019@meta.data %>% group_by("cyto") %>% pairwise_wilcox_test(HALLMARK_INTERFERON_GAMMA_RESPONSE~cyto, p.adjust.method = "BH") %>% add_significance("p.adj") %>% add_xy_position(x="cell_type", fun="max", step.increase = 0.1)
png("Research/manuscript/cell_of_origin/Figure1_suppl/vangalenPre_IFNGscore_box.png",res=200,width=9, height=10,units="in")
p+annotate(geom="text",x=1.5,y=0.3,label="Wilcoxon test, p = 1.49e-50",cex=8)
dev.off()

p=VlnPlot(vangalen2019,features = "BST2",group.by = "cyto1",pt.size=0) + 
  geom_boxplot(width=0.15, color="grey", alpha=0.2) + ggtitle("") + labs(y="BST2") +
  theme_classic() + scale_fill_manual(values=c("indianred2","turquoise4")) + xlab("BST2")+
  theme(legend.position = "none", axis.title = element_text(size = 26),
        axis.text.x = element_text(size = 30,angle = 0,hjust = 0.5),axis.title.x = element_blank(),
        axis.text.y = element_text(size = 32),axis.ticks = element_blank()) 
###vanGalen post treatment
vangalen.alldata=readRDS("Research/AML/info_files/VanGalen_Cell2019/RDS_and_MD/vanGalen_Allsample.rds")
vangalen.post=subset(vangalen.alldata,subset = orig.ident %in% c("AML314-D31","AML371-D34","AML475-D29","AML997-D35","AML556--D31","AML707B-D41","AML707B-D97","AML328-D29","AML328-D113","AML328-D171"))  #only choose del7 and diloid
vangalen.post$cyto=ifelse(vangalen.post$orig.ident %in% c("AML314-D31","AML371-D34","AML475-D29","AML997-D35","AML556--D31","AML707B-D41","AML707B-D97"),"Diploid","Del7.7q")
vangalen.post=subset(vangalen.post,subset = CellType %in% c("cDC-like","GMP-like","HSC-like","Mono-like","Prog-like","ProMono-like")) #only keep AML cells
AML328=subset(vangalen.alldata,subset = orig.ident %in% c("AML328-D0","AML328-D29","AML328-D113","AML328-D171")) 


#Wu et al., Journal of Hematology and Oncology 40 patients
#wu2020=readRDS("/Users/bwang8/Research/AML/info_files/public_scAML/Wu2020_JHO.rds")
petti2022=readRDS("Research/AML/info_files/public_scAML/petti_processed.rds")






#------------------------------------------------------------------------------
############################################################################
#Figure : TCGA,BEATAML, Abbas deconvolution and mutation analysis 
#         
############################################################################
#------------------------------------------------------------------------------
bulk.deconv=read.csv("Research/manuscript/cell_of_origin/Figure4_new/deconvolution/Bulk3sets_deconv_class3.csv")
colnames(bulk.deconv)[1]="Sample"
bulk.deconv$Sample[762:934]=substr(bulk.deconv$Sample[762:934],1,12)
bulk.deconv=merge(bulk.deconv[,1:10],bulk.meta,by="Sample")
fraction=bulk.deconv %>% group_by(CG_group,celltype)  %>% summarise(count=n()) %>% mutate(Percent = round(count/sum(count),4))
#fraction$celltype=factor(fraction$celltype,levels = c("HSC-like","CMP.LMPP-like","CLP-like","GMP-like","lympho-like", "Mono-like", "Baso-like","Erythroid-like","DC-like"))
fraction$celltype=factor(fraction$celltype,levels = c("Primitive-like","GMP-like","Lympho-like","Committed-like", "Erythroid-like"))
fraction$CG_group=factor(fraction$CG_group,levels = c("Diploid","Del5.5q","Del7.7q","Double del","Other"))
#fraction$celltype=factor(fraction$celltype,levels = c("LSPC-Quiescent","LSPC-Primed","LSPC-Cycle","GMP-like","ProMono-like", "Mono-like","cDC-like"))
p=ggplot(fraction, aes_string(fill="celltype",y="Percent",x="CG_group")) +
  theme_classic()+xlab("")+ylab("Percentage") + geom_bar(position="stack",stat="identity") + 
  scale_fill_manual(values=brewer.pal(n = 10, name = 'Set3')[c(1,4,5,6,8)])+
  #scale_x_discrete(labels=c("Diploid","Del5.5q","Del7.7q","Double del","Other")) +
  theme(axis.text.x=element_text(angle=30, vjust = 0.8,hjust=0.8, size=20), #labs(fill='Cell type') +
        axis.title.y = element_text(size = 28),axis.text.y = element_text(size = 22),
        legend.text=element_text(size=22),legend.title=element_blank()) 

temp=as.data.frame(reshape(as.data.frame(bulk.deconv %>% group_by(CG_group,celltype)  %>% summarise(count=n())),timevar = "celltype",idvar = "CG_group",direction = "wide"))
temp[is.na(temp)]=0;chisq.test(temp[,2:6])
png("Research/manuscript/cell_of_origin/Figure4_new/deconvolution/composition_bulk_cyto3_class3.png",res=200,width=9, height=10,units="in")
p+annotate(geom="text",x=3,y=1.05,label="Chi-squar, p = 6.891e-05",cex=8) 
dev.off()

#####stacked barplot show relative abundance of each cell type in each patient, as Andy's figure 2
df=bulk.deconv[,c(1:6,10,533)]  #need to order the CG_group of bulk_deconv
df=df[df$CG_group=="Other",]
df=reshape(df, direction = "long", varying = list(names(df)[2:6]),v.names = "Percentage",
           idvar = c("Sample","CG_group"),timevar = "celltype",times = colnames(df)[2:6])
df$celltype=str_replace(df$celltype, "[.]", "-")
#df=df[order(match(df$Sample,bulk.deconv$Sample)),]
df$celltype=factor(df$celltype,levels = c("Primitive-like","GMP-like","Lympho-like","Committed-like", "Erythroid-like"))
png("Research/manuscript/cell_of_origin/Figure4_new/deconvolution/composition_bulk_other_class3_patient.png",res=200,width=12, height=5,units="in")
ggplot(df, aes_string(fill="celltype",y="Percentage",x="Sample")) + 
  theme_classic()+xlab("")+ylab("Percentage") + geom_bar(position="stack",stat="identity") + 
  scale_fill_manual(values=brewer.pal(n = 10, name = 'Set3')[c(1,4,5,6,8)])+ theme_bw()+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.15))+ ylab('')+
  theme(axis.text.x=element_blank(), axis.title.y = element_text(size = 28), axis.ticks = element_blank(),
        axis.text.y = element_blank(),legend.position = "none",
        panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + 
  annotate(ymin = 1.005, ymax = 1.15,xmin = c(-Inf),xmax = c(Inf),geom = "rect",
           fill = c("turquoise4")) #, "indianred2","darkseagreen3", "gold2","deepskyblue",
dev.off()

rm(fraction,p,temp,df)

#####analyze clinical characteristics related to celltypes
trial2=bulk.deconv # TCGA: trial2=bulk.deconv[520:672,]
#colnames(trial2)[17]="OS_days";colnames(trial2)[22]="Age"
trial2[trial2$label=="","label"]=NA;trial2[trial2==""]=NA
trial2 <- trial2 %>% dplyr::select(celltype,Age,OS_days,WBC_Count,ELN2017_incomplete, label,FLT3_mut,FLT3ITD_mut,NPM1_mut,RUNX1_mut,IDH1_mut,IDH2_mut,TP53_mut)
colnames(trial2)[5:6]=c("ELN2017","Karyotype")
trial2$ELN2017[trial2$ELN2017 %in% c("IntermediateOrAdverse","NotInitial","FavorableOrIntermediate","NonAML")]="Other"
fisher.test.simulate.p.values <- function(data, variable, by, ...) {
  result <- list()
  test_results <- stats::fisher.test(data[[variable]], data[[by]], simulate.p.value = TRUE)
  result$p <- test_results$p.value
  result$test <- test_results$method
  result
}
trial2 %>% tbl_summary(by = celltype,
                       statistic = list(all_continuous() ~ "{mean} ({sd})",all_categorical() ~ "{n} / {N} ({p}%)"),
                       digits = all_continuous() ~ 2,
                       #label = c(Cytogenetic.risk	~ "Cytogenetic risk",IsDeNovo.~"DeNovo",Not_Complex.~"Phenotype"),
                       missing_text = "(Missing)") %>% add_p(test = list(all_categorical() ~ "fisher.test.simulate.p.values")) # %>% add_p()

#trial2 <- trial2 %>% select(celltype,IDH2_mut,TP53_mut,NPM1_mut,RUNX1_mut)#TCGA
#select(celltype,CEBPA_mut,FLT3_mut,FLT3ITD_mut,NPM1_mut,JAK2_mut,RUNX1_mut)#beataml
#select(celltype,IDH1_mut,RUNX1_mut)#MDA

trial2=trial2[,c(10,73:111)]
trial2[trial2==""]=NA
#Multi-variate model
trial2 <- trial2 %>% dplyr::select(celltype,Age,OS_days,Dead_Alive,ELN2017_incomplete,CEBPA_mut,DNMT3A_mut,FLT3_mut,FLT3ITD_mut,NPM1_mut,NRAS_mut,RUNX1_mut,IDH1_mut,IDH2_mut,TP53_mut)
trial2$Age_group <- ifelse(trial2$Age > 60, "Greater", "Less")
trial2$Age_group <- factor(trial2$Age_group, levels=c("Less", "Greater"))
trial2$Dead_Alive <- ifelse(trial2$Dead_Alive == "Dead", 1, 0)
colnames(trial2)[5]="ELN2017"
trial2$ELN2017=factor(trial2$ELN2017,levels = c("Favorable","Intermediate","Adverse","Other",NA ))
covariate_names <- c(`Age_group:Greater` ="Age > 60", `ELN2017:Adverse`="ELN Adverse", 
                     `ELN2017:Favorable`="ELN Favorable",`ELN2017:Intermediate`="ELN Intermediate",
                     `FLT3ITD_mut:Positive` ="FLT3ITD mutation",`IDH1_mut:Positive` ="IDH1 mutation",`NPM1_mut:Positive`="NPM1 mutation",`FLT3_mut:Positive` ="FLT3 mutation",`RUNX1_mut:Positive` ="RUNX1 mutation",`IDH2_mut:Positive` ="IDH2 mutation",`TP53_mut:Positive` ="TP53 mutation")

trial2 %>% analyse_multivariate(vars(OS_days, Dead_Alive),
                                covariates = vars(Age_group, ELN2017,IDH1_mut,IDH2_mut,FLT3_mut,FLT3ITD_mut,NPM1_mut,RUNX1_mut,TP53_mut),covariate_name_dict = covariate_names) -> result


forest_plot(result,factor_labeller = covariate_names,
            endpoint_labeller = c(time="OS_days"),orderer = ~order(HR),
            labels_displayed = c("factor"),
            ggtheme = ggplot2::theme_bw(base_size = 20),
            relative_widths = c(0.8, 2, 0.8),
            HR_x_breaks = c(0, 0.5, 1, 2, 3))
dev.off()
###Multiple Univariate Analyses
map(vars(Age_group, Cytogenetic.risk, WBC_Count,IDH2_mut,NPM1_mut,TP53_mut,CEBPA_mut,FLT3_mut), 
    function(by){analyse_multivariate(trial2,vars(OS_days, Dead_Alive.),
                                      covariates = list(by), # covariates expects a list
                                      covariate_name_dict = covariate_names)
    }) %>% forest_plot(factor_labeller = covariate_names,
                       endpoint_labeller = c(time="OS_days"),orderer = ~order(HR),labels_displayed = c("factor"),
                       ggtheme = ggplot2::theme_bw(base_size = 20),relative_widths = c(0.6, 2, 0.8),)


rm(trial2,covariate_names,result)

##########################
#heatmap of RUNX1 positive and negative patients, showing other mutations by celltype
runx1.pos=bulk.deconv[bulk.deconv$RUNX1_mut=="Positive",];runx1.neg=bulk.deconv[bulk.deconv$RUNX1_mut=="Negative",]
runx1.pos=runx1.pos[,c(1,10,13,73:111)];runx1.neg=runx1.neg[,c(1,10,13,73:111)]
runx1.pos=runx1.pos[order(match(runx1.pos$celltype,c("Primitive-like","GMP-like","Lympho-like", "Committed-like", "Erythroid-like"))),]
runx1.neg=runx1.neg[order(match(runx1.neg$celltype,c("Primitive-like","GMP-like","Lympho-like", "Committed-like", "Erythroid-like"))),]
column_ha=HeatmapAnnotation(df=data.frame(
  cell_type=factor(runx1.pos[,2],levels = c("Primitive-like","GMP-like","Lympho-like", "Committed-like", "Erythroid-like"))),annotation_name_side = "left",annotation_name_gp= gpar(fontsize = 10),show_legend = c(TRUE),
  col = list(cell_type=c("Primitive-like"="#8DD3C7","GMP-like"="#FB8072","Lympho-like"="#80B1D3", "Committed-like"="#FDB462", "Erythroid-like"="#FCCDE5")),height=unit(600, "mm"))
pdf("Research/manuscript/intratumor_heterogeneity/ITH_figure1/deconvolution/RUNX1_pos_mut.pdf",width=10, height=10)
Heatmap(as.matrix(t(runx1.pos[,4:41])),name="Mutation",top_annotation = column_ha,
        col=structure(c('firebrick1',"steelblue1","white"), names = c("Positive","Negative","NA")),
        cluster_columns=FALSE, cluster_rows = FALSE, column_names_rot = 90, #width = ncol(gm)*unit(20, "mm"),
        column_title=NULL,show_heatmap_legend = T,row_names_side = "left", show_column_names = FALSE,
        #heatmap_legend_param = list(title = ""),
        height = unit(150, "mm"), width = unit(200, "mm"),
        row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 11))
dev.off()















