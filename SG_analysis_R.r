library(Seurat)
library(dplyr)
library(harmony)
library(scDblFinder)
library(SeuratDisk)
library(BiocParallel)
library(ggplot2)

library(future)
plan(strategy = "multicore", workers = 4)

XNCX_1 = Read10X_h5('/home/qukun/xuhao/workspace/Sjs/cellbender/outs/post_cellbender_XNCX_1_FPR_0.01_filtered.h5')
XNCX_2 = Read10X_h5('/home/qukun/xuhao/workspace/Sjs/cellbender/outs/post_cellbender_XNCX_2_FPR_0.01_filtered.h5')
XNCX_3 = Read10X_h5('/home/qukun/xuhao/workspace/Sjs/cellbender/outs/post_cellbender_XNCX_3_FPR_0.01_filtered.h5')
XNCX_4 = Read10X_h5('/home/qukun/xuhao/workspace/Sjs/cellbender/outs/post_cellbender_XNCX_4_FPR_0.01_filtered.h5')
XNCX_5 = Read10X_h5('/home/qukun/xuhao/workspace/Sjs/cellbender/outs/post_cellbender_XNCX_5_FPR_0.01_filtered.h5')
XNCX_6 = Read10X_h5('/home/qukun/xuhao/workspace/Sjs/cellbender/outs/post_cellbender_XNCX_6_FPR_0.01_filtered.h5')
XNCX_7 = Read10X_h5('/home/qukun/xuhao/workspace/Sjs/cellbender/outs/post_cellbender_XNCX_7_FPR_0.01_filtered.h5')
XNCX_8 = Read10X_h5('/home/qukun/xuhao/workspace/Sjs/cellbender/outs/post_cellbender_XNCX_8_FPR_0.01_filtered.h5')
XNCX_9 = Read10X_h5('/home/qukun/xuhao/workspace/Sjs/cellbender/outs/post_cellbender_XNCX_9_FPR_0.01_filtered.h5')
XNCX_10 = Read10X_h5('/home/qukun/xuhao/workspace/Sjs/cellbender/outs/post_cellbender_XNCX_10_FPR_0.01_filtered.h5')
XNCX_11 = Read10X_h5('/home/qukun/xuhao/workspace/Sjs/cellbender/outs/post_cellbender_XNCX_11_FPR_0.01_filtered.h5')
HCCX_1 = Read10X_h5('/home/qukun/xuhao/workspace/Sjs/cellbender/outs/post_cellbender_HCCX_1_FPR_0.01_filtered.h5')
HCCX_2 = Read10X_h5('/home/qukun/xuhao/workspace/Sjs/cellbender/outs/post_cellbender_HCCX_2_FPR_0.01_filtered.h5')
HCCX_3 = Read10X_h5('/home/qukun/xuhao/workspace/Sjs/cellbender/outs/post_cellbender_HCCX_3_FPR_0.01_filtered.h5')
HCCX_4 = Read10X_h5('/home/qukun/xuhao/workspace/Sjs/cellbender/outs/post_cellbender_HCCX_4_FPR_0.01_filtered.h5')
HCCX_5 = Read10X_h5('/home/qukun/xuhao/workspace/Sjs/cellbender/outs/post_cellbender_HCCX_5_FPR_0.01_filtered.h5')

XNCX_1 = CreateSeuratObject(XNCX_1,project = 'XNCX_1')
XNCX_2 = CreateSeuratObject(XNCX_2,project = 'XNCX_2')
XNCX_3 = CreateSeuratObject(XNCX_3,project = 'XNCX_3')
XNCX_4 = CreateSeuratObject(XNCX_4,project = 'XNCX_4')
XNCX_5 = CreateSeuratObject(XNCX_5,project = 'XNCX_5')
XNCX_6 = CreateSeuratObject(XNCX_6,project = 'XNCX_6')
XNCX_7 = CreateSeuratObject(XNCX_7,project = 'XNCX_7')
XNCX_8 = CreateSeuratObject(XNCX_8,project = 'XNCX_8')
XNCX_9 = CreateSeuratObject(XNCX_9,project = 'XNCX_9')
XNCX_10 = CreateSeuratObject(XNCX_10,project = 'XNCX_10')
XNCX_11 = CreateSeuratObject(XNCX_11,project = 'XNCX_11')
HCCX_1 = CreateSeuratObject(HCCX_1,project = 'HCCX_1')
HCCX_2 = CreateSeuratObject(HCCX_2,project = 'HCCX_2')
HCCX_3 = CreateSeuratObject(HCCX_3,project = 'HCCX_3')
HCCX_4 = CreateSeuratObject(HCCX_4,project = 'HCCX_4')
HCCX_5 = CreateSeuratObject(HCCX_5,project = 'HCCX_5')

sample.list <- c(XNCX_1,XNCX_2,XNCX_3,XNCX_4,XNCX_5,XNCX_6,XNCX_7,XNCX_8,XNCX_9,XNCX_10,XNCX_11,HCCX_1,HCCX_2,HCCX_3,HCCX_4,HCCX_5)

for (i in seq_along(sample.list)){
    sample.list[[i]] <- NormalizeData(sample.list[[i]]) %>% FindVariableFeatures(nfeatures = 3000)
}

hvg_counts <- data.frame(row.names = rownames(HCCX_5))
hvg_counts$counts <- 0
for (i in seq_along(sample.list)){
   hvg_counts[sample.list[[i]]@assays$RNA@var.features, 'counts'] <- hvg_counts[sample.list[[i]]@assays$RNA@var.features, 'counts'] + 1
}
length(rownames(hvg_counts)[hvg_counts$counts > 5])

CX <- merge(XNCX_1,c(XNCX_2,XNCX_3,XNCX_4,XNCX_5,XNCX_6,XNCX_7,XNCX_8,XNCX_9,XNCX_10,XNCX_11,HCCX_1,HCCX_2,HCCX_3,HCCX_4,HCCX_5))

options(future.globals.maxSize= 3*1024*1024^2)

CX <- NormalizeData(CX)
CX@assays$RNA@var.features <- rownames(hvg_counts)[hvg_counts$counts > 5]
CX <- ScaleData(CX) %>% RunPCA(npcs=30)
CX[['percent.mt']] <- PercentageFeatureSet(CX,pattern = '^MT-')
CX <- subset(CX,subset=nCount_RNA > 200 & nFeature_RNA > 100)

CX.sce <- as.SingleCellExperiment(CX)
CX.sce <- scDblFinder(CX.sce, samples = 'orig.ident', clusters=FALSE, BPPARAM=MulticoreParam(4))
CX <- as.Seurat(CX.sce)
CX[['sample']] <- CX.sce[['orig.ident']]
CX <- CX[,CX@meta.data[['scDblFinder.class']] == 'singlet']
CX <- subset(CX,subset=nFeature_RNA > 500 & percent.mt < 25)

VlnPlot(CX, c('nFeature_RNA','nCount_RNA','percent.mt'),pt.size=0,group.by='sample')

sample.list <- SplitObject(CX,split.by='sample')

for (i in seq_along(sample.list)){
    sample.list[[i]] <- NormalizeData(sample.list[[i]]) %>% FindVariableFeatures(nfeatures = 3000)
}

hvg_counts <- data.frame(row.names = rownames(HCCX_5))
hvg_counts$counts <- 0
for (i in seq_along(sample.list)){
   hvg_counts[sample.list[[i]]@assays$RNA@var.features, 'counts'] <- hvg_counts[sample.list[[i]]@assays$RNA@var.features, 'counts'] + 1
}
length(rownames(hvg_counts)[hvg_counts$counts > 5])

CX <- NormalizeData(CX)
CX@assays$RNA@var.features <- rownames(hvg_counts)[hvg_counts$counts > 5]

CX <- ScaleData(CX) %>% RunPCA(npcs=30,verbose=F) %>% RunUMAP(dims=1:30)
CX <- RunHarmony(CX, group.by.vars='sample')
CX <- RunUMAP(CX,dims=1:30,reduction='harmony')
CX <- FindNeighbors(CX, reduction='harmony', dims=1:30)
CX <- FindClusters(CX, resolution=2)

DimPlot(CX,group.by='sample')

DimPlot(CX,label=T)

allmarkers <- FindAllMarkers(CX)

allmarkers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC) -> top5

options(repr.plot.width = 25, repr.plot.height = 10, repr.plot.res = 150)
DotPlot(CX,feature=unique(top5$gene), group.by='seurat_clusters') + RotatedAxis()

celltype_submarkers_dict=c(
    '12'='Serous_acini_ANKRD36C_SRRM2_GOLGB1',
    '14'='Serous_acini_JUN_FOSB_ELF3',
    '11'='Serous_acini_SCGB3A1_WFDC2_S100A6_KRT19',
    '32'='Serous_acini_TAGLN_CST3_TPM2_ACTA2',
    '21'='Serous_acini_SPARCL1_LRRC26',
    '7'='Serous_acini_MUC7_CRISP3_PIGR',
    '10'='Serous_acini_CLDN10_PRDX4_SEC11C',
    '5'='Serous_acini_PRR4_SMR3B_LYZ',
    '36'='Mucous_acini_SMR3B_PI3_ZG16B',
    '17'='Mucous_acini_FCGBP_MUC5B_SCGB3A1',
    '29'='Mucous_acini_LTF_KRT17_LCN2',
    '24'='Ductal_RARRES1_ALDH1A3_CCL28',
    '38'='Ductal_LYZ_PI3_PRR4',
    '30'='Ductal_IGFBP7_MIR205HG_FST',
    '35'='Ionocytes',
    '13'='Endothelial_RGCC_CA4_PIP',
    '22'='Endothelial_ACKR1_SELE_CCL14',
    '34'='Endothelial_RGS5_NDUFA4L2_IGFBP3',
    '2'='Myoepithelial',
    '33'='Melanocytes',
    '1'='Fibroblast_DNAJB1_SOD2_HSP90AB1_NFKBIA',
    '3'='Fibroblast_THBS4_SPON2_F2R',
    '9'='Fibroblast_PIP_APOE_LYZ',
    '26'='Fibroblast_MFAP5_PCOLCE2_FBN1',
    '31'='Fibroblast_CCL19_CXCL10_CXCL9',
    '16'='Myofibroblast',
    '27'='Pericytes',
    '37'='Unknown',
    '8'='Plasma_IGKV2D-28_IGKC_JCHAIN',
    '15'='Plasma_IGHG2_HSPA1A_HSPE1',
    '25'='Plasma_IER3_SAT1_HSPB1',
    '40'='Plasma_IGHV7-4-1_IGHA2_TXNDC5',
    '23'='Mac/mo_C3_FCGR3A_CLDN1_IL1B',
    '28'='Mac/mo_RNASE1_SELENOP_CCL18',
    '4'='CD8_T_ZG16B_HSPA6_LYZ',
    '18'='CD8_T_GIMAP7_GIMAP4_ACTB_GZMA',
    '20'='CD8_T_CCL4L2_CCL4_XCL2_GZMB',
    '6'='CD4_T_CXCL13_CREM_ICOS',
    '19'='CD4_T_CCL4_CCL5_CCL4L2',
    '0'='B',
    '39'='Mast'
)

Idents(CX) <- "seurat_clusters"
CX <- RenameIdents(CX, celltype_submarkers_dict)
CX$celltype_submarkers <- Idents(CX)

options(repr.plot.width = 25, repr.plot.height = 10, repr.plot.res = 150)
DotPlot(CX,feature=unique(top10$gene), group.by='celltype_submarkers') + RotatedAxis()

# Serous acini minor markers
options(repr.plot.width = 12, repr.plot.height = 9, repr.plot.res = 150)
DotPlot(CX,feature=c('ANKRD36C','SRRM2','GOLGB1','JUN','FOSB','ELF3','SCGB3A1','WFDC2','S100A6','KRT19','TAGLN','CST3','TPM2','ACTA2','SPARCL1','LRRC26','MUC7','CRISP3','PIGR','CLDN10','PRDX4','SEC11C','PRR4','SMR3B','LYZ'),group.by='celltype_submarkers') + RotatedAxis()

# Mucous acini minor markers
options(repr.plot.width = 8, repr.plot.height = 9, repr.plot.res = 150)
DotPlot(CX,feature=c('SMR3B','PI3','ZG16B','FCGBP','MUC5B','SCGB3A1','LTF','LCN2'),group.by='celltype_submarkers') + RotatedAxis()

# Ductal minor markers
options(repr.plot.width = 8, repr.plot.height = 9, repr.plot.res = 150)
DotPlot(CX,feature=c('RARRES1','ALDH1A3','CCL28','LYZ','PI3','PRR4','IGFBP7','MIR205HG','FST'),group.by='celltype_submarkers') + RotatedAxis()

# Endothelial minor markers
options(repr.plot.width = 8, repr.plot.height = 9, repr.plot.res = 150)
DotPlot(CX,feature=c('RGCC','CA4','PIP','ACKR1','SELE','CCL14','RGS5','NDUFA4L2','IGFBP3'),group.by='celltype_submarkers') + RotatedAxis()

# Fibroblast minor markers
options(repr.plot.width = 9, repr.plot.height = 9, repr.plot.res = 150)
DotPlot(CX,feature=c('DNAJB1','SOD2','HSP90AB1','NFKBIA','THBS4','SPON2','F2R','PIP','APOE','LYZ','MFAP5','PCOLCE2','FBN1','CCL19','CXCL10','CXCL9'),group.by='celltype_submarkers') + RotatedAxis()

# Plasma minor markers
options(repr.plot.width = 9, repr.plot.height = 9, repr.plot.res = 150)
DotPlot(CX,feature=c('IGHG1','IGKV2D-28','IGKC','JCHAIN','IGHG2','HSPA1A','HSPE1','IER3','SAT1','HSPB1','IGHV7-4-1','IGHA2','TXNDC5'),group.by='celltype_submarkers') + RotatedAxis()

# Mac/Mo minor markers
options(repr.plot.width = 9, repr.plot.height = 9, repr.plot.res = 150)
DotPlot(CX,feature=c('C3','FCGR3A','CLDN1','IL1B','RNASE1','SELENOP','CCL18'),group.by='celltype_submarkers') + RotatedAxis()

# T minor markers
options(repr.plot.width = 10, repr.plot.height = 9, repr.plot.res = 150)
DotPlot(CX,feature=c('ZG16B','HSPA6','LYZ','GIMAP7','GIMAP4','ACTB','GZMA','CCL4L2','CCL4','XCL2','GZMB','CXCL13','CREM','ICOS','CCL5','CD69','CXCR4','GPR183','CD2','SLFN5','MYL12A','HSPH1','HSPA1A'),group.by='celltype_submarkers') + RotatedAxis()

celltype_minor_dict=c(
    '12'='Serous_acini_ANKRD36C',
    '14'='Serous_acini_ELF3',
    '11'='Serous_acini_SCGB3A1',
    '32'='Serous_acini_TAGLN',
    '21'='Serous_acini_LRRC26',
    '7' ='Serous_acini_MUC7',
    '10'='Serous_acini_CLDN10',
    '5' ='Serous_acini_PRR4',
    '36'='Mucous_acini_PI3',
    '17'='Mucous_acini_FCGBP',
    '29'='Mucous_acini_LTF',
    '24'='Ductal_ALDH1A3_CCL28',
    '38'='Ductal_PI3',
    '30'='Ductal_MIR205HG',
    '35'='Ionocytes',
    '13'='Endothelial_CA4',
    '22'='Endothelial_ACKR1_CCL14',
    '34'='Endothelial_RGS5',
    '2' ='Myoepithelial',
    '33'='Melanocytes',
    '1' ='Fibroblast_DNAJB1',
    '3' ='Fibroblast_THBS4',
    '9' ='Fibroblast_APOE',
    '26'='Fibroblast_PCOLCE2',
    '31'='Fibroblast_CXCL9',
    '16'='Myofibroblast',
    '27'='Pericytes',
    '37'='Unknown',
    '8' ='Plasma_IGKV2D-28',
    '15'='Plasma_IGHG2',
    '25'='Plasma_HSPB1',
    '40'='Plasma_IGHA2',
    '23'='Mac/mo_IL1B',
    '28'='Mac/mo_CCL18',
    '4' ='CD8_T_CREM_HSPA1A',
    '18'='CD8_T_GZMA',
    '20'='CD8_T_GZMB',
    '6' ='CD4_T_ICOS',
    '19'='CD4_T_CCL5',
    '0' ='B',
    '39'='Mast'
)

Idents(CX) <- "seurat_clusters"
CX <- RenameIdents(CX, celltype_minor_dict)
CX$celltype_minor <- Idents(CX)

options(repr.plot.width = 12, repr.plot.height = 9, repr.plot.res = 150)
DotPlot(CX,feature=unique(c('ANKRD36C','ELF3','SCGB3A1','TAGLN','LRRC26','MUC7','CLDN10','PRR4','PI3','FCGBP','LTF','ALDH1A3','CCL28','PI3','MIR205HG','CA4','ACKR1','CCL14','RGS5','DNAJB1','THBS4','APOE','PCOLCE2','CXCL9','IGKV2D-28','IGHG2','HSPB1','IGHA2','IL1B','CCL18','CREM','HSPA1A','GZMA','GZMB','ICOS','CCL5')),group.by='celltype_minor') + RotatedAxis()

options(repr.plot.width = 16,repr.plot.height = 10, repr.plot.res = 100)
DimPlot(CX,group.by='celltype_minor',label=T)

cx = sc.read_h5ad('s02_cx_annotated.h5ad')

def get_cluster_proportions(adata,
                            cluster_key="cluster_final",
                            sample_key="replicate",
                            drop_values=None):
    """
    Input
    =====
    adata : AnnData object
    cluster_key : key of `adata.obs` storing cluster info
    sample_key : key of `adata.obs` storing sample/replicate info
    drop_values : list/iterable of possible values of `sample_key` that you don't want
    
    Returns
    =======
    pd.DataFrame with samples as the index and clusters as the columns and 0-100 floats
    as values
    """
    
    adata_tmp = adata.copy()
    sizes = adata_tmp.obs.groupby([cluster_key, sample_key]).size()
    props = sizes.groupby(level=1).apply(lambda x: 100 * x / x.sum()).reset_index() 
    props = props.pivot(columns=sample_key, index=cluster_key).T
    props.index = props.index.droplevel(0)
    props.fillna(0, inplace=True)
    
    if drop_values is not None:
        for drop_value in drop_values:
            props.drop(drop_value, axis=0, inplace=True)
    return props


def plot_cluster_proportions(cluster_props, 
                             cluster_palette=None,
                             xlabel_rotation=0): 
    fig, ax = plt.subplots(dpi=300)
    fig.patch.set_facecolor("white")
    
    cmap = None
    if cluster_palette is not None:
        cmap = sns.palettes.blend_palette(
            cluster_palette, 
            n_colors=len(cluster_palette), 
            as_cmap=True)
   
    cluster_props.plot(
        kind="bar", 
        stacked=True, 
        ax=ax, 
        legend=None, 
        colormap=cmap
    )
    
    ax.legend(bbox_to_anchor=(1.01, 1), frameon=False, title="Cluster")
    sns.despine(fig, ax)
    ax.tick_params(axis="x", rotation=xlabel_rotation)
    ax.set_xlabel(cluster_props.index.name.capitalize())
    ax.set_ylabel("Proportion")
    fig.tight_layout()
    
    return fig

class_p_df = get_cluster_proportions(cx,cluster_key='cell_class',sample_key='sample')

class_p_df_count = pd.DataFrame()
count = 0
for i in class_p_df.index:
    for c in class_p_df.columns:
        class_p_df_count.loc[count,['cell_class','sample','proportion']] = [c,i,class_p_df.loc[i,c]]
        count += 1
class_p_df_count.loc[class_p_df_count['sample'].isin(['HCCX_1','HCCX_2','HCCX_3','HCCX_4','HCCX_5']),'cli_state'] = 'HC'
class_p_df_count.loc[class_p_df_count['sample'].isin(['XNCX_1','XNCX_2','XNCX_3','XNCX_4','XNCX_5','XNCX_6','XNCX_7','XNCX_8','XNCX_9','XNCX_10','XNCX_11']),'cli_state'] = 'SjS'
class_p_df_count

class_p_df_count = pd.concat([class_p_df_count[class_p_df_count['cell_class'] == c] for c in order],axis=0)

fig,ax=plt.subplots(figsize=(3,2.5))
sns.boxplot(data=class_p_df_count,y='proportion',x='cell_class',hue='cli_state',ax=ax,fliersize=0,linewidth=0.5)
sns.stripplot(data=class_p_df_count,y='proportion',x='cell_class',hue='cli_state',ax=ax,dodge=True,size=3,linewidth=0.5)
plt.xticks(rotation=90)
# box lines
for i, box in enumerate(ax.artists):
    color = box.get_facecolor()
    box.set_edgecolor(color)
    box.set_facecolor('white')

    # iterate over whiskers and median lines
    for j in range(6*i, 6*(i+1)):
        ax.lines[j].set_color(color)

# scatter points
for coll in ax.collections:
    color = coll.get_facecolor()
    coll.set_edgecolor(color)
    coll.set_facecolor('white')
plt.savefig('figures/boxplot_cx_cell_class_proportion.pdf',bbox_inches='tight')
plt.show()

celltype_major_major_map = pd.DataFrame(np.array(['Epithelial','Epithelial','Epithelial','Plasma','Myeloid','Endothelial','T','Ionocytes','T','SMC','Pericytes','B','Epithelial','Fibroblast','Melanocytes','Unknown','Myeloid']),index=np.array(cx.obs.celltype_major.unique()))
cx.obs.loc[:,'celltype_major_major'] = pd.Categorical(celltype_major_major_map.loc[cx.obs.celltype_major,:][0].values)
cx.obs.celltype_major_major = pd.Categorical(cx.obs.celltype_major_major, categories=['Epithelial','Ionocytes','Fibroblast','SMC','Endothelial','Pericytes','Melanocytes','T','B','Plasma','Myeloid'])

cx.obs.celltype_major_major = np.array(cx.obs.celltype_major_major)

cx.obs.loc[cx.obs.celltype_minor == 'Endothelial_ACKR1_CCL14','celltype_major_major'] = 'Endothelial_ACKR1'
cx.obs.loc[cx.obs.celltype_minor == 'Endothelial_RGS5','celltype_major_major'] = 'Endothelial_RGS5'
cx.obs.loc[cx.obs.celltype_minor == 'Endothelial_CA4','celltype_major_major'] = 'Endothelial_CA4'

cx.obs.loc[cx.obs.celltype_major == 'Serous_acini','celltype_major_major'] = 'Serous_acini'
cx.obs.loc[cx.obs.celltype_major == 'Mucous_acini','celltype_major_major'] = 'Mucous_acini'
cx.obs.loc[cx.obs.celltype_major == 'Ductal','celltype_major_major'] = 'Ductal'
cx.obs.loc[cx.obs.celltype_major == 'Myoepithelial','celltype_major_major'] = 'Myoepithelial'

celltype_major_major_map.index = [c.split('_')[0] for c in celltype_major_major_map.index]

cx_cli_state_df = []
for c in cx.obs.celltype_minor.unique():
    cx_ctmp = cx[cx.obs.celltype_minor==c].copy()
    if len(np.unique(cx_ctmp.obs.cli_state)) > 1:
        sc.tl.rank_genes_groups(cx_ctmp,groupby='cli_state',key_added=f'{c}_cli',method='wilcoxon',pts=True,use_raw=False)
        cx_ctmp_df = sc.get.rank_genes_groups_df(cx_ctmp,group='SjS',key=f'{c}_cli')
        cx_ctmp_df['celltype_minor'] = c
        cx_cli_state_df.append(cx_ctmp_df)
cx_cli_state_df = pd.concat(cx_cli_state_df,axis=0)

cx_cli_state_logfc = pd.DataFrame(index=cx_cli_state_df.celltype_minor.unique(),columns=cx_cli_state_df.names.unique())
cx_cli_state_pval = pd.DataFrame(index=cx_cli_state_df.celltype_minor.unique(),columns=cx_cli_state_df.names.unique())
cx_cli_state_pts = pd.DataFrame(index=cx_cli_state_df.celltype_minor.unique(),columns=cx_cli_state_df.names.unique())
cx_cli_state_scores = pd.DataFrame(index=cx_cli_state_df.celltype_minor.unique(),columns=cx_cli_state_df.names.unique())
for c in cx_cli_state_df.celltype_minor.unique():
    cx_cli_state_df_tmp = cx_cli_state_df[cx_cli_state_df.celltype_minor == c]
    cx_cli_state_logfc.loc[c,cx_cli_state_df_tmp.names] = cx_cli_state_df_tmp.logfoldchanges.values
    cx_cli_state_pval.loc[c,cx_cli_state_df_tmp.names] = cx_cli_state_df_tmp.pvals_adj.values
    cx_cli_state_pts.loc[c,cx_cli_state_df_tmp.names] = cx_cli_state_df_tmp.pct_nz_group.values
    cx_cli_state_scores.loc[c,cx_cli_state_df_tmp.names] = cx_cli_state_df_tmp.scores.values
cx_cli_state_logfc = cx_cli_state_logfc.fillna(0)
cx_cli_state_pval = cx_cli_state_pval.fillna(0)
cx_cli_state_pts = cx_cli_state_pts.fillna(0)
cx_cli_state_scores = cx_cli_state_scores.fillna(0)

cx_cli_state_up_logfc = cx_cli_state_logfc.loc[:,(((cx_cli_state_pts.values > 0.1)&(cx_cli_state_logfc.values > 0.5)&(cx_cli_state_pval.values < 0.01)).sum(0) != 0)]
cx_cli_state_down_logfc = cx_cli_state_logfc.loc[:,(((cx_cli_state_pts.values > 0.1)&(cx_cli_state_logfc.values < -0.5)&(cx_cli_state_pval.values < 0.01)).sum(0) != 0)]

cx_cli_state_up_logfc = cx_cli_state_up_logfc*(cx_cli_state_pval.loc[:,cx_cli_state_up_logfc.columns].values < 0.01)
cx_cli_state_down_logfc = cx_cli_state_down_logfc*(cx_cli_state_pval.loc[:,cx_cli_state_down_logfc.columns].values < 0.01)

major_up_deg_counts = pd.DataFrame([((cx_cli_state_up_logfc[(celltype_major_major_map.loc[[c.split('_')[0] for c in cx_cli_state_up_logfc.index],0] == c).values]>0.5).sum(0)!=0).sum() for c in np.unique(celltype_major_major_map[0].values)],index=np.unique(celltype_major_major_map[0].values))
major_down_deg_counts = pd.DataFrame([((cx_cli_state_down_logfc[(celltype_major_major_map.loc[[c.split('_')[0] for c in cx_cli_state_down_logfc.index],0] == c).values]<-0.5).sum(0)!=0).sum() for c in np.unique(celltype_major_major_map[0].values)],index=np.unique(celltype_major_major_map[0].values))

major_up_deg_counts.to_csv('cx_celltype_major_major_up_deg_counts.csv')

order = (major_up_deg_counts+major_down_deg_counts).sort_values(by=0,ascending=False).index

sns.set_style('ticks')
plt.subplots(figsize=(3,2.5))
sns.barplot(x=np.arange(len(major_up_deg_counts)),y=major_up_deg_counts.loc[order,0])
sns.barplot(x=np.arange(len(major_down_deg_counts)),y=-major_down_deg_counts.loc[order,0])
plt.xticks(np.arange(len(major_down_deg_counts)),order,rotation=90)
plt.ylabel('DEGs number')
# plt.savefig('figures/barplot_cx_cell_class_degs_number.pdf',bbox_inches='tight',dpi=300)
plt.show()

sc.pl.stacked_violin(cx,['EPCAM','KRT7','ASCL3','FOXI1', 'DCN', 'COL1A1', 'MUSTN1', 'MYH11', 'PECAM1', 'SOX18', 'RGS5', 'CCL21', 'GFRA3', 'FOXD3', 'CD3D', 'CD3E', 'MS4A1', 'CD79A', 'MZB1', 'JCHAIN', 'CD68', 'AIF1'],groupby='celltype_major_major',cmap='Reds',use_raw=False,vmax=4,edgecolor=None,figsize=(7,7),save='cx_celltype_major_major_makers.pdf')
