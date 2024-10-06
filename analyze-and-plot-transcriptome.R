## Cale Seymour, MS
## University of Nevada, Las Vegas
## 2023
## This script plots the transcriptome data and iterates over the different treatment comparisons to run stats.
library('DESeq2')
library('data.table')
library('magrittr')
library('ggplot2')
library('ClusterR')
library('vegan')
library('ggpubr')
library('pheatmap')
library('cowplot')
library('svglite')
library('plotly')
library('htmlwidgets')


palette = fread('palette.txt', sep='\t', header=TRUE)[,.(HEX, Label = Substrate, Column = Substrate %>% make.names())]
source('include/read.gff.R')
source('include/themes.R')

gffname = 'FSac-sealed-genome/GCF_030520105.1.gff3'
gff = read.gff(gffname)

cazymeannots = fread('CAZymes.tsv.xls', sep='\t', header=TRUE, check.names=TRUE) %>% setkey(Gene.ID)
cazymes = cazymeannots[, Gene.ID]
GHs = cazymeannots[grepl('GH', Consensus), Gene.ID]
allprots = gff[,protein_id] %>% na.omit()

annotations = gff[!is.na(protein_id), .(protein_id, product, go_function, go_process, dbCAN = as.character(NA))] %>% setkey(protein_id)
annotations[cazymeannots, dbCAN := Consensus]
annotations[dbCAN == 'N', dbCAN := NA]

## Makes a volcano plot
volcano_plot = function(r)
{
    r$sig = ifelse(r$padj<0.05, "q < 0.05", "N.S.")
    
    p = ggplot(r, aes(log2FoldChange, -log10(pvalue))) +
        geom_point(aes(color = sig)) +
        scale_color_manual(values = c("red", "black")) +
        ggtitle("Volcano Plot of DESeq2 analysis") +
        theme_pubclean()
      
    return(p)
}

dir.create('output-transcriptome')

## Read in metadata, do some filtering.
metadata = fread('read_mappings.txt.xls', header=TRUE, sep='\t', check.names=TRUE)
metadata = metadata[X6.conditionNumber != 0,]
metadata$X7.groupName = make.names(metadata$X7.groupName)
colData = metadata %>% as.data.frame()
rownames(colData) = metadata[,X1.libraryName]

## Read in counts.
x = fread('FSac-sealed-genome_transcriptome/feature-counts.tsv.xls', header=TRUE, sep='\t')
counts = x[,-c(1:6)] %>% as.data.frame()
rownames(counts) = x[,Geneid]
colnames(counts) = colnames(counts) %>% gsub('\\.Aligned\\.out\\.bam','', .) %>% gsub('.*\\/','',.)

countData = counts[,metadata[,X1.libraryName]]

sink('count-summary.txt')
    summary(counts)
sink()

insituactivity = fread('metaproteome-activity.tsv.xls', sep='\t', header=TRUE)

## Expressed activity
expressed = c('WP_259102106.1', 'WP_259097438.1', 'WP_259092593.1', 'WP_259093347.1', 'WP_259101925.1', 'WP_259101730.1')
protein_names = c('Fsa01040Cel', 'Fsa02490Xyn', 'Fsa15405Xyn', 'Fsa16295Glu' , 'Fsa01670Gal', 'Fsa11540Glu') %>% setNames(expressed)

## Build DESeq object.
dds = DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~ X7.groupName)
dds = dds[rowSums(counts(dds)) >= 10,]

## Normalize counts.
dds = estimateSizeFactors(dds)
rld = rlog(dds, blind = FALSE)
logcpm = assay(rld)

if(!exists('ddsLRT')) ddsLRT = DESeq(dds, test = "LRT", reduced = ~ 1)
lrt_res = results(ddsLRT) %>% as.data.table(keep.rownames=TRUE)
outfname = paste0('output-transcriptome/de/raw/diffexp-lrt.txt.xls')
de_genes = lrt_res[rn %in% paste0('cds-',allprots),]
de_genes[,rn := gsub('cds-','',rn)]
setkey(de_genes, rn)
de_genes_annotations = annotations[de_genes,][!duplicated(protein_id),]
fwrite(de_genes_annotations[order(-abs(log2FoldChange)),], outfname, sep='\t')

logmeandt = lrt_res[,.(rn = gsub('^cds-','', rn), baseMean)] %>% setkey(rn)

normcounts = counts(dds, normalized=TRUE) %>% as.data.table(keep.rownames=TRUE)
mnormcounts = normcounts %>% melt(id.var = 'rn') %>% setkey(variable)
setkey(metadata, X1.libraryName)

mnormedmeans = metadata[mnormcounts,][,.(groupedMean = mean(value)),by=.(rn, X7.groupName)]
normedmeans = mnormedmeans %>% dcast(rn ~ X7.groupName, value.var = 'groupedMean')

fwrite(normedmeans, 'output-transcriptome/normalized-group-means.tsv.xls', sep='\t')

for (prots_i in c('cazymes', 'GHs', 'allprots'))
{
    genes_i = get(prots_i) %>% paste0('cds-',.)
    logcpm_i = logcpm[intersect(rownames(logcpm), genes_i),]
    
    # Scaled log-counts, with the minimum value added back to the counts in order
    # to make all count values positive after scaling.
    logcpm_scaled = logcpm_i %>% t() %>% scale()
    logcpm_positive = logcpm_scaled - logcpm_scaled %>% min()
    
    ## NMDS and adonis
    mds = metaMDS(logcpm_positive, distance = 'euclidean', trymax=1000, wascores = TRUE)
    dist_matrix = logcpm_positive %>% vegdist('euclidean')
    
    adonis = adonis2(dist_matrix ~ X7.groupName, data = metadata)
    pval = adonis$`Pr(>F)`[1]
    
    adoniscombns = metadata[,X7.groupName] %>% unique() %>% combn(2)
    
    onevsoneadonis = apply(adoniscombns, 2,
        function(adoniscombn)
        {
            a = adoniscombn[1]
            b = adoniscombn[2]
            combnmetadata = metadata[X7.groupName %in% c(a, b),]
            combnsamples = combnmetadata[,X1.libraryName]
            combndist =  (dist_matrix %>% as.matrix)[combnsamples,combnsamples] %>% as.dist()
            
            combnres = adonis2(combndist ~ X7.groupName, data = combnmetadata)
            data.table(a, b, F=combnres$`F`[1], 'Pr(>F)'=combnres$`Pr(>F)`[1])
        }
    ) %>% do.call('rbind',.)
    
    onevsalladonis = lapply(metadata[,X7.groupName] %>% unique(),
        function(adoniscombn)
        {
            a = adoniscombn
            combnmetadata = copy(metadata)
            combnmetadata[,grouping := 'Not']
            combnmetadata[X7.groupName == a, grouping := 'In']
            combnsamples = combnmetadata[,X1.libraryName]
            combndist =  (dist_matrix %>% as.matrix)[combnsamples,combnsamples] %>% as.dist()
            
            combnres = adonis2(combndist ~ grouping, data = combnmetadata)
            data.table(group = a, F=combnres$`F`[1])
        }
    ) %>% do.call('rbind',.)


    ord = cbind(colnames(logcpm_i) %>% as.data.table(), mds$points) %>% setnames(c('sample','x','y'))
    setkey(ord, sample)
    stress = mds$stress

    sitepoints = ord[metadata,]
    
    spp = mds$species %>% as.data.table(keep.rownames=TRUE) %>% setnames(c('species','x','y'))
    spp[,species := gsub('^cds-', '', species)]
    setkey(spp, species)
    
    speciespoints = annotations[spp,]
    speciespoints[,isCAZyme := protein_id %in% cazymes]
    speciespoints[,isGHs := protein_id %in% GHs]
        
    speciespoints = speciespoints[order(isGHs + isCAZyme),][isGHs == TRUE,]
    
    speciespoints[grepl('GH', dbCAN), GHs := 'Other GH']
    speciespoints[grepl('GH109', dbCAN), GHs := 'GH109']
    speciespoints[grepl('GH177', dbCAN), GHs := 'GH177']
    speciespoints[grepl('GH179', dbCAN), GHs := 'GH179']
    
    speciespoints$color = GHPal[speciespoints$GH]

    xlimits = c(min(sitepoints$x, speciespoints$x), max(sitepoints$x, speciespoints$x))
    xroundlims = ceiling(abs(xlimits)/5)*5 * sign(xlimits)
    
    ylimits = c(min(sitepoints$y, speciespoints$y), max(sitepoints$y, speciespoints$y))
    yroundlims = ceiling(abs(ylimits)/5)*5 * sign(ylimits)
    
    statistics_label = paste0('stress = ', round(stress, 3), '\nadonis Pr(>F) = ', pval)
    
    ordplot = ggplot() +
        geom_polygon(data = sitepoints, aes(x = x, y = y, color=X7.groupName, fill=X7.groupName, group = X7.groupName), alpha=0.2, linewidth = 0.5) +
        geom_point(data = speciespoints, aes(x = x, y = y, protein_id = protein_id, sample = ''), color = speciespoints$color, size = 1) +
        geom_point(data = sitepoints, aes(x = x, y = y, color = gsub('.*_','', X7.groupName), protein_id = '', sample = sample), size = 2) +
        annotate(geom='text', x = min(xroundlims), y = max(yroundlims), label = statistics_label, hjust=0, vjust=1, size = 8 * (5/14), lineheight = 0.75, fontfamily='Arial') +
        scale_fill_manual(values = palette[,HEX], breaks = palette[,Column], labels = palette[,Label], guide='none') +
        scale_color_manual(values = palette[,HEX], breaks = palette[,Column], labels = palette[,Label]) +
        guides(color=guide_legend('Substrate', override.aes = list(linetype = 0, fill = 'transparent', size = 2), nrow = 4, title.position = 'top')) +
        xlab('MDS1') +
        ylab('MDS2') +
        theme_pubclean() +
        theme_nmds
        
    fwrite(sitepoints, paste0('output-transcriptome/', prots_i, '-sites.txt'), sep='\t')
    fwrite(speciespoints, paste0('output-transcriptome/', prots_i, '-spp.txt'), sep='\t')

    sink(paste0('output-transcriptome/', prots_i, '-adonis.txt'))
        adonis %>% print()
        cat('\nOne-vs-all\n')
        print(onevsalladonis)
        cat('\nOne-vs-One\n')
        print(onevsoneadonis)
    sink()


    pdf(paste0('output-transcriptome/', prots_i, '-ordination.pdf'), width = 3.5, height = 3.5)
        plot(ordplot)
    dev.off()
    
    fig_nmds = ggplotly(ordplot, tooltip=c('protein_id', 'sample')) %>% layout(title = paste0(prots_i, 'Interactive NMDS'), width = 800, height = 800)
    saveWidget(fig_nmds, paste0('output-transcriptome/', prots_i, '-ordination.html'), selfcontained = TRUE)
    
    assign(paste0(prots_i, '_ordplot'), ordplot)

    ## Run DESeq.
    dds = DESeq(dds)

    mod_mat = model.matrix(design(dds), colData(dds))

    ## Iterate over comparisons to generate lists of DE genes for each treatment.
    dir.create('output-transcriptome/de')
    dir.create('output-transcriptome/de/raw')
    dir.create(paste0('output-transcriptome/de/raw/', prots_i))
    dir.create('output-transcriptome/de/summary_plots')
    dir.create(paste0('output-transcriptome/de/summary_plots/', prots_i))

    comparisons = dds$X7.groupName %>% unique() %>% combn(2) %>% t()
    for(i in 1:nrow(comparisons))
    {
        ## Use a DESeq2 contrast to extract this comparison.
        a = comparisons[i,1]
        b = comparisons[i,2]
        
        a_means = colMeans(mod_mat[dds$X7.groupName == a,])
        b_means = colMeans(mod_mat[dds$X7.groupName == b,])
        
        ab_res = results(dds, contrast = a_means - b_means) %>% as.data.table(keep.rownames=TRUE)
        
        ## Create volcano plot
        outfname = paste0('output-transcriptome/de/summary_plots/volcano_plot.', a, '_vs_', b, '.pdf')
        vp = volcano_plot(ab_res)
        pdf(outfname)
            plot(vp)
        dev.off()
        
        ## Write de genes to a file
        outfname = paste0('output-transcriptome/de/raw/', prots_i, '/', prots_i, '-differential_expression.', a, '_vs_', b, '.txt.xls')
        de_genes = ab_res[rn %in% genes_i,][padj < 0.05 & abs(log2FoldChange) >= 1,]
        de_genes[,rn := gsub('cds-','',rn)]
        setkey(de_genes, rn)
        de_genes_annotations = annotations[de_genes,][!duplicated(protein_id),]
        fwrite(de_genes_annotations[order(-abs(log2FoldChange)),], outfname, sep='\t')
        
        de_genes_annotations[is.na(dbCAN), dbCAN_annot := 'Not in dbCAN']
        de_genes_annotations[grepl('GH', dbCAN), dbCAN_annot := 'GH']
        de_genes_annotations[!grepl('GH', dbCAN) & !is.na(dbCAN), dbCAN_annot := 'Other CAZymes']
        
        de_genes_annotations[,protein_name := paste0(product, ' (', protein_id, ')')]
        de_genes_annotations[dbCAN_annot != 'Not in dbCAN',protein_name := paste0(dbCAN, ' (', protein_id, ')')]
        de_genes_annotations[, protein_factor := factor(protein_id, levels=protein_id[order(log2FoldChange)] %>% unique())]
        
        label_x = paste0("Log2-fold change,\n(positive = ", a, ")\n","(negative = ", b, ")")
        p = ggplot(de_genes_annotations, aes(log2FoldChange, protein_factor)) + 
            geom_col(aes(fill = dbCAN_annot), color='black', width = 0.8) +
            geom_text(aes(label = protein_name), x = 0 - sign(de_genes_annotations$log2FoldChange) * 0.05, hjust = ifelse(de_genes_annotations$log2FoldChange > 0, 1, 0), size = 8 * (5/14), fontfamily='Arial') +
            labs(x = label_x, y = NULL, fill = "Annotation") +
            #scale_fill_manual(breaks = c('ko', 'wildtype')) +
            scale_x_continuous(expand = c(0, 0), breaks=c(-10:-1, 0:10), limits = c(-10, 10)) +
            theme_bw() + theme_lefse
        pdf(paste0('output-transcriptome/de/summary_plots/', prots_i, '/', prots_i, '-differential_expression.', a, '_vs_', b, '.pdf'), height = 1.5+nrow(de_genes_annotations)*0.2, width = 10)
            print(p)
        dev.off()

    }
    
    # Normalize.
    dds_norm = rlog(dds)
    
    # Filter the data choosing only genes whose variances are in the upper quartile
    df_by_var = data.frame(assay(dds_norm)) #%>% dplyr::filter(variances > 0)
    df_by_var = df_by_var[rownames(df_by_var) %in% genes_i,]
    annotation_df = attributes(dds_norm)$colData %>% as.data.frame()
    annotation_df = annotation_df['X7.groupName']
    annotation_df$X7.groupName = factor(annotation_df$X7.groupName, levels = palette[,Column])
    
    annotation_colors = list(X7.groupName=palette[,setNames(HEX, Column)])
    
    heatmap_annotated = pheatmap(df_by_var,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_rownames = FALSE,
        annotation_col = annotation_df,
        annotation_colors = annotation_colors,
        main = "Treatment Clustering",
        colorRampPalette(c("deepskyblue", "black", "yellow"))(25),
        scale = "row"
    )
  
  
    pdf(paste0('output-transcriptome/', prots_i, '-heatmap.pdf'))
        print(heatmap_annotated)
    dev.off()
}

source('plot-transcriptome.R')