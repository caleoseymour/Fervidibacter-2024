## Cale Seymour, MS
## University of Nevada, Las Vegas
## 2023
## This script plots the proteome data and iterates over the different treatment comparisons to run stats.
library('proDA')
library('data.table')
library('magrittr')
library('ggplot2')
library('vegan')
library('ggpubr')
library('pheatmap')
library('dplyr')
library('cowplot')
library('plotly')
library('htmlwidgets')

palette = fread('palette.txt', sep='\t', header=TRUE)[,.(HEX, Label = Substrate, Column = Substrate %>% make.names())]
source('include/read.gff.R')
source('include/themes.R')

## Makes a volcano plot
volcano_plot = function(r)
{
    r$sig = ifelse(r$adj_pval<0.05, "q < 0.05", "N.S.")
    
    p = ggplot(r, aes(diff, -log10(pval))) +
        geom_point(aes(color = sig)) +
        scale_color_manual(values = c("red", "black")) +
        ggtitle("Volcano Plot of proDA analysis") +
        theme_pubclean()
      
    return(p)
}

gffname = 'GCF_030520105.1.gff3'
gff = read.gff(gffname)

cazymeannots = fread('CAZymes.tsv.xls', sep='\t', header=TRUE, check.names=TRUE) %>% setkey(Gene.ID)
cazymes = cazymeannots[, Gene.ID]
GHs = cazymeannots[grepl('GH', Consensus), Gene.ID]
allprots = gff[,protein_id] %>% na.omit()

annotations = gff[!is.na(protein_id), .(protein_id, product, go_function, go_process, dbCAN = as.character(NA))] %>% setkey(protein_id)
annotations[cazymeannots, dbCAN := Consensus]
annotations[dbCAN == 'N', dbCAN := NA]

dir.create('output-proteome/de/raw', recursive = TRUE)

normalized_counts = fread('normalized-counts.txt', sep='\t', header=TRUE, check.names=TRUE)
metadata = fread('proteome-mappings.txt.xls', sep='\t', header=TRUE, check.names=TRUE) %>% setkey(X1.libraryName)
metadata[,X7.groupName := make.names(X7.groupName)]

count_matrix = normalized_counts[, mget(metadata[,X1.libraryName])] %>% as.matrix()
rownames(count_matrix) = normalized_counts[,PG.ProteinGroups]

logtrans_matrix = count_matrix %>% apply(2, log2)
sample_info_df = metadata[,.(name = X1.libraryName, condition = make.names(X7.groupName))] %>% tibble()
if(!exists('proDALRT')) fit = proDA(logtrans_matrix, design = ~ 0 + condition, col_data = sample_info_df)

logcounts = logtrans_matrix
logcounts[is.nan(logcounts)] = 0
logcounts = logcounts[rowSums(logcounts) != 0,]

insituactivity = fread('metaproteome-activity.tsv.xls', sep='\t', header=TRUE)

## Expressed activity
expressed = c('WP_259102106.1', 'WP_259097438.1', 'WP_259092593.1', 'WP_259093347.1', 'WP_259101925.1', 'WP_259101730.1')
protein_names = c('Fsa01040Cel', 'Fsa02490Xyn', 'Fsa15405Xyn', 'Fsa16295Glu' , 'Fsa11540Glu', 'Fsa01670Gal') %>% setNames(expressed)

## Likelihood ratio test
if(!exists('proDALRT')) proDALRT = test_diff(fit, reduced_model = ~ 1)
lrt_res = proDALRT %>% as.data.table(keep.rownames=TRUE)
outfname = paste0('output-proteome/de/raw/diffexp-lrt.txt.xls')
de_genes = lrt_res[name %in% allprots,]
setkey(de_genes, name)
de_genes_annotations = annotations[de_genes,][!duplicated(protein_id),]
fwrite(de_genes_annotations[order(-abs(f_statistic)),], outfname, sep='\t')

logmeandt = lrt_res[,.(rn = name, avg_abundance = 2^avg_abundance)] %>% setkey(rn)


normcounts = normalized_counts[,mget(c('PG.ProteinGroups', metadata[,X1.libraryName]))]
mnormcounts = normcounts %>% melt(id.var = 'PG.ProteinGroups') %>% setkey(variable)
setkey(metadata, X1.libraryName)

mnormedmeans = metadata[mnormcounts,][,.(groupedMean = mean(value)),by=.(PG.ProteinGroups, X7.groupName)]
normedmeans = mnormedmeans %>% dcast(PG.ProteinGroups ~ X7.groupName, value.var = 'groupedMean')

fwrite(normedmeans, 'output-proteome/normalized-group-means.tsv.xls', sep='\t')

for (prots_i in c('cazymes', 'GHs', 'allprots'))
{
    genes_i = get(prots_i)
    logcounts_i = logcounts[intersect(rownames(logcounts), genes_i),]
    
    
    logcpm_scaled = logcounts_i %>% t() %>% scale()
    logcpm_positive = logcpm_scaled - (logcpm_scaled %>% min())
    
    
    ## NMDS and adonis
    mds = metaMDS(logcpm_positive, distance = 'euclidean', trymax=1000)
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

    
    ord = cbind(colnames(logcounts_i) %>% as.data.table(), mds$points) %>% setnames(c('sample','x','y'))
    setkey(ord, sample)
    stress = mds$stress
    
    spp = cbind(rownames(logcounts_i) %>% as.data.table(), mds$species) %>% as.data.table(keep.rownames=TRUE) %>% setnames(c('species','x','y'))
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

    sitepoints = ord[metadata,]
    
    xlimits = c(min(sitepoints$x, speciespoints$x), max(sitepoints$x, speciespoints$x))
    xroundlims = ceiling(abs(xlimits)/5)*5 * sign(xlimits)
    
    ylimits = c(min(sitepoints$y, speciespoints$y), max(sitepoints$y, speciespoints$y))
    yroundlims = ceiling(abs(ylimits)/5)*5 * sign(ylimits)
    
    statistics_label = paste0('stress = ', round(stress, 3), '\nadonis Pr(>F) = ', pval)
    
    ordplot = ggplot() +
        geom_polygon(data = sitepoints, aes(x = x, y = y, color=X7.groupName, fill=X7.groupName, group = X7.groupName), alpha=0.2, linewidth = 0.5) +
        geom_point(data = speciespoints, aes(x = x, y = y, protein_id = protein_id, sample = ''), color = speciespoints$color, size = 1) +
        geom_point(data = sitepoints, aes(x = x, y = y, color = gsub('.*_','', X7.groupName), protein_id = '', sample = sample), size = 2) +
        annotate(geom='text', x =  min(xroundlims), y = max(yroundlims), label = statistics_label, hjust=0, vjust=1, size = 8 * (5/14), lineheight = 0.75) +
        scale_fill_manual(values = palette[,HEX], breaks = palette[,Column], labels = palette[,Label], guide='none') +
        scale_color_manual(values = palette[,HEX], breaks = palette[,Column], labels = palette[,Label]) +
        guides(color=guide_legend('Substrate', override.aes = list(linetype = 0, fill = 'transparent', size = 2), nrow = 2, title.position = 'top')) +
        xlab('MDS1') +
        ylab('MDS2') +
        theme_pubclean() +
        theme_nmds
    
    fwrite(sitepoints, paste0('output-proteome/', prots_i, '-sites.txt'), sep='\t')
    fwrite(speciespoints, paste0('output-proteome/', prots_i, '-spp.txt'), sep='\t')
    
    sink(paste0('output-proteome/', prots_i, '-adonis.txt'))
        adonis %>% print()
        cat('\nOne-vs-all\n')
        print(onevsalladonis)
        cat('\nOne-vs-One\n')
        print(onevsoneadonis)
    sink()

    pdf(paste0('output-proteome/', prots_i, '-ordination.pdf'), width = 3.5, height = 3.5)
        plot(ordplot)
    dev.off()
    
    fig_nmds = ggplotly(ordplot, tooltip=c('protein_id', 'sample')) %>% layout(title = paste0(prots_i, 'Interactive NMDS'), width = 800, height = 800)
    saveWidget(fig_nmds, paste0('output-proteome/', prots_i, '-ordination.html'), selfcontained = TRUE)
    
    assign(paste0(prots_i, '_ordplot'), ordplot)

    ## Iterate over comparisons to generate lists of DE genes for each treatment.
    dir.create('output-proteome/de')
    dir.create('output-proteome/de/raw')
    dir.create(paste0('output-proteome/de/raw/', prots_i))
    dir.create('output-proteome/de/summary_plots')
    dir.create(paste0('output-proteome/de/summary_plots/', prots_i))

    comparisons = metadata$X7.groupName %>% unique() %>% combn(2) %>% t()
    for(i in 1:nrow(comparisons))
    {
        ## Use a DESeq2 contrast to extract this comparison.
        a = comparisons[i,1]
        b = comparisons[i,2]
        
        ab_res =  test_diff(fit, eval(paste0('(condition', a, ' - condition', b,')'))) %>% as.data.table(keep.rownames=TRUE)
        
        ## Create volcano plot
        outfname = paste0('output-proteome/de/summary_plots/volcano_plot.', a, '_vs_', b, '.pdf')
        vp = volcano_plot(ab_res)
        pdf(outfname)
            plot(vp)
        dev.off()
        
        ## Write de genes to a file
        outfname = paste0('output-proteome/de/raw/', prots_i, '/', prots_i, '-differenential_peptide_abundance.', a, '_vs_', b, '.txt.xls')
        de_genes = ab_res[name %in% genes_i,][adj_pval < 0.05 & abs(diff) >= 1,]
        setkey(de_genes, name)
        de_genes_annotations = annotations[de_genes,][!duplicated(protein_id),]
        fwrite(de_genes_annotations[order(-abs(diff)),], outfname, sep='\t')
        
        de_genes_annotations[is.na(dbCAN), dbCAN_annot := 'Not in dbCAN']
        de_genes_annotations[grepl('GH', dbCAN), dbCAN_annot := 'GH']
        de_genes_annotations[!grepl('GH', dbCAN) & !is.na(dbCAN), dbCAN_annot := 'Other CAZymes']
        
        de_genes_annotations[,protein_name := paste0(product, ' (', protein_id, ')')]
        de_genes_annotations[dbCAN_annot != 'Not in dbCAN',protein_name := paste0(dbCAN, ' (', protein_id, ')')]
        de_genes_annotations[, protein_factor := factor(protein_id, levels=protein_id[order(diff)] %>% unique())]
        
        label_x = paste0("Log2-fold change,\n(positive = ", a, ")\n","(negative = ", b, ")")
        p = ggplot(de_genes_annotations, aes(diff, protein_factor)) + 
            geom_col(aes(fill = dbCAN_annot), color='black', width = 0.8) +
            geom_text(aes(label = protein_name), x = 0 - sign(de_genes_annotations$diff) * 0.05, hjust = ifelse(de_genes_annotations$diff > 0, 1, 0), size = 8 * (5/14)) +
            labs(x = label_x, y = NULL, fill = "Annotation") +
            scale_x_continuous(expand = c(0, 0), breaks=c(-10:-1, 0:10), limits = c(-10, 10)) +
            theme_bw() + theme_lefse
        pdf(paste0('output-proteome/de/summary_plots/', prots_i, '/', prots_i, '-differenential_peptide_abundance.', a, '_vs_', b, '.pdf'), height = 1.5+nrow(de_genes_annotations)*0.2, width = 10)
            print(p)
        dev.off()

    }
    
    if (any(rowSums(logcounts_i) == 0))
    {
        df_by_var = data.frame(logcounts_i)[-which(rowSums(logcounts_i) == 0),]
    } else {
        df_by_var = data.frame(logcounts_i)
    }
    df_by_var = df_by_var[rownames(df_by_var) %in% genes_i,]
    annotation_df = metadata %>% as.data.frame()
    annotation_df = annotation_df['X7.groupName']
    rownames(annotation_df) = metadata[,X1.libraryName]
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
  
  
    pdf(paste0('output-proteome/', prots_i, '-heatmap.pdf'))
        print(heatmap_annotated)
    dev.off()
}

source('plot-proteome.R')