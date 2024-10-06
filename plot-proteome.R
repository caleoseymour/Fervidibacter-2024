library('purrr')
library('tibble')
library('stringr')
library('dplyr')
library('grid')
library('gridExtra')
library('ggpubr')
library('ggrepel')
library('svglite')
library('data.table')

source('include/themes.R')

gffname = 'GCF_030520105.1.gff3'
gff = read.gff(gffname)

cazymeannots = fread('CAZymes.tsv.xls', sep='\t', header=TRUE, check.names=TRUE) %>% setkey(Gene.ID)
cazymes = cazymeannots[, Gene.ID]
GHs = cazymeannots[grepl('GH', Consensus), Gene.ID]
allprots = gff[,protein_id] %>% na.omit()

annotations = gff[!is.na(protein_id), .(protein_id, product, go_function, go_process, dbCAN = as.character(NA))] %>% setkey(protein_id)
annotations[cazymeannots, dbCAN := Consensus]
annotations[dbCAN == 'N', dbCAN := NA]

add_dodged_ticks <- function(plot, axis = 'bottom') {
  plot_grob <- ggplotGrob(plot)
  plot_grob$layout %>% as.tibble() %>%
    rownames_to_column('index') %>%
    filter(str_detect(name, paste0('axis-',axis[0]))) %>% 
    pull(index) %>% as.numeric() %>% 
    map(function(x) {
      ticks <- plot_grob$grobs[[x]]$children[[2]]$grobs[[1]]$y
      fourseq <- ifelse(length(ticks) < 2, 3, length(ticks)) 
      for (i in 1:length(ticks)) {
        if (i %in% seq(3, fourseq, 4)) {
          plot_grob$grobs[[x]]$children[[2]]$grobs[[1]]$y[[i]] <<- unit(1,'npc') - unit(7,'pt')
        }
      }
    })
  return(plot_grob)
}

genenamekey = fread('geneid-to-name-key.tsv.xls', sep='\t', header=TRUE)[,setNames(GeneName, GeneID)]

logmeanGHs = annotations[logmeandt,][protein_id %in% GHs,][order(-avg_abundance),]
logmeanGHs[,protein_id := factor(protein_id, levels = protein_id)]
logmeanGHs[grepl('GH', dbCAN), GHs := 'Other GH']
logmeanGHs[grepl('GH109', dbCAN), GHs := 'GH109']
logmeanGHs[grepl('GH177', dbCAN), GHs := 'GH177']
logmeanGHs[grepl('GH179', dbCAN), GHs := 'GH179']

logmeanGHs[,label := '']
logmeanGHs[,axlabel := '']
logmeanGHs[,multidomain := FALSE]
logmeanGHs[,insitu := FALSE]
logmeanGHs[,biochem := FALSE]

logmeanGHs[grepl('\\+', dbCAN), multidomain := TRUE]
logmeanGHs[protein_id %in% insituactivity[,protein_id], insitu := TRUE]
logmeanGHs[protein_id %in% expressed, biochem := TRUE]

logmeanGHs[insitu | biochem, label := .I]
logmeanGHs[, axlabel := paste0(ifelse(.I %% 2 == 0,'\n',''), label)]

expressed_here = protein_names[intersect(expressed, logmeanGHs[,as.character(protein_id)])]
GHPal = c('GH109' = '#58111A', 'GH177' = '#BF4F51', 'GH179' = '#E03C31', 'Other GH' = 'grey90')

ymax = (logmeanGHs[,avg_abundance] * 1.05) %>% max()
xmax = logmeanGHs[,protein_id] %>% as.numeric() %>% max()

multidomain_legend = paste0(logmeanGHs[label != '',][multidomain == TRUE, paste0(label, '. ', genenamekey[as.character(protein_id)])] %>% as.character(), collapse='\n')
biochem_legend = paste0(logmeanGHs[label != '',][biochem == TRUE, paste0(label, '. ', genenamekey[as.character(protein_id)])] %>% as.character(), collapse='\n')

number_rows = 7
number_columns = ceiling((logmeanGHs[insitu == TRUE, .N]) / number_rows)

insitu_legends = paste0(logmeanGHs[label != '',][insitu == TRUE, paste0(label, '. ', genenamekey[as.character(protein_id)])] %>% as.character(), collapse='\n')

coord_plane = expand.grid(1:number_rows, 1:number_columns) %>% as.data.table()
colnames(coord_plane) = c('row','col')
coord_plane[,ii := .I]
coord_plane[,label := logmeanGHs[label != '',][insitu == TRUE, paste0(label, '. ', genenamekey[as.character(protein_id)])][ii]]
insitu_leg_columns = columns = coord_plane[!is.na(label),paste0(label, collapse='\n'), by=col][,V1]

logmeanGHs[, angle := 45 ]
logmeanGHs[as.numeric(protein_id) > 100, angle := 90 ]


logmeanGHs[,lab := paste0(ifelse(label != '', dbCAN, ''), ifelse(multidomain, '🞶', ''))]

barplots = ggplot() + 
    geom_bar(data = logmeanGHs, aes(x=protein_id, y=avg_abundance, fill = GHs), stat = "identity", width=0.9, color='black', linewidth=0.25) +
    geom_text(data = logmeanGHs, aes(x = protein_id, y = avg_abundance + (ymax * 0.01), label = lab), angle = logmeanGHs[,angle], size = 6 * (5/14), hjust = 0, vjust = 0.5, lineheight = 0.5, family='Arial') +
    #annotate(geom='text', x = xmax*0.725, y = ymax, label = multidomain_legend, hjust=0, vjust=1, size = 2, lineheight = 0.75, family='Arial') +
    #annotate(geom='text', x = xmax*0.725, y = ymax*1.05, label = 'Multidomain enzymes', hjust=0, vjust=1, size = 2, lineheight = 0.75, fontface='bold', family='Arial') +
    annotate(geom='text', x = xmax* (0.85 + 0.03), y = ymax*1.05, label = 'Expressed', hjust=0, vjust=1, size = 6 * (5/14), lineheight = 0.75, fontface='bold', family='Arial') +
    annotate(geom='text', x = xmax* (0.85 + 0.03), y = ymax, label = biochem_legend, hjust=0, vjust=1, size = 6 * (5/14), lineheight = 0.75, family='Arial')

min_x_leg = (0.89 - 0.1) - (number_columns - 1) * 0.135
barplots = barplots + annotate(geom='text', x = xmax * (min_x_leg), y = ymax*1.05, label = 'Metaproteome', hjust=0, vjust=1, size = 2, lineheight = 0.75, fontface='bold', family='Arial')
for (ii in number_columns:1)
{
    barplots = barplots + annotate(geom='text', x = xmax * (min_x_leg + (0.11 * (ii - 1))), y = ymax, label = insitu_leg_columns[ii], hjust=0, vjust=1, size = 2, lineheight = 0.75, family='Arial')
}

barplots = barplots + scale_fill_manual(values = GHPal, breaks = names(GHPal), labels = names(GHPal)) +
scale_x_discrete(labels = logmeanGHs[,axlabel], expand = expansion(mult = c(0.01, 0.01))) +
scale_y_continuous(expand = expansion(mult = c(0, 0.025))) +
xlab('Glycoside hydrolase') +
ylab('Average Normalized Abundance (LRT)') +
guides(fill=guide_legend('GH type', nrow = 2, title.position = 'top')) +
theme_pubclean() +
theme_nmds +
theme(axis.text.x = element_text(size=6, vjust=1, hjust = 0.5, angle=0, family='Arial', lineheight=0.7),
      axis.ticks.length=unit(3, "pt"),
      axis.ticks.x = element_line(linewidth=0.25))

barplots2 = add_dodged_ticks(barplots + theme(legend.position = 'none'), axis = 'bottom') %>% as_ggplot()

pdf(paste0('output-proteome/barplot.pdf'), width = 7, height = 3, family='Arial')
    plot(barplots2)
dev.off()

svg(paste0('output-proteome/barplot.svg'), width = 7, height = 3, family='Arial')
    plot(barplots2)
dev.off()


## Plot!
figure = plot_grid(
    barplots2,
    plot_grid(ggplot() + theme_void(), get_legend(GHs_ordplot), get_legend(barplots), rel_widths = c(0.2, 1.5, 1), nrow = 1),
    plot_grid(
        allprots_ordplot + theme(legend.position = 'none'),
        cazymes_ordplot + theme(legend.position = 'none'),
        label_size = max_text_size,
        labels = c('b','c'), label_fontfamily = 'Arial'
        )
    , nrow = 3
    , rel_heights = c(1 + 1/6, 0.3, 1 + 1/6)
    , label_size = max_text_size
    , labels = c('a','',''), label_fontfamily = 'Arial'
)

svglite(paste0('output-proteome/proteome-figure-', format(Sys.Date(), "%m-%d-%Y"), '.svg'), width = 6.05*8/7, height = 5.09*8/7)
    plot(figure)
dev.off()

pdf(paste0('output-proteome/proteome-figure-', format(Sys.Date(), "%m-%d-%Y"), '.pdf'), width = 6.05*8/7, height = 5.09*8/7, family='Arial')
    plot(figure)
dev.off()