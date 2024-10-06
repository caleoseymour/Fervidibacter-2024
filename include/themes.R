library('extrafont')
loadfonts(device="win")

max_text_size = 8
min_text_size = 5 * 8/7


GHPal = c('GH109' = '#58111A', 'GH177' = '#BF4F51', 'GH179' = '#E03C31', 'Other GH' = 'grey90')


theme_nmds = theme(
    text = element_text(color='black', size = 8, family = 'Arial'),
    panel.background = element_rect(fill='white', color = 'black'),
    panel.border = element_rect(fill='transparent', color = 'black'),
    #plot.background = element_rect(fill='transparent', color = 'black'),
    axis.text.x = element_text(family='Arial', size=8, vjust=1, hjust = 0.5),
    axis.ticks.y = element_line(linewidth=0.25),
    axis.ticks.x = element_line(linewidth=0.25),
    panel.grid.major.y = element_blank(), #element_line(linetype = "dotted", color = "grey"),
    panel.grid.major.x = element_blank(), #element_line(linetype = "dotted", color = "grey"),
    axis.text.y = element_text(family='Arial', size=8, vjust=0.5, hjust = 1),
    axis.title = element_blank(),
    axis.title.x = element_text(family='Arial', size=8, face='bold'),
    axis.title.y = element_text(family='Arial', size=8, face='bold'),
    legend.text = element_text(family='Arial', size=8),
    legend.title = element_text(family='Arial', size=8, face='bold'),
    legend.spacing.x = unit(0.1, 'cm'),
    legend.spacing.y = unit(0, 'cm'),
    #panel.border = ggplot2::element_rect(color='black'),
    #panel.spacing = unit(0, 'cm'),
    #strip.placement = 'outside',
    legend.key.size = unit(0.25,"cm"),
    legend.key.height = unit(0.25, 'cm'), #change legend key height
    legend.key.width = unit(0.25, 'cm'),
    legend.justification = c(0,1),
    #legend.background=element_blank(),
    legend.box = "vertical",
    legend.position = 'top'
)

theme_lefse = theme(
    legend.position=c(0.01,0.99),
    panel.background = element_blank(),
    plot.background = element_rect(fill='white'),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(family='Arial', size=8, vjust=1, hjust = 0.5),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.title.x = element_text(family='Arial', size=8, face='bold'),
    axis.title.y = element_text(family='Arial', size=8, face='bold'),
    legend.text = element_text(family='Arial', size=8),
    legend.title = element_text(family='Arial', size=8, face='bold'),
    legend.spacing.x = unit(0, 'cm'),
    legend.spacing.y = unit(0, 'cm'),
    panel.border = ggplot2::element_rect(color='black'),
    panel.spacing = unit(0, 'cm'),
    strip.background = ggplot2::element_rect(color='black',fill='white'),
    strip.background.x = ggplot2::element_rect(color='black',fill='white'),
    strip.background.y = ggplot2::element_rect(color='black',fill='white'),
    strip.placement = 'outside',
    legend.key.size = unit(0.25,"cm"),
    legend.key.height = unit(0.25, 'cm'), #change legend key height
    legend.key.width = unit(0.25, 'cm'),
    legend.justification = c(0,1),
    #legend.background=element_blank(),
    legend.box = "vertical",
    plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
)
