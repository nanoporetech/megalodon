library(ggplot2)
library(dplyr)
dat <- read.table('llrs.txt')
colnames(dat) <- c('llhr','pos','mean_llhr','truth','mega_call')
dat$truth_mega <- paste(dat$truth, dat$mega_call)
pos_dat <- dat %>% group_by(pos) %>% summarize(mean_llhr=first(mean_llhr))
dat$pos <- factor(dat$pos, ordered=TRUE, levels=pos_dat$pos[order(pos_dat$mean_llhr)])

pdf('llrs.pdf', width=15)
ggplot(dat) + geom_point(aes(x=pos, y=llhr, color=truth_mega), alpha=0.05) +
    guides(color=guide_legend(override.aes=list(alpha=1))) +
    theme_bw() + theme(axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       panel.grid.major=element_blank(),
                       panel.grid.minor=element_blank()) +
    geom_hline(yintercept=0) +
    facet_grid(. ~ truth_mega, scale='free_x', space='free_x')
foo <- dev.off()
