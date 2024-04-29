# setwd('/home/timaeus/Desktop/Projects/Plants/UMN/BioCon_spatdiv')

library(dplyr)
library(ggplot2)
library(patchwork)

bcdat <- read.csv('../hi_res_SR_biocon.csv')

# remove bare ground
bcspec <- bcdat[bcdat$Species != 'Bargr',]
bcspec$Species <- as.factor(bcspec$Species)
bcspec$num.spec <- as.numeric(bcspec$Species)
uspec <- unique(bcspec$Species)

# order plots
plots <- unique(bcdat$Plot)
plots <- plots[order(plots)]

plist <- vector("list",length=31)

# get lists of species richness across all sampling squares
idx <- 1
for (p in seq(1,length(plots))) {
    # set up array to store cover fractions
    sp_cov_arr <- array(dim=c(18,18,16),
                        dimnames=(list(row=1:18,column=1:18,
                                  species=unique(bcspec$Species)))) 
    pl <- which(bcspec$Plot==plots[p])
    bc <- bcspec[pl,]
    # only BioCON
    if (bc$CountOfSpecies[1]==16) {
        # only count rooted plants?
        #r <- bc$Root == "R"
        #bc <- bc[r,]
        for (row in seq(1,18)) {
            for (col in seq(1,18)) {
                r <- match(bc$Row,row)
                c <- match(bc$Column,col)
                # species cover
                s_cover <- bc[!is.na(r) & !is.na(c),c('Species','Coverage')]
                sp <- match(s_cover$Species,uspec)
                sp_cov_arr[col,row,sp] <- s_cover$Coverage
            }
        }
        #print(idx)
        plist[[idx]]$plot <- plots[p]
        plist[[idx]]$treat <- paste(bc$CO2Treatment[1],bc$Ntreatment[1],sep='-')  
        plist[[idx]]$array <- sp_cov_arr
        idx <- idx + 1
    }
}

spnames <- as.character(unique(bcspec$Species))

sp10 <- vector(length=31)
sp20 <- vector(length=31)
sp30 <- vector(length=31)
sp40 <- vector(length=31)
sp50 <- vector(length=31)
sp60 <- vector(length=31)
sp70 <- vector(length=31)
sp80 <- vector(length=31)
sp90 <- vector(length=31)

for (i in seq(1,length(plist))){
    sp10[i] <- length(spnames[!is.nan(apply(plist[[i]]$array,3,mean,na.rm=T))])
    sp20[i] <- length(spnames[!is.nan(apply(plist[[i]]$array[2:17,2:17,],3,mean,na.rm=T))])
    sp30[i] <- length(spnames[!is.nan(apply(plist[[i]]$array[3:16,3:16,],3,mean,na.rm=T))])
    sp40[i] <- length(spnames[!is.nan(apply(plist[[i]]$array[4:15,4:15,],3,mean,na.rm=T))])
    sp50[i] <- length(spnames[!is.nan(apply(plist[[i]]$array[5:14,5:14,],3,mean,na.rm=T))])
    sp60[i] <- length(spnames[!is.nan(apply(plist[[i]]$array[6:13,6:13,],3,mean,na.rm=T))])
    sp70[i] <- length(spnames[!is.nan(apply(plist[[i]]$array[7:12,7:12,],3,mean,na.rm=T))])
    sp80[i] <- length(spnames[!is.nan(apply(plist[[i]]$array[8:11,8:11,],3,mean,na.rm=T))])
    sp90[i] <- length(spnames[!is.nan(apply(plist[[i]]$array[9:10,9:10,],3,mean,na.rm=T))])
}

plots <- as.data.frame(lapply(plist,"[","plot"))
treats <- as.data.frame(lapply(plist,"[","treat"))

tmp <- cbind(t(plots),t(treats),sp10,sp20,sp30,sp40,sp50,sp60,sp70,sp80,sp90)
row.names(tmp) <- 1:31
colnames(tmp)[1] <- "plot"
colnames(tmp)[2] <- "treat"
tmp <- as.data.frame(tmp)
tmp$plot <- as.numeric(tmp$plot)
tmp$sp10 <- as.numeric(tmp$sp10) # 3.24
tmp$sp20 <- as.numeric(tmp$sp20) # 2.56
tmp$sp30 <- as.numeric(tmp$sp30) # 1.96
tmp$sp40 <- as.numeric(tmp$sp40) # 1.44
tmp$sp50 <- as.numeric(tmp$sp50) # 1
tmp$sp60 <- as.numeric(tmp$sp60) # 0.64
tmp$sp70 <- as.numeric(tmp$sp70) # 0.36
tmp$sp80 <- as.numeric(tmp$sp80) # 0.16
tmp$sp90 <- as.numeric(tmp$sp90) # 0.04

patch.sizes <- c(0.04,0.16,0.36,0.64,1,1.44,1.96,2.56)

trt_vec <- rep(tmp$treat,9)
area_vec <- c(rep(3.24,31),rep(2.56,31),rep(1.96,31),rep(1.44,31),rep(1,31),rep(0.64,31),
              rep(0.36,31),rep(0.16,31),rep(0.04,31))
sp_vec <- c(tmp$sp10,tmp$sp20,tmp$sp30,tmp$sp40,tmp$sp50,tmp$sp60,tmp$sp70,tmp$sp80,tmp$sp90)

edgeDF <- data.frame('treat'= trt_vec, 'patch.size'= area_vec, 'sr'=sp_vec, 'aggregate'='edge')

sr_bc_trt <- read.csv('../sr_df.csv')
sr_bc_trt$aggregate <- 'full_random'

out <- rbind(edgeDF,sr_bc_trt[,c('treat','patch.size','sr','aggregate')])
out$aggregate <- factor(out$aggregate)

agg_comp <- out[out$patch.size %in% patch.sizes,]

acan <- ggplot() + geom_boxplot(aes(x=factor(patch.size),y=sr,fill='black',color=aggregate),
                                 data=agg_comp[agg_comp$treat=='AC-AN',]) +
                   scale_fill_manual(values='black',guide="none")+
                   scale_color_manual(values=c('red','blue'))+
                   labs(x=expression('Patch Size ['*m^2*']'),y='SR',title='AC-AN') +
                   ylim(0,12.5) +
                   theme(text = element_text(size=12),
                         axis.text.x = element_text(size=8),
                         panel.border = element_blank(),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.line = element_line(colour = "black"))

ecan <- ggplot() + geom_boxplot(aes(x=factor(patch.size),y=sr,fill='grey',color=aggregate),
                                 data=agg_comp[agg_comp$treat=='EC-AN',]) +
                   scale_fill_manual(values='grey',guide="none")+
                   scale_color_manual(values=c('red','blue'))+
                   labs(x=expression('Patch Size ['*m^2*']'),y='SR',title='EC-AN') +
                   ylim(0,12.5) +
                   theme(text = element_text(size=12),
                         axis.text.x = element_text(size=8),
                         panel.border = element_blank(),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.line = element_line(colour = "black"))

acen <- ggplot() + geom_boxplot(aes(x=factor(patch.size),y=sr,fill='green',color=aggregate),
                                 data=agg_comp[agg_comp$treat=='AC-EN',]) +
                   scale_fill_manual(values='green',guide="none")+
                   scale_color_manual(values=c('red','blue'))+
                   labs(x=expression('Patch Size ['*m^2*']'),y='SR',title='AC-EN') +
                   ylim(0,12.5) +
                   theme(text = element_text(size=12),
                         axis.text.x = element_text(size=8),
                         panel.border = element_blank(),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.line = element_line(colour = "black"))

ecen <- ggplot() + geom_boxplot(aes(x=factor(patch.size),y=sr,
                                fill='light green',color=aggregate),
                                data=agg_comp[agg_comp$treat=='EC-EN',]) +
                   scale_fill_manual(values='light green',guide="none")+
                   scale_color_manual(values=c('red','blue'))+
                   labs(x=expression('Patch Size ['*m^2*']'),y='SR',title='EC-EN') +
                   ylim(0,12.5) +
                   theme(text = element_text(size=12),
                         axis.text.x = element_text(size=8),
                         panel.border = element_blank(),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.line = element_line(colour = "black"))

fedge <- acan + acen + ecan + ecen + plot_annotation(tag_levels="a",tag_suffix=")") +
                                     plot_layout(guides="collect")

ggsave('figureSX_edge.pdf',fedge,width=7.5,height=7.5,units="in")

agg_msd <- agg_comp %>% group_by(patch.size,treat,aggregate) %>% summarize(msr=mean(sr),sdsr=sd(sr))

acanm <- agg_msd[agg_msd$treat=='AC-AN','msr']
acan_diff <- acanm[seq(2,16,2),] - acanm[seq(1,15,2),]

acenm <- agg_msd[agg_msd$treat=='AC-EN','msr']
acen_diff <- acenm[seq(2,16,2),] - acenm[seq(1,15,2),]

ecanm <- agg_msd[agg_msd$treat=='EC-AN','msr']
ecan_diff <- ecanm[seq(2,16,2),] - ecanm[seq(1,15,2),]

ecenm <- agg_msd[agg_msd$treat=='EC-EN','msr']
ecen_diff <- ecenm[seq(2,16,2),] - ecenm[seq(1,15,2),]

tr_ed_diff <- round(cbind(acan_diff,acen_diff,ecan_diff,ecen_diff),2)
rownames(tr_ed_diff) <- c(0.04,0.16,0.36,0.64,1,1.44,1.96,2.56)
names(tr_ed_diff) <- c('acan','acen','ecan','ecen')

