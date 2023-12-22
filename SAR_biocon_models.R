library('plyr')
library('dplyr') 
library('minpack.lm')
library('ggplot2')
library('patchwork')

bcdat <- read.csv('hi_res_SR16_biocon.csv')
sr_bc_trt <- read.csv('sr_df.csv')
load('perm_boot.Rdata')
plots <- unique(bcdat$Plot)
plots <- plots[order(plots)]

trt_vec <- paste(bcdat$CO2Treatment,bcdat$Ntreatment,sep='-')
shrt_trt <- unique(trt_vec)

set.seed(23)

##### begin modeling and plotting code
base_mod <- list()
boot_mod <- list()
base_modH <- list()
boot_modH <- list()
base_modD <- list()
boot_modD <- list()
base_modBP <- list()
boot_modBP <- list()

for (i in seq(1,4)){
    mod <- data.frame()
    modD <- data.frame()
    modH <- data.frame()
    modBP <- data.frame()
    pl <- sr_bc_trt$treat==shrt_trt[i]
    # richness
    mod.base <- nlsLM("sr~c*patch.size^z",data=sr_bc_trt[pl,],start=list(c=5,z=0.2))
    base_mod[[i]] <- coefficients(mod.base)
    pred <- base_mod[[i]]['c']*sr_bc_trt[pl,]$patch.size^base_mod[[i]]['z']
    base_mod[[i]]['r2'] <- cor(sr_bc_trt[pl,]$sr,pred)^2 
    # shannon
    mod.baseH <- nlsLM("H~c*patch.size^z",data=sr_bc_trt[pl,],start=list(c=3,z=0.1))
    base_modH[[i]] <- coefficients(mod.baseH)
    predH <- base_modH[[i]]['c']*sr_bc_trt[pl,]$patch.size^base_modH[[i]]['z']
    base_modH[[i]]['r2'] <- cor(sr_bc_trt[pl,]$H,predH)^2 
    # simpson
    mod.baseD <- nlsLM("D~c*patch.size^z",data=sr_bc_trt[pl,],start=list(c=3,z=0.05))
    base_modD[[i]] <- coefficients(mod.baseD)
    predD <- base_modD[[i]]['c']*sr_bc_trt[pl,]$patch.size^base_modD[[i]]['z']
    base_modD[[i]]['r2'] <- cor(sr_bc_trt[pl,]$D,predD,use='c')^2 
    # berger-parker
    mod.baseBP <- nlsLM("BP~c*patch.size^z",data=sr_bc_trt[pl,],start=list(c=2,z=0.05))
    base_modBP[[i]] <- coefficients(mod.baseBP)
    predBP <- base_modBP[[i]]['c']*sr_bc_trt[pl,]$patch.size^base_modBP[[i]]['z']
    base_modBP[[i]]['r2'] <- cor(sr_bc_trt[pl,]$BP,predBP,use='c')^2 
    for (b in seq(1,1000)){
        bootdat <- slice_sample(sr_bc_trt[pl,],prop=1,replace=TRUE)
        mod.out <- nlsLM("sr~c*patch.size^z",data=bootdat,start=list(c=5,z=0.2))
        mod <- rbind(mod,coefficients(mod.out))
        mod.out <- nlsLM("D~c*patch.size^z",data=bootdat,start=list(c=4,z=0.1))
        modD <- rbind(modD,coefficients(mod.out))
        mod.out <- nlsLM("H~c*patch.size^z",data=bootdat,start=list(c=3,z=0.05))
        modH <- rbind(modH,coefficients(mod.out))
        mod.out <- nlsLM("BP~c*patch.size^z",data=bootdat,start=list(c=0.3,z=0))
        modBP <- rbind(modBP,coefficients(mod.out))
    }
    colnames(mod) <- c('c','z')
    boot_mod[[i]] <- mod
    colnames(modD) <- c('c','z')
    boot_modD[[i]] <- modD
    colnames(modH) <- c('c','z')
    boot_modH[[i]] <- modH
    colnames(modBP) <- c('c','z')
    boot_modBP[[i]] <- modBP
}

moddf <- data.frame(matrix(unlist(base_mod),nrow=4,byrow=TRUE))
colnames(moddf) <- c('c','z','r2')
modHdf <- data.frame(matrix(unlist(base_modH),nrow=4,byrow=TRUE))
colnames(modHdf) <- c('c','z','r2')
modDdf <- data.frame(matrix(unlist(base_modD),nrow=4,byrow=TRUE))
colnames(modDdf) <- c('c','z','r2')
modBPdf <- data.frame(matrix(unlist(base_modBP),nrow=4,byrow=TRUE))
colnames(modBPdf) <- c('c','z','r2')

# modeled richness dataframe
x <- seq(0.01,3.24,.01) 
xf <- data.frame(treat=NULL,x=NULL,y=NULL)
boot_conf <- list()
xfq <- data.frame(treatment=NULL,x=NULL,yl=NULL,yh=NULL)
# modeled Shannon Diversity
xfH <- data.frame(treat=NULL,x=NULL,y=NULL)
boot_confH <- list()
xfqH <- data.frame(treatment=NULL,x=NULL,yl=NULL,yh=NULL)
# modeled Simpson Diversity
xfD <- data.frame(treat=NULL,x=NULL,y=NULL)
boot_confD <- list()
xfqD <- data.frame(treatment=NULL,x=NULL,yl=NULL,yh=NULL)
# modeled Berger-Parker Diversity
xfBP <- data.frame(treat=NULL,x=NULL,y=NULL)
boot_confBP <- list()
xfqBP <- data.frame(treatment=NULL,x=NULL,yl=NULL,yh=NULL)
for (i in seq(1,4)) {
    mc <- base_mod[[i]]
    y <- mc[1]*x^mc[2]
    treat <- rep(shrt_trt[i],length(x))
    df <- data.frame(treatment=treat,x=x,y=y)
    xf <- rbind(xf,df)
    xfb <- matrix(nrow=1000,ncol=324)
    # Shannon Diversity
    mcH <- base_modH[[i]]
    yH <- mcH[1]*x^mcH[2]
    treatH <- rep(shrt_trt[i],length(x))
    dfH <- data.frame(treatment=treatH,x=x,y=yH)
    xfH <- rbind(xfH,dfH)
    xfbH <- matrix(nrow=1000,ncol=324)
    # Simpson Diversity
    mcD <- base_modD[[i]]
    yD <- mcD[1]*x^mcD[2]
    treatD <- rep(shrt_trt[i],length(x))
    dfD <- data.frame(treatment=treatD,x=x,y=yD)
    xfD <- rbind(xfD,dfD)
    xfbD <- matrix(nrow=1000,ncol=324)
    # Berger-Parker Diversity
    mcBP <- base_modBP[[i]]
    yBP <- mcBP[1]*x^mcBP[2]
    treatBP <- rep(shrt_trt[i],length(x))
    dfBP <- data.frame(treatment=treatBP,x=x,y=yBP)
    xfBP <- rbind(xfBP,dfBP)
    xfbBP <- matrix(nrow=1000,ncol=324)
    for (j in seq(1,1000)){
        mc <- boot_mod[[i]][j,]
        y <- mc[1,1]*x^mc[1,2]
        xfb[j,] <- y
        mcD <- boot_modD[[i]][j,]        
        yD <- mcD[1,1]*x^mcD[1,2]
        xfbD[j,] <- yD
        mcH <- boot_modH[[i]][j,]
        yH <- mcH[1,1]*x^mcH[1,2]
        xfbH[j,] <- yH
        mcBP <- boot_modBP[[i]][j,]
        yBP <- mcBP[1,1]*x^mcBP[1,2]
        xfbBP[j,] <- yBP
    }
    boot_conf[[i]] <- xfb
    tmp <- apply(boot_conf[[i]],2,quantile,probs=c(0.0275,0.975))
    df <- data.frame(treatment=treat,x=x,yl=tmp[1,],yh=tmp[2,])
    xfq <- rbind(xfq,df)
    boot_confD[[i]] <- xfbD
    tmp <- apply(boot_confD[[i]],2,quantile,probs=c(0.0275,0.975))
    dfD <- data.frame(treatment=treatD,x=x,yl=tmp[1,],yh=tmp[2,])
    xfqD <- rbind(xfqD,dfD)
    boot_confH[[i]] <- xfbH
    tmp <- apply(boot_confH[[i]],2,quantile,probs=c(0.0275,0.975))
    dfH <- data.frame(treatment=treatH,x=x,yl=tmp[1,],yh=tmp[2,])
    xfqH <- rbind(xfqH,dfH)
    boot_confBP[[i]] <- xfbBP
    tmp <- apply(boot_confBP[[i]],2,quantile,probs=c(0.0275,0.975))
    dfBP <- data.frame(treatment=treatBP,x=x,yl=tmp[1,],yh=tmp[2,])
    xfqBP <- rbind(xfqBP,dfBP)
}    

## plotting functions
# model plot, log-log
ln_plot <- function(line_model,conf_model,ylow,yhigh,y_label,leg.pos='none'){
    p <- ggplot() + geom_line(aes(x=x,y=y,color=treatment,linetype=treatment),
                                data=line_model,size=1) + 
                    scale_color_manual(values=c('black','green','black','green'))+
                    scale_linetype_manual(values=c('solid','solid','dashed','dashed'))+
                    geom_ribbon(aes(x=x,ymin=yl,ymax=yh,color=treatment,linetype=treatment),
                                data=conf_model,alpha=0.1) +
                    labs(x=expression('Patch Size ['*m^2*']'),
                         y=y_label) +
                   scale_x_log10()+scale_y_log10(limits=c(ylow,yhigh))+
                    theme(text = element_text(size=16),
                          axis.text.x = element_text(size=12),
                          axis.text.y = element_text(size=12),
                          panel.border = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank(),
                          axis.line = element_line(colour = "black"),
                          legend.position=leg.pos)
}

# model plot, linear
lnn_plot <- function(line_model,conf_model,y_label,leg.pos='none'){
    p <- ggplot() + geom_line(aes(x=x,y=y,color=treatment,linetype=treatment),
                                data=line_model,size=1) +
                    scale_color_manual(values=c('black','green','black','green'))+
                    scale_linetype_manual(values=c('solid','solid','dashed','dashed'))+
                    geom_ribbon(aes(x=x,ymin=yl,ymax=yh,color=treatment,linetype=treatment),
                                data=conf_model,alpha=0.1) +
                    labs(x=expression('Patch Size ['*m^2*']'),
                         y=y_label) +
                    theme(text = element_text(size=16),
                          axis.text.x = element_text(size=12),
                          axis.text.y = element_text(size=12),
                          panel.border = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank(),
                          axis.line = element_line(colour = "black"),
                          legend.position=leg.pos)
}

# data box plot
bx_plot <- function(rich_met,indata,y_label){
    rich_met <- sym(rich_met)
    p <- ggplot() + geom_boxplot(aes(x=factor(patch.size),y=!!rich_met,fill=factor(treat)),
                                data=indata) + 
                   scale_fill_manual(values=c('black','green','grey','light green'))+
                   labs(x=expression('Patch Size ['*m^2*']'),
                        y=y_label) +
                   theme(text = element_text(size=12),
                         axis.text.x = element_text(size=8),
                         panel.border = element_blank(),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.line = element_line(colour = "black"))
}

# data and linear model plot
bx_ln_plot <- function(ln_model,conf_model,rich_met,indata,y_label){
    rich_met <- sym(rich_met)
    p <- ggplot() + geom_line(aes(x=x,y=y,color=treatment,linetype=treatment),
                                data=ln_model,size=1) +
                    scale_color_manual(values=c('black','green','black','green'))+
                    scale_linetype_manual(values=c('solid','solid','dashed','dashed'))+
                    geom_ribbon(aes(x=x,ymin=yl,ymax=yh,color=treatment,linetype=treatment),
                                data=conf_model,alpha=0.1) +
                    geom_boxplot(inherit.aes=FALSE,width=0.1,aes(x=patch.size,y=!!rich_met,
                                 group=interaction(factor(patch.size),factor(treat)),fill=factor(treat)),
                                data=indata) + 
                   scale_fill_manual(values=c('black','green','grey','light green'))+
                   labs(x=expression('Patch Size ['*m^2*']'),
                        y=y_label) +
                   theme(text = element_text(size=12),
                         axis.text.x = element_text(size=8),
                         panel.border = element_blank(),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.line = element_line(colour = "black"))
}

ndat <- bcdat %>% filter(Plot==163 & Species!="Bargr")
f1a <-  ggplot(ndat, aes(Row, Column))+geom_point(aes(col=Species,size=Coverage),shape=1)+theme_classic()+ 
        scale_x_continuous(limits=c(0,18))+ scale_y_continuous(limits=c(0,18))#+
        scale_color_viridis_d(option="magma") #(for a color blind friendly palette)
ggsave('figure2a.pdf',f1a)

# useful stats reported in manuscript
# average difference between control and a treatment
# mean((xf[xf$treatment=="AC-AN","y"]-xf[xf$treatment=="EC-EN","y"])/xf[xf$treatment=="AC-AN","y"])

# Species Richness
f2a <- lnn_plot(xf,xfq,paste0('Species Richness [',pres_cov,']'),c(0.75,0.25))
f2b <- ln_plot(xf,xfq,1.5,10,'')
f2 <- f2a + f2b + plot_annotation(tag_levels="a",tag_suffix=")")
ggsave(paste0('figure2_',pres_cov,'.pdf'),f2,width=6.5,height=4,units="in")

fs1a <- bx_plot('sr',sr_bc,'SR')
fs2a <- bx_ln_plot(xf,xfq,'sr',sr_bc,'SR')

## Cover
pres_cov <- 'cover'
# Shannon Diversity
f3ca <- ln_plot(xfH,xfqH,1,4.5,paste0("Shannon Div."),c(0.75,0.275))
# Simpson Index
f3cb <- ln_plot(xfD,xfqD,1,4.5,paste0("Inv. Simpson"))
# Berger-Parker Index
f3cc <- ln_plot(xfBP,xfqBP,1,4.5,paste0("Inv. Berger-Parker"))

f3 <- f3ca + f3cb + f3cc + plot_annotation(tag_levels="a",tag_suffix=")") 

ggsave('figure5.pdf',f3,width=8,height=4,units="in")

# Shannon Diversity
fs1b <- bx_plot('H',sr_bc,"H")
fs2b <- bx_ln_plot(xfH,xfqH,'H',sr_bc,"H")

# Simpson Index
fs1c <- bx_plot('D',sr_bc,"1/D")
fs2c <- bx_ln_plot(xfD,xfqD,'D',sr_bc,"1/D")

# Berger-Parker Index
fs1d <- bx_plot('BP',sr_bc,"1/BP")
fs2d <- bx_ln_plot(xfBP,xfqBP,'BP',sr_bc,"1/BP")

# assemble plots
f2 <- f2a + f2b + plot_annotation(tag_levels="a",tag_suffix=")")+plot_layout(guides="collect")
fs1 <- fs1a / fs1b / fs1c / fs1d + 
        plot_annotation(title=paste0("Diversity metrics using species ",pres_cov),tag_levels="a",tag_suffix=")")  
        plot_layout(guides="collect")
fs2 <- fs2a / fs2b / fs2c / fs2d + 
        plot_annotation(title=paste0("Diversity metrics using species ",pres_cov),tag_levels="a",tag_suffix=")") + 
        plot_layout(guides="collect")

# supplemental figures
ggsave(paste0('figure_s1_',pres_cov,'.pdf'),fs1,width=6.5,height=6.5,units="in")
ggsave(paste0('figure_s2_',pres_cov,'.pdf'),fs2,width=6.5,height=6.5,units="in")

# Richness compared
spc_comp <- function(treat,col,y_label){
    p <- ggplot() + geom_line(aes(x=x,y=y),
                                data=xf[xf$treatment %in% treat,],
                                color=col,linetype='solid',size=1) +
                    geom_line(aes(x=x,y=y),                             
                                data=xfbq[xfbq$treatment %in% treat,],size=1,
                                color=col,linetype='dotted',alpha=0.5) +
                    scale_x_log10()+scale_y_log10()+
                    labs(title=treat,
                        x=expression('Patch Size ['*m^2*']'),
                         y=y_label) +
                    theme(text = element_text(size=16),
                          axis.text.x = element_text(size=12),
                          axis.text.y = element_text(size=12),
                          panel.border = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank(),
                          axis.line = element_line(colour = "black"))
}

cols <- c('green','black','green','black')
typs <- c('solid','solid','dashed','dashed')

f5a <- spc_comp(shrt_trt[2],cols[2],paste0("Species Richness"))
f5b <- spc_comp(shrt_trt[1],cols[1],"")
f5c <- spc_comp(shrt_trt[4],cols[4],paste0("Species Richness"))
f5d <- spc_comp(shrt_trt[3],cols[3],"")

f5 <- f5a + f5b + f5c + f5d + plot_annotation(tag_levels="a",tag_suffix=")") + plot_layout(guides="collect")

ggsave('figureSX.pdf',f5,width=6.5,height=6.5,units="in")

sr_sp_diff <- list()
x <- seq(0.01,3.24,.01)
sp_diff <- data.frame()
# spatial sorting difference figures
for (i in 1:4) {
    s1 <- sample(1:1000,replace=TRUE)
    s2 <- sample(1:1000,replace=TRUE)
    sr_sp_diff[[i]] <- big_boot_conf[[i]][s1,] - boot_conf[[i]][s2,]
    tmp <- apply(sr_sp_diff[[i]],2,quantile,probs=c(0.0275,0.5,0.975))
    treat <- rep(shrt_trt[i],length(x))
    df <- data.frame(treatment=treat,x=x,y=tmp[2,],yl=tmp[1,],yh=tmp[3,])
    sp_diff <- rbind(sp_diff,df)
}

sp_nosp_p <- ggplot() + geom_line(aes(x=x,y=y,color=treatment,linetype=treatment),
                            data=xfbq,linewidth=1) +
                scale_color_manual(values=c('black','green','black','green'))+
                scale_linetype_manual(values=c('solid','solid','dashed','dashed'))+
                labs(x=expression('Patch Size ['*m^2*']'),
                     y='Richness, No Spatial Sorting') +
                scale_x_log10()+scale_y_log10()+
                theme(text = element_text(size=16),
                      axis.text.x = element_text(size=12),
                      axis.text.y = element_text(size=12),
                      panel.border = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black"),
                      legend.position=c(0.75,0.275))

sp_diff_p <- ggplot() + geom_line(aes(x=x,y=y,color=treatment,linetype=treatment),
                            data=sp_diff,linewidth=1) +
                scale_color_manual(values=c('black','green','black','green'))+
                scale_linetype_manual(values=c('solid','solid','dashed','dashed'))+
                geom_ribbon(aes(x=x,ymin=yl,ymax=yh,color=treatment),data=sp_diff,alpha=0.1) +
                labs(x=expression('Patch Size ['*m^2*']'),
                     y='Species Gain, No Spatial Sorting') +
                scale_x_log10()+scale_y_log10()+
                theme(text = element_text(size=16),
                      axis.text.x = element_text(size=12),
                      axis.text.y = element_text(size=12),
                      panel.border = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black"),
                      legend.position='none')

f3ns <- sp_nosp_p + sp_diff_p + plot_annotation(tag_levels="a",tag_suffix=")") 

ggsave('figure3.pdf',f3ns,width=7.5,height=4,units="in")

### append permuted plots and extend estimated sr to aggregate
agg_plots_trt <- data.frame('patch.size'=NULL,'sr'=NULL,'treat'=NULL)
for (tr in 1:4){
    agg_plots <- data.frame('patch.size'=NULL,'sr'=NULL,'treat'=NULL)
    for (p in 1:dim(permL[[tr]])[2]) {
        df <- data.frame('patch.size'=rep(3.24*p,dim(permL[[tr]])[1]),
                         'sr'=permL[[tr]][,p],'treat'=shrt_trt[tr])
        agg_plots <- rbind(agg_plots,df)
    }
    agg_plots_trt <- rbind(agg_plots_trt,agg_plots)
}

# modeled richness dataframe 
xfL <- data.frame(treat=NULL,x=NULL,y=NULL)
boot_confL <- list()
xfqL <- data.frame(treatment=NULL,x=NULL,yl=NULL,yh=NULL)
for (i in seq(1,4)) {
    x <- seq(3.24,3.24*dim(permL[[i]])[2],length.out=100)
    mc <- base_mod[[i]]
    y <- mc[1]*x^mc[2]
    treat <- rep(shrt_trt[i],length(x))
    df <- data.frame(treatment=treat,x=x,y=y)
    xfL <- rbind(xfL,df)
    xfb <- matrix(nrow=1000,ncol=100)
    for (j in seq(1,1000)){
        mc <- boot_mod[[i]][j,]
        y <- mc[1,1]*x^mc[1,2]
        xfb[j,] <- y
    }
    boot_confL[[i]] <- xfb
    tmp <- apply(boot_confL[[i]],2,quantile,probs=c(0.0275,0.975))
    df <- data.frame(treatment=treat,x=x,yl=tmp[1,],yh=tmp[2,])
    xfqL <- rbind(xfqL,df)
}
 
f4 <- ggplot() + geom_line(aes(x=x,y=y,color=treatment,linetype=treatment),
                            data=xf,size=1) +
                scale_color_manual(values=c('black','green','black','green'))+
                scale_linetype_manual(values=c('solid','solid','dashed','dashed'))+
                geom_ribbon(aes(x=x,ymin=yl,ymax=yh,color=treatment,linetype=treatment),
                            data=xfq,alpha=0.1) +
                ### extrapolated estimates               
                geom_line(aes(x=x,y=y,color=treatment,linetype=treatment),
                            data=xfL,size=1,show.legend=FALSE) +
                scale_color_manual(values=c('black','green','black','green'))+
                scale_linetype_manual(values=c('solid','solid','dashed','dashed'))+
                geom_ribbon(aes(x=x,ymin=yl,ymax=yh,color=treatment,linetype=treatment),
                            data=xfqL,alpha=0.1,show.legend=FALSE) +
                geom_boxplot(inherit.aes=FALSE,width=1,aes(x=patch.size,y=sr,
                             group=interaction(factor(patch.size),factor(treat)),fill=factor(treat),alpha=0.3),
                            data=agg_plots_trt,show.legend=FALSE) + 
                scale_fill_manual(values=c('black','green','grey','light green')) +
               labs(x=expression('Patch Size ['*m^2*']'),
                    y=paste0('Species Richness [',pres_cov,']')) +
               theme(text = element_text(size=12),
                     axis.text.x = element_text(size=8),
                     panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     axis.line = element_line(colour = "black"))

ggsave(paste0('figure6_',pres_cov,'.pdf'),f4,width=6.5,height=4,units="in")
