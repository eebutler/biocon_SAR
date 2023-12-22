library('dplyr')
library('combinat')
library('minpack.lm')

# data processing code
bcdat <- read.csv('hi_res_SR16_biocon.csv')
plots <- unique(bcdat$Plot)
plots <- plots[order(plots)]
sr_plot <- vector("list",length=length(plots))
spec_plot <- vector("list",length=length(plots))

# flag for whether presence or cover should be used
pres_cov <- 'cover'

set.seed(23)

# Diversity and evenness from Chao and Ricotta, added in Berger-Parker for Infinite q
qD <- function(p,q){
  p <- p[p>0]
  if(q==1){
    exp(-sum(p*log(p)))
  }else if(is.finite(q)){
    (sum(p^q))^(1/(1-q))
  }else{
    1/max(p)
  }
}

# get lists of species richness across all sampling squares
for (p in seq(1,length(plots))) {
    pl <- which(bcdat$Plot==plots[p])
    bc <- bcdat[pl,]
    # only count rooted plants?
    #r <- bc$Root == "R"
    #bc <- bc[r,]
    sr_list <- vector("list",length=18)
    for (i in seq(1,18)) {
        total_sq <- (19-i)^2
        sq_area <- i^2*0.01
        sr <- vector(length=total_sq,mode="numeric")
        D <- vector(length=total_sq,mode="numeric")
        H <- vector(length=total_sq,mode="numeric")
        BP <- vector(length=total_sq,mode="numeric")
        np <- vector(length=total_sq,mode="numeric")
        idx <- 1 
        for (row in seq(1,19-i)) {
            for (col in seq(1,19-i)) {
                r <- match(bc$Row,seq(row,row+i-1))
                c <- match(bc$Column,seq(col,col+i-1))
                # table of species to tabulate richness
                s_tab <- table(bc[!is.na(r) & !is.na(c),'Species'])
                spec <- names(s_tab)!="Bargr"
                ste <- s_tab[spec]
                # probability of finding a species within a quadrat
                p_s <- ste/sum(ste)
                num_plants <- sum(ste)
                ## alternately species cover within each area
                s_cover <- bc[!is.na(r) & !is.na(c),c('Species','Coverage')]
                s_cover$area <- s_cover$Coverage * 0.01^2 # convert to percent then m^2
                tot_area <- s_cover %>% group_by(Species) %>% summarize(total_area = sum(area))
                tot_area <- tot_area[tot_area$Species !="Bargr",]
                p_c <- tot_area$total_area/sum(tot_area$total_area)
                names(p_c) <- names(ste)
                # get species richness from Hill Numbers function above
                if (pres_cov=='presence'){
                    sr[idx] <- qD(p_s,0)
                    H[idx] <- qD(p_s,1)
                    D[idx] <- qD(p_s,2)
                    BP[idx] <- qD(p_s,Inf)
                } else {
                    sr[idx] <- qD(p_c,0)               
                    H[idx] <- qD(p_c,1)
                    D[idx] <- qD(p_c,2)
                    BP[idx] <- qD(p_c,Inf)
                }
                np[idx] <- num_plants
                # set infinite values to NA
                D[D==Inf] <- NA
                BP[BP==-Inf] <- NA  
                if (i==1) {
                    spec_plot[[p]][[idx]] <- names(ste)
                }
                idx <- idx + 1
            }
        }
        sr_list[[i]]$sr <- sr
        sr_list[[i]]$D <- D
        sr_list[[i]]$H <- H
        sr_list[[i]]$BP <- BP
        sr_list[[i]]$p_c <- p_c
        sr_list[[i]]$np <- np
    }
    print(p)
    sr_plot[[p]] <- sr_list
}

ring1 <- c(7,11,12,13,16,17,21,23,26,33,34,43,45,48,50,51,52)
ring2 <- c(64,65,66,67,69,83,84,87,97,104,106,107,108,110,117,118,119)
ring3 <- c(128,130,131,133,147,149,151,161,163,165,166,170,171,173,175,177,179,180)
ring4 <- c(184,185,186,187,188,194,201,205,213,222,224,226,230,232,233,234,236,243)

rp <- data.frame('ring'=c(rep(1,17),rep(2,17),rep(3,18),rep(4,18)),'plot'=c(ring1,ring2,ring3,ring4))

sr_plots_df <- data.frame('ring'=NULL,'plot'=NULL,'patch.size'=NULL,'sr'=NULL,'H'=NULL,'D'=NULL,'BP'=NULL)
# convert list to dataframe with ring and plot information
for (p in seq(1,length(sr_plot))){
    for (ps in seq(1,18)){
        sr <- sr_plot[[p]][[ps]]$sr
        D <- sr_plot[[p]][[ps]]$D
        H <- sr_plot[[p]][[ps]]$H
        BP <- sr_plot[[p]][[ps]]$BP
        np <- sr_plot[[p]][[ps]]$np
        patch_size <- rep(ps^2*.01,length(sr))
        plt <- rep(plots[p],length(sr))
        ring <- rep(rp$ring[which(rp$plot==plt[1])],length(sr)) 
        df <- data.frame('ring'=ring,'plot'=plt,'patch.size'=patch_size,'sr'=sr,'D'=D,'H'=H,'BP'=BP)
        sr_plots_df <- rbind(sr_plots_df,df)
    }
}

# treatment vector
trt_vec <- paste(bcdat$CO2Treatment,bcdat$Ntreatment,sep='-')
shrt_trt <- unique(trt_vec)

# organize treatments and plots in BioCON
pt16 <- list()
for (i in seq(1,4)){
    pt16[[i]] <- unique(bcdat$Plot[trt_vec==shrt_trt[i]])
}
sr_plots_df$treat <- NA
for (i in 1:4){
    sr_plots_df$treat[sr_plots_df$plot %in% pt16[[i]]] <- shrt_trt[i]
}

# remove NAs and group by treatment and patch size
sr_bc <- sr_plots_df[!is.na(sr_plots_df$treat),]
sr_bc_trt <- sr_bc %>% group_by(treat,patch.size)

# save data frame for use in modeling and plotting code
write.csv(sr_bc_trt,'sr_df.csv')

# extending to larger areas
permL <- list()
for (tr in seq(1,4)){
    pl <- which(trt_vec==shrt_trt[tr])
    tplots <- unique(bcdat$Plot[pl]) 
    tpl <- which(plots %in% tplots)
    
    # aggregate species lists
    sbplot <- list() 
    for (t in 1:length(tpl)){
        sbplot[[t]] <- names(sr_plot[[tpl[t]]][[18]]$p_c)
    }
    num_p <- length(tpl)
    ordL <- permn(seq(1,num_p),num_p,fun=NULL)
    perm_mat <- matrix(nrow=length(ordL),ncol=num_p)
    for (p in 1:length(ordL)){
        ord <- ordL[[p]]
        for (i in 1:num_p){
            perm_mat[p,i] <- length(unique(unlist(sbplot[ord[1:i]])))
        }
    }
   permL[[tr]] <- perm_mat 
}

### note the analysis below is very slow and can take many hours to run
### on a desktop machine
# ignoring spatial sorting - much more rapid increase in SR
big_boot <- list()
for (tr in seq(1,4)){
    pl <- which(trt_vec==shrt_trt[tr])
    tplots <- unique(bcdat$Plot[pl]) 
    tpl <- which(plots %in% tplots)
    
    # bootstrap generate larger plot sizes, and uncertainties
    tmp <- spec_plot[tpl]
    bplot <- list() # aggregate all 0.01 patches for each treatment into a single list
    for (t in seq(1,length(tmp))){
        bplot <- c(bplot,tmp[[t]])
    }
    tr_tab <- sr_bc_trt$treat==shrt_trt[tr]
    rep_vec <- as.vector(table(sr_bc_trt[tr_tab,]$patch.size))

    ## plot sizes, up to sqrt(length(bplot))
    SR_mod_coef <- data.frame()
    for (i in seq(1,1000)){
        if(i%%10==0){print(i)}
        SR_nospace <- data.frame()
        for (j in seq(1,18)){ 
            for (k in seq(1,rep_vec[j])){
                out <- sample(bplot,j^2,replace=TRUE)
                bptab <- table(unlist(out))
                spec <- names(bptab)!="Bargr"
                # Number of species "richness"
                sr <- length(bptab[spec])
                # table of species richness
                ste <- bptab[spec]
                # probability of finding a species
                p_s <- ste/sum(ste)
                # Diversity Metrics
                sr <- qD(p_s,0)
                bp_dat <- data.frame('sr'=sr,patch.size=j^2*.01)
                SR_nospace <- rbind(SR_nospace,bp_dat)
            }
        }
    SR.mod <- nlsLM("sr~c*patch.size^z",data=SR_nospace,start=list(c=5,z=0.2))
    SR_mod_coef <- rbind(SR_mod_coef,coefficients(SR.mod))
    }
    colnames(SR_mod_coef) <- c('c','z')
    big_boot[[tr]] <- SR_mod_coef
}

save(list=c('big_boot','permL'),file='perm_boot.RData')
