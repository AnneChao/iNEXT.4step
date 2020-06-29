library(Rcpp)
library(dplyr)
library(ggplot2)
library(iNEXT)
library(ggpubr)
source("Taxonomic_profile.r")
source("EstimateD_general.r")
myggeven <- function(out){
  ans <- ggplot(data = out)+geom_line(aes(x = q, y = Evenness, col = Community),size = 1.1)+theme_bw()+
    geom_point(data = filter(out,q %in% c(0,1,2)),aes(x = q, y = Evenness, col = Community),show.legend = FALSE, size = 3)+
    labs(x="Order q")+ggtitle("(e) Evenness profiles")+
    theme(text=element_text(size=14),legend.position = "bottom",legend.title=element_blank(),
          legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,-10,-5,-10),
          plot.margin = margin(9.5, 5.5, 5.5, 5.5, "pt"),
          legend.key.width = unit(2,"line"))
  return(ans)
}


nms <- c("Fossil","Spider","Vegetation","Stony coral")
dts <- lapply(nms,function(x) read.table(paste0("E:/YanHan/Japan_meeting(incl 4 steps)/",x," data.txt")))
names(dts) <- nms

q = seq(0,2,0.1)
B <- 10
#=====Fig_1=====
x <- dts$Fossil %>% .[rowSums(.)>0,]

out1 <- apply(x,2,function(i) sc_profile(freq = i,datatype = "abundance",q = q,B = B,conf = 0.95)) %>% 
  do.call(rbind,.)
out1$Community <- rep(colnames(x),each=length(q))
out1_p <- plot_sc_profile(data = out1) + theme_bw() + ylab("Sample completeness")+
  #geom_text(data = filter(out1,Order.q %in% c(0,1,2)),aes(x = Order.q, y = Estimate, label=round(Estimate,2),col = Community),nudge_x = 0.05,nudge_y = -0.05,show.legend = FALSE)+
  geom_point(data = filter(out1,Order.q %in% c(0,1,2)),aes(x = Order.q, y = Estimate, col = Community),show.legend = FALSE, size = 3)+
  theme(text = element_text(size=14), legend.position = "bottom",legend.title = element_blank()) + ggtitle("(a) Sample completeness profiles")
out3_est <- apply(x,2,function(i) MakeTable_Proposeprofile(data = i,B = B,q =q,conf = 0.95)) %>% 
  do.call(rbind,.) %>% mutate(Community = rep(colnames(x),each=length(q)), method = "Estimate") %>% 
  rename(Diversity = Estimate)
out3_mle <- apply(x,2,function(i) MakeTable_Empericalprofile(data = i,B = B,q =q,conf = 0.95)) %>% 
  do.call(rbind,.) %>% mutate(Community = rep(colnames(x),each=length(q)), method = "Empirical") %>% 
  rename(Diversity = Emperical)
out3_mle$LCL <- out3_mle$UCL <- out3_mle$Diversity
out3 <- rbind(out3_est,out3_mle)
out3$method[out3$method=="Estimate"] <- 'Asymptotic'
out3$method <- factor(out3$method,levels = unique(out3$method))
out3_p <- plot_diversity_profile_yhc(output = out3) + theme_bw() + ggtitle("(c) Asymptotic and empirical diversity profiles")+
  geom_point(data = filter(out3,Order.q %in% c(0,1,2)),aes(x = Order.q, y = Diversity, col = Community),show.legend = FALSE, size = 3)+
  theme(legend.position = "bottom",legend.key.width = unit(2,"line"),text = element_text(size=14),
        legend.title=element_blank(),
        legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,-10,-5,-10)) + ylab("Species diversity")
# + coord_cartesian(ylim = c(6, 70)) 
out24 <- iNEXT::iNEXT(x = x,q = c(0,1,2),datatype = "abundance",nboot = B,endpoint = 4800)
out24$iNextEst$Miocene <- out24$iNextEst$Miocene %>% filter(m<=1500)
out24$iNextEst <- lapply(out24$iNextEst,function(out){
  out$method[out$method =='interpolated'] <- 'Interpolated'
  out$method[out$method =='extrapolated'] <- 'Extrapolated'
  out
})
out2_p <- ggiNEXT.iNEXT(x = out24,type = 1,facet.var = "order")+ ggtitle("(b) Size-based rarefaction/extrapolation")+
  theme(text = element_text(size=14),legend.position = "none") + coord_cartesian(xlim = c(0 , 5000))+
  scale_x_continuous(breaks =c(0,1500,3000,4500))
out4_p <- ggiNEXT.iNEXT(x = out24,type = 3,facet.var = 'order')+ggtitle("(d) Coverage-based rarefaction/extrapolation")+
  theme(text = element_text(size=14),legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,-10,-5,-10))+
  #annotate("text", x=0.5, y=100, label=c("C[max] == 98.9*'%' "), parse=TRUE) + 
  scale_x_continuous(breaks = c(0.25,0.5,0.75,1),labels = c("0.25","0.50","0.75","1.00"))
SC <- min(sapply(x, function(x) iNEXT:::Chat.Ind(x, 2*sum(x))))
out5 <- estimateDyhc(x = x,q = q,datatype = "abundance",base = "coverage",level = SC,conf = NULL)
out5 <- out5 %>% group_by(Community) %>% mutate(Evenness = (Diversity-1)/(Diversity[q==0]-1))
out5_p <- myggeven(out = out5)

out_p <- ggarrange(out1_p,out2_p,out3_p,out4_p,out5_p,ncol = 2, nrow = 3)
ggsave(filename = "Figure_1.pdf",plot = out_p,device = "pdf",width = 36,height = 25,
       units = "cm",dpi = 320)

#=====Fig_2=====
x <- dts$Spider %>% .[rowSums(.)>0,]

out1 <- apply(x,2,function(i) sc_profile(freq = i,datatype = "abundance",q = q,B = B,conf = 0.95)) %>% 
  do.call(rbind,.)
out1$Community <- rep(colnames(x),each=length(q))
out1_p <- plot_sc_profile(data = out1) + theme_bw() + ylab("Sample completeness")+
  #geom_text(data = filter(out1,Order.q %in% c(0,1,2)),aes(x = Order.q, y = Estimate, label=round(Estimate,2),col = Community),nudge_x = 0.05,nudge_y = -0.05,show.legend = FALSE)+
  geom_point(data = filter(out1,Order.q %in% c(0,1,2)),aes(x = Order.q, y = Estimate, col = Community),show.legend = FALSE, size = 3)+
  theme(text = element_text(size=14), legend.position = "bottom",legend.title = element_blank()) + ggtitle("(a) Sample completeness profiles")
out3_est <- apply(x,2,function(i) MakeTable_Proposeprofile(data = i,B = B,q =q,conf = 0.95)) %>% 
  do.call(rbind,.) %>% mutate(Community = rep(colnames(x),each=length(q)), method = "Estimate") %>% 
  rename(Diversity = Estimate)
out3_mle <- apply(x,2,function(i) MakeTable_Empericalprofile(data = i,B = B,q =q,conf = 0.95)) %>% 
  do.call(rbind,.) %>% mutate(Community = rep(colnames(x),each=length(q)), method = "Empirical") %>% 
  rename(Diversity = Emperical)
out3_mle$LCL <- out3_mle$UCL <- out3_mle$Diversity
out3 <- rbind(out3_est,out3_mle)
out3$method[out3$method=="Estimate"] <- 'Asymptotic'
out3$method <- factor(out3$method,levels = unique(out3$method))
out3_p <- plot_diversity_profile_yhc(output = out3) + theme_bw() + ggtitle("(c) Asymptotic and empirical diversity profiles")+
  geom_point(data = filter(out3,Order.q %in% c(0,1,2)),aes(x = Order.q, y = Diversity, col = Community),show.legend = FALSE, size = 3)+
  theme(legend.position = "bottom",legend.key.width = unit(2,"line"),text = element_text(size=14),
        legend.title=element_blank(),
        legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,-10,-5,-10)) + ylab("Species diversity")
# + coord_cartesian(ylim = c(6, 70)) 
out24 <- iNEXT::iNEXT(x = x,q = c(0,1,2),datatype = "abundance",nboot = B)
out24$iNextEst <- lapply(out24$iNextEst,function(out){
  out$method[out$method =='interpolated'] <- 'Interpolated'
  out$method[out$method =='extrapolated'] <- 'Extrapolated'
  out
})
out2_p <- ggiNEXT.iNEXT(x = out24,type = 1,facet.var = "order")+ ggtitle("(b) Size-based rarefaction/extrapolation")+
  theme(text = element_text(size=14),legend.position = "none")
out4_p <- ggiNEXT.iNEXT(x = out24,type = 3,facet.var = 'order')+ggtitle("(d) Coverage-based rarefaction/extrapolation")+
  theme(text = element_text(size=14),legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,-10,-5,-10))+
  #annotate("text", x=0.5, y=100, label=c("C[max] == 98.9*'%' "), parse=TRUE) + 
  scale_x_continuous(breaks = c(0.25,0.5,0.75,1),labels = c("0.25","0.50","0.75","1.00"))
SC <- min(sapply(x, function(x) iNEXT:::Chat.Ind(x, 2*sum(x))))
out5 <- estimateDyhc(x = x,q = q,datatype = "abundance",base = "coverage",level = SC,conf = NULL)
out5 <- out5 %>% group_by(Community) %>% mutate(Evenness = (Diversity-1)/(Diversity[q==0]-1))
out5_p <- myggeven(out = out5)

out_p <- ggarrange(out1_p,out2_p,out3_p,out4_p,out5_p,ncol = 2, nrow = 3)
ggsave(filename = "Figure_2.pdf",plot = out_p,device = "pdf",width = 36,height = 25,
       units = "cm",dpi = 320)
#=====Fig_3=====
x <- dts$Vegetation %>% .[rowSums(.)>0,]

out1 <- apply(x,2,function(i) sc_profile(freq = i,datatype = "incidence",q = q,B = B,conf = 0.95)) %>% 
  do.call(rbind,.)
out1$Community <- rep(colnames(x),each=length(q))
out1_p <- plot_sc_profile(data = out1) + theme_bw() + ylab("Sample completeness")+
  #geom_text(data = filter(out1,Order.q %in% c(0,1,2)),aes(x = Order.q, y = Estimate, label=round(Estimate,2),col = Community),nudge_x = 0.05,nudge_y = -0.05,show.legend = FALSE)+
  geom_point(data = filter(out1,Order.q %in% c(0,1,2)),aes(x = Order.q, y = Estimate, col = Community),show.legend = FALSE, size = 3)+
  theme(text = element_text(size=14), legend.position = "bottom",legend.title = element_blank()) + ggtitle("(a) Sample completeness profiles")
out3_est <- apply(x,2,function(i) MakeTable_Proposeprofile.inc(data = i,B = B,q =q,conf = 0.95)) %>% 
  do.call(rbind,.) %>% mutate(Community = rep(colnames(x),each=length(q)), method = "Estimate") %>% 
  rename(Diversity = Estimate)
out3_mle <- apply(x,2,function(i) MakeTable_EmpericalDiversityprofile.inc(data = i,B = B,q =q,conf = 0.95)) %>% 
  do.call(rbind,.) %>% mutate(Community = rep(colnames(x),each=length(q)), method = "Empirical") %>% 
  rename(Diversity = Emperical)
out3_mle$LCL <- out3_mle$UCL <- out3_mle$Diversity
out3 <- rbind(out3_est,out3_mle)
out3$method[out3$method=="Estimate"] <- 'Asymptotic'
out3$method <- factor(out3$method,levels = unique(out3$method))
out3_p <- plot_diversity_profile_yhc(output = out3) + theme_bw() + ggtitle("(c) Asymptotic and empirical diversity profiles")+
  geom_point(data = filter(out3,Order.q %in% c(0,1,2)),aes(x = Order.q, y = Diversity, col = Community),show.legend = FALSE, size = 3)+
  theme(legend.position = "bottom",legend.key.width = unit(2,"line"),text = element_text(size=14),
        legend.title=element_blank(),
        legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,-10,-5,-10)) + ylab("Species diversity")
# + coord_cartesian(ylim = c(6, 70)) 
out24 <- iNEXT::iNEXT(x = x,q = c(0,1,2),datatype = "incidence_freq",nboot = B)
out24$iNextEst <- lapply(out24$iNextEst,function(out){
  out$method[out$method =='interpolated'] <- 'Interpolated'
  out$method[out$method =='extrapolated'] <- 'Extrapolated'
  out
})
out2_p <- ggiNEXT.iNEXT(x = out24,type = 1,facet.var = "order")+ ggtitle("(b) Size-based rarefaction/extrapolation")+
  theme(text = element_text(size=14),legend.position = "none")+
  scale_x_continuous(breaks =c(0,150,300))
out4_p <- ggiNEXT.iNEXT(x = out24,type = 3,facet.var = 'order')+ggtitle("(d) Coverage-based rarefaction/extrapolation")+
  theme(text = element_text(size=14),legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,-10,-5,-10))+
  #annotate("text", x=0.5, y=100, label=c("C[max] == 98.9*'%' "), parse=TRUE) + 
  scale_x_continuous(breaks = c(0.25,0.5,0.75,1),labels = c("0.25","0.50","0.75","1.00"))
SC <- min(sapply(x, function(x) iNEXT:::Chat.Sam(x, 2*x[1])))
out5 <- estimateDyhc(x = x,q = q,datatype = "incidence_freq",base = "coverage",level = SC,conf = NULL)
out5 <- out5 %>% group_by(Community) %>% mutate(Evenness = (Diversity-1)/(Diversity[q==0]-1))
out5_p <- myggeven(out = out5)

out_p <- ggarrange(out1_p,out2_p,out3_p,out4_p,out5_p,ncol = 2, nrow = 3)
ggsave(filename = "Figure_3.pdf",plot = out_p,device = "pdf",width = 36,height = 25,
       units = "cm",dpi = 320)

#=====Fig_5=====
x <- dts$`Stony coral` %>% .[rowSums(.)>0,]

out1 <- apply(x,2,function(i) sc_profile(freq = i,datatype = "incidence",q = q,B = B,conf = 0.95)) %>% 
  do.call(rbind,.)
out1$Community <- rep(colnames(x),each=length(q))
out1_p <- plot_sc_profile(data = out1) + theme_bw() + ylab("Sample completeness")+
  #geom_text(data = filter(out1,Order.q %in% c(0,1,2)),aes(x = Order.q, y = Estimate, label=round(Estimate,2),col = Community),nudge_x = 0.05,nudge_y = -0.05,show.legend = FALSE)+
  geom_point(data = filter(out1,Order.q %in% c(0,1,2)),aes(x = Order.q, y = Estimate, col = Community),show.legend = FALSE, size = 3)+
  theme(text = element_text(size=14), legend.position = "bottom",legend.title = element_blank()) + ggtitle("(a) Sample completeness profiles")
out3_est <- apply(x,2,function(i) MakeTable_Proposeprofile.inc(data = i,B = B,q =q,conf = 0.95)) %>% 
  do.call(rbind,.) %>% mutate(Community = rep(colnames(x),each=length(q)), method = "Estimate") %>% 
  rename(Diversity = Estimate)
out3_mle <- apply(x,2,function(i) MakeTable_EmpericalDiversityprofile.inc(data = i,B = B,q =q,conf = 0.95)) %>% 
  do.call(rbind,.) %>% mutate(Community = rep(colnames(x),each=length(q)), method = "Empirical") %>% 
  rename(Diversity = Emperical)
out3_mle$LCL <- out3_mle$UCL <- out3_mle$Diversity
out3 <- rbind(out3_est,out3_mle)
out3$method[out3$method=="Estimate"] <- 'Asymptotic'
out3$method <- factor(out3$method,levels = unique(out3$method))
out3_p <- plot_diversity_profile_yhc(output = out3) + theme_bw() + ggtitle("(c) Asymptotic and empirical diversity profiles")+
  geom_point(data = filter(out3,Order.q %in% c(0,1,2)),aes(x = Order.q, y = Diversity, col = Community),show.legend = FALSE, size = 3)+
  theme(legend.position = "bottom",legend.key.width = unit(2,"line"),text = element_text(size=14),
        legend.title=element_blank(),
        legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,-10,-5,-10)) + ylab("Species diversity")
# + coord_cartesian(ylim = c(6, 70)) 
out24 <- iNEXT::iNEXT(x = x,q = c(0,1,2),datatype = "incidence_freq",nboot = B)
out24$iNextEst <- lapply(out24$iNextEst,function(out){
  out$method[out$method =='interpolated'] <- 'Interpolated'
  out$method[out$method =='extrapolated'] <- 'Extrapolated'
  out
})
out2_p <- ggiNEXT.iNEXT(x = out24,type = 1,facet.var = "order")+ ggtitle("(b) Size-based rarefaction/extrapolation")+
  theme(text = element_text(size=14),legend.position = "none")
out4_p <- ggiNEXT.iNEXT(x = out24,type = 3,facet.var = 'order')+ggtitle("(d) Coverage-based rarefaction/extrapolation")+
  theme(text = element_text(size=14),legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,-10,-5,-10))+
  #annotate("text", x=0.5, y=100, label=c("C[max] == 98.9*'%' "), parse=TRUE) + 
  scale_x_continuous(breaks = c(0.25,0.5,0.75,1),labels = c("0.25","0.50","0.75","1.00"))
SC <- min(sapply(x, function(x) iNEXT:::Chat.Sam(x, 2*x[1])))
out5 <- estimateDyhc(x = x,q = q,datatype = "incidence_freq",base = "coverage",level = SC,conf = NULL)
out5 <- out5 %>% group_by(Community) %>% mutate(Evenness = (Diversity-1)/(Diversity[q==0]-1))
out5_p <- myggeven(out = out5)

out_p <- ggarrange(out1_p,out2_p,out3_p,out4_p,out5_p,ncol = 2, nrow = 3)
ggsave(filename = "Figure_5.pdf",plot = out_p,device = "pdf",width = 36,height = 25,
       units = "cm",dpi = 320)