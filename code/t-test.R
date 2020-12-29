library(tidyverse)
library(here)
library(broom)

#load slope data
d <- read_csv(here("data/HeatDecaySlopes.csv"))

# calculate difference in resistance between infected and non infected spores

d.goe2 <- 
  d%>%
    filter(phage=="Goe2"|phage=="noPhage")%>%
  #resistance  = 1+slope
  mutate(resistance=1+slope)%>%
  # make each spore type a column, for diff
  pivot_wider(id_cols=c("host", "colony"), 
              names_from=spore,
              values_from=resistance)%>%
  #calc. diff for CFU spores
  mutate(diff.cfu=`infected host spore`-`non-infected host spore`)%>%
  #calc. diff for PFU spores
  mutate(diff.pfu=`viral spore`-`non-infected host spore`)%>%
  #return to long format with diffs
  select(-`non-infected host spore`,-`infected host spore`,-`viral spore`)%>%
  pivot_longer(cols=c("diff.cfu", "diff.pfu"),names_to = "spore.type", values_to="slope.diff")



# one sample t.test
  # test if infection with ssp phage caused a significant shift in spore resistance
  # that would be manifested as the calculated difference deviating from 0
  # we hypothesized that infection by ssp-carrying phages will 
    # restore heat resistance of dSspAB => diff > 0
      # OR
    # compromise resistance of the WT host => diff < 0 
  # so using one0sided tests with appropriate alternatives.
d.ttest <- tibble()

for (h in unique(d.goe2$host)){
  for(sp.type in unique(d.goe2$spore.type)){
    t.alt <- if_else(h=="wt", "l","g")
    
    t.res <- d.goe2%>%
      filter(host==h)%>%
      filter(spore.type==sp.type)%>%
      pull(slope.diff)%>%
      t.test(., alternative = t.alt)
   
    d.ttest <- 
      tibble(host=h,phage="Goe2",spore.type=sp.type,glance(t.res))%>%
          bind_rows(d.ttest)

  }
}


##############
# same thing for SPO1
d.spo1 <- 
  d%>%
  filter(phage=="SPO1"|phage=="noPhage")%>%
  #resistance  = 1+slope
  mutate(resistance=1+slope)%>%
  # make each spore type a column, for diff
  pivot_wider(id_cols=c("host", "colony"), 
              names_from=spore,
              values_from=resistance)%>%
  #calc. diff for CFU spores
  mutate(diff.cfu=`infected host spore`-`non-infected host spore`)%>%
  #calc. diff for PFU spores
  mutate(diff.pfu=`viral spore`-`non-infected host spore`)%>%
  #return to long format with diffs
  select(-`non-infected host spore`,-`infected host spore`,-`viral spore`)%>%
  pivot_longer(cols=c("diff.cfu", "diff.pfu"),names_to = "spore.type", values_to="slope.diff")



# one sample t.test
# test if infection with ssp phage caused a significant shift in spore resistance
# that would be manifested as the calculated difference deviating from 0
for (h in unique(d.spo1$host)){
  for(sp.type in unique(d.spo1$spore.type)){
    t.alt <- if_else(h=="wt", "l","g")
    
    t.res <- d.spo1%>%
      filter(host==h)%>%
      filter(spore.type==sp.type)%>%
      pull(slope.diff)%>%
      t.test(., alternative = t.alt)
    
    d.ttest <- 
      tibble(host=h,phage="SPO1",spore.type=sp.type,glance(t.res))%>%
      bind_rows(d.ttest)
    
  }
}

d.ttest$p.adj <- p.adjust(d.ttest$p.value, method = "BH")

write_csv(d.ttest,here("data/diff-ttest.csv"))

d.ttest%>%
  mutate(lab=paste("P=",p.adj%>%signif(3)))%>%
  ggplot(aes(x=spore.type, y=estimate))+
  geom_hline(yintercept = 0)+
  geom_pointrange(aes(ymin=conf.low, ymax=conf.high))+
  geom_label(aes(label=lab), y=-0.03)+
  ylim(-0.03,0.03)+
  facet_grid(phage~host)+
  theme_bw()+
  ggsave(here("figures/diff-ttest.png"), width = 7, height = 7)
