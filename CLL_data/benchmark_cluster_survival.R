#-------------------------------------------------生存分析-----------------------------------------------------------
library("ggplot2")
library("survival")
library("survminer")


# get the chisq value
check.survival <- function(surv.data) {
  return(survdiff(Surv(T5, treatedAfter ) ~ cluster, data=surv.data))
}

get.logrank.pvalue <- function(survdiff.res) {
  1 - pchisq(survdiff.res$chisq, length(survdiff.res$n) - 1)  
}

get.empirical.surv <- function(surv.data) {
  # surv fucn2(input: clustering, subtype; output: empirical.surv)
  set.seed(521)
  surv.ret = check.survival(surv.data)
  orig.chisq = surv.ret$chisq
  orig.pvalue = get.logrank.pvalue(surv.ret)
  
  # The initial number of permutations to run
  num.perms = round(min(max(10 / orig.pvalue, 1000), 1e6))
  should.continue = T
  
  total.num.perms = 0
  total.num.extreme.chisq = 0
  cur.surv.data = surv.data
  
  while (should.continue) {
    print('Another iteration in empirical survival calculation')
    print(num.perms)
    # could not find function "mclapply", replace mclapplay by lapply
    perm.chisq = as.numeric(lapply(1:num.perms, function(i) {
      cur.clustering = sample(surv.data$cluster)
      cur.surv.data$cluster = cur.clustering
      cur.chisq = check.survival(cur.surv.data)$chisq
      return(cur.chisq)
    }))
    
    total.num.perms = total.num.perms + num.perms
    total.num.extreme.chisq = total.num.extreme.chisq + sum(perm.chisq >= orig.chisq)
    
    binom.ret = binom.test(total.num.extreme.chisq, total.num.perms)
    cur.pvalue = binom.ret$estimate
    cur.conf.int = binom.ret$conf.int
    
    print(c(total.num.extreme.chisq, total.num.perms))
    print(cur.pvalue)
    print(cur.conf.int)
    
    sig.threshold = 0.05
    is.conf.small = ((cur.conf.int[2] - cur.pvalue) < min(cur.pvalue / 10, 0.01)) & ((cur.pvalue - cur.conf.int[1]) < min(cur.pvalue / 10, 0.01))
    is.threshold.in.conf = cur.conf.int[1] < sig.threshold & cur.conf.int[2] > sig.threshold
    if ((is.conf.small & !is.threshold.in.conf) | (total.num.perms > 1e6)) {
      #if (is.conf.small) {
      should.continue = F
    } else {
      print(paste0("is.conf.small: ", is.conf.small))
      print(paste0("!is.threshold.in.conf: ", !is.threshold.in.conf))
      num.perms = 1e5
    }
  }
  
  print(paste0("orig.pvalue: ", orig.pvalue))
  print(paste0("empirical.pvalue: ", cur.pvalue))
  
  return(list(pvalue = cur.pvalue, conf.int = cur.conf.int, total.num.perms=total.num.perms, 
              total.num.extreme.chisq=total.num.extreme.chisq))
}



# fused cluster surviral test & clinical parameters enrichment test
# first cluster choise
all.surv.pvalues = list()
ALGORITHM.NAMES = c('kmeans', 'spectral', 'lracluster', 'pins', 'snf', 'mcca', 'iCluster')

print("--------------------final-----------------")
all.clustering = read.csv("RESULTS_DIR_PATH//clusterings//clusterings_CLL_all.csv", header = T, row.names = 1)
patmeta = read.csv("CLL_patmeta.csv", header = T, row.names = 1)
idx = intersect(rownames(all.clustering), rownames(patmeta))

for (method in ALGORITHM.NAMES)  {
  print(method)
  surv.data = patmeta[idx, c('T5', 'treatedAfter')]
  surv.data$cluster = all.clustering[idx, method]
  surv.data = surv.data[complete.cases(surv.data), ]
  
  survie <- Surv(as.numeric(surv.data$T5),as.numeric(surv.data$treatedAfter))
  fit <- survfit(survie~as.factor(surv.data$cluster))
  
  surv.ret = check.survival(surv.data = surv.data)
  raw.pvalue = get.logrank.pvalue(surv.ret)
  
  print(surv.ret)
  
  
  surv.empirical.ret = get.empirical.surv(surv.data)
  all.surv.pvalues[[method]] = c(get.logrank.pvalue(surv.ret), surv.empirical.ret$pvalue)

}
df = as.data.frame(all.surv.pvalues)
rownames(df) = c('origin', 'empirical')
file_name = paste0("all.surv.pvalues", ".csv")
# df = cbind(df, -log10(df))
write.csv(df, file = file_name)
print(df)
















# 
# png('DLSF-survival-first.png', height=600,width=800)
# ggsurvplot(fit, data = surv.data,
#            conf.int = FALSE,
#            pval = TRUE, # 添加P值
#            legend = c(0.8,0.75), # 指定图例位置
#            legend.title = "", # 设置图例标题
#            title = 'DLSF-survival',
#            ggtheme = theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
# )
# dev.off()
# 





# ggsurvplot(fit,data = surv.data,
#            title="DLSF",
#            legend.title = "", # legend.labs = c("cluster1","cluster2","cluster3"),
#            linetype = "solid",size =0.8,palette = "jco",xlab = "Time" , pval=T,# pval = paste0("p=",surv.empirical.ret$pvalue),
#            ggtheme = theme(plot.title = element_text(hjust = 0.5),axis.ticks.length = unit(0.2, 'cm'),axis.line = element_line(colour = 'black'),
#                            axis.text = element_text(size = 12),axis.title = element_text(size = 16),
#                            panel.grid = element_blank(),panel.border = element_blank(),panel.background = element_blank())
# )

