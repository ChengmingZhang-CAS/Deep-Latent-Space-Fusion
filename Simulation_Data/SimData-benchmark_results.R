# simulation data test
load.libraries <- function() {
  library('PMA')
  library('R.matlab')
  library('SNFtool')
  library('PINSPlus')
  library('LRAcluster')
  library('kernlab')
  library('survival')
  library('NMF')
  library('iClusterPlus')
  
  # # bioconductor packages
  # source("https://bioconductor.org/biocLite.R")
  # biocLite("impute")
  # # biocLite("iClusterPlus")
}
load.libraries()


SIMULATION.DATA = list(
  list(name='SimData1', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='SimData1'),
  list(name='SimData2', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='SimData2'),
  list(name='SimData3', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='SimData3'),
  list(name='SimData4', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='SimData4'),
  list(name='SimData5', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='SimData5'),
  list(name='SimData6', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='SimData6'),
  list(name='SimData7', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='SimData7'),
  list(name='SimData8', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='SimData8'),
  list(name='SimData9', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='SimData9'),
  list(name='SimData10', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='SimData10'))



MAX.NUM.CLUSTERS = 15

OMIC.SUBSETS = list('multi_omics', 'exp', 'methy', 'mirna')
names(OMIC.SUBSETS) = c('all', '1', '2', '3')

get.clustering.results.dir.path <- function() {
  if (!dir.exists('RESULTS_DIR_PATH')) {
    dir.create('RESULTS_DIR_PATH')
  }
  return('RESULTS_DIR_PATH')
}

get.plots.dir.path <- function() {
  results.dir.path = get.clustering.results.dir.path()
  return(file.path(results.dir.path, 'plots'))
}

get.tables.dir.path <- function() {
  results.dir.path = get.clustering.results.dir.path()
  return(file.path(results.dir.path, 'tables'))
}


subtype.to.display.name <- function(subtype) {
  for (i in 1:length(SIMULATION.DATA)) {
    if (SIMULATION.DATA[[i]]$name == subtype) {
      return(SIMULATION.DATA[[i]]$display.name)
    }
  }
}

set.omics.list.attr <- function(subtype.raw.data, subtype.data) {
  attr(subtype.raw.data[[1]], 'is.seq') = subtype.data$is.rna.seq
  attr(subtype.raw.data[[2]], 'is.seq') = F
  attr(subtype.raw.data[[3]], 'is.seq') = subtype.data$is.mirna.seq
  return(subtype.raw.data)
}

# ALGORITHM.NAMES = c('kmeans', 'spectral', 'lracluster', 'pins', 'snf', 'mkl', 'mcca', 'nmf', 'iCluster')
# ALGORITHM.DISPLAY.NAMES = as.list(c('K-means', 'Spectral', 'LRAcluster', 'PINS', 'SNF', 'rMKL-LPP', 'MCCA', 'MultiNMF', 'iClusterBayes'))

ALGORITHM.NAMES = c('kmeans', 'spectral', 'lracluster', 'pins', 'snf', 'mcca', 'iCluster')
ALGORITHM.DISPLAY.NAMES = as.list(c('K-means', 'Spectral', 'LRAcluster', 'PINS', 'SNF', 'MCCA', 'iClusterBayes'))

names(ALGORITHM.DISPLAY.NAMES) = ALGORITHM.NAMES
			

run.benchmark <- function() {
  # run all algorithm in all datasets(run and save algorithm.ret)
  for (i in 1:length(SIMULATION.DATA)) {
    current.subtype.data = SIMULATION.DATA[[i]]
    subtype = current.subtype.data$name
    subtype.raw.data = get.raw.data(subtype, 
                                    only.primary=current.subtype.data$only.primary)
    
    subtype.raw.data = set.omics.list.attr(subtype.raw.data, 
                                           current.subtype.data)
    
    for (algorithm.name in ALGORITHM.NAMES) {
      for (j in c('all', '1', '2', '3')) {
        set.seed(42)
        print(paste('subtype', subtype, 'running algorithm', algorithm.name, j))
        clustering.path = file.path(get.clustering.results.dir.path(),
                                    paste(subtype, algorithm.name, j, sep='_'))
        timing.path = file.path(get.clustering.results.dir.path(),
                                paste(subtype, algorithm.name, j, 'timing', sep='_'))
        
        
        if (!file.exists(clustering.path)) {
          algorithm.func.name = paste0('run.', algorithm.name)
          algorithm.func = get(algorithm.func.name)
          if (j == 'all') {
            cur.iteration.data = subtype.raw.data
          } else {
            cur.iteration.data = subtype.raw.data[as.numeric(j)]
          }
          algorithm.ret = algorithm.func(cur.iteration.data, current.subtype.data)
          clustering = algorithm.ret$clustering
          timing = algorithm.ret$timing
          print('before saving')
          save(clustering, file = clustering.path)
          save(timing, file = timing.path)
        }
      }
    }
  }
}


analyze.benchmark <- function() {
  # get all the algorithm clusterings results(load algorithm.ret and get benchmark.ret)
  all.clusterings = list()
  all.timings = list()
  for (i in 1:length(SIMULATION.DATA)) {
    current.subtype.data = SIMULATION.DATA[[i]]
    subtype = current.subtype.data$name
    subtype.raw.data = get.raw.data(subtype, 
                                    only.primary=current.subtype.data$only.primary)
    
    all.clusterings[[subtype]] = list()
    all.timings[[subtype]] = list()
    
    for (algorithm.name in ALGORITHM.NAMES) {
      all.clusterings[[subtype]][[algorithm.name]] = list()
      all.timings[[subtype]][[algorithm.name]] = list()
      for (j in c('all', '1', '2', '3')) {
        clustering.path = file.path(get.clustering.results.dir.path(),
                                    paste(subtype, algorithm.name, j, sep='_'))
        timing.path = file.path(get.clustering.results.dir.path(),
                                paste(subtype, algorithm.name, j, 'timing', sep='_'))
        load(clustering.path)
        load(timing.path)
        if (!any(is.na(clustering))) {
          names(clustering) = colnames(subtype.raw.data[[1]])
        }
        
        all.clusterings[[subtype]][[algorithm.name]][[j]] = clustering
        all.timings[[subtype]][[algorithm.name]][[j]] = timing
      }
    }
  }
  return(list(all.clusterings=all.clusterings, all.timings=all.timings))
}


write.clusterings.to.csv <- function(benchmark.ret) {
  all.clusterings = benchmark.ret$all.clusterings
  all.timings = benchmark.ret$all.timings
  folder.path = file.path(get.clustering.results.dir.path(), 'clusterings')
  if (!dir.exists(folder.path)) {
    dir.create(folder.path)
  }
  
  # all.clusterings$SimData$algorithm$all
  for (i in 1:length(SIMULATION.DATA)) {
    subtype = SIMULATION.DATA[[i]]$name
    print(subtype)
    for (omic in c('all', '1', '2', '3')) {
      cur.clusterings = list()
      for (algorithm.name in ALGORITHM.NAMES) {
        cur.clusterings[[algorithm.name]] = all.clusterings[[subtype]][[algorithm.name]][[omic]]
      }
      # colnames(cur.clusterings) = ALGORITHM.NAMES
      cur.clusterings.path = file.path(folder.path, 
                                       paste0('clusterings', '_', subtype, '_', omic, '.csv'))
      cur.clusterings = as.data.frame(cur.clusterings)
      rownames(cur.clusterings) = matrix(1:length(cur.clusterings[[1]]))
      write.csv(cur.clusterings, file = cur.clusterings.path)
    }
  }
}


benchmark.omics.num.clusters <- function(benchmark.ret) {
  num.clusters = matrix(1, nrow=length(OMIC.SUBSETS), ncol=length(ALGORITHM.NAMES))
  rownames(num.clusters) = sapply(OMIC.SUBSETS, function(x) x)
  colnames(num.clusters) = ALGORITHM.NAMES
  all.clusterings = benchmark.ret$all.clusterings
  folder.path = file.path(get.clustering.results.dir.path(), 'clusterings')
  if (!dir.exists(folder.path)) {
    dir.create(folder.path)
  }
  omics = c('all', '1', '2', '3')
  for (k in 1:length(SIMULATION.DATA)) {
    subtype = SIMULATION.DATA[[k]]$name
    print(subtype)
    for (i in 1:length(OMIC.SUBSETS)) {
      for (j in 1:length(ALGORITHM.NAMES)) {
        clustering = all.clusterings[[subtype]][[ALGORITHM.NAMES[j]]][[omics[i]]]
        num.clusters[i, j] = max(clustering)
      }
    }
    cur.num.clusters.path = file.path(folder.path, 
                                     paste0('clusterings', '_', subtype, '_', 'num_cluster.csv'))
    num.clusters = as.data.frame(num.clusters)
    write.csv(num.clusters, file = cur.num.clusters.path)
  }
}


log.and.normalize <- function(omics.data, subtype.data, normalize=T,
                              filter.var=F) {
  # filter features with no variance at all
  for (i in 1:length(omics.data)) {
    omics.data[[i]] = omics.data[[i]][apply(omics.data[[i]], 1, var) > 0,]
  }
  
  for (i in 1:length(omics.data)) {
    if (attr(omics.data[[i]], 'is.seq')) {
      omics.data[[i]] = log(1+omics.data[[i]])
    }
  }
  
  if (filter.var) {
    omics.data = lapply(omics.data, keep.high.var.features)
  }
  
  if (normalize) {
    omics.data = lapply(omics.data, normalize.matrix)    
  }
  
  return(omics.data)
}

normalize.matrix <- function(data.matrix) {
  temp = data.matrix - rowMeans(data.matrix)
  should.keep = (apply(temp, 1, sd) != 0)
  return ((temp / apply(temp, 1, sd))[should.keep, ])
}


get.raw.data <- function(subtype, only.primary=NA) {
  # omics.dir = file.path(datasets.path, subtype)
  omics.files = c("gene_expression.csv", "DNA_methylation.csv", "mirna_expression.csv")
  omics.dir = paste0(subtype, "_", omics.files)
  raw.data = list()
  for (i in 1:length(omics.files)){
    print(paste("omics: ", i))
    raw.data[[i]] = read.csv(file.path(subtype, omics.dir[i]), header = TRUE, row.names = 1)
  }
  
  return(raw.data)
}

get.elbow <- function(values, is.max) {
  second.derivatives = c()
  for (i in 2:(length(values) - 1)) {
    second.derivative = values[i + 1] + values[i - 1] - 2 * values[i]
    second.derivatives = c(second.derivatives, second.derivative)
  }
  print(second.derivatives)
  if (is.max) {
    return(which.max(second.derivatives) + 1)
  } else {
    return(which.min(second.derivatives) + 1)
  }
}

# Does not support a single omic dataset
run.mcca <- function(omics.list, subtype.data) {
  if (length(omics.list) == 1) {
    return(list(clustering=rep(NA, ncol(omics.list[[1]])), timing=1))
  }
  start = Sys.time()
  omics.list = log.and.normalize(omics.list, subtype.data, 
                                 normalize = T,
                                 filter.var = T)
  
  subtype = subtype.data$name
  omics.transposed = lapply(omics.list, t)
  cca.ret = PMA::MultiCCA(omics.transposed, 
                          ncomponents = MAX.NUM.CLUSTERS)
  sample.rep = omics.transposed[[1]] %*% cca.ret$ws[[1]]
  
  explained.vars = sapply(1:MAX.NUM.CLUSTERS, 
                          function(i) sum(unlist(apply(sample.rep[1:i,,drop=F], 2, var))))
  
  dimension = get.elbow(explained.vars, is.max=F)
  print(dimension)
  sample.rep = sample.rep[,1:dimension]
  sils = c()
  clustering.per.num.clusters = list()
  for (num.clusters in 2:MAX.NUM.CLUSTERS) {
    cur.clustering = kmeans(sample.rep, num.clusters, iter.max=100, nstart=30)$cluster  
    sil = get.clustering.silhouette(list(t(sample.rep)), cur.clustering)
    sils = c(sils, sil)
    clustering.per.num.clusters[[num.clusters - 1]] = cur.clustering
  }
  # NOTE: the next line contains an error. We mistakenly selected the minimal rather maximal silhouette.
  # See more details in: http://acgt.cs.tau.ac.il/multi_omic_benchmark/download.html.
  # Now we corrected the error and selcted the maximal silhouette again!
  cca.clustering = clustering.per.num.clusters[[which.max(sils)]]
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=cca.clustering, timing=time.taken))
}

run.snf <- function(omics.list, subtype.data) {
  start = Sys.time()
  omics.list = log.and.normalize(omics.list, subtype.data)
  subtype = subtype.data$name
  alpha=0.5
  T.val=30
  num.neighbors = round(ncol(omics.list[[1]]) / 10)
  similarity.data = lapply(omics.list, function(x) {SNFtool::affinityMatrix(SNFtool::dist2(as.matrix(t(x)),as.matrix(t(x))), 
                                                                            num.neighbors, alpha)})
  if (length(similarity.data) == 1) {
    W = similarity.data[[1]]
  } else {
    W = SNFtool::SNF(similarity.data, num.neighbors, T.val)  
  }
  
  num.clusters = SNFtool::estimateNumberOfClustersGivenGraph(W, 2:MAX.NUM.CLUSTERS)[[3]]  
  clustering = SNFtool::spectralClustering(W, num.clusters)
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=clustering, timing=time.taken))
}

run.iCluster <- function(omics.list, subtype.data) {
  omics.list = log.and.normalize(omics.list, subtype.data, normalize = F)
  
  start = Sys.time()
  subtype = subtype.data$name
  dev.ratios = c()
  icluster.rets = list()
  
  if (length(omics.list) == 1) {
    # icluster.ret = iClusterPlus::tune.iClusterBayes(cpus=(MAX.NUM.CLUSTERS - 1), t(omics.list[[1]]), 
                                                    # K=1:(MAX.NUM.CLUSTERS - 1), type=c('gaussian'))$fit
    icluster.ret = iClusterPlus::tune.iClusterBayes(cpus=1, t(omics.list[[1]]), 
                                                    K=1:(MAX.NUM.CLUSTERS - 1), type=c('gaussian'))$fit
  } else {
    icluster.ret = iClusterPlus::tune.iClusterBayes(cpus=1, t(omics.list[[1]]), 
                                                    t(omics.list[[2]]), 
                                                    t(omics.list[[3]]), 
                                                    K=1:(MAX.NUM.CLUSTERS - 1), type=rep('gaussian', 3))$fit
  }
  dev.ratios = lapply(1:(MAX.NUM.CLUSTERS - 1), function(i) icluster.ret[[i]]$dev.ratio)
  
  print('dev.ratios are:')
  print(dev.ratios)
  
  optimal.solution = icluster.ret[[which.max(dev.ratios)]]
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=optimal.solution$clusters, 
              timing=time.taken))
}

get.mkl.binary.path = function() {
  return('MKL_BINARY_PATH')
}

get.mkl.arguments.path = function() {
  return('MKL_ARGS_PATH')
}


run.mkl <- function(omics.list, subtype.data) {
  start = Sys.time()
  omics.list = log.and.normalize(omics.list, subtype.data)
  subtype = subtype.data$name
  omics.list = lapply(omics.list, normalize.matrix)
  time.taken = as.numeric(Sys.time() - start, units='secs')
  export.subtype.to.mkl(omics.list, subtype)
  
  start = Sys.time()
  bin.path = get.mkl.binary.path()
  subtype.dir = paste0(get.mkl.arguments.path(), subtype, '\\')
  paste0(subtype.dir, 'kernels')
  command = paste(bin.path, paste0(subtype.dir, 'kernels'),
                  paste0(subtype.dir, 'ids'),
                  paste0(subtype.dir, 'output'), 
                  '9', '5')
  command.return = system(command)
  stopifnot(command.return == 0)
  time.taken2 = as.numeric(Sys.time() - start, units='secs')
  clustering = get.mkl.clustering(subtype)
  return(list(clustering=clustering, 
              timing=time.taken + time.taken2))
}

run.nmf <- function(omics.list, subtype.data) {
  total.time.taken = 0
  start = Sys.time()
  omics.list = log.and.normalize(omics.list, subtype.data,
                                 filter.var = T, normalize = F)
  time.taken = as.numeric(Sys.time() - start, units='secs')
  total.time.taken = total.time.taken + time.taken
  
  save.subtype.matlab.format(omics.list)
  subtype = subtype.data$name
  if (length(omics.list) > 1) {
    command.ret = system('MULTI_NMF_COMMAND')
    stopifnot(command.ret == 0)
    nmf.timing = read.csv('SOME_TEMP_PATH_TIMING', header=F)[1, 1]
    total.time.taken = total.time.taken + nmf.timing
  } else {
    for (k in 1:MAX.NUM.CLUSTERS) {
      start = Sys.time()
      file.name = paste0('SOME_TEMP_PATH/', k, '_consensus')
      nmf.ret = nmf(omics.list[[1]], k, method='lee')
      coef.mat = t(coef(nmf.ret))
      time.taken = as.numeric(Sys.time() - start, units='secs')
      total.time.taken = total.time.taken + time.taken
      write.table(coef.mat, file=file.name, quote=F, row.names=F, col.names=F, sep=',')
    }
  }
  
  explained.vars = c()
  clustering.per.num.clusters = list()
  for (k in 1:MAX.NUM.CLUSTERS) {
    file.name = paste0('SOME_TEMP_PATH/', k, '_consensus')
    consensus.mat = read.csv(file.name, header=F)
    
    start = Sys.time()
    cur.clustering = apply(consensus.mat, 1, which.max)
    explained.var = sum(unlist(apply(consensus.mat, 2, var)))
    explained.vars = c(explained.vars, explained.var)
    clustering.per.num.clusters[[k]] = cur.clustering
    time.taken = as.numeric(Sys.time() - start, units='secs')
    total.time.taken = total.time.taken + time.taken
  }
  
  dimension = get.elbow(explained.vars, is.max=F)
  nmf.clustering = clustering.per.num.clusters[[dimension]]
  return(list(clustering=nmf.clustering, timing=total.time.taken))  
}

run.pins <- function(omics.list, subtype.data) {
  start = Sys.time()
  omics.list = log.and.normalize(omics.list, subtype.data, normalize = F)
  subtype = subtype.data$name
  omics.transposed = lapply(omics.list, t)
  if (length(omics.list) == 1) {
    pins.ret = PINSPlus::PerturbationClustering(data=omics.transposed[[1]],
                                                kMax = MAX.NUM.CLUSTERS)
    clustering = pins.ret$cluster
    
  } else {
    pins.ret = PINSPlus::SubtypingOmicsData(dataList=omics.transposed,
                                            kMax = MAX.NUM.CLUSTERS)
    clustering = pins.ret$cluster2
  }
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=clustering, timing=time.taken))
}

run.lracluster <- function(omics.list, subtype.data) {
  omics.list = log.and.normalize(omics.list, subtype.data, normalize = F)
  
  subtype = subtype.data$name
  start = Sys.time()
  
  dim.range = 1:MAX.NUM.CLUSTERS
  all.clustering.results = list()
  
  omics.matrix.list = lapply(omics.list, as.matrix)
  for (dimension in dim.range) {
    print(paste('running lra cluster for dimension', dimension))
    data.names = c('gene expression', 'methylation', 'miRNA expression')
    clustering.results = LRAcluster::LRAcluster(omics.matrix.list, 
                                                rep('gaussian', length(omics.list)), 
                                                dimension=dimension, data.names)
    all.clustering.results[[dimension]] = clustering.results
  }
  explained.var = sapply(all.clustering.results, function(x) x$potential)
  print(explained.var)
  dimension = get.elbow(explained.var, is.max=F)
  print(dimension)
  solution = all.clustering.results[[dimension]]$coordinate
  
  sils = c()
  clustering.per.num.clusters = list()
  for (num.clusters in 2:MAX.NUM.CLUSTERS) {
    print(paste('running kmeans in lra cluster for num clusters', num.clusters))
    cur.clustering = kmeans(t(solution), num.clusters, iter.max=100, nstart=60)$cluster
    sil = get.clustering.silhouette(list(solution), cur.clustering)
    sils = c(sils, sil)
    clustering.per.num.clusters[[num.clusters - 1]] = cur.clustering
  }
  print(sils)
  # NOTE: the next line contains an error. We mistakenly selected the minimal rather maximal silhouette.
  # See more details in: http://acgt.cs.tau.ac.il/multi_omic_benchmark/download.html.
  # Now we corrected the error and selcted the maximal silhouette again!
  chosen.clustering = clustering.per.num.clusters[[which.max(sils)]]
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=chosen.clustering, timing=time.taken))
}

run.kmeans <- function(omics.list, subtype.data) {
  start = Sys.time()
  omics.list = log.and.normalize(omics.list, subtype.data, 
                                 filter.var = T)
  
  subtype = subtype.data$name
  all.withinss = c()
  all.clusterings = list()
  k.range = 1:MAX.NUM.CLUSTERS
  for (k in k.range) {
    concat.omics = do.call(rbind, omics.list)
    kmeans.ret = kmeans(t(concat.omics), k, iter.max=100, nstart=60)
    all.withinss = c(all.withinss, kmeans.ret$tot.withinss)
    all.clusterings[[k]] = kmeans.ret$cluster
  }
  
  best.k = get.elbow(all.withinss, is.max=T)
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=all.clusterings[[best.k]], 
              timing=time.taken))
}

run.spectral <- function(omics.list, subtype.data) {
  start = Sys.time()
  omics.list = log.and.normalize(omics.list, subtype.data, 
                                 filter.var = T)
  subtype = subtype.data$name
  concat.omics = do.call(rbind, omics.list)
  
  similarity.data = SNFtool::affinityMatrix(SNFtool::dist2(as.matrix(t(concat.omics)),
                                                           as.matrix(t(concat.omics))), 20, 0.5)
  
  num.clusters = SNFtool::estimateNumberOfClustersGivenGraph(similarity.data, 
                                                             2:MAX.NUM.CLUSTERS)[[3]]  
  clustering = SNFtool::spectralClustering(similarity.data, num.clusters)
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=clustering, timing=time.taken))
}

load.libraries <- function() {
  library('PMA')
  library('R.matlab')
  library('SNFtool')
  library('PINSPlus')
  library('LRAcluster')
  library('kernlab')
  library('survival')
  library('NMF')
  library('iClusterPlus')
  
  # # bioconductor packages
  # source("https://bioconductor.org/biocLite.R")
  # biocLite("impute")
  # # biocLite("iClusterPlus")
}



########################################
###############   MKL    ###############
########################################

radial.basis <- function(mat, gamma) {
  if (missing(gamma)) {
    gamma = 1 / (2*nrow(mat)**2)
  }
  npatients = ncol(mat)
  output.mat = matrix(0, ncol=npatients, nrow=npatients)
  for (i in 1:npatients) {
    for (j in 1:npatients) {
      output.mat[i, j] = exp(-norm(as.matrix(mat[,i] - mat[,j]), type = 'F')**2 * gamma)
    }
  }
  
  D = apply(output.mat, 2, sum) / npatients
  E = sum(D) / npatients
  J = matrix(1, nrow=npatients, ncol=1) %*% D
  ret = output.mat - J - t(J) + E * matrix(1, ncol=npatients, nrow=npatients)
  ret = diag(1/sqrt(diag(ret))) %*% ret %*% diag(1/sqrt(diag(ret)))
  return(ret)
}

clear.dir <- function(dir.path) {
  files.in.dir = list.files(dir.path)
  for (file.in.dir in files.in.dir) {
    full.file.path = file.path(dir.path, file.in.dir)
    file.remove(full.file.path)
  }
}

export.subtype.to.mkl <- function(omics.list, dir.name) {
  
  folder.path = file.path(get.mkl.arguments.path(), dir.name)
  if (!dir.exists(folder.path)) {
    dir.create(folder.path)
  }
  
  kernels.path = file.path(folder.path, 'kernels')
  
  if (!dir.exists(kernels.path)) {
    dir.create(kernels.path)
  }
  clear.dir(kernels.path)
  
  gammas = 10 ** seq(-6, 6, by=3)
  for (i in 1:length(omics.list)) {
    for (j in 1:length(gammas)) {
      datum = omics.list[[i]]
      gamma = gammas[[j]] / (2*nrow(datum)**2)
      mat = radial.basis(datum, gamma)
      R.matlab::writeMat(file.path(kernels.path, paste(i, '_', j, sep='')), mat=mat)
    }
  }
  
  output.path = file.path(folder.path, 'output')
  if (!dir.exists(output.path)) {
    dir.create(output.path)
  }
  clear.dir(output.path)
  
  write.table(colnames(omics.list[[1]]), file=file.path(folder.path, 'ids'),
              quote=F, row.names = F, col.names = F)
  
}

get.mkl.clustering <- function(dir.name) {
  folder.path = file.path(get.mkl.arguments.path(), dir.name)
  output.path = file.path(folder.path, 'output')
  output.files = list.files(output.path)
  clustering = read.csv(file.path(output.path, output.files[grep('clusters', output.files)]))[,2]
  return(clustering)
}

check.survival <- function(groups, subtype, survival.file.path) {
  # package.surv(input: clustering, subtype; output: )
  if (missing(survival.file.path)) {
    survival.file.path = get.subtype.survival.path(subtype)
  }
  survival.data = read.table(survival.file.path, header = TRUE)
  patient.names = names(groups)
  patient.names.in.file = as.character(survival.data[, 1])
  patient.names.in.file = toupper(gsub('-', '\\.', substring(patient.names.in.file, 1, 12)))
  
  stopifnot(all(patient.names %in% patient.names.in.file))
  
  indices = match(patient.names, patient.names.in.file)
  ordered.survival.data = survival.data[indices,]
  ordered.survival.data["cluster"] <- groups
  rownames(ordered.survival.data) = patient.names
  # ordered.survival.data$Survival[is.na(ordered.survival.data$Survival)] = 0
  # ordered.survival.data$Death[is.na(ordered.survival.data$Death)] = 0
  ordered.survival.data = ordered.survival.data[!is.na(ordered.survival.data$Survival), ]
  return(survdiff(Surv(Survival, Death) ~ cluster, data=ordered.survival.data))
  
}


get.clustering.silhouette <- function(raw.data, clustering) {
  sils = c()
  for (i in 1:length(raw.data)) {
    x = raw.data[[i]]
    distmatrix = SNFtool::dist2(as.matrix(t(x)),as.matrix(t(x)))
    sil = cluster::silhouette(clustering, dmatrix = distmatrix)[,3]
    sils = c(sils, mean(sil))
  }
  return(mean(sils))
}


keep.high.var.features <- function(omic, num.features=2000) {
  if (nrow(omic) < num.features) {
    return(omic)
  } else {
    feature.vars = apply(omic, 1, var)
    threshold = feature.vars[order(feature.vars, decreasing = T)][num.features]
    return(omic[feature.vars >= threshold,])    
  }
}

run.benchmark()
benchmark.ret = analyze.benchmark()
write.clusterings.to.csv(benchmark.ret)
benchmark.omics.num.clusters(benchmark.ret)


