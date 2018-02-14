library(mcMST)
library(BBmisc)

# setup
nruns = 100L
graph.files = list.files("../tests/instances", pattern = ".graph$", full.names = TRUE)
graph.files = graph.files[-4]

runtimes = data.frame()

for (graph.file in graph.files) {
  catf("%s", basename(graph.file))

  # import of file
  st = proc.time()
  g = grapherator::readGP(graph.file)
  time.import = (proc.time() - st)[3L]

  # now run mcMSTPrim x times
  st = proc.time()
  res = mcMST::mcMSTPrim(g, lambdas = rep(0.5, nruns))
  runtime.total = (proc.time() - st)[3L]

  runtime.mean = runtime.total / nruns

  runtimes = rbind(runtimes,
    data.frame(
      prob = basename(graph.file),
      time.import = time.import,
      runtime.mean   = runtime.mean,
      runtime.total = runtime.total,
      time.total  = time.import + runtime.total
    ))
}

write.table(runtimes, file = "r.csv", quote = FALSE, row.names = FALSE, col.names = TRUE)


