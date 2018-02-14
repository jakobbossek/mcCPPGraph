library(BBmisc)
library(ggplot2)
library(tidyverse)

rres = read_delim("r.csv", delim = " ")
cres = read_delim("cpp.csv", delim = " ")

rres$Language = "R"
cres$Language = "C++"

res = rbind(rres, cres)

# remove n=1000 instances (no chance for R version)
#res = res[!grepl("E499500", res$prob), , drop = FALSE]

res = reshape2::melt(res, id.vars = c("prob", "Language"), value.name = "Value", variable.name = "Time")

pl = ggplot(res, aes(x = prob, y = Value, fill = Language))
pl = pl + geom_bar(stat = "identity", position = "dodge")
#pl = pl + geom_point(aes(colour = prob, shape = Language))
pl = pl + labs(
  title = "Runtime comparisson",
  x = "Problem instance"
)
pl = pl + facet_wrap(~ Time, nrow = 1L, scale = "free_y")
pl = pl + theme(legend.position = "right", axis.text.x = element_text(hjust = 1, angle = 45))
pl = pl + viridis::scale_fill_viridis(discrete = TRUE, end = 0.6)
print(pl)
