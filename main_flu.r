library(here)
source(here("code","setup01.R"))  # 
clc()

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### where is the data ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### MEM package ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Package ‘mem’
# Moving Epidemic Method: <doi:10.1111/j.1750-2659.2012.00422.x>, https://onlinelibrary.wiley.com/doi/10.1111/irv.12330
#install.packages("mem")
library("mem")
help("mem")

# Castilla y Leon Influenza Rates data
data(flucyl)
# simulate my own data
# understand the outcome
# get data from ECDC for 1 single setting

# Data of the last season
cur <- flucyl[8]
# The model
epi <- memmodel(flucyl[1:7])
# Epidemic thresholds
e.thr <- epi$epidemic.thresholds # based on past curves
# Intensity threhsolds
i.thr <- epi$intensity.thresholds # 

# The graph, default values
# uncomment to execute
m1 <- memsurveillance(cur, e.thr, i.thr,
   i.graph.file = F,
   i.graph.file.name = "graph 1"
)
# No intensity levels
m2 <- memsurveillance(cur, e.thr, i.thr,
   i.graph.file = F,
   i.graph.file.name = "graph 2", i.no.intensity = TRUE
)
# No start/end tickmarks
m3 <- memsurveillance(cur, e.thr, i.thr,
   i.graph.file = F,
   i.graph.file.name = "graph 3", i.start.end.marks = FALSE
)
# Post-epidemic threshold
m4 <- memsurveillance(cur, e.thr, i.thr,
   i.graph.file = F,
   i.graph.file.name = "graph 4", i.pos.epidemic = TRUE
)
# Report for week 2, instead of all data
m5 <- memsurveillance(cur, e.thr, i.thr,
   i.graph.file = F,
   i.graph.file.name = "graph 5", i.week.report = 2
)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### New section ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
