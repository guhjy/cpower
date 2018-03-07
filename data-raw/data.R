anchor <-read.csv('data-raw/anchor.csv')
anchor$anchor<-factor(anchor$anchor)
devtools::use_data(anchor,overwrite = T)
