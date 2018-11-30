

Dat<-read.csv('Great_tit_AdultsOnly.csv')


str(Dat)

Tab<-table(Dat$code, Dat$winteryear)

str(Tab)
head(Tab)

Tab<-apply(Tab, c(1,2), sign)

# remove last year first marking

remove.last<-function(CH) {
  Index<-apply(CH, 1, FUN=function(x) ifelse(which(x==1)[1] == length(x), FALSE, TRUE))
  return(CH[Index,])
}

Tab<-remove.last(Tab)

