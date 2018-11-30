

Dat<-read.csv('File_name_csv')


str(Dat)

Tab<-table(Dat$code, Dat$winteryear)

str(Tab)
head(Tab)

Tab<-apply(Tab, c(1,2), sign)


