require(lubridate)
require(openxlsx)
require(plyr)
require(googlesheets4)

members<-read.xlsx("Qbic_members.xlsx",colNames = F)
head(members)
#members$X1<-members[!members$X1 %in% "Luis K.",]
wed<-ymd("2020-08-05")

### manual, the +1 is to add an extra week in case there is a holiday. if a list has 3 holiday,
### you should add +3
qbicDays<-wed + weeks(x = 0:(length(unique(members$X1))-1+1))

### manual  
Holidays = c(ymd("2020-09-16"),ymd("2020-01-01"),ymd("2018-12-26"),ymd("2019-01-02"))


qbicDays<-qbicDays[!qbicDays %in% Holidays]

ff<-data.frame(Day=day(x = qbicDays),month= as.vector(month(x = qbicDays,label = T)),members=sample(members$X1),Talk="X")

ff$talkday<- paste(ff$Day,ff$month,sep = " ")
head(ff)
ff1<- ff
ff

ff<-df_defactorise(ff)

ff.matrix<-daply(ff,.(members,talkday),function(x)x$Talk)

colnames(ff.matrix)

### dummy dataframe to sort matrix column 
x<-data.frame(Day=day(x = qbicDays),month= as.vector(month(x = qbicDays,label = T)),members=sample(members$X1),Talk="X")
x$talkday<- paste(x$Day,x$month,sep = " ")

ff.matrix<-ff.matrix[,order(match(colnames(ff.matrix),x$talkday))] # sort column based on date
ff.matrix[is.na(ff.matrix)]<- ""

write.xlsx(x = ff.matrix,file = "~/ownCloud_Qbic/R/Utility_scripts/Wednesday_talks_scheduler/wednesdaytalk_R2_2020.xlsx",sheetName = "R1",col.names = T,row.names = T,showNA = F)

#### send the xlsx file to qbic staff  and wait for feedback


##### modify slots after feedback
ff.matrix<-read.xlsx("~/ownCloud_Qbic/R/Utility_scripts/Wednesday_talks_scheduler/wednesdaytalk_R2_2020.xlsx",colNames = T)
ff.matrix

dd<-melt(ff.matrix,id.vars = "")
dd$variable<-gsub(pattern = "X",replacement = "",dd$variable)
dd$value<-as.vector(dd$value)
dd[!dd$value %in% "X","value"]<-"-"
gsub(pattern = ".*\\.","",dd$variable)
#dd$variable<-paste(as.numeric(gsub(pattern = "\\..*","",dd$variable))+2,gsub(pattern = ".*\\.","",dd$variable),sep = ".")

dd$variable<-gsub("31.Jun","1.Jul",dd$variable)

##### to switch places 

dd[dd$variable %in% "25.Sep" & dd$NA. %in% "Chris","value"]<-"-"
dd[dd$variable %in% "2.Oct" & dd$NA. %in% "Chris","value"]<-"X"

dd[dd$variable %in% "25.Sep" & dd$NA. %in% "Lukas H.","value"]<-"X"
dd[dd$variable %in% "2.Oct" & dd$NA. %in% "Lukas H.","value"]<-"-"


#dd[dd$variable %in% "23.Oct" & dd$NA. %in% "Friederike","value"]<-"-"
#dd[dd$variable %in% "4.Dec" & dd$NA. %in% "Friederike","value"]<-"X"

#dd[dd$variable %in% "23.Oct" & dd$NA. %in% "Julian","value"]<-"X"
#dd[dd$variable %in% "4.Dec" & dd$NA. %in% "Julian","value"]<-"-"


dd$variable<-factor(dd$variable,unique(dd$variable))
colnames(dd)[1]<-"NA."
df1<-as.data.frame(daply(dd,.(NA.,variable),function(x)x$value))
#df1$`31.Jul`<-""
#dim(df1)
#df1[20,]<- c(rep(x = "",19),"X")
#rownames(df1)[20]<-"Tobias"


#df1

ff.matrix<-as.matrix(df1)
ff.matrix[ff.matrix == "-"]<-""
#as.Date(x = "10Apr2019","%d%B%Y")


ff.matrix
write.xlsx(x = ff.matrix,file = "~/ownCloud_Qbic/R/Utility_scripts/Wednesday_talks_scheduler/wednesdaytalk_R2_2020.xlsx",sheetName = "R1",col.names = T,row.names = T,showNA = F)


library(grid)
library(gridExtra)
library(gtable)
pdf("~/Desktop/wednesdaytalk_R2_2020.pdf",encoding = "ISOLatin9.enc",width = 16,height = 8,pointsize = 12)
title <- textGrob(label = "QBiC wednesday talks R2 2020",gp=gpar(fontsize=30))
padding <- unit(5,"mm")

table <- gtable_add_rows(
  tableGrob(ff.matrix), 
  heights = grobHeight(title) + padding,
  pos = 0)
table <- gtable_add_grob(
  table, 
  title, 
  1, 1, 1, ncol(table))
#grid.newpage()
grid.draw(table)
#grid.table(ff.matrix)
dev.off()




### automatic calender event ####

ff.matrix<-read.xlsx("~/ownCloud_Qbic/R/Utility_scripts/Wednesday_talks_scheduler/wednesdaytalk_R2_2020.xlsx",colNames = T)
ff.matrix


dd<-melt(ff.matrix,id.vars = "")
head(dd)

new.dd<-dd[dd$value %in% "X",]

wed<-dmy_hms("05/08/2020 13:30:00",tz = "Europe/Berlin")
qbicDays<-wed + weeks(x = 0:21)
Holidays = c(ymd("2020-09-16"),ymd("2020-01-01"),ymd("2018-12-26"),ymd("2019-01-02"))
#qbicDays<-qbicDays[!qbicDays %in% Holidays]
qbicDays<-grep(paste(Holidays,collapse = "|"),qbicDays,value = T,invert = T)
#df.dm<- data.frame(mon=unique(as.vector(month(qbicDays,label = TRUE))),
#             monn=unique(month(qbicDays)))
#dd$variable

endTime<-dmy_hms("05/08/2020 14:00:00",tz = "Europe/Berlin")
endTime<-endTime + weeks(x = 0:21)
#qbicDays<-qbicDays[!qbicDays %in% Holidays]
endTime<-grep(paste(Holidays,collapse = "|"),endTime,value = T,invert = T)

#df.dm<-data.frame(day=day(qbicDays),mon=month(qbicDays,label = F))

#new.dd$startTime<-paste(df.dm$mon,df.dm$day,"2020 13:00:00",sep = "/")
#new.dd$endTime<-paste(df.dm$mon,df.dm$day,"2020 14:00:00",sep = "/")
new.dd$startTime<- qbicDays
new.dd$endTime<-endTime

colnames(new.dd)[1]<-"name"



gmail<-read.xlsx("~/ownCloud_Qbic/R/Utility_scripts/Wednesday_talks_scheduler/QBiC - Gmail Accounts and Phone numbers.xlsx",colNames = T)

#gmail$name[24]<-"Silvia"
new.dd$name
gmail$name[c(19,20)]<- c("Steffen.L","Steffen.G")


new.dd$name %in% gmail$name


new.dd$name[12]<-"Sven.N"
new.dd$name[17]<-"Sven.F"

gmail$name[4]<-"Sven.F"
gmail$name[16]<-"Sven.N"


ngf<-merge(new.dd,gmail,by="name",all.x=T)
nggf<-ngf[,-c(8,9)]

nggf$main<-"qbic.tue@gmail.com"

gs4_auth()
TFtalk<-gs4_create(name = "TFTalksCaleventTest",sheets = "sheet 1",
                   locale = "de_DE",
                   timeZone = "Europe/Berlin")


nggf %>% sheet_write(TFtalk,sheet = "sheet 1")
#write.xlsx(nggf,file = "ownCloud_Qbic/R/Utility_scripts/Wednesday_talks_scheduler/wednesdaytalk_R2_2020_calenderEvent.xlsx")

#get the addon on google sheets "Get addons" Mail Merge with Attachments
#have a look at QBiC_R12020_Mondayremainder
#to create the events QBiC_talk_R12020_calenderEvent (TFTalksCalevent)


