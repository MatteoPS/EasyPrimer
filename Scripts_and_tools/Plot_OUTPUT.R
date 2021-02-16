#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

genename = args[1]
out_folder = args[2]
cwd = args[3]
snp = args[4]
prim_thr = args[5]
amp_thr = args [6]
jobname = args[7]
conslimit = args[8]
minpri = args[9]
maxpri = args[10]
minamp = args[11]
maxamp = args[12]
tmp_folder = args[13]

#adding the complete paths
tmp_folder = paste(cwd,tmp_folder,"/", sep="")
out_folder = paste(cwd,out_folder,"/", sep="")

today = format(Sys.Date(), "%d-%m-%Y")

logfile = paste(tmp_folder, 'log/LOG_Plotting_',today,'.txt', sep = "")


cat(paste(format(Sys.time()),"===================START on ",genename,"=================\n", sep = ""), file = logfile, append = TRUE)
# this allows the use of a prefix in the name of the output as well as in the tmp files
if (jobname == 'empty')
{

  tab<-read.csv(paste(tmp_folder, 'tmp/', genename,'_TAB.csv', sep = ""))
  score<-read.csv(paste(tmp_folder, 'tmp/', genename,'_SCORES.csv', sep = ""))
  cons_thr_tab <-read.table(paste(tmp_folder, 'tmp/',genename,'_consensus.txt', sep = ''), sep = '\t')
  cons_thr <- read.table(paste(tmp_folder, 'tmp/',genename,'_consensus_threshold.txt', sep = ''), sep = '\t')$V2
}else
{

  cat(paste(format(Sys.time()),"trying importing: ",paste(tmp_folder, 'tmp/',jobname, genename,'_TAB.csv', sep = ""), sep = ''), file = logfile, append = TRUE)
  tab<-read.csv(paste(tmp_folder, 'tmp/',jobname, genename,'_TAB.csv', sep = ""))
  score<-read.csv(paste(tmp_folder, 'tmp/',jobname, genename,'_SCORES.csv', sep = ""))
  cons_thr_tab <-read.table(paste(tmp_folder, 'tmp/',jobname,genename,'_consensus.txt', sep = ''), sep = '\t')
  cons_thr <- read.table(paste(tmp_folder, 'tmp/',jobname,genename,'_consensus_threshold.txt', sep = ''), sep = '\t')$V2
}
prim_thr = as.numeric(prim_thr)
amp_thr = as.numeric(amp_thr)
conslimit = as.numeric(conslimit)

genelen <- length(tab$Position)
num_seq <- sum(tab[1,2:6])
q_prim <-  unname(quantile(score$Score_MaxPrimer, prim_thr))
q_amp <-  unname(quantile(score$Score_amplicon, amp_thr))

cat(paste(format(Sys.time())," CSV table files corretly imported\n", sep = ''), file = logfile, append = TRUE)

# selection from the score Dataframe by the thresholds given
cat(paste(format(Sys.time())," Selection from the score Dataframe by the thresholds given\n", sep = ''), file = logfile, append = TRUE)
select<-subset(score, (score$Score_amplicon >= q_amp & score$Score_MaxPrimer <= q_prim))
if (nrow(select) == 0)
{
  print("###ERROR### PROGRAM HAS STOPPED---> one or both the treshold parameters (-prithr, -ampthr) are too stringent, NO data is selected\n")
  cat(paste(format(Sys.time()),"###ERROR### PROGRAM HAS STOPPED---> one or both the treshold parameters (-prithr, -ampthr) are too stringent, NO data is selected\n", sep = ''), file = logfile, append = TRUE)
  quit()
}


cat(paste(format(Sys.time())," Selection completed\n", sep = ''), file = logfile, append = TRUE)



prim_pos_list=c()

for (i in 1:length(row.names(select)))
{
  prim_pos_list = c(prim_pos_list,(select[i,"p1"]:select[i,"p2"]),(select[i,"p3"]:select[i,"p4"]))

}

amp_pos_list=c()
for (i in 1:length(row.names(select)))
{
  amp_pos_list = c(amp_pos_list,select[i,"p2"]:select[i,"p3"])
}

all_pos_list<-c(prim_pos_list, amp_pos_list)

# to create the lines represented in the final graph
h_primer<-hist(prim_pos_list, breaks=1:max(all_pos_list), plot=F)
h_amplicon<-hist(amp_pos_list, breaks=1:max(all_pos_list), plot=F)

y_max=max(c(h_primer$counts, h_amplicon$counts))
x_max=max(score$p4)

# consenus bases written above their position in the gene
consensus = tab$consensus
cons_str = as.character(consensus)
marks = seq(0, x_max, by=5)

# corrects the legend and the filename on the basis of the analysis (HRM detectable snps or ALL the snps)
if (snp == 'HRM')
{
  snp_legend = "HRM SNPs"
  name_sufx = "_HRM_analysis"
  snp_col = "gray"
}else if (snp == 'ALL')
{
  snp_legend = ""
  name_sufx = "_ALL_SNPs_analysis"
  snp_col = "white"
}else
{
  print('###ERROR### PROGRAM HAS STOPPED---> -snp parameter provided is not allowed: it must be either «HRM»(default) or «ALL»')
  cat(paste(format(Sys.time())," ###ERROR### PROGRAM HAS STOPPED---> -snp parameter provided is not allowed: it must be either «HRM»(default) or «ALL»\n", sep = ''), file = logfile, append = TRUE)
  quit()
}

if (jobname == 'empty')
{
  title1 = paste(genename, "Analysis",sep=" ")
  outplotname = paste(out_folder, genename, name_sufx, ".pdf", sep = "")

}else
{
  outplotname = paste(out_folder, jobname, genename, name_sufx, ".pdf", sep = "")
  title1 = paste(jobname, genename, "Analysis",sep=" ")
}

col_picchi = 'gray'

cat(paste(format(Sys.time()),paste(genename, name_sufx, ".pdf", sep = ""), " will be the output PDF file\n", sep = ' '), file = logfile, append = TRUE)
pdf(outplotname, height = 5, width = (genelen*0.05))

#display peaks only  when fHRM > conslimit
#tab$fHRM[tab$fHRM < conslimit] <- 0

plot(tab$fSNP , type="l",xlim=c(1,x_max),ylim=c(0,1), xaxt= 'n',
     main=title1, las=2, col=col_picchi, xlab = 'Gene positions', ylab = '', bty="n")
cat(paste(format(Sys.time()),"Plotting...\n", sep = ' '), file = logfile, append = TRUE)

# add gray dotted lines on ambiguous residues HRM SNPs
if (snp == 'HRM')
{
  for (k in 1:genelen)
  {
    if (toString(tab$consensus[k]) == 'R')
    {
      abline(v=k, col="gray",lwd=0.6, lty = 2)
    }
    else if (toString(tab$consensus[k]) == 'M')
    {
      abline(v=k, col="gray",lwd=0.6, lty = 2)
    }
    else if (toString(tab$consensus[k]) == 'K')
    {
      abline(v=k, col="gray",lwd=0.6, lty = 2)
    }
    else if (toString(tab$consensus[k]) == 'Y')
    {
      abline(v=k, col="gray",lwd=0.6, lty = 2)
    }
  }
}

lines(h_primer$mids, h_primer$counts/y_max,lwd=2,col="dodgerblue4")
lines(h_amplicon$mids, h_amplicon$counts/y_max, type="l",col="red",lwd=2)
axis(1, (1:x_max), labels = FALSE ,cex.axis=0.2 ,tcl=0.5,lwd.ticks=1, col.ticks='gray',mgp=c(0,0,0.5))
axis(1, marks, las=2, cex.axis=0.7 ,tcl=0.6,lwd.ticks=1, mgp=c(0,0.2,0.5))


#red colored gap means is greater than thr and less than 1-thr of the consensus
vec <- rep(FALSE, length(tab$fgap))

for (k in 1:length(tab$fgap))
{
  if (tab$fgap[k] >= cons_thr & tab$fgap[k] <= 1-cons_thr)
  {
    vec[k] <- TRUE
  }

}
colors <- ifelse(vec, "red", "black")

text(1:x_max, y=0.02 , cons_str, cex = 0.3, pos = 1, font=2, col = colors)

legend("topleft", legend = c("Suggested Primers Areas","Suggested Amplicon Areas", "SNPs frequency", snp_legend, paste("- Primer range:", minpri,"-", maxpri, "nt", sep=" "),
                             paste("- Amplicon range:", minamp,"-", maxamp, "nt", sep=" "),paste("- Consensus limit:", (conslimit*100),"%", sep=" "),paste("- Sequences:", num_seq, sep=" ")),
       col=c("dodgerblue4","red", col_picchi, snp_col, 'white', 'white', 'white'),lty = c(0,0,0,2,0,0,0),pch = c(22, 22, 22,NA,22,22,22),pt.bg = c("dodgerblue4","red", col_picchi,NA, NA, NA, NA),pt.cex = 2, bty = "n")

cat(paste(format(Sys.time()),paste(genename, name_sufx, ".pdf", sep = ""), "is avilable\n" , sep = ' '), file = logfile, append = TRUE)
cat(paste(format(Sys.time()),"===================DONE! on ",genename,"=================\n\n", sep = ''), file = logfile, append = TRUE)

#####plot to see where primer thr is into the score distribution

#pdf(tmp_plot_name, height = 5, width = (genelen*0.05))
#plot(density(score$Score_MaxPrimer),main=title2, xlab= 'Sannon index', bty='n', lwd=1.5 )
#abline(v=q_prim, col="red",lwd=2, lty = 2)
#legend("topright", legend=c("Score distribution", "Primer threshold"), col=c("black", "red"), lty=1:2, lwd = 1.5:2, bty='n')
#cat(paste(format(Sys.time()),tmp_plot_name, "contains score distribution\n", sep = ' '), file = logfile, append = TRUE)
#cat(paste(format(Sys.time()),"=====================R-FINISH====================\n\n", sep = ' '), file = logfile, append = TRUE)
