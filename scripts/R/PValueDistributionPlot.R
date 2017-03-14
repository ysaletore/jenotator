pvals = sort(read.table("Sample_ActD1.bwa.001-raw-windows.bed", colClasses=c("NULL", "NULL", "NULL", "NULL", "numeric", "NULL"))$V5);
adjusted = sort(read.table("Sample_ActD1.bwa.003-filtered-windows-sorted.bed", colClasses=c("NULL", "NULL", "NULL", "NULL", "numeric", "NULL"))$V5);
plot(pvals, type='l', col='blue')
lines(adjusted, type='l', col='red')