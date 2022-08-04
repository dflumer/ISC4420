dir.create("../Results")
pdf(paste("../Results/MergedReadlengths.pdf",sep=""))
for(i in 1:100000){
	if(!file.exists(paste("../I",i,"/I",i,"_lengths.txt",sep=""))){next}
	x=scan(paste("../I",i,"/I",i,"_lengths.txt",sep=""));
	maxlen=length(x)
	plot(x[1:(maxlen-1)],type="l",main=paste("I",i,sep=""),xlab="Length of Merged Read",ylab="Number of Merged Reads")
}
dev.off()
