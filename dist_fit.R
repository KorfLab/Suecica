library(MASS)	#fitdistr

files = as.matrix(read.csv("HMMInputfiles.txt", header=FALSE, sep="\n"))
for (file in files) {
	print(file)#
	data <- read.csv(file,header=FALSE, col.names = c("Coverage","Frequency"), sep="\t", skip=1)
	cov <- unlist(data[1])
	freq <- unlist(data[2])
	actual_cov <- rep(cov,freq)
	
	pois_fit <- fitdistr(actual_cov, "poisson")
	
	pois_lambda <- unlist(pois_fit$estimate[1])
	pois_negtwologlik <- unlist(pois_fit$loglik[1]) * -2
	
	variables <- c("Poisson_lambda", "Poisson_-2loglik")
	var_data <- c(pois_lambda, pois_negtwologlik)
	outtable <- data.frame(variables, var_data)
	outfile <- paste(gsub(".txt","",file),"_distfits.txt", sep="")
	write.matrix(outtable,outfile, sep="\t")
}
