GRM <- function(Z, effect=c("A", "D", "AA", "AD", "DD"))
{
# Note: Z is the genotype matrix coded as 0, 1, 2

	fast_cp <- function(x){
		r.open <- !inherits(try(Revo.version, silent=TRUE),"try-error")
		if(is.null(x)){
			stop("Please assign the matrix!")
		}
		if(r.open){
			cpd <- crossprod(x)
		}else{
			cpd <- crossprodcpp(x)
		}
		return(cpd)
	}
	if(!is.matrix(Z))	Z <- as.matrix(Z)
	switch(
		match.arg(priority),
		"A" = {
			SNPmean <- rowMeans(Z)
			Z <- Z - SNPmean
			K <- fast_cp(Z)
		},
		"D" = {
			Z[which(Z == 2, arr.ind=TRUE)] <- 0
			SNPmean <- rowMeans(Z)
			Z <- Z - SNPmean
			K <- fast_cp(Z)
		},
		"AA" = {
			SNPmean <- rowMeans(Z)
			Z <- Z - SNPmean
			K <- fast_cp(Z) * fast_cp(Z) - fast_cp(Z * Z)
		},
		"AD" = {
			W <- Z
			W[which(W == 2, arr.ind=TRUE)] <- 0
			SNPmean <- rowMeans(Z)
			Z <- Z - SNPmean
			SNPmean <- rowMeans(W)
			W <- W - SNPmean
			K <- fast_cp(Z) * fast_cp(W) - fast_cp(Z * W)
		},
		"DD" = {
			Z[which(Z == 2, arr.ind=TRUE)] <- 0
			SNPmean <- rowMeans(Z)
			Z <- Z - SNPmean
			K <- fast_cp(Z) * fast_cp(Z) - fast_cp(Z * Z)
		}
	)
	return(K/mean(diag(K)))
}
