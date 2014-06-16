hypergeom.test <-
function( x, alternative='two.sided', pval.method="minimum.likelihood", nthreads=2, ... ){
	res <- .hypergeax.test( x=x, nthreads=nthreads, ... )

	p.value <- 1;

	if(alternative %in% c("less", "greater")){
		p.value <- res$prob[alternative]
		pval.method <- "exact"
	}else{
		p.value <- switch(pval.method, minimum.likelihood = res$prob['two.sided_ml'], double = res$prob['two.sided_d'])
	}

	dm <- dim(x)
	method <- "Exact hypergeometric test for a "
	if( length(dm) == 3 ){ method <- paste(method, "2x2x2 table", sep="") }
	if( length(dm) == 2 ){ I <- dm[1]; J <- dm[2]; method <- paste(method, I, "x", J, " table", sep="") }

	obj <- list(p.value=as.numeric(p.value), pval.method=as.character(pval.method), alternative=alternative, statistic=setNames(res$count.obs, "a"), estimate=setNames(res$odds.ratio, "odds ratio"), prob.obs=res$prob.obs, method=method
	, data.name=paste(c(deparse(substitute(x)), res$call), collapse=", "), conf.int=res$conf.int)
	class(obj) <- "htest"

	if( !is.null(list(...)$r) ){ obj <- list(obj, res) }
return(obj)
}
