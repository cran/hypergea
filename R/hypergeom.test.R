hypergeom.test <-
function( x, alternative='two.sided', pval.method="fisheR", nthreads=2, ... ){
	res <- .hypergeax.test( x=x, nthreads=nthreads, ... )
	dm <- dim(x)

	if( !( alternative %in% c("less", "greater", "two.sided") ) ){
        alternative <- switch( alternative, left="less", right="greater", both="two.sided" )
        message("Changed 'alternative' to valid value.")
    }
    
	PVAL.METHOD <- c("fisheR", "minimum.likelihood", "double")
	PVAL.METHODC <- c("fisheR", "ml", "d")
	p.value <- 1;
	P.value <- NULL;
	if(alternative %in% c("less", "greater")){
		p.value <- res$prob[alternative]
		pval.method <- "exact"
	}else{
		if( length(pval.method)==1 ){
			p.value <- switch(pval.method, fisheR=res$prob['two.sided_fisheR'], minimum.likelihood = res$prob['two.sided_ml'], double = res$prob['two.sided_d'])
			if( p.value == 0.0 && pval.method != 'two.sided_fisheR' ){ # p.value equal to zero can occur when, e.g. x[1,0]==0 and x[0,0] small (<<10)
				warning("P-value exactly zero. Returning P-value for minimum-likelihood approach.")
				p.value <- res$prob['two.sided_ml']
				if( p.value == 0.0 ){ 
					warning("P-value still exactly zero. Returning P-value for double approach.")
					p.value <- res$prob['two.sided_d'] 
				} 
			}
			P.value <- setNames(p.value, pval.method)
		}else{
            if( is.null(pval.method) ){ pval.method <- PVAL.METHOD }

			pval.method0 <- sapply(pval.method, function( i ){ PVAL.METHODC[ which(PVAL.METHOD==i) ] } )
			P.value <- setNames(res$prob[ paste0("two.sided_", pval.method0) ], pval.method)
			p.value <- P.value[1]
		}
	}
	p.value <- min(p.value, 1) # double approach can result in P-value equal 2, e.g. x <- cbind(c(0,0,0),c(5,0,1))
	P.value <- sapply( P.value, function( v ){ min(v, 1) } )

	
	METHOD <- "Exact hypergeometric test for a "
	if( length(dm) == 3 ){ METHOD <- paste(METHOD, "2x2x2 table", sep="") }
	if( length(dm) == 2 ){ I <- dm[1]; J <- dm[2]; METHOD <- paste(METHOD, I, "x", J, " table", sep="") }

	obj <- list(p.value=as.numeric(p.value), P.value=P.value, alternative=alternative, statistic=setNames(res$count.obs, "a"), estimate=setNames(res$odds.ratio, "odds ratio"), prob.obs=res$prob.obs, method=METHOD, conf.int=res$conf.int
	, data.name=paste(c(deparse(substitute(x)), res$call), collapse=", "))
	class(obj) <- "htest"

	if( !is.null(list(...)$r) ){ obj <- list(obj, res) }
return(obj)
}
