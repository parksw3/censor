library(tidyverse)

set.seed(2341)

runInd <- function(tsymp, tinf
	, gbar, gk, ibar, ik, Rbar, Rk
){
	numInf <- 1 + rnbinom(n=1, mu=Rbar, size=Rk)
	if (numInf==0) return(
		data.frame(pinf=NULL, psymp=NULL, inf=NULL, symp=NULL)
	)
	return(map_dfr(1:numInf, function(x){
		data.frame(pinf=tinf, psymp=tsymp
			, inf = tinf + rnbinom(n=1, mu=gbar, size=gk)
			, symp = tinf + rnbinom(n=1, mu=ibar, size=ik)
		)
	}))
}

lineMax <- 1e4

line = 0
l <- data.frame(pinf=NA, psymp=NA, inf=0, symp=0)

while((nrow(l) <= lineMax) && (line<=nrow(l))){
	line = line+1
	with(as.list(l[line, ]), {
		nline <- runInd(tsymp=symp, tinf=inf
			, gbar=2, gk=2, ibar=2, ik=2, Rbar=2, Rk=10
		)
		l <<- bind_rows(l, nline)
	})
}

print(l)
