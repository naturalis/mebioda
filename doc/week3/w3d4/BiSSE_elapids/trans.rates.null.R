trans.rates =TransMatMaker(hidden.states=TRUE)
trans.rates.nodual =ParDrop(trans.rates,c(3,5,8,10))
trans.rates.nodual.allequal =ParEqual(trans.rates.nodual,c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal
