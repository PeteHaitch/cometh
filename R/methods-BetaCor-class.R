## TODO: A show method for BetaCor objects

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###
showBetaCor <- function(x){
  

  if (print.seqlengths) {

  }
}

setMethod("show", "BetaCor",function(object){
  cat("class:", class(object), "\n")
  cat("dim:", dim(object), "\n")
  cat("methylation type:", sort(unique(object@methylation_type)), "\n")
  if (!is.na(object@NIL[1])){
    cat("NIL:", object@NIL, "\n")
  }
  if (!is.na(object@IPD[1])){
    cat("min. IPD:", min(object@IPD), "\n")
    cat("max. IPD:", max(object@IPD), "\n")
  }
  cat("feature name:", object@feature_name, "\n")
  cat("correlation method:", object@method, "\n")
  cat("genome:", unique(genome(object@seqinfo)), "\n")
  
  callNextMethod()
})