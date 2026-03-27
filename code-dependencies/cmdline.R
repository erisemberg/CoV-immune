pkg.globals <- new.env()
pkg.globals$cmdArgs <- commandArgs(trailingOnly=TRUE)

cmdline.setCommandArgs <- function(ca){
  pkg.globals$cmdArgs <- ca
}

cmdline.getCommandArgs <- function(){
  pkg.globals$cmdArgs
}

cmdline.integer <- function(key, ...)
  # return the value part of a command line option as an integer
{
  cmdline.integers(key, ..., howmany=1)
}

cmdline.integers <- function(key, ...){
  s <- cmdline.strings(key, ...)
  if (!is.null(s)) return (as.integer(s))
}

cmdline.flag <- function(name)
  # return TRUE iff command line flag was present
  # Note: --arg=value are not considered flags
{
  0 != length(grep(paste("^--",name,"$", sep=""), commandArgs(trailingOnly=TRUE)))
}

cmdline.logical <- function(key, ...)
  # return value part of command line option as numeric
{
  cmdline.logicals(key, ..., howmany=1)
}

cmdline.logicals <- function(key, ...)
  # return value part of command line option as numeric
{
  s <- cmdline.strings(key,...)
  if (!is.null(s)){           
    ints=integer(length(s))
    for (i in 1:length(s)){
      ints[i] = as.logical(as.integer(switch(s[i], "T"=1, "F"=0, "TRUE"=1, "FALSE"=0, s[i])))
    }
    return (as.logical(ints))
  }
}

cmdline.numeric <- function(key, ...)
  # return value part of command line option as numeric
{
  cmdline.numerics(key, ..., howmany=1)
}

cmdline.numerics <- function(key, ...)
  # return value part of command line option as numeric
{
  s <- cmdline.strings(key, ...)
  if (!is.null(s)) return (as.numeric(s))
}


cmdline.has.option <- function(key)
  # return true if option was specified
{
  !is.null(cmdline.option(key, allow.omit=TRUE))
}

cmdline.option <- function(key, default=NULL,
                           stop.on.fail=TRUE,
                           allow.omit=!stop.on.fail,
                           allowed.values=NULL)
  # return the value part of a command line option
{
  ca <- grep("=", grep("^--", pkg.globals$cmdArgs, value=TRUE), value=TRUE)

  keys <- sub(pattern="=.*", replacement="", ca)
  keys <- sub(keys, pattern="--", replacement="")
  i <- match(key, keys)
  if (is.na(i))
  {
    if (is.null(default))
    {
      if (!allow.omit) stop("Could not find key ", key, "\n")
      return (NULL)
    }
    else
    {
      return (default)
    }
  }
  values <- sub(pattern=".*=", replacement="", ca)

  if (!is.null(allowed.values))
  {
    ok <- values[i] %in% allowed.values
    if (!all(ok))
    {
      stop("Illegal values for key ", key, ": ",
           paste(sep=" ,", values[i][!ok]), "\n")
    }
  }

  return (values[i])
}

# for consistency
cmdline.string <- function(key, ...)
{
  cmdline.option(key, ...)
}

cmdline.strings <- function(key, howmany=c(0,Inf), default=NULL, ...)
  # return comma separated values
{
  string <- cmdline.option(key, default=default, ...)
  if (identical(string, default)){
    return (default)
  }
  if (is.null(string)) return (NULL)
  strings <- unlist(strsplit(string, split=",", perl=TRUE))
  k <- length(strings)
  
  howmany=rep(howmany, length.out=2)
  if (k < howmany[1] | k > howmany[2]){
    stop("Argument ",key," requires ",howmany[1],"-",howmany[2]," values, but got ",k,": \"",paste(strings,collapse="\",\""),"\"")
  }
  strings
}
