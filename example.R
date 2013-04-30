#------------------------------------------------------------------------------+
# FILE: example.R                                                              |
# DATE: 2012-02-28, 16:50:05                                                   |
#------------------------------------------------------------------------------;

## Import every R file in the directory 'R'
for (filename in dir(path = './R', pattern = "[^example].*\\.R")) {
    source(paste('./R', filename, sep = '/'))
}

## Compile the library if it is missing or out of date, then load it.
mylib <- paste('./src/rpbounds', .Platform$dynlib.ext, sep = '')
mysrc <- paste('./src/rpbounds', '.c', sep = '')
tryCatch(expr  = dyn.unload(mylib),
         error = function(.e) invisible(NULL))
if (!file.exists(mylib) || file.info(mylib)$mtime < file.info(mysrc)$mtime) {
    system(paste('R CMD SHLIB -c', mysrc))
}
dyn.load(mylib)

## Define a 'riskproc' object containing a compound Poisson risk process
## with diffusion with hypo-exponentially(1, 10) distributed claims,
## a premium rate of 2, a claim frequency of 1 and a diffusion
## volatility of sqrt(0.4).  Note that all objects necessary for the
## subsequent calculations are stored within this 'riskproc' object.
myprocess <- riskproc(claims   = claiminfo(hypoexp = list(rates = c(1.0, 10.0))),
                      premium  = 2.0,
                      freq     = 1.0,
                      variance = 0.4)

## Uncomment this in order to force 'boundsRuinprob()' to compute
## a rough approximation to the integrated tail area distribution
## function.
# myprocess$claims$cdf.tailarea <- NULL

## In case of hypo-exponentially distributed claim amounts, we are able
## to compute the theoretical probability of ruin.  So this is the real
## deal!
psi <- ruinprob(process = myprocess,
                method  = 'hypoexp')

## Use the other methods to compute an approximation to the probability
## of ruin.

psi.fft <- ruinprob(process    = myprocess,
                    method     = 'fft',
                    interval   = 0.01,
                    maxreserve = 40.96)

psi.s <- ruinprob(process   = myprocess,
                  method    = 'saddlepoint',
                  jensen    = FALSE,
                  normalize = FALSE)

psi.sn <- ruinprob(process   = myprocess,
                   method    = 'saddlepoint',
                   jensen    = FALSE,
                   normalize = TRUE)

psi.s.star <- ruinprob(process = myprocess,
                       method  = 'saddlepoint',
                       jensen  = TRUE)

psi.b <- ruinprob(process     = myprocess,
                  method      = 'bounds',
                  richardson  = TRUE,
                  interval    = 0.01,
                  maxreserve  = 10.0,
                  use.splines = FALSE)

## Unload the library again (so that it can be recompiled in the
## meantime) and do some clean-up.
tryCatch(expr  = dyn.unload(mylib),
         error = function(.e) invisible(NULL))
rm(mylib, mysrc, filename)
