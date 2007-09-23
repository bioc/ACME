setClass("ACME",
         contains="ExpressionSet")

setClass("ACMECalc",
         representation(windowsize="integer",
                        threshold="numeric"),
         prototype=list(windowsize=integer(),
           threshold=numeric()),
         contains="ACME")

