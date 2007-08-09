setClass("ACME",
         contains="ExpressionSet")

setClass("ACMECalc",
         representation(windowsize="integer",
                        threshold="numeric"),
         contains="ExpressionSet")

