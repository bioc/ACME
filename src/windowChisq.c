#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>


double chi_square_calc( int a, int b, int c, int d ) {
  double chi_square;
  double da = (double)a;
  double db = (double)b;
  double dc = (double)c;
  double dd = (double)d;

  chi_square = (((da * dd) - (db * dc)) * ((da * dd) - (db * dc)) * (da + db + dc + dd))/((da + db) * (dc + dd) * (db + dd) * (da + dc));
  return(chi_square);
}

SEXP windowChisq(SEXP locations,SEXP values,SEXP windowSize, SEXP totprobes, SEXP posprobes)
{
  int i,j,na,firstprobe,lastprobe;
  double *xa,*xlocations,*xchivals;
  int *xwindowSize,*xnprobes,*xtotprobes,*xposprobes,*xvalues,*xsums;
  double halfWindowSize;
  SEXP sums,nprobes,list,chivals;
  
  PROTECT(locations = AS_NUMERIC(locations));
  xlocations = NUMERIC_POINTER(locations);
  na = LENGTH(locations);

  PROTECT(values = AS_INTEGER(values));
  xvalues = INTEGER_POINTER(values);

  PROTECT(windowSize = AS_INTEGER(windowSize));
  xwindowSize = INTEGER_POINTER(windowSize);
  halfWindowSize = xwindowSize[0]/2;

  PROTECT(totprobes = AS_INTEGER(totprobes));
  xtotprobes = INTEGER_POINTER(totprobes);

  PROTECT(posprobes = AS_INTEGER(posprobes));
  xposprobes = INTEGER_POINTER(posprobes);

  PROTECT(sums = NEW_INTEGER(na));
  xsums = INTEGER_POINTER(sums);

  PROTECT(chivals = NEW_NUMERIC(na));
  xchivals = NUMERIC_POINTER(chivals);

  PROTECT(nprobes = NEW_INTEGER(na));
  xnprobes = INTEGER_POINTER(nprobes);

  for(i = 0; i < na; i++) {
    // Find window extremes for each point
    firstprobe=i;
    lastprobe=i;
    // look for leftmost probe in window
    while((xlocations[i]-xlocations[firstprobe-1] < (halfWindowSize)) &
	  firstprobe>0) {
      firstprobe += -1;
    }
    // look for rightmost probe in window
    while((xlocations[lastprobe+1]-xlocations[i] < (halfWindowSize)) &
	  lastprobe<(na-1)) {
      lastprobe += 1;
    }
    xnprobes[i] = lastprobe-firstprobe+1;
    // get sum of the number of probes in the window
    xsums[i] = 0;
    for(j = firstprobe; j <= lastprobe; j++) {
      xsums[i] += xvalues[j];
    }
    // compute chi_square VALUE, not p-value
    xchivals[i] = chi_square_calc(xposprobes[0],xtotprobes[0]-xposprobes[0],xsums[i],xnprobes[i]-xsums[i]);
  }
  PROTECT(list = allocVector(VECSXP,4));
  SET_VECTOR_ELT(list, 0, sums);
  SET_VECTOR_ELT(list, 1, nprobes);
  SET_VECTOR_ELT(list, 2, chivals);
  SET_VECTOR_ELT(list, 3, values);
  
  UNPROTECT(9);
  return(list);
}

R_CallMethodDef callMethods[] = {
    {"windowChisq",&windowChisq,5},
    {NULL,NULL,0}
};

void R_init_ACME(DllInfo *info)
{
    R_registerRoutines(info,NULL,callMethods,NULL,NULL);
}
