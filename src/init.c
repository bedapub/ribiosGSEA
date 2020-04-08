#include <R_ext/Rdynload.h>
#include "ribiosGSEA.h"

static const R_CallMethodDef callMethods[] = {
  CALLMETHOD_DEF("cpp_geneSetPerm", 3),
  CALLMETHOD_DEF("list2mat", 1),
  {NULL, NULL, 0}
};

void R_init_ribiosGSA(DllInfo *info) {
 R_registerRoutines(info, NULL, callMethods, NULL, NULL);
 /* the line below says that the DLL is not to be searched
  * for entry points specified by character strings so
  * .C etc calls will only find registered symbols
  */
 R_useDynamicSymbols(info, FALSE);
 /* R_forceSymbols call only allows .C etc calls which
  * specify entry points by R objects such as C_routineName
  * (and not by character strings)
  */
  R_forceSymbols(info, TRUE);
}
