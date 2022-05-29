#include "WolframLibrary.h"

EXTERN_C DLLEXPORT int Initialize_squared(WolframLibraryData libData);

EXTERN_C DLLEXPORT void Uninitialize_squared(WolframLibraryData libData);

EXTERN_C DLLEXPORT int squared(WolframLibraryData libData, mreal A1, mreal *Res);

