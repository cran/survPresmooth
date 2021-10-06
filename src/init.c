#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern void alphaintegrand(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void c1integrand1(void *, void *, void *, void *, void *, void *, void *);
extern void c1integrand2(void *, void *, void *, void *, void *, void *, void *);
extern void dintegrand(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void funplugin(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void isevect(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void lscv(void *, void *, void *, void *, void *, void *, void *);
extern void nadarayawatson(void *, void *, void *, void *, void *, void *, void *, void *);
extern void pilot2forhintegrand(void *, void *, void *, void *, void *, void *);
extern void presmdensfast(void *, void *, void *, void *, void *, void *, void *, void *);
extern void presmestim(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void presmtwfast(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void simpson(void *, void *, void *);
extern void termsmise(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void termsmisenopresmooth(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"alphaintegrand",       (DL_FUNC) &alphaintegrand,        9},
    {"c1integrand1",         (DL_FUNC) &c1integrand1,          7},
    {"c1integrand2",         (DL_FUNC) &c1integrand2,          7},
    {"dintegrand",           (DL_FUNC) &dintegrand,            9},
    {"funplugin",            (DL_FUNC) &funplugin,            13},
    {"isevect",              (DL_FUNC) &isevect,              17},
    {"lscv",                 (DL_FUNC) &lscv,                  7},
    {"nadarayawatson",       (DL_FUNC) &nadarayawatson,        8},
    {"pilot2forhintegrand",  (DL_FUNC) &pilot2forhintegrand,   6},
    {"presmdensfast",        (DL_FUNC) &presmdensfast,         8},
    {"presmestim",           (DL_FUNC) &presmestim,           11},
    {"presmtwfast",          (DL_FUNC) &presmtwfast,           9},
    {"simpson",              (DL_FUNC) &simpson,               3},
    {"termsmise",            (DL_FUNC) &termsmise,            16},
    {"termsmisenopresmooth", (DL_FUNC) &termsmisenopresmooth, 12},
    {NULL, NULL, 0}
};

void R_init_survPresmooth(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
