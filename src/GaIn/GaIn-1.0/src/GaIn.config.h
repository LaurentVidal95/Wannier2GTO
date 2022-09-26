/* src/GaIn.config.h.  Generated from GaIn.config.h.in by configure.  */
/* src/GaIn.config.h.in.  Generated from configure.ac by autoheader.  */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef FC_DUMMY_MAIN */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define FC_FUNC(name,NAME) name ## _

/* As FC_FUNC, but for C identifiers containing underscores. */
#define FC_FUNC_(name,NAME) name ## _

/* Name of package */
#define PACKAGE "GaIn"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "ivan.duchemin@cea.fr"

/* Define to the full name of this package. */
#define PACKAGE_NAME "LibGaIn - A library for Gaussian Integrals"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "LibGaIn - A library for Gaussian Integrals 1.0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "GaIn"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.0"

/* Version number of package */
#define VERSION "1.0"
