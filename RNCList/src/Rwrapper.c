/* R can only get info back from C via its arguments, so we 
   wrap the NCList functionality into code which makes this 
   possible and tell R about the available functions */


#include "R_ext/Rdynload.h"


/* Register functions with R */
void R_init_libnclist(){

}
