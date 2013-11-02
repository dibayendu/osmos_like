#include "utility.h"

/* displayed the message as command line statement and exits the application */
void panic(const char *m) 
{
  	printf("%s", m);
  	exit(1);
}

