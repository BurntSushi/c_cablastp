#ifndef __CABLASTP_UTIL_H__
#define __CABLASTP_UTIL_H__

#include <stdio.h>

char *
trim(char *s, const char *totrim);

char *
trim_space(char *s);

int
readline(FILE *f, char **line);

#endif
