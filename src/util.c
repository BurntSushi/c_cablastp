#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "util.h"

char *
trim_space(char *s)
{
    return trim(s, " \r\n\t");
}

char *
trim(char *s, const char *totrim)
{
    int i, j;
    int start, end;
    int slen, totrimlen, newlen;
    bool trimmed;
    char *news;

    slen = strlen(s);
    totrimlen = strlen(totrim);

    start = 0;
    for (i = 0; i < slen; i++) {
        trimmed = false;
        for (j = 0; j < totrimlen; j++)
            if (totrim[j] == s[i]) {
                start++;
                trimmed = true;
                break;
            }
        if (!trimmed)
            break;
    }

    end = slen - start;
    for (i = slen - 1; i >= 0; i--) {
        trimmed = false;
        for (j = 0; j < totrimlen; j++)
            if (totrim[j] == s[i]) {
                end--;
                trimmed = true;
            }
        if (!trimmed)
            break;
    }

    newlen = (end >= start) ? (end - start) : 0;
    news = malloc((newlen + 1) * sizeof(*news));
    assert(news);

    strncpy(news, s + start, newlen);
    news[newlen] = '\0';
    free(s);

    return news;
}

int
readline(FILE *f, char **line)
{
    char buf[1024];
    int allocated;

    allocated = 1; /* for \0 */
    *line = malloc(allocated * sizeof(**line));
    assert(line);
    (*line)[0] = '\0';
    while (NULL != fgets(buf, 1024, f)) {
        allocated += strlen(buf);
        *line = realloc(*line, allocated * sizeof(**line));
        assert(*line);
        strcat(*line, buf);

        /* if we have found a new line, quit */
        if ((*line)[allocated - 2] == '\n')
            break;
    }
    return allocated - 1;
}
