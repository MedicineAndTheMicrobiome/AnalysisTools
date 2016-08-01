/* 
  Input parameters are:
    - file with list of ids to extract
    - <fasta/fastq> for mate
    - keep/remove
    - fasta/fastq
    - mate1/mate2/single
  or mate file from stdin and other three.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXLINELEN 1000000
#define WHITE_SPACE " \n\r\t\v"
#define UNKNOWN -2
#define NIL '\0'

#define FIRST_MATE  1
#define SECOND_MATE  2
#define SINGLE  3

/* TRUE and FALSE may be already defined, */
/* let's not get in the way of other definitions. */
 
#if     defined(TRUE)
#undef TRUE 
#endif
 
#if     defined(FALSE)
#undef FALSE
#endif
 
typedef enum {FALSE = 0, TRUE = 1} bool; /* boolean type */

int token_list_count;
char **token_list;

void sort_token_lists(int, int);
void swap_c(int, int);
int find_index(char *);
void trim_id(char *);

bool KEEP, FASTQ;
int READ_TYPE;

const int MinArgs = 4;
const char ArgInfo[] =
"Syntax: extract_fullseq <id_list> <-keep|-remove> <-fasta|-fastq> <-single|-mate1|-mate2> [<largefile>]\n";

int main(int argc, char *argv[])
{ 
    int index, line_number;
    FILE *infile = NULL;
    char descriptions[MAXLINELEN+1];
    char copy_line[MAXLINELEN+1];
    char id[MAXLINELEN+1];
    char *token;
    int BLOCK = 100000;
    bool PRINT;

    if (((argc != MinArgs + 1) && (argc != MinArgs+2)) ||
	((strcmp(argv[2],"-keep") != 0) &&
         (strcmp(argv[2],"-remove") != 0)) ||
	((strcmp(argv[3],"-fasta") != 0) &&
         (strcmp(argv[3],"-fastq") != 0)) ||
	((strcmp(argv[4],"-mate1") != 0) &&
         (strcmp(argv[4],"-mate2") != 0) &&
         (strcmp(argv[4],"-single") != 0)))
    {
       fprintf(stderr,ArgInfo);
       exit(EXIT_FAILURE);
    }

    if (strcmp(argv[2],"-keep") == 0)
        KEEP = TRUE;
    else
        KEEP = FALSE;

    if (strcmp(argv[3],"-fastq") == 0)
        FASTQ = TRUE;
    else
        FASTQ = FALSE;

    if (strcmp(argv[4],"-single") == 0)
        READ_TYPE = SINGLE;
    else
    {
        if (strcmp(argv[4],"-mate1") == 0)
            READ_TYPE = FIRST_MATE;
        else
            READ_TYPE = SECOND_MATE;
    }

    infile = fopen(argv[1], "r"); 
    if (infile == NULL)
    { 
       fprintf(stderr, "Cannot open %s\n", argv[1]); 
       exit(EXIT_FAILURE);
    }
    
    token_list_count = 0;
    while(fgets(descriptions, MAXLINELEN , infile) != NULL) 
    {
        if ((token_list_count%BLOCK) == 0)
        {
            token_list = (char **) realloc(token_list, (token_list_count+BLOCK) * sizeof(char *));
            if (token_list == NULL)
            {
                fprintf(stderr,"Cannot realloc for %d ids\n", token_list_count+BLOCK);
                exit(EXIT_FAILURE);
            }
        }
    
        token = strtok(descriptions, WHITE_SPACE);
        if (token == NULL)
        {
            fprintf(stderr,"No id on line %d\n", token_list_count+1);
            exit(EXIT_FAILURE);
        }
        token_list[token_list_count] = (char *) strdup(token);
        token_list_count++;
    }
    fclose(infile);
    sort_token_lists(0, token_list_count-1);

    if (argc == MinArgs + 1)
        infile = stdin;
    else
    {
        infile = fopen(argv[5], "r"); 
        if (infile == NULL)
        { 
           fprintf(stderr, "Cannot open %s\n", argv[5]); 
           exit(EXIT_FAILURE);
        }
    }
    
    line_number = 0;
    while(fgets(descriptions, MAXLINELEN , infile) != NULL) 
    {
        line_number++;

        strcpy(copy_line, descriptions);
        if (((FASTQ == TRUE) && ((line_number%4) == 1 )) || //descriptions[0] == '@')) ||
            ((FASTQ == FALSE) && (descriptions[0] == '>')))
        {
            token = strtok(descriptions, WHITE_SPACE);
            strcpy(id, &(token[1]));
            trim_id(id);
            index = find_index(id);
            if (((KEEP == TRUE) && (index != UNKNOWN)) ||
                ((KEEP == FALSE) && (index == UNKNOWN)))
                PRINT = TRUE;
            else
                PRINT = FALSE;
        }

        if (PRINT == TRUE)
            printf("%s", copy_line);
    }
    if (argc != MinArgs + 1)
        fclose(infile);
}

void sort_token_lists(left, right)
int left;
int right;
{
  int i, last;

  if (left >= right) return;
  swap_c(left, (left+right)/2);
  last =left;
  for (i=left+1; i<=right; i++){
    if (strcmp(token_list[left], token_list[i]) > 0)
      swap_c(++last, i);
  }
  swap_c(left,last);

  sort_token_lists(left,last-1);
  sort_token_lists(last+1, right);
}

void swap_c(index,current)
int index;
int current;
{
    char *name;

    name = token_list[current];
    token_list[current] = token_list[index];
    token_list[index] = name;
}

int find_index(first)
char *first;
{
    int index, left, right;

    left = 0;
    right = token_list_count - 1;

    while (left <= right)
    {
        index = (int)((left+right)/2);
        if (strcmp(first, token_list[index]) == 0)
            return(index);
        if (strcmp(first, token_list[index]) < 0)
            right = index-1;
        else
            left = index+1;
    }
    return(UNKNOWN);
}

/* handle case where even though mates are present, the read id may not have suffix */
void trim_id(given)
char *given;
{
    int length, mate;

    switch (READ_TYPE) {
        case SINGLE: return;
        case FIRST_MATE: mate = 1; break;
        case SECOND_MATE: mate = 2; break;
        default: fprintf(stderr,"Unrecognized mate type: %d\n", READ_TYPE);
                 exit(EXIT_FAILURE);
    }

    length = strlen(given);
    if ((length > 2) && (given[length-1] == ('0'+mate)) &&
        ((given[length-2] == '.') ||
         (given[length-2] == '/') ||
         (given[length-2] == '_')))
        given[length-2] = NIL;
}
