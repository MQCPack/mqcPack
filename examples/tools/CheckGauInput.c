/* C Example */
#include <math.h>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
/*
#define DEBUG
*/
typedef struct InputLine {
  char *next;
  char *last;
  int Blank;
  int Keyword;
  int Link1;
  char *Name;
  char *Lower_Name;
} InputLine;

void mqc_LOWER( char * );
void mqc_Remove_Init_WhiteSpace( char * );
int mqc_FINDstrINline( char *, char * );
void mqc_error( char * );
char *mqc_GetMatrixFile( char *, char * );
void mqc_Global_init( void );
char *mqc_DupString( char *, char * );
InputLine *Alloc_InputLine( char *, InputLine * );
void mqc_Write_InputLines( char *, InputLine * );
void mqc_Remove_Return( char * );
char *mqc_Replace_String( char *, char *, char * );

char CarriageReturn[8];
char Tab[8];
int int_tolower;
FILE *GAU_File=NULL;
FILE *CheckGAU_File=NULL;
FILE *ProfGAU_File=NULL;
InputLine *Last_Keyword = (InputLine *)NULL;
InputLine *Last_NonBlank = (InputLine *)NULL;

int main(int argc, char *argv[])
{
  char *MatFileName;
  char charBuf[1024];

#ifdef DEBUG
  printf( "Start of Program\n" );
  fflush(stdout);
#endif
  if ( argc-1 != 2 ) {
    printf( "Error in number of argments to this program.  Arguments should be:\n     - Full path to a Gaussian file\n         - This can either be a MatrixFile\n         - This can be a Gaussian input file\n     - The name of the binary for the program that \n       generates the MatrixFile, most likely \"g16\".  \n       The environment should be setup to run the program.\n" );
    if ( argc-1 == 0 ) {
      mqc_error( "No arguments to this program\n"); 
    } else {
      sprintf( charBuf, "%d arguments to this program.  The first 3 are:\n  %s\n  %s\n  %s\n", argc-1, argv[0],argv[1],argv[3]);
      mqc_error( charBuf );
    }
  }

  mqc_Global_init( );
  MatFileName = mqc_GetMatrixFile( argv[1], argv[2] );

}

char *mqc_GetMatrixFile( char *FileName, char *Program )
{
  char *MatFileName;
  char charBuf[1024];
  char lowerBuf[1024];
  char tmp_CheckGAU_File_Name[1024];
  char tmp_CheckGAU_File_NameA[1024];
  char GAUversion[1024];
  char GAUroot[1024];
  char JobName[1024];
  char InputFile[1024];
  char *external_proc;
  char *scrdir;
  char *rtn_fgets;
  int fclose_return;
  int PID;
  int Total_Jobs;
  int MatrixFile_Jobs;
  int Blank_line;
  int rtn_fputs;
  int Restart;
  int last_job;
  int Job_Number;
  int Current_job;
  int system_failure;
  int i, j, ilen;
  int GAU_dir_level;
  int found;
  int NoModNeeded;
  int external;
  InputLine *Start_Line;
  InputLine *Current_Line;
  InputLine *Last_Line;
  InputLine *Next_Line;

  rtn_fgets = (char *)NULL;
  Total_Jobs=0;
  MatrixFile_Jobs=0;
  last_job=-1;
  Job_Number=0;

#ifdef DEBUG
  printf( "Open file to store temporaries.  This should be temporary\n" );
  fflush(stdout);
#endif
  PID = (int)getpid();
  scrdir=getenv( "GAUSS_SCRDIR" );
  if ( scrdir == NULL ) {
    sprintf( tmp_CheckGAU_File_Name, "./mqc_pass_File_%d", PID);
  } else {
    sprintf( tmp_CheckGAU_File_Name, "%s/mqc_pass_File_%d", scrdir, PID);
  }
#ifdef DEBUG
  printf( " tmp_CheckGAU_File_Name %s\n", tmp_CheckGAU_File_Name );
  fflush(stdout);
#endif

  CheckGAU_File = fopen( tmp_CheckGAU_File_Name,"w");
  if ( CheckGAU_File == NULL ) {
    sprintf( charBuf, "Could not open CheckGAU file %s\n", tmp_CheckGAU_File_Name);
    mqc_error( charBuf );
  }
  fclose_return = fclose(CheckGAU_File);
  CheckGAU_File = NULL;
  if ( fclose_return != 0 ) {
    sprintf( charBuf, "Could not close CheckGAU File %s\n", tmp_CheckGAU_File_Name);
    mqc_error( charBuf );
  }
  unlink(tmp_CheckGAU_File_Name);
  CheckGAU_File = NULL;

#ifdef DEBUG
  printf( "Check if the first argument is the name of a binary MatrixFile or if it is an ASCII Gaussian input file\n" );
  fflush(stdout);
#endif

  sprintf( charBuf, "file %s > %s\n", 
	   FileName, tmp_CheckGAU_File_Name );
#ifdef DEBUG
  printf( "execute %s\n", charBuf );
  fflush(stdout);
#endif
  system_failure = system(charBuf);

  if ( system_failure == -1 ) {
    /* the file on the First input failed.  Check if the file exists. */
    GAU_File = fopen( FileName,"r");
    if ( GAU_File == NULL ) {
      sprintf( charBuf, "Could not open Gaussian File specified by the first argument %s\n", FileName);
    } else {
      sprintf( charBuf, "file command on %s failed.  Tried to do this in the Bash shell\n", FileName);
    }
    CheckGAU_File = fopen( tmp_CheckGAU_File_Name,"w");
    unlink(tmp_CheckGAU_File_Name);
    mqc_error( charBuf );
  }

  CheckGAU_File = fopen( tmp_CheckGAU_File_Name,"r");
  if ( CheckGAU_File == NULL ) {
    sprintf( charBuf, "Could not open CheckGAU file %s, should contain the output of a file command\n", tmp_CheckGAU_File_Name);
    mqc_error( charBuf );
  }
  rtn_fgets=fgets(charBuf,1024,CheckGAU_File);
  fclose_return = fclose(CheckGAU_File);
  CheckGAU_File = NULL;
  if ( -1 == mqc_FINDstrINline( " ASCII ", charBuf ) ) {
    /* The file is Binary, not ACSII. It is a MatrixFile */
    printf( "%s is binary, so a MatrixFile\n",FileName);
    unlink(tmp_CheckGAU_File_Name);
    MatFileName = mqc_DupString(FileName, (char *)NULL);
    return( MatFileName );	     
  }

/* Input is a Gaussian Input file */
#ifdef DEBUG
  printf( "Check if the second argument is the name of an executable which is in the path\n" );
  fflush(stdout);
#endif

  sprintf( tmp_CheckGAU_File_NameA, "%sa%", tmp_CheckGAU_File_Name );
  sprintf( charBuf, "#! /bin/bash -f\nif command -v %s > /dev/null 2>&1; then\n echo 1 &> %s\nelse\n echo 2 &> %s\nfi\n", Program, tmp_CheckGAU_File_NameA, tmp_CheckGAU_File_NameA );

  CheckGAU_File = fopen( tmp_CheckGAU_File_Name,"w");
  rtn_fputs=fputs(charBuf,CheckGAU_File);
  fclose_return = fclose(CheckGAU_File);
  CheckGAU_File = NULL;
  sprintf( charBuf, "chmod 777 %s;%s\n", tmp_CheckGAU_File_Name, tmp_CheckGAU_File_Name );

#ifdef DEBUG
  printf( "execute %s\n", charBuf );
  fflush(stdout);
#endif
  system_failure = system(charBuf);
  unlink(tmp_CheckGAU_File_Name);
  if ( system_failure == -1 ) {
    sprintf(charBuf, "Failure in Script to determine if \"%s\" is an executable\n", Program );
    unlink(tmp_CheckGAU_File_Name);
    mqc_error( charBuf );
  }

  CheckGAU_File = fopen( tmp_CheckGAU_File_NameA,"r");
  rtn_fgets=fgets(charBuf,1024,CheckGAU_File);
  fclose_return = fclose(CheckGAU_File);
  CheckGAU_File = NULL;
  unlink(tmp_CheckGAU_File_NameA);

  if ( strncmp( charBuf, "1", 1) != 0 ) {
    sprintf( charBuf, "Did not find %s in the path\n", Program );
    mqc_error( charBuf );
  }

#ifdef DEBUG
  printf( "Open Gaussian input file\n" );
  fflush(stdout);
#endif
  GAU_File = fopen( FileName,"r");
  if ( GAU_File == NULL ) {
    sprintf( charBuf, "Could not open Gaussian Input File %s\n", FileName);
    mqc_error( charBuf );
  }
  last_job=-1;
  rtn_fgets=fgets(charBuf,1024,GAU_File);
  if ( rtn_fgets!= NULL ) {
#ifdef DEBUG
    printf( " Read from input: %s\n", charBuf );
    fflush(stdout);
#endif
    Current_Line = Alloc_InputLine( charBuf, (InputLine *)NULL );
    Start_Line = Current_Line;
  }

  while( rtn_fgets!= NULL ) {
    if ( Current_Line->Keyword == 1 ) {
      Restart=0;
      Job_Number++;
      Current_job=0;
      Total_Jobs++;
      while ( Current_Line->Blank == 0 && rtn_fgets!= NULL ) {
	if (mqc_FINDstrINline( "restart", Current_Line->Lower_Name ) != -1 ) {
	  Restart=1;
	}
	found = 0;
	i = mqc_FINDstrINline( "output=", Current_Line->Lower_Name );
#ifdef DEBUG
	if ( i != -1 ) {
	  printf( "output= found\n" );
	  fflush(stdout);
	} else {
	  printf( "output= not found\n" );
	  fflush(stdout);
	}
#endif
	if ( i != -1 ) {
	  strcpy( lowerBuf, Current_Line->Lower_Name );
	  strcpy( charBuf, &lowerBuf[i] );
	  i = mqc_FINDstrINline( " ", charBuf );
	  if ( i != -1 ) {
	    strcpy( &charBuf[i], "\0" );
	  }
#ifdef DEBUG
	  printf( "output instructions - \"%s\"\n Mat detection %d\n Raw detection %d\n", charBuf, mqc_FINDstrINline( "mat", charBuf ), mqc_FINDstrINline( "raw", charBuf ));
	  fflush(stdout);
#endif
	  Current_Line->Keyword = 1;
	  if ( mqc_FINDstrINline( "mat", Current_Line->Lower_Name ) != -1 || 
	       mqc_FINDstrINline( "raw", Current_Line->Lower_Name ) != -1 ) {
	    found = 1;
	  }
	}

	external_proc = mqc_DupString( Current_Line->Lower_Name, (char *)NULL );
	external_proc = mqc_Replace_String( ",external", ",", external_proc );
	external_proc = mqc_Replace_String( "(external", "(", external_proc );
	if ( mqc_FINDstrINline( "external", external_proc ) != -1 ){
	  external = 1;
	} else {
	  external = 0;
	}
#ifdef DEBUG
	printf( "found is - %d\n", found );
	printf( "external is - %d %s\n", external, external_proc );
	fflush(stdout);
#endif
	free( external_proc );

	if ( found == 1 || external == 1 ){
	  if ( Current_job == 0 ) {
	    MatrixFile_Jobs++;
	  } else {
	    printf( "WARNING: Detected multiple triggers to write MatrixFile in job %d\n", Job_Number );
	  }
	  Current_job=1;
	}
	/* Keep reading keywords until reach a blank line */
	rtn_fgets=fgets(charBuf,1024,GAU_File);
	if ( rtn_fgets!= NULL ) {
#ifdef DEBUG
	  printf( " Read from input: %s\n", charBuf );
	  fflush(stdout);
#endif
	  Last_Line = Current_Line;
	  Current_Line = Alloc_InputLine( charBuf, Last_Line );
	}
	/* Find the name of the MartixFile in the keywords */
      }
      if ( Current_job == 0 && Restart == 1 ) {
	if ( last_job == -1 ) {
	  /* This is a restart job, and the keywords for the job that has */
	  /* been restarted my trigger a write of a MatrixFile.  It's not */
	  /* possible to tell. */
	  printf( "WARNING: Whether to write the MatrixFile in job %d is determined by the keywords in the job that has been restarted.\n", Job_Number);

	} else {
	  /* This is a restart job, and the keywords from the last job in */
	  /* the iinput will determine whether we write the MatrixFile */
	  /* or not. */
	  if ( last_job != 0 ) {
	    MatrixFile_Jobs++;
	  }
	  Current_job=last_job;
	}
      } else {
	if ( last_job != -1 || Current_job != 0 ) {
	  last_job=Current_job;
	}
      }
    }
    if ( rtn_fgets!= NULL ) {
      rtn_fgets=fgets(charBuf,1024,GAU_File);
      if ( rtn_fgets!= NULL ) {
#ifdef DEBUG
	printf( " Read from input: %s\n", charBuf );
	fflush(stdout);
#endif
	Last_Line = Current_Line;
	Current_Line = Alloc_InputLine( charBuf, Last_Line );
      }
    }
  }
  fclose_return = fclose(GAU_File);
  GAU_File = NULL;
  if ( fclose_return != 0 ) {
    sprintf( charBuf, "Could not close Gaussian Input File %s\n", FileName);
    mqc_error( charBuf );
  }
#ifdef DEBUG
  printf( "%d inputs of which %d generate Matrix files\n", 
	   Total_Jobs, MatrixFile_Jobs );
  printf( " Current_job %d, last_job %d\n", Current_job, last_job );
#endif
  NoModNeeded = 1;
  if ( Current_job == 0 && last_job == -1 ) {
    if ( Restart != 1 ) {
#ifdef DEBUG
      printf( "WARNING: This input does not write a MartrixFile.\n");
#endif
      NoModNeeded = 0;
    }
  } else if ( last_job == 0 ) {
#ifdef DEBUG
    printf( "WARNING: This input writes a MartrixFile, but the last job does not write the MatrixFile.\n");
#endif
    NoModNeeded = 0;
  }
#ifdef DEBUG
  printf( " NoModNeeded pass 1 %d\n", NoModNeeded );
  fflush(stdout);
#endif
  if ( NoModNeeded == 0 ) {
    mqc_Remove_Return( Last_Keyword->Name );
    sprintf( charBuf, "%s output=mat\n", Last_Keyword->Name );
    Last_Keyword->Name = mqc_DupString(charBuf, Last_Keyword->Name);
  }

  strcpy( JobName, FileName );
  ilen = strlen( JobName );
  if ( strncmp( &JobName[ilen-4], ".com", 4 ) == 0 ) {
    strcpy( &JobName[ilen-4], "\0" );
  }
  Last_Line = (InputLine *)Last_NonBlank->last;
  if ( Last_Line->Blank == 0
       || mqc_FINDstrINline( " ", Last_Line->Name ) != -1 ){
    NoModNeeded = 0;
  }
#ifdef DEBUG
  printf( "before nonblank \"%s\"\n",   Last_Line->Name );
  printf( "Blank %d\n", Last_Line->Blank );
  printf( "blank \"%s\"\n", Last_NonBlank->Name );

  printf( " NoModNeeded pass 2 %d\n", NoModNeeded );
  fflush(stdout);
#endif
  if ( NoModNeeded == 0 ) {
    /* Append a blank line, MatrixFile Name, blank line behind the */
    /* last non-blank line */
    sprintf( charBuf, "\n%s.mat\n\n", JobName );
    Current_Line = Alloc_InputLine( charBuf, Last_NonBlank );
    sprintf( charBuf, "%s.mat", JobName );
    MatFileName = mqc_DupString( charBuf, (char *)NULL);
    sprintf( InputFile, "%s_mod", FileName );
    mqc_Write_InputLines( InputFile, Start_Line );
  } else {
    sprintf( InputFile, FileName );
    MatFileName = mqc_DupString( Last_NonBlank->Name, (char *)NULL);
  }

  sprintf( charBuf, "%s < %s >%s.log", Program, InputFile, JobName );

#ifdef DEBUG
  printf( "execute %s\n", charBuf );
  fflush(stdout);
#endif
  system_failure = system(charBuf);
  if ( system_failure == -1 ) {
    printf( " Job failed. Command was:\n%s\n", charBuf );
  }

  mqc_Remove_Return (MatFileName);
  sprintf( charBuf, "file %s > %s\n", 
	   MatFileName, tmp_CheckGAU_File_Name );
#ifdef DEBUG
  printf( "Execute %s\n", charBuf );
  fflush(stdout);
#endif
  system_failure = system(charBuf);

  if ( system_failure == -1 ) {
    /* the file on the First input failed.  Check if the file exists. */
    sprintf( charBuf, "Problem with MatrixFile.  It may not exist \"file %s > %s\n failed\n", 
	   MatFileName, tmp_CheckGAU_File_Name );
    mqc_error( charBuf );
  }

  CheckGAU_File = fopen( tmp_CheckGAU_File_Name,"r");
  if ( CheckGAU_File == NULL ) {
    sprintf( charBuf, "Could not open CheckGAU file %s, should contain the output of a file command\n", tmp_CheckGAU_File_Name);
    mqc_error( charBuf );
  }
  rtn_fgets=fgets(charBuf,1024,CheckGAU_File);
  fclose_return = fclose(CheckGAU_File);
  unlink(tmp_CheckGAU_File_Name);
  CheckGAU_File = NULL;
  if ( -1 == mqc_FINDstrINline( " data", charBuf ) ) {
    /* The file is Binary, not ACSII. It is a MatrixFile */
#ifdef DEBUG
    printf( "result of file command %s\n", charBuf );
    fflush(stdout);
#endif
    sprintf( charBuf, "%s is not binary, so not a MatrixFile\n",MatFileName);
    mqc_error( charBuf );
  }

  return( MatFileName );
}

void mqc_error( char *Message )
{
  printf( "ERROR: %s\n", Message );
  fflush(stdout);
  if (CheckGAU_File != NULL ) fclose(CheckGAU_File);
  if (GAU_File != NULL ) fclose(GAU_File);
  if (ProfGAU_File != NULL ) fclose(ProfGAU_File);
  exit(EXIT_FAILURE);
}

void mqc_Remove_Init_WhiteSpace( char *charBuf )
{
  while( strncmp( &charBuf[0], " ", 1) == 0 ||
	 strncmp( &charBuf[0], Tab, 1) == 0 ) 
	 strcpy( charBuf, &charBuf[1] );
  while( strncmp( &charBuf[0], "\n", 1) == 0 ) 
    strcpy( charBuf, &charBuf[1] );
  while( strncmp( charBuf, CarriageReturn, 1) == 0 ) 
    strcpy( charBuf, &charBuf[1] );

  return;
}

void mqc_Remove_Return( char *charBuf )
{
  int ilen;

  ilen = strlen(charBuf);
  while( strncmp( &charBuf[ilen-1], "\n", 1) == 0 ||
	 strncmp( &charBuf[ilen-1], CarriageReturn, 1) == 0 ) {
    strcpy( &charBuf[ilen-1], "\0" );
    ilen = strlen(charBuf);
  }

  return;
}

void mqc_LOWER( char *charBuf )
{
  int i, len;
  int ucase;

  len = strlen( charBuf );
  for ( i=0; i < len; i++ ) {
    if( isupper( (int)charBuf[i] ) ) {
      ucase = (int)charBuf[i];
      charBuf[i] = (char) (ucase+int_tolower);
    }
  }
  return;
}

int mqc_FINDstrINline( char *string, char *line )
{
  int i,imax;
  int stringlen;

  stringlen = (int) strlen( string );
  if ( stringlen == 0 ) return(-1);
  imax = (int) strlen( line ) - stringlen + 1;
  for( i=0 ;  i < imax ; i++ ) {
    if ( 0==strncmp ( &line[i], string, (size_t)stringlen ) ) {
      return(i);
    }
  }
  return(-1);
}

void mqc_Global_init(void)
{
  char charBuf[1024];

  CarriageReturn[0]=(char)13;
  CarriageReturn[1]=(char)0;
  strcpy( charBuf, "aA");
  int_tolower = (int)charBuf[0] - (int)charBuf[1];
}

char *mqc_DupString( char *old_str, char *free_this )
{
  char *new_str=NULL;
  int ilen;

  ilen = strlen( old_str );
  new_str = malloc( ilen+8 );
  strcpy( new_str, old_str );
  if (free_this != (char *)NULL) {
    free(free_this);
  }
  return( new_str );
}

InputLine *Alloc_InputLine( char *InputLine_name, InputLine *LAST_InputLine_struct )
{
  static Blank_count = 0;

  InputLine *new_struct = (InputLine *)NULL;

  new_struct = (InputLine *)malloc ( (size_t)sizeof( InputLine )+128);
  if ( new_struct == (InputLine *)NULL ) {
    printf( "Error: Not able to allocate memory in Alloc_InputLine.");
    exit( 0 );
  }
  if ( LAST_InputLine_struct == (InputLine *)NULL ) {
    new_struct->last = (char *)NULL;
    new_struct->next = (char *)NULL;
  } else {
    new_struct->next = LAST_InputLine_struct->next;
    LAST_InputLine_struct->next = (char *)new_struct;
    new_struct->last = (char *)LAST_InputLine_struct;
  }
  new_struct->Name = mqc_DupString(InputLine_name, (char *)NULL) ;
  mqc_Remove_Init_WhiteSpace( InputLine_name );
  new_struct->Lower_Name = mqc_DupString(InputLine_name, (char *)NULL) ;
  mqc_LOWER( new_struct->Lower_Name );

  new_struct->Keyword = 0;
  if ( strlen( new_struct->Lower_Name ) == 0 ) {
    /* when this is really smart, we will know what the input is by the */
    /* location between blank lines and the keywords. */
    Blank_count++;
    new_struct->Blank = Blank_count;
  } else {
    new_struct->Blank = 0;
    new_struct->Link1 = 0;
    Last_NonBlank = new_struct;
    /* The very last line in the input is the name of the MartixFile, */
    /* if we are writing one.  If we aren't writing one, need to append a */
    /* blank line and then the MartixFile name afterwards. */
    if ( strncmp( new_struct->Lower_Name, "#", 1 ) == 0 ) {
      new_struct->Keyword = 1;
      Last_Keyword = new_struct;
      /* We may need to append a keyword to cause a write of a Martrix file */
    } else if ( mqc_FINDstrINline( "--link1--", new_struct->Lower_Name ) != -1 ) {
      /* New job */
      Blank_count=0;
      new_struct->Link1 = 1;
    }
  }

  return(new_struct);
}

void mqc_Write_InputLines( char *File_Name, InputLine *Start_InputLine )
{
  FILE *MY_GAU_File=NULL;
  int fclose_return;
  int rtn_fputs;

  InputLine *Current_InputLine = (InputLine *)NULL;
  InputLine *Next_InputLine = (InputLine *)NULL;

  MY_GAU_File = fopen( File_Name,"w");
  Current_InputLine = Start_InputLine;
  while (Current_InputLine != (InputLine *)NULL ) {
    Next_InputLine = (InputLine *)Current_InputLine->next;
    rtn_fputs=fputs( Current_InputLine->Name,MY_GAU_File);
    Current_InputLine = Next_InputLine;
  }
  fclose_return = fclose(MY_GAU_File);
}

char *mqc_Replace_String( char *Find, char *Insert, char *Line )
{
  int Flen, Ilen;
  int i;
  char *Return_Line;

  i = mqc_FINDstrINline( Find, Line );
  if ( i != -1 ) {
    Ilen = strlen( Insert );
    Return_Line = malloc( strlen( Line ) + Ilen + 80 );
    strncpy( Return_Line, Line, i );
    strncpy( &Return_Line[i], Insert, Ilen );
    strcpy( &Return_Line[i+Ilen], &Line[i+strlen(Find)] );
    free( Line );
    return( mqc_Replace_String( Find, Insert, Return_Line ) );
  }
  return( Line );
}
