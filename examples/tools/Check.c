/* C Example */
#include <math.h>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
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
char *mqc_GetMatrixFile( char *, char *, int );
void mqc_Global_init( void );
char *mqc_DupString( char *, char * );
InputLine *mqc_Alloc_InputLine( char *, InputLine *, int );
void mqc_Write_InputLines( char *, InputLine * );
void mqc_Free_InputLines( InputLine * );
void mqc_Remove_Return( char * );
char *mqc_Replace_String( char *, char *, char * );
void mqc_get_MatrixFile_Names(char *, char *, int);
void mqc_error_i_c2f_0( char *, int * ); 
void print_line_c2f ( char *, int * ); 
void flush_c2f( int * );
char CarriageReturn[8];
char Tab[8];
int int_tolower;
FILE *GAU_File=NULL;
FILE *CheckGAU_File=NULL;
FILE *ProfGAU_File=NULL;
InputLine *Last_Keyword = (InputLine *)NULL;
InputLine *Last_NonBlank = (InputLine *)NULL;

void mqc_get_MatrixFile_Names(char *FileName, char *Program, int iout)
{
  static int mqc_first_call=0;
  char *MatFileName;

  if ( mqc_first_call == 0 ){
    mqc_Global_init( );
    mqc_first_call=1;
  }
  MatFileName = mqc_GetMatrixFile( FileName, Program, iout);
}

char *mqc_GetMatrixFile( char *FileName, char *Program, int iout )
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
  char prnt_str[2048];
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
  print_line_c2f ( "Open file to store temporaries.  This should be temporary\0", &iout ); 
  flush_c2f( &iout );
#endif
  PID = (int)getpid();
  scrdir=getenv( "GAUSS_SCRDIR" );
  if ( scrdir == NULL ) {
    sprintf( tmp_CheckGAU_File_Name, "./mqc_pass_File_%d", PID);
  } else {
    sprintf( tmp_CheckGAU_File_Name, "%s/mqc_pass_File_%d", scrdir, PID);
  }
#ifdef DEBUG
  sprintf( prnt_str, " tmp_CheckGAU_File_Name %s\0", tmp_CheckGAU_File_Name );
  print_line_c2f ( prnt_str, &iout ); 
  flush_c2f( &iout );
#endif
  CheckGAU_File = fopen( tmp_CheckGAU_File_Name,"w");
  if ( CheckGAU_File == NULL ) {
    sprintf( charBuf, "Could not open CheckGAU file %s\0", tmp_CheckGAU_File_Name);
    mqc_error_i_c2f_0( charBuf, &iout); 
  }
  fclose_return = fclose(CheckGAU_File);
  CheckGAU_File = NULL;
  if ( fclose_return != 0 ) {
    sprintf( charBuf, "Could not close CheckGAU File %s\0", tmp_CheckGAU_File_Name);
    mqc_error_i_c2f_0( charBuf, &iout); 
  }
  unlink(tmp_CheckGAU_File_Name);
  CheckGAU_File = NULL;
#ifdef DEBUG
  print_line_c2f ( "Check if the first argument is the name of a binary MatrixFile or if it is an ASCII Gaussian input file\0", &iout ); 
  flush_c2f( &iout );
#endif

  sprintf( charBuf, "file %s > %s\n", 
	   FileName, tmp_CheckGAU_File_Name );
#ifdef DEBUG
  sprintf( prnt_str, "execute %s\0", charBuf );
  print_line_c2f ( prnt_str, &iout ); 
  flush_c2f( &iout );
#endif
  system_failure = system(charBuf);

  if ( system_failure == -1 ) {
    /* the file on the First input failed.  Check if the file exists. */
    GAU_File = fopen( FileName,"r");
    if ( GAU_File == NULL ) {
      sprintf( charBuf, "Could not open Gaussian File specified by the first argument %s\0", FileName);
    } else {
      sprintf( charBuf, "file command on %s failed.  Tried to do this in the Bash shell\0", FileName);
    }
    CheckGAU_File = fopen( tmp_CheckGAU_File_Name,"w");
    unlink(tmp_CheckGAU_File_Name);
    mqc_error_i_c2f_0( charBuf, &iout); 
  }

  CheckGAU_File = fopen( tmp_CheckGAU_File_Name,"r");
  if ( CheckGAU_File == NULL ) {
    sprintf( charBuf, "Could not open CheckGAU file %s, should contain the output of a file command\0", tmp_CheckGAU_File_Name);
    mqc_error_i_c2f_0( charBuf, &iout); 
  }
  rtn_fgets=fgets(charBuf,1024,CheckGAU_File);
  fclose_return = fclose(CheckGAU_File);
  CheckGAU_File = NULL;
  if ( -1 == mqc_FINDstrINline( " ASCII ", charBuf ) ) {
    /* The file is Binary, not ACSII. It is a MatrixFile */
#ifdef DEBUG
    sprintf( prnt_str, "%s is binary, so a MatrixFile\0",FileName);
    print_line_c2f ( prnt_str, &iout ); 
    flush_c2f( &iout );
#endif
    unlink(tmp_CheckGAU_File_Name);
    MatFileName = mqc_DupString(FileName, (char *)NULL);
    return( MatFileName );	     
  }

/* Input is a Gaussian Input file */
#ifdef DEBUG
  print_line_c2f ( "Check if the second argument is the name of an executable which is in the path\0", &iout ); 
  flush_c2f( &iout );
#endif

  sprintf( tmp_CheckGAU_File_NameA, "%sa", tmp_CheckGAU_File_Name );
  sprintf( charBuf, "#! /bin/bash -f\nif command -v %s > mqc_tmp_file 2>&1; then\n echo 1 > %s\nelse\n echo 2 > %s\nfi\n", Program, tmp_CheckGAU_File_NameA, tmp_CheckGAU_File_NameA );

  CheckGAU_File = fopen( tmp_CheckGAU_File_Name,"w");
  rtn_fputs=fputs(charBuf,CheckGAU_File);
  fclose_return = fclose(CheckGAU_File);
  CheckGAU_File = NULL;
  sprintf( charBuf, "chmod 755 %s;%s\0", tmp_CheckGAU_File_Name, tmp_CheckGAU_File_Name );

#ifdef DEBUG
  sprintf( prnt_str, "execute %s\0", charBuf );
  print_line_c2f ( prnt_str, &iout ); 
  flush_c2f( &iout );
#endif
  system_failure = system(charBuf);
  unlink("mqc_tmp_file");
  unlink(tmp_CheckGAU_File_Name);
  if ( system_failure == -1 ) {
    sprintf(charBuf, "Failure in Script to determine if \"%s\" is an executable\0", Program );
    mqc_error_i_c2f_0( charBuf, &iout); 
  }

  CheckGAU_File = fopen( tmp_CheckGAU_File_NameA,"r");
  rtn_fgets=fgets(charBuf,1024,CheckGAU_File);
  fclose_return = fclose(CheckGAU_File);
  CheckGAU_File = NULL;
  unlink(tmp_CheckGAU_File_NameA);
  if ( strncmp( charBuf, "1", 1) != 0 ) {
    sprintf( charBuf, "Did not find %s in the path\0", Program );

printf( "%s\n", charBuf );

    mqc_error_i_c2f_0( charBuf, &iout); 
  }
#ifdef DEBUG
  print_line_c2f ( "Open Gaussian input file\0", &iout ); 
  flush_c2f( &iout );
#endif
  GAU_File = fopen( FileName,"r");
  if ( GAU_File == NULL ) {
    sprintf( charBuf, "Could not open Gaussian Input File %s\0", FileName);
    mqc_error_i_c2f_0( charBuf, &iout);
  }
  last_job=-1;
  rtn_fgets=fgets(charBuf,1024,GAU_File);
  if ( rtn_fgets!= NULL ) {
#ifdef DEBUG
    sprintf( prnt_str, " Read from input: %s\0", charBuf );
    print_line_c2f ( prnt_str, &iout ); 
    flush_c2f( &iout );
#endif
    Current_Line = mqc_Alloc_InputLine( charBuf, (InputLine *)NULL, iout );
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
	  print_line_c2f ( "output= found\0", &iout ); 
	} else {
	  print_line_c2f ( "output= not found\0", &iout ); 
	}
	flush_c2f( &iout );
#endif
	if ( i != -1 ) {
	  strcpy( lowerBuf, Current_Line->Lower_Name );
	  strcpy( charBuf, &lowerBuf[i] );
	  i = mqc_FINDstrINline( " ", charBuf );
	  if ( i != -1 ) {
	    strcpy( &charBuf[i], "\0" );
	  }
#ifdef DEBUG
	  sprintf( prnt_str, "output instructions - \"%s\"\n Mat detection %d\n Raw detection %d\0", charBuf, mqc_FINDstrINline( "mat", charBuf ), mqc_FINDstrINline( "raw", charBuf ));
	  print_line_c2f ( prnt_str, &iout ); 
	  flush_c2f( &iout );
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
	sprintf( prnt_str, "found is - %d\0", found );
	print_line_c2f ( prnt_str, &iout ); 
	sprintf( prnt_str, "external is - %d %s\0", external, external_proc );
	print_line_c2f ( prnt_str, &iout ); 
	flush_c2f( &iout );
#endif
	free( external_proc );

	if ( found == 1 || external == 1 ){
	  if ( Current_job == 0 ) {
	    MatrixFile_Jobs++;
	  } else {
	    sprintf( prnt_str, "WARNING: Detected multiple triggers to write MatrixFile in job %d\0", Job_Number );
	    print_line_c2f ( prnt_str, &iout ); 
	  }
	  Current_job=1;
	}
	/* Keep reading keywords until reach a blank line */
	rtn_fgets=fgets(charBuf,1024,GAU_File);
	if ( rtn_fgets!= NULL ) {
#ifdef DEBUG
	  sprintf( prnt_str, " Read from input: %s\0", charBuf );
	  print_line_c2f ( prnt_str, &iout ); 
	  flush_c2f( &iout );
#endif
	  Last_Line = Current_Line;
	  Current_Line = mqc_Alloc_InputLine( charBuf, Last_Line, iout );
	}
	/* Find the name of the MartixFile in the keywords */
      }
      if ( Current_job == 0 && Restart == 1 ) {
	if ( last_job == -1 ) {
	  /* This is a restart job, and the keywords for the job that has */
	  /* been restarted my trigger a write of a MatrixFile.  It's not */
	  /* possible to tell. */

	  sprintf( prnt_str, "WARNING: Whether to write the MatrixFile in job %d is determined by the keywords in the job that has been restarted.\0", Job_Number);
	  print_line_c2f ( prnt_str, &iout ); 

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
	sprintf( prnt_str, " Read from input: %s\0", charBuf );
	print_line_c2f ( prnt_str, &iout ); 
	flush_c2f( &iout );
#endif
	Last_Line = Current_Line;
	Current_Line = mqc_Alloc_InputLine( charBuf, Last_Line, iout );
      }
    }
  }
  fclose_return = fclose(GAU_File);
  GAU_File = NULL;
  if ( fclose_return != 0 ) {
    sprintf( charBuf, "Could not close Gaussian Input File %s\0", FileName);
    mqc_error_i_c2f_0( charBuf, &iout);
  }
#ifdef DEBUG
  sprintf( prnt_str, "%d inputs of which %d generate Matrix files\0", 
	   Total_Jobs, MatrixFile_Jobs );
  print_line_c2f ( prnt_str, &iout ); 
  sprintf( prnt_str, " Current_job %d, last_job %d\0", Current_job, last_job );
  print_line_c2f ( prnt_str, &iout ); 
#endif
  NoModNeeded = 1;
  if ( Current_job == 0 && last_job == -1 ) {
    if ( Restart != 1 ) {
#ifdef DEBUG
      print_line_c2f ( "WARNING: This input does not write a MartrixFile.\0", &iout ); 
#endif
      NoModNeeded = 0;
    }
  } else if ( last_job == 0 ) {
#ifdef DEBUG
    print_line_c2f ( "WARNING: This input writes a MartrixFile, but the last job does not write the MatrixFile.\0", &iout );
#endif
    NoModNeeded = 0;
  }
#ifdef DEBUG
  sprintf( prnt_str, " NoModNeeded pass 1 %d\0", NoModNeeded );
  print_line_c2f ( prnt_str, &iout ); 
  flush_c2f( &iout );
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
  sprintf( prnt_str, "before nonblank \"%s\"\0",   Last_Line->Name );
  print_line_c2f ( prnt_str, &iout ); 
  sprintf( prnt_str, "Blank %d\0", Last_Line->Blank );
  print_line_c2f ( prnt_str, &iout ); 
  sprintf( prnt_str, "blank \"%s\"\0", Last_NonBlank->Name );
  print_line_c2f ( prnt_str, &iout ); 
  sprintf( prnt_str, " NoModNeeded pass 2 %d\0", NoModNeeded );
  print_line_c2f ( prnt_str, &iout ); 
  flush_c2f( &iout );
#endif
  if ( NoModNeeded == 0 ) {
    /* Append a blank line, MatrixFile Name, blank line behind the */
    /* last non-blank line */
    sprintf( charBuf, "\n%s.mat\n\n", JobName );
    Current_Line = mqc_Alloc_InputLine( charBuf, Last_NonBlank, iout );
    sprintf( charBuf, "%s.mat", JobName );
    MatFileName = mqc_DupString( charBuf, (char *)NULL);
    sprintf( InputFile, "%s.mod", JobName );
    mqc_Write_InputLines( InputFile, Start_Line );
  } else {
    strcpy( InputFile, FileName );
    MatFileName = mqc_DupString( Last_NonBlank->Name, (char *)NULL);
  }
  mqc_Free_InputLines( Start_Line );

  sprintf( charBuf, "%s < %s >%s.log", Program, InputFile, JobName );

#ifdef DEBUG
  sprintf( prnt_str, "execute %s\0", charBuf );
  print_line_c2f ( prnt_str, &iout ); 
  flush_c2f( &iout );
#endif
  system_failure = system(charBuf);
  if ( system_failure == -1 ) {
    sprintf( prnt_str, " Job failed. Command was:\n%s\0", charBuf );
    print_line_c2f ( prnt_str, &iout ); 
  }

  mqc_Remove_Return (MatFileName);
  sprintf( charBuf, "file %s > %s\n", 
	   MatFileName, tmp_CheckGAU_File_Name );
#ifdef DEBUG
  sprintf( prnt_str, "Execute %s\0", charBuf );
  print_line_c2f ( prnt_str, &iout ); 
  flush_c2f( &iout );
#endif
  system_failure = system(charBuf);

  if ( system_failure == -1 ) {
    /* the file on the First input failed.  Check if the file exists. */
    sprintf( charBuf, "Problem with MatrixFile.  It may not exist \"file %s > %s\n failed\0", 
	   MatFileName, tmp_CheckGAU_File_Name );
    mqc_error_i_c2f_0( charBuf, &iout);
  }

  CheckGAU_File = fopen( tmp_CheckGAU_File_Name,"r");
  if ( CheckGAU_File == NULL ) {
    sprintf( charBuf, "Could not open CheckGAU file %s, should contain the output of a file command\0", tmp_CheckGAU_File_Name);
    mqc_error_i_c2f_0( charBuf, &iout);
  }
  rtn_fgets=fgets(charBuf,1024,CheckGAU_File);
  fclose_return = fclose(CheckGAU_File);
  unlink(tmp_CheckGAU_File_Name);
  CheckGAU_File = NULL;
  if ( -1 == mqc_FINDstrINline( " data", charBuf ) ) {
    /* The file is Binary, not ACSII. It is a MatrixFile */
#ifdef DEBUG
    sprintf( prnt_str, "result of file command %s\0", charBuf );
    print_line_c2f ( prnt_str, &iout ); 
    flush_c2f( &iout );
#endif
    sprintf( charBuf, "%s is not binary, so not a MatrixFile\0",MatFileName);
    mqc_error_i_c2f_0( charBuf, &iout);
  }

  return( MatFileName );
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

InputLine *mqc_Alloc_InputLine( char *InputLine_name, InputLine *LAST_InputLine_struct, int iout )
{
  static int Blank_count = 0;

  InputLine *new_struct = (InputLine *)NULL;

  new_struct = (InputLine *)malloc ( (size_t)sizeof( InputLine )+128);
  if ( new_struct == (InputLine *)NULL ) {
    mqc_error_i_c2f_0("Error: Not able to allocate memory in mqc_Alloc_InputLine.", &iout);
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

void mqc_Free_InputLines( InputLine *Start_InputLine )
{
  InputLine *Current_InputLine = (InputLine *)NULL;
  InputLine *Next_InputLine = (InputLine *)NULL;

  Current_InputLine = Start_InputLine;
  while (Current_InputLine != (InputLine *)NULL ) {
    Next_InputLine = (InputLine *)Current_InputLine->next;
    /* Free the strings */
    free(Current_InputLine->Name);
    free(Current_InputLine->Lower_Name);
    /* Set values to NULL, just in case we accidently to use a */
    /* freed structure */
    Current_InputLine->Name = (char *)NULL;
    Current_InputLine->Lower_Name = (char *)NULL;
    Current_InputLine->next=NULL;
    Current_InputLine->last=NULL;
    /* Free the structure */
    free(Current_InputLine);
    Current_InputLine = Next_InputLine;
  }
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
