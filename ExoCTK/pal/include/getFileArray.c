/* This file is part of Exo_Transmit.

    Exo_Transmit is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Exo_Transmit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Exo_Transmit.  If not, see <http://www.gnu.org/licenses/>.
*/


/*-------------------------- getFileArray.c ----------------------

Author: Patrick Slough (sloughpa@grinnell.edu)

------------------------------------------------------------------ */

/* This routine reads in data from the various input files and 
   stores the data in a relevant data structure which will be 
   proccessed in other parts of the program.

------------------------------------------------------------------ */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

char** getFileArray(){
  
  
  /* variables*/
  char pwd[256], extension[35][256], **fileArray;
  
  int i, j, lineCount; // array and line counters
  
  /* filenames */
  char input_filename [30] = "userInput.in"; // Input file Name for the program
  char opac_filename [30] = "otherInput.in"; // Input file with opacity data
  
  FILE *input; //File pointer varible
  
  
  /* Get data from input files */
  
  int contentLines[3] = {6,8,10};  // 6,8,10 are data containing lines in 'userInput.in'
  
  input = fopen(input_filename, "r"); //open stream
  
  if(input){
    
    /* initialize counters */
    i = 0; 
    lineCount = 1;
    
    char line[256]; // stores line data
    
    /* read in data line by line until end of file */
    while(fgets(line, sizeof(line), input) != NULL){
      /* check for the line with 'pwd' and store 'pwd' as string in 'line' */
      if(lineCount == 4){
        line[strlen(line) - 1] = '\0';
        strcpy(pwd, line);
        lineCount++;
      } /* if line is not pwd and contains data*/ 
      else if ((i < 3) && (lineCount == contentLines[i])) {
        line[strlen(line) - 1] = '\0';
        strcpy(extension[i], line);
        i++;
        lineCount++;
      } /* if line is not pwd and does not contain data*/
      else{
        lineCount++;
      }
    }
  } /* if input file does not exist */
  else {
    printf("Error reading input file\n");
    exit(1);
  }
  
  fclose(input);  // close stream
  
  /* Get data from 'otherInput.in' */
  
  int opac_contentLines[31] = {4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35};  // data containing fields in 'otherInput.in'
  
  input = fopen(opac_filename, "r");   //open stream
  
  if(input){
    
    // initialize counters
    lineCount = 1;
    int k = 0; // contentLine counter
    
    char line[256]; // stores line data
    
    /*read in data line by line till end of file*/
    while(fgets(line, sizeof(line), input) != NULL){
      /* if line contains data */
      if ((k < 32) && (lineCount == opac_contentLines[k])) {
        line[strlen(line) - 1] = '\0';
        strcpy(extension[i], line);
        i++;
        lineCount++;
        k++;
      } /* if line is not pwd and does not contain data*/
      else{
        lineCount++;
      }
    }
  } /* if input file does not exist */
  else {
    printf("Error reading input file\n");
    exit(1);
  }

  fclose(input);  // close stream
  
  fileArray = malloc(35*sizeof(char*)); //allocate the 35 string "slots"
  for ( j = 0; j < 35; j++){
    fileArray[j] = (char*)malloc(strlen(pwd)+strlen(extension[j])+1);
    strcpy(fileArray[j], pwd);
    strcat(fileArray[j], extension[j]);
    printf("%s\n", fileArray[j]);
  }

  return fileArray; //return the array 
  
}


/* -------- end ----------------- GetFileArray ------------------- */
