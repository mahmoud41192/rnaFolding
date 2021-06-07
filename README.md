# rnaFolding

## program description
this program simulates the folding of an RNA molecule as nucleotides are added along the process of transcription. 
the user inputs a file that contains the initial parameters of the ensemble of RNA molecules, and determines the 
number and type of nucleotides to be added with the option of ligand binding as well, then the program does the calculation
outputing a file containing all the information about the simulation. 

## program design
there are the main libraries which define the structure of RNA, and computational tools. <br/>
rna.h: <br/>
* this contains three classes: <br/>
  * the rnaMolecule class which contains the basic information about an RNA molecule and its structure, as well as 
** some of the basic operations that can be done to modify such molecule. <br/>
>>(2)the rnaList class which mediates how to handle a list of rnaMolecule objects, as well as how to read and fetch 
>>the structures from an input file. <br/>
>>(3)the rnaReactionSpace class which handles the interactions and transitions between different RNA confomers. <br/>
structure.h: <br/>
>contains basic functions that are used to find the possible new RNA confomers after adding a new nucleotide. <br/>
masterEquation.h: <br/>
>contains functions that are used to calculate the updated populations after adding a new nucleotides. <br/>
lignad.h: <br/>
>this contains two classes: <br/>
>>(1) ligand which includes the basic information about such ligand and its binding properties. <br/>
>>(2) ligandList which handles an input list of ligands and checks for a particular RNA confomer 
>>if a ligand in the list can bind to it. <br/>
energy.h: <br/>
>this contains functions that uses a precompiled program efn2.exe to calculate the energy for each 
>RNA confomer in a certain RNA list. <br/>
utility.h <br/>
>this contains a function that writes the output of the calculation to a file. <br/>
main.cpp <br/>
>this is the main source code for the program containing the main steps of the calculation. 
