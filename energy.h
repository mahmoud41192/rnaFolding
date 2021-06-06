#include <fstream>
#include <stdio.h>

std::string structsFile = "structs.dat";
std::string energiesFile = "energies.dat";

/*THIS FUNCTION WRITES A LIST OF SECONDARY STRUCTURES TO A 
 FILE THAT WILL BE USED TO CALCULATE THE FREE ENEERGY FOR 
 THEM USING DAVE MATHEWS PROGRAM*/
void writeStructToFile(rnaMolecule molecule) {
    std::ofstream file; 
    file.open(structsFile,std::ofstream::out | std::ofstream::trunc); 
    //std::cout << "entered function" << std::endl; 
    file << ">RNA-1" << std::endl; 
    file << molecule.sequence << std::endl;
    //std::cout << molecule.structure << std::endl;
    file << molecule.structure << std::endl;

}


/*COMPUTES THE FREE ENERGY FROM THE FILE AND GETS*/
void computeFreeEnergy(rnaMolecule& molecule) { 
     
    std::string command = "./efn2 " + structsFile + " " + energiesFile; 
    
    //USE THE PROGRAM TO COMPUTE FREE ENERGY
    std::system(command.c_str()); 
    
    /*WE USE C I/O RATHER THAT C++ INTERFACE BECAUSE THE 
     C STYLE STREAMS PROVIDE MUCH BETTER INTERFACING FUNCTIONS*/
    FILE *file = fopen(energiesFile.c_str(),"r"); 
    char buffer[20]; 
    int buffer2;
    int index = 0; 
    double energy;
    while(!feof(file)) {
        fscanf(file,"%s %d %s %s %lf %s %s\n",buffer,&buffer2,buffer,buffer,&energy
        ,buffer,buffer);
    }
    //std::cout << energy << std::endl; 
    molecule.freeEnergy = energy;
    
    fclose(file); //CLOSE FILE 
}

