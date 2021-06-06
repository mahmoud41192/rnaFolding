#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cstdlib>
#include "rna.h"
#include "ligand.h"
#include "energy.h"
#include "structure.h"
#include "masterEquation.h"
#include "utility.h"


int main(int argc, char* argv[]) {
    
    /*========================================================================*/
                   /*READING INPUT AND PREPARE FOR THE COMPUTATION*/
    /*========================================================================*/
    
    std::string ligandsFile = "null";
    std::string outputFile = "output.dat";
    std::string inputFile = "null"; 
    std::string reactionsMapFile = "null";
    
    ligandList ligands;
    rnaList list;
    rnaReactionSpace reactions;
    
    double transcriptionalSpeed = 80;  // nt/sec
    double baseAdditionTime = 1/transcriptionalSpeed;
    
    if(argc < 5) {
        std::cout << "ERROR USAGE: ./compute -i [input file] -l [ligands file] -o [output file] -m [reactions matrix]" 
                << std::endl;
        return -1; 
    }
    
    for(int i = 1; i < argc; i++) {
        std::string option = argv[i]; 
        
        if(option == "-i") {
            inputFile = argv[++i]; 
            //std::cout << inputFile << std::endl;
        }
        else if(option == "-l") {
            ligandsFile = argv[++i];
            //std::cout << ligandsFile << std::endl;
        }
        else if(option == "-o") {
            outputFile = argv[++i];
            //std::cout << outputFile << std::endl;
        }
        else if(option == "-m") {
            reactionsMapFile = argv[++i];
        }
    }
    
    if(inputFile == "null") {
        std::cout << "ERROR: NO INPUT FILE" << std::endl;
        return -1; 
    }
    
    if(ligandsFile != "null") {
        std::cout << "ligands file " << ligandsFile << std::endl;
        ligands.fetchFromFile(ligandsFile);
    }
    
    //std::cout << ligands.size << std::endl;
     
    list.fetchFromFile(inputFile); 
    //list.showList(); 
    
    std::cout << "ligands: " << ligands.size << std::endl;
    std::cout << "initial molecules: " << list.size << std::endl; 
     
    reactions.numOfConfomers = list.size;
    //std::cout << list.size << std::endl;
    reactions.fetchFromFile(reactionsMapFile);
    //reactions.showReactions(); 
    
    
    /*WE HAVE TO ADD LIGANDS TO INPUT STRUCTURES*/
    //THIS IS FOR THE FIRST TIME TO PUT LIGANDS, PLEASE REMEMBER THAT
    for(int i = 0; i < list.size; i++) {
        rnaMolecule buffer = list.ptr[i];
        if(buffer.ligand == -1) {
            int ligandIndex = ligands.checkForLigandBinding(buffer);
            if(ligandIndex != -1) {
                std::cout << "found ligand" << std::endl;
                //WE PUSH ANOTHER ELEMENT WITH THE LIGAND BOUND
                rnaMolecule buffer2(buffer.sequence,buffer.structure);
                buffer2.ligand = ligandIndex;
                //WE SHOULD ADD THE BONUS FREE ENERGY DUE TO LIGAND BINDING
                buffer2.freeEnergy = buffer.freeEnergy+ligands.list[ligandIndex].bonusFreeEnergy;
                std::cout << buffer2.freeEnergy << std::endl;
                list.addConfomer(buffer2);
                reactions.addNewReaction(reactions.numOfConfomers-1,ligands.list[ligandIndex].bindingRateConstant);
            }
        }
    }
    
    //list.showList();
    
    //READ THE SEQUENCE WHICH WE COMPUTE FOR FROM THE USER
    std::string parsingSequence;
    std::cout << "ENTER THE SEQUENCE TO COMPUTE FOLDING FOR : " << std::endl;
    std::cin >> parsingSequence;
     
    
    //std::cout << parsingSequence << std::endl;
    
    
    
    
    /*========================================================================*/
                                /*COMPUTATION*/
    /*========================================================================*/
    
    
    for(int i = 0; i < parsingSequence.length(); i++) {
        
        //WE NEED TO SAVE ALL DATA FROM ALL TRANSCRIPTION
        //STAGES TO AN OUTPUT FILE
        //writeListToFile("data.dat",list);  
        
        //FIRST: APPEND THE NEW NUCLEOTIDE
        list.appendNucleotideToList(parsingSequence[i]);
        //list.showList()
        
        //SECOND: UPDATE THE FREE ENERGIES
        for(int j = 0; j < list.size; j++) {
            writeStructToFile(*(list.ptr+j)); 
            computeFreeEnergy(*(list.ptr+j));
            if(list.ptr[j].ligand != -1) {
                list.ptr[j].freeEnergy += ligands.list[list.ptr[j].ligand].bonusFreeEnergy;
                //std::cout << list.ptr[i].freeEnergy << std::endl;
            }
        }
        
        //THIRD: PARSE FOR NEW STRUCTURES
        while(findNewStructures(list,reactions,ligands)) {}
        
        Eigen::MatrixXd reactionsMatrix = getTransferMatrix(list,reactions);
        Eigen::MatrixXcd populations = getPopulations(list);
        
        //FOURTH: COMPUTE THE FOLDING DYNAMICS
        updatePopulation(reactionsMatrix,populations,baseAdditionTime);
        
        //FIFTH: WE PASTE THE CORRECTIONS TO POPULATIONS TO THE LIST
        for(int i = 0; i < list.size; i++) {
            (list.ptr+i)->population = populations(i,0).real(); 
        }
        
        std::cout << populations.sum() << std::endl; 
        std::cout << "step " << i << std::endl; 
        writeListToFile(outputFile,list,reactions);
    }
    
    
    
    return 0; 
}
