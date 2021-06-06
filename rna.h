#include <stdio.h>
#include<string>
#include "Eigen/Dense"


/*THIS IS THE MAXIMUM AMOUT WE CAN 
 *ALLOCATE ON THE CACHE OF THE PROCESSOR 
 *THIS MAKES A PROBLEM IN RUNNING OF 
 * THE PROGRAM*/
#define MAXCONFOMERS 1000
#define NUMPERADDITION 1
#define MAXMOLELENGTH 1000


/*NOW WE SHOULD MAKE A CLASSS TO HOLD AN RNA MOLECULE*/
class rnaMolecule {
public:
    std::string sequence;         
    std::string structure; 
    double freeEnergy;
    double population;
    bool notParsed; 
    //int rateLimitingReactions;
    int ligand; //HOLDS THE LIGAND INDEX IN THE LIGAND LIST AND -1 IF LIGAND IS NOT BOUND
    
    rnaMolecule(std::string seq,std::string dot) : sequence(seq), 
    structure(dot), population(0), 
    notParsed(true), ligand(-1), freeEnergy(0) {}
    
    rnaMolecule(): population(0), notParsed(true), ligand(-1), freeEnergy(0) {}
    
    void appendNucleotide(char nucleotide) {
        //WHEN YOU APPEND A NEW NUCLEOTIDE YOU FRESH THE MOLECULE FOR PARSING
        sequence.push_back(nucleotide); 
        structure.push_back('.');
        notParsed = true;
        
    } 
    
};


class rnaList {
public:
    /*WE HAVE ALLOCATED THE ARRAY FOR THE LIST IN THE HEAP WE CAN 
      ALTERNATIVELY ALLOCATE IT IN THE PROGRAM SPACE FOR SIMPLICITY 
      OF HANDLING THE ARRAY*/
    rnaMolecule *ptr; //THIS SHOULD HOLD THE ADDRESS OF AN ARRAY IN THE HEAP
    int size; //THIS IS THE SIZE OF THE ARRAY
    double minFreeEnergy; 
    
    rnaList(): size(0), minFreeEnergy(0) {
        ptr = new rnaMolecule[MAXCONFOMERS];
    } 
    
    ~rnaList() {
        delete[] ptr; 
    }
    
    int addConfomer(rnaMolecule molecule) {
        if(size < MAXCONFOMERS) {
            *(ptr+size) = molecule;
            //std::cout << (ptr+size)->freeEnergy << std::endl;
            std::cout << molecule.freeEnergy << std::endl;
            size++;
            if(molecule.freeEnergy < this->minFreeEnergy) {
                this->minFreeEnergy = molecule.freeEnergy; 
            }
            return 1; 
        }
        else {
            return 0;
        }
    }
    
    void showList() {
        std::cout << ptr->sequence << std::endl;
        for(int i = 0; i < size; i++) {
            std::cout << (ptr+i)->structure << "\t" 
                    << "\t" << (ptr+i)->ligand << "\t" <<
                    (ptr+i)->freeEnergy << "\t" << (ptr+i)->population 
                     << std::endl; 
        }
    }
    
    void appendNucleotideToList(char newBase) {
        for(int i = 0; i < this->size; i++) {
            (ptr+i)->appendNucleotide(newBase); 
        }
    }
    
    void fetchFromFile(std::string file) {
        
        FILE *fptr = fopen(file.c_str(), "r");
        
        char sequence[MAXMOLELENGTH];
        char structure[MAXMOLELENGTH]; 
        int lgnd; 
        double freeEnergy; 
        double population; 
        
        fscanf(fptr,"%s \n",sequence);
        
        while(!feof(fptr)) { 
            
            fscanf(fptr,"%s %d %lf %lf \n", structure,&lgnd,&freeEnergy,&population); 
            
            rnaMolecule molecule(sequence,structure); 
            molecule.freeEnergy = freeEnergy; 
            molecule.ligand = lgnd; 
            molecule.population = population; 
            
            this->addConfomer(molecule); 
        }
    }
};



/*WE WANT FIRST TO MAKE THE NETWORK MATRIX*/
class rnaReactionSpace {
public:
    
    /*THIS MATRIX HOLDS THE RATE CONSTANTS FOR THE REACTION WE EXPECT 
      TO TAKE PLACE, AND IT WOULD HAVE A VALUE OF 0 IF REACTION
      IS TO OCCUR BETWEEN TWO CONFOMERS*/
    Eigen::MatrixXd reactionsMatrix;
    
     int numOfConfomers;
    
    rnaReactionSpace(): numOfConfomers(1) {
        reactionsMatrix = Eigen::MatrixXd::Zero(MAXCONFOMERS,MAXCONFOMERS);  
    }
    
    void addReactionToElements(int firstReactant,int secondReactant,double rateConstant) {
        reactionsMatrix(firstReactant,secondReactant) = rateConstant; 
        reactionsMatrix(secondReactant,firstReactant) = rateConstant; 
    }
    
    void addNewReaction(int oldMolecule,double rateConstant) { 
        reactionsMatrix(oldMolecule,numOfConfomers) = rateConstant; 
        reactionsMatrix(numOfConfomers,oldMolecule) = rateConstant;
        numOfConfomers++;
    }
    
    void showReactions() {
        for(int i = 0; i < numOfConfomers; i++) {
            for(int j = 0; j < numOfConfomers; j++) {
                std::cout << reactionsMatrix(i,j) << "\t";
            }
            std::cout << std::endl;
        }
    }
    
    void fetchFromFile(std::string file) {
        
        FILE *fptr = fopen(file.c_str(),"r"); 
        
        for(int i = 0; i < numOfConfomers; i++) {
            for(int j = 0; j < numOfConfomers; j++) {
                fscanf(fptr,"%lf",&reactionsMatrix(i,j));
            }
            fscanf(fptr,"\n"); 
        }
    }
    
};
 
