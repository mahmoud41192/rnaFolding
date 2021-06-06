#include <iostream>
#include <stdio.h>

#define MAXLIGANDS 100
#define MAXLIGANDLENGTH 1000

/*WE ASSUME THAT USUALLY THE LIGAND WILL BIND
  TO A ONE SPECIFIC PART OF THE MOLECULE*/


class ligand {
public: 
    std::string bindingSequence; 
    std::string bindingStructure;
    
    double bonusFreeEnergy; 
    double bindingRateConstant; 
    
    ligand(std::string seq, std::string structure): bindingSequence(seq), 
    bindingStructure(structure) {}
    
    ligand() {}
    
};



class ligandList {
public: 
    ligand list[MAXLIGANDS]; 
    int size; 
    
    ligandList(): size(0) {}
    
    void addLigand(ligand lg) {
        list[size] = lg; 
        size++; 
    }
    
    void fetchFromFile(std::string fileName); 
    
    int checkForLigandBinding(rnaMolecule molecule); //TODO
};

void ligandList::fetchFromFile(std::string fileName) {
    
    FILE *file = fopen(fileName.c_str(),"r"); 
    
    char sequence[MAXLIGANDLENGTH]; 
    char structure[MAXLIGANDLENGTH]; 
    double bonusEnergy; 
    double rateConstant; 
    
    while(!feof(file)) {
        fscanf(file,"%s %s %lf %lf \n",sequence,structure,&bonusEnergy,&rateConstant);
        
        //NOW WE MAKE THE LIGAND AND ADD IT TO THE LIST
        ligand buffer(sequence,structure);
        buffer.bonusFreeEnergy = bonusEnergy; 
        buffer.bindingRateConstant = rateConstant; 
        
        //std::cout << buffer.bindingSequence << std::endl;
        
        this->addLigand(buffer); 
    }
   
}

int ligandList::checkForLigandBinding(rnaMolecule molecule) {
    
    if(molecule.ligand == -1) {
        for(int i = 0; i < size; i++) {
            int ligandSize = list[i].bindingSequence.length();
            //WE SEARCH OVER THE WHOLE LENGTH OF THE MOLECULE FOR THE LIGAND MOTIF
            if(molecule.sequence.length() < ligandSize) {
                return -1;
            }
            else{
                for(int j = 0; j < molecule.sequence.length()-ligandSize; j++) {
                    //WE CHOP THE MOLECULE
                    std::string seqBuffer = molecule.sequence.substr(j,ligandSize);
                    std::string strucBuffer = molecule.structure.substr(j,ligandSize);
            
                    //std::cout << seqBuffer << "\t" << strucBuffer << std::endl;
            
                    //CHECK FOR LIGAND BINDING
                    if(seqBuffer == list[i].bindingSequence && strucBuffer == list[i].bindingStructure) {
                        return i;
                    }
                }
            }
        
        }
    }
    
    return -1; 
}

