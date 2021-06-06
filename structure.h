#include <math.h>

/*THIS HOLDS THE DIFFERENT RATE CONSTANT*/
double rateConstants[] = {1e9,1e9};


/*THIS FUNCTION CHECKS IF THE REACTION 
 *IS RATE LIMITING*/
/*THE PROBLEM IS THAT WE ONLY CHECK WHEN 
 *WE HAVE THE REACTION BEING ADDED 
 *THROUHG THE MAKEONEMOVE FUNCTION*/
bool isRateLimiting(rnaMolecule molecule,int first,int second) {
    if(first > 1 && first < molecule.sequence.length()-1 && second > 1 && 
            second < molecule.sequence.length()-1) {
        if(molecule.structure[first-1] == '(' && molecule.structure[first-2] == '.') 
            return true; 
        if(molecule.structure[second-1] == ')' && molecule.structure[second-2] == '.')
            return true; 
        if(molecule.structure[first+1] == '(' && molecule.structure[first+2] == '.')
            return true;
        if(molecule.structure[second+1] == ')' && molecule.structure[second+2] == '.')
            return true;
    }
    
    return false; 
}


/*THIS CHECKS IF A PAIR OF BASES IS A BASEPAIR
  WE MAKE IT INT TO KNOW WHICH BASEPAIR IS 
  BEING FORMED IN THE NEW MOLECULE AND THAT 
  IS SUCH WE CAN ASSIGN THE CORRECT RATE CONSTANT*/

int isBasePair(char first,char second) {
    if(first == 'G' && second == 'C')
        return 1; 
    else if(first == 'C' && second == 'G')
        return 1; 
    else if(first == 'A' && second == 'U')
        return 2; 
    else if(first == 'U' && second == 'A')
        return 2; 
    else 
        return -1; 
}

/*THIS SEES IF AN ELEMENT IS ALREADY 
 *IN THE LIST OR NOT*/
/*ROBUST*/
bool isNotInList(rnaMolecule& molecule,rnaList& list) {
    for(int i = 0; i < list.size; i++) {
        if((list.ptr+i)->structure == molecule.structure)
            return false;
    }
    return true; 
}

/*THIS FUNCTION GETS THE INDEX OF A CONFOMER 
 *IN THE LIST AND RETURNS -1 IF NOT IN THE LIST*/
int findInList(rnaMolecule& molecule,rnaList& list) {
    for(int i = 0; i < list.size; i++) {
        if((list.ptr+i)->structure == molecule.structure) {
            return i; 
        }
    }
    return -1; 
}

/*THIS FUNCTION CHECKS IF THE REACTION IS MUCH SLOW COMPARED TO THE 
 *THE MIN FREE ENERGY PATH*/
/*ROBUST*/
bool notSlow(double freeEnergy,double min) {
    double slowRate = 6;  //kcal/mole 
    if(freeEnergy-min > slowRate) {
        return 0;
    }
    else 
        return true; 
}


/*THIS FUNCTION RETURNS AN INDEX */

/*THIS FUNCTION MAKES ONE ADDITION MOVE 
  FOR THE ARGUMENT CONFOMER AND SAVES IT 
  IN THE ARGUMENT LIST*/
/*ROBUST*/
void makeOneMove(rnaMolecule molecule, rnaList& list, rnaReactionSpace& reactions,int index, ligandList ligands) {
    //FIRST WE ITERATE OVER ALL POSSIBLE COMBINATIONS
    int lastNucleotide = molecule.sequence.length()-1;
    int minLoop = 3; 
    
    int helices = 0; 
    
    std::cout << molecule.structure << std::endl;
    
    for(int i = 0; i < molecule.sequence.length()-minLoop; i++) {
        //std::cout << i << "\t" << lastNucleotide << std::endl;  
        //CHECK FOR UNBINDED ELEMENTS
        
        /*WE HAVE TO AVOID HAVING CROSS PAIRING*/
        if(molecule.structure[i] == '(') {
            helices++;
        }
        else if(molecule.structure[i] == ')') {
            helices--; 
        }
        
        
        if(molecule.structure[i] == '.' && molecule.structure[lastNucleotide] == '.' && helices == 0) {
            
            int pair = isBasePair(molecule.sequence[i],molecule.sequence[lastNucleotide]);
            
            if(pair != -1) {
                //std::cout << "sequence " << molecule.sequence[i] << std::endl;
                rnaMolecule buffer(molecule.sequence,molecule.structure); 
                buffer.structure[i] = '('; 
                buffer.structure[lastNucleotide] = ')';
                std::cout << buffer.structure << std::endl;
                /*WE NEED TO CHECK IF IT IS RATE LIMITING*/
                /*if(isRateLimiting(buffer,i,lastNucleotide)) {
                    buffer.rateLimitingReactions++;
                    std::cout << "rate limiting reactions" << std::endl
                          << buffer.structure << std::endl;
                }*/
                
                    
                /*A NEW VERSION OF THE IN LIST WHICH ADDS 
                  THE REACTION*/
                int var = findInList(buffer,list);
                //std::cout << var << std::endl;
                if(var == -1) {
                    //A NEW CONFOMER
                    writeStructToFile(buffer); 
                   
                    //COMPUTE THE FREE ENERGY
                    computeFreeEnergy(buffer);
                    
                    //ADD LIGAND FREE ENERGY
                    if(buffer.ligand != -1) {
                        buffer.freeEnergy = buffer.freeEnergy+ligands.list[buffer.ligand].bonusFreeEnergy;
                        std::cout << buffer.freeEnergy << std::endl;
                    }
                    
                    
                    //std::cout << buffer.freeEnergy << std::endl; 
                    
                    //WE SHOULD CHECK IF IT IS FAST OR SLOW
                    if(notSlow(buffer.freeEnergy,list.minFreeEnergy)) {
                        
                        list.addConfomer(buffer);
                        reactions.addNewReaction(index,rateConstants[pair-1]);
                    
                        //HERE IS THE PROBLEM
                        int ligandIndex = ligands.checkForLigandBinding(buffer);
                        std::cout << "found ligand" << std::endl;
                        if(ligandIndex != -1 && buffer.ligand == -1) {
                            //WE PUSH ANOTHER ELEMENT WITH THE LIGAND BOUND
                            rnaMolecule buffer2(buffer.sequence,buffer.structure); 
                            buffer2.ligand = ligandIndex;
                            //WE SHOULD ADD THE BONUS FREE ENERGY DUE TO LIGAND BINDING
                            buffer2.freeEnergy = buffer.freeEnergy+ligands.list[ligandIndex].bonusFreeEnergy; 
                            list.addConfomer(buffer2); 
                            reactions.addNewReaction(reactions.numOfConfomers-1,ligands.list[ligandIndex].bindingRateConstant);
                            
                        }
                        
                    }
                }
                else {  
                    //IS A REACTION BETWEEN EXISTING ELEMENTS
                    //IT CAN NEVER BE A LIGAND BINDING REACTION
                    reactions.addReactionToElements(var,index,rateConstants[pair-1]);
                }
            }
        }
    }
    
    molecule.notParsed = false; 
}


/*THIS MAKES A ONE SCAN AND PARSE LOOP*/

int findNewStructures(rnaList& list,rnaReactionSpace& space, ligandList ligands) {
    
    int initialSize = list.size; //TO KEEP TRACK OF ELEMENTS
    
    for(int i = 0; i < initialSize; i++) {
        if((list.ptr+i)->notParsed) {
            makeOneMove(*(list.ptr+i),list,space,i,ligands);
        }
    }
    
    return list.size - initialSize; 
}


