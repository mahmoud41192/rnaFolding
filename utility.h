

void writeListToFile(std::string fileName, rnaList& list,rnaReactionSpace reactions) {
    
    std::fstream file;
    file.open(fileName,std::ios::out|std::ios::app); 
    
    file << list.ptr->sequence << std::endl;
    
    for(int i = 0; i < list.size; i++) {
        file << (list.ptr+i)->structure << "\t";
        
        file << "\t" << (list.ptr+i)->ligand << "\t" << (list.ptr+i)->freeEnergy 
                << "\t" << (list.ptr+i)->population << std::endl;
    }
    
    file << std::endl;
    
    //WE WRITE THE REACTION MATRIX AS WELL
    for(int i = 0; i < reactions.numOfConfomers; i++) {
        for(int j = 0; j < reactions.numOfConfomers; j++) {
            file << reactions.reactionsMatrix(i,j) << "\t"; 
        }
        
        file << std::endl;
    }
    
}

