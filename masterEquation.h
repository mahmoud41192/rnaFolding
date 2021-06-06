/*NOW WE HAVE TO SOLVE THE NUMERICAL STABILITY PROBLEM*/
/*SOLVED !!!!!*/
#include "Eigen/Dense"
#include <vector>
#include <math.h>


double T = 300.0; // IN Kelvins
double R = 0.001987;  //IN Kcal/mole.K

int getMinReaction(rnaList& list,rnaReactionSpace space,int index) {
    
    int minIndex = index; 
    double minFreeEnergy = (list.ptr+index)->freeEnergy; 
    for(int i = 0; i < list.size; i++) {
        if(space.reactionsMatrix(index,i)) {
            if((list.ptr+i)->freeEnergy < minFreeEnergy) {
                minFreeEnergy = (list.ptr+i)->freeEnergy; 
                minIndex = i; 
            }
        }
    }
    
    return minIndex; 
}

/*NOW WE ARE READY TO RUN*/
Eigen::MatrixXd getTransferMatrix(rnaList& list,rnaReactionSpace space) {
    
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(space.numOfConfomers,space.numOfConfomers); 
    
    //ITERATE OVER ALL POSSIBLE REACTIONS
    //THE REACTION READS "PROBABILITY TO GO FROM I TO J"
    for(int i = 0; i < space.numOfConfomers; i++) {
        
        /*WE SHOULD CALCULATE THE DIAGONAL ELEMENTS 
          WHICH NAMES THE REACTIONS THAT CONSUME THE 
          THE MOLECULE OF THAT ROW*/
        /*WE SHOULD SUM THE COLUMNS RATHER THAN THE ROWS */ 
        
        for(int j = 0; j < space.numOfConfomers; j++) {
            if(space.reactionsMatrix(i,j)) {
                //std::cout << i << "\t" << j << std::endl;
                //WE HAVE TO HANDLE IT IF IT IS INFINITE
                //IT OCCURS VERY OFTEN
                double Gi = (list.ptr+i)->freeEnergy; 
                double Gj = (list.ptr+j)->freeEnergy;
                /*THE FREE ENERGY TO GO FROM STATE j TO STATE i*/
                double exponent = -1*(Gi-Gj)/(2*R*T);
                //std::cout << exponent << std::endl;
                //CHECK IF EXPOENET IS TOO LARGE
                //IT IS ONLY A PROBLEM WHEN IT IS TOO BIG
                if(fabs(exponent) > 5) {
                    /*THIS MOVES THE POPULATION TO THE STATE OF MIN FREE ENERGY*/
                    int min = getMinReaction(list,space,i); 
                    double pop = (list.ptr+i)->population;
                    (list.ptr+i)->population = 0; 
                    (list.ptr+min)->population += pop; 
                }
                else {
                    //THE TRANSITION FROM CONFOMER j TO CONFOMER i
                    matrix(i,j) = space.reactionsMatrix(i,j)*exp(exponent); 
                    //std::cout << i << "\t" << j << std::endl;
                    //WE ADD THIS TO THE REVERSE REACTIONS 
                    //std::cout << diagonalReactions << std::endl;
                }
            }
        }
        
    }

    //NOW WE SHOULD ADD THE DIAGONAL ELEMENTS
    for(int i = 0; i < space.numOfConfomers; i++) {
        matrix(i,i) = -1*matrix.col(i).sum(); 
    }
    
    return matrix; 
}

Eigen::MatrixXd getPopulations(rnaList& list) {
    
    Eigen::MatrixXd populations = Eigen::MatrixXd::Zero(list.size,1); 
    
    for(int i = 0; i < list.size; i++) {
        populations(i,0) = (list.ptr+i)->population; 
    }
    
    return populations; 
}


/*THIS WORKS JUST FINE, NOT VERY ACCURATE BUT WORKABLE INDEED
 *LATER ON WE CAN USE A MORE ACCURATE METHOD THAN EULER'S METHOD*/

//WE COMMENT THIS FUNCTION BECAUSE IT SLOWS DOWN COMPILATION

void updatePopulation(Eigen::MatrixXd reactions,Eigen::MatrixXcd& populations,double timeInterval) {  
    
    //FIRST: WE DIAGONALIZE THE REACTIONS MATRIX
    Eigen::EigenSolver<Eigen::MatrixXd> eigen(reactions); 
    Eigen::MatrixXcd eigenVectors = eigen.eigenvectors(); 
    Eigen::VectorXcd eigenValues = eigen.eigenvalues(); 
    
    //std::cout << eigenVectors << std::endl << eigenValues << std::endl;
    
    int size = eigenValues.rows(); 
    //std::cout << size << std::endl; 
    
    //SECOND: WE FIND THE COOEFICIENTS FOR THE THE EIGEN VECTORS 
    Eigen::VectorXcd coefficients = eigenVectors.colPivHouseholderQr().solve(populations); 
    
    //THIRD: WE COMPUTE THE FINAL POPULATION
    Eigen::VectorXcd solution = Eigen::VectorXcd::Zero(size); 
    
    for(int i = 0; i < size; i++) {
        solution += coefficients(i)*eigenVectors.col(i)*exp(eigenValues(i)*timeInterval); 
    }
    
    populations = solution;
}






















