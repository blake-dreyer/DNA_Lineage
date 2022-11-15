// lineage.cpp
// Purpose: This program prompts the user for a file and then for multiple queries
// related to martian genetic date stored in the file. Queries are answered via
// recursive functions in some cases.
// Written By: Blake Dreyer (bdreye02) for CS 11 at Tufts University
// Date: 09 Nov 2022
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

struct Gene;

struct Mutation {
    int cost;
    Gene *target;
};

struct Gene {
    string sq;
    bool seen; //Used in recursive functions to avoid loops
    Mutation *muts;
};

//Hard coded available queries
const string EVOLVE = "e";
const string E_STEPS = "es";
const string ENE_EVOLVE = "ene";
const string PATH = "path";
const string QUIT = "q";


int  create_graph(Gene **);
void populate_nodes(Gene *dna, ifstream *inf, int num_nodes);
void populate_links(Gene *dna, ifstream *inf, int num_nodes);
int  index_of(Gene *, string, int);
void init(Gene *, int);
bool can_evolve(Gene *dna, Gene *src, Gene *tgt);
int can_evolve_steps(Gene *dna, Gene *src, Gene *tgt);
bool energetic_evolution(Gene *dna, Gene *src, Gene *tgt, int energy);
bool evolutionary_path(Gene *dna, Gene *src, Gene *tgt, string *ps);

int main(){

    Gene *dna;

    //After this line executes, dna is a pointer to the array of gene
    //information read in from a given data file, and num_nodes contains the
    //number of genes stored in that array. 
    int num_nodes = create_graph(&dna);

    //USE dna AND num_nodes AFTER THIS LINE TO CONTINUE THE COMPUTATION

    //Take user queries
    string query;
    cout << "Enter a query: ";
     while (cin >> query) { 
        string source;
        string target;
        int energy;
        cin.ignore();
    
        if (query == EVOLVE){
            cin >> source;
            cin >> target;
            int source_index = index_of(dna, source, num_nodes);
            int target_index = index_of(dna, target, num_nodes);
            bool possible = can_evolve(dna, &dna[source_index], &dna[target_index]);
            if(possible){
                cout << source << " can evolve into " << target << endl;
            }else{
                cout << source << " cannot evolve into " << target << endl;
            }
        } else if (query == E_STEPS){
            cin >> source;
            cin >> target;
            int source_index = index_of(dna, source, num_nodes);
            int target_index = index_of(dna, target, num_nodes);
            int n = can_evolve_steps(dna, &dna[source_index], &dna[target_index]);
            cout << "It will take " << n << " evolutionary steps to get from " << source << " to " << target << endl;
        }else if (query == ENE_EVOLVE){
            cin >> source;
            cin >> target;
            cin >> energy;
            int source_index = index_of(dna, source, num_nodes);
            int target_index = index_of(dna, target, num_nodes);
            bool possible = energetic_evolution(dna, &dna[source_index], &dna[target_index], energy);
            if(possible){
                cout << source << " can evolve into " << target << " with at most " << energy << " evolutionary cost" << endl;
            }else{
                cout << source << " cannot evolve into " << target << " with at most " << energy << " evolutionary cost" << endl;
            }
        }else if (query == PATH){
            cin >> source;
            cin >> target;
            string path_steps = "";
            int source_index = index_of(dna, source, num_nodes);
            int target_index = index_of(dna, target, num_nodes);
            bool poss = evolutionary_path(dna, &dna[source_index], &dna[target_index], &path_steps);
            if(!poss){
                cout << "There is no path from " << source << " to " << target << endl;
            }else{
                cout << source << path_steps << " -> " << target << endl;
            }
            
        } else if (query == QUIT){
            break;
        } else {
            cout << query << " not recognized." << endl;
        }
        cout << endl << "Enter a query: ";
    }

    for (int i = 0; i < num_nodes; i++) {
        delete dna[i].muts;
    }
    delete [] dna;
    return 0;
}



// Creates graph and stores it inside of array at passed location.
// create_graph 
// Input: An uninitialized pointer to a Gene pointer.
// Description: Read in a file provided the user and use the data to populate an
//              array of Genes.
// Output: Populates the array pointed to by *dna_p with the contents of a data
//         file, and returns the number of array elements
int create_graph(Gene **dna_p){
    string filename;
    cout << "Enter data file name: ";
    cin >> filename;

    ifstream inf(filename);
    if (inf.fail()) {
        cerr << "ERROR OPENING FILE: Exiting Program" << endl;
        exit(EXIT_FAILURE);
    }
    
    int num_nodes;
    inf >> num_nodes;
    inf.ignore(); 
    
    // Create and populate the nodes in the array
    *dna_p = new Gene[num_nodes];
    init(*dna_p, num_nodes);
    populate_nodes(*dna_p, &inf, num_nodes);

    // Reset the file to read the links
    inf.close();
    inf.open(filename);
    inf >> num_nodes;
    inf.ignore();

    populate_links(*dna_p, &inf, num_nodes);
    
    return num_nodes;
}

// init 
// Input: A pointer to an uninitialized array of genes and the size of 
//        that array.
// Description: Initializes the array of genes with non-garbage information.
// Output: Initializes the array pointed to by 'dna'.
void init(Gene *dna, int num_nodes){
    for (int i = 0; i < num_nodes; i++) {
        dna[i].sq = "";
        dna[i].seen = false;
        dna[i].muts = nullptr;
    }
}

// populate_nodes
// Input: A pointer to an array of genes, a pointer to a file, 
//        and the size of the array.
// Description: Read the file and populate the genes in the array with their
//              sequence information.
// Output: Populates the sequence information for the gene array pointed to by
//         'dna'. Moves the file pointer further into the file.
void populate_nodes(Gene *dna, ifstream *inf, int num_nodes){
    string line;
    stringstream sstr;
    for (int i = 0; i < num_nodes; i++) {
        getline(*inf, line);
        sstr.str(line);
        sstr >> dna[i].sq;
        sstr.clear();
    }
}

// populate_links
// Input: A pointer to an array of genes, a pointer to a file, 
//        and the size of the array.
// Description: Read the file and populate the mutation information between
//              genes in the array. We assume the data file only contains
//              well-formed input.
// Output: Populates the mutation information for the gene array pointed to by
//         'dna'. Moves the file pointer further into the file.
void populate_links(Gene *dna, ifstream *inf, int num_nodes){
    string line, seq, mut_seq;
    stringstream sstr;
    int cost;
    for (int i = 0; i < num_nodes; i++) {
        getline(*inf, line);

        sstr.str(line);
        sstr >> seq;//Skip over first gene name on each line
        
        // Add a mutation if it exists
        sstr >> mut_seq >> cost;
        if (!sstr.fail()) {
            int mut_index = index_of(dna, mut_seq, num_nodes);
            dna[i].muts = new Mutation;
            dna[i].muts->cost = cost;
            dna[i].muts->target = &dna[mut_index];
        }
        //Prepare to read another line from this string stream
        sstr.clear();
    }
}

// index_of
// Input: An array of genes, a gene sequence, and the size of the array.
// Description: Find the index of the given sequence in the array.
// Output: The index of the sequence in the array, or -1 if the sequence is not
//         in the array.
int index_of(Gene *dna, string seq, int num_nodes){
  for (int i = 0; i < num_nodes; i++) {
    if (dna[i].sq == seq) {
      return i;
    }
  }
  return -1;
}

// can_evolve
// Input: The genes array and two pointers to genes in the Genes array, one pointed at the target gene
// and one pointed at the source gene
// Description: Checks whether or not the given source gene can mutate into the 
// given target gene
// Output: Returns a boolean that represents whether the proposed evolution is 
// possible.
bool can_evolve(Gene *dna, Gene *src, Gene *tgt){
    if (src->muts == nullptr){
        return false;
    }
    if (src->muts->target->sq == tgt->sq){
        return true;
    }
    if (src->seen == true){
        return false;
    }
    src->seen = true;
    Gene *new_src = src->muts->target;
    return can_evolve(dna, new_src, tgt);
}

// energetic_evolution
// Input: The genes array and two pointers to genes in the Genes array, one pointed at the target gene
// and one pointed at the source gene. An int representing the user's amount of specified energy
// Description: Checks whether or not the given source gene can mutate into the 
// given target gene with the amount of specified energy.
// Output: Returns a boolean that represents whether the proposed evolution is 
// possible.
bool energetic_evolution(Gene *dna, Gene *src, Gene *tgt, int energy){
    if (src->muts == nullptr){
        return false;
    }
    if (src->muts->target->sq == tgt->sq && energy >= src->muts->cost){
        return true;
    }
    if (src->seen == true){
        return false;
    }
    src->seen = true;
    Gene *new_src = src->muts->target;
    int price = src->muts->cost;
    energy -= price;
    return energetic_evolution(dna, new_src, tgt, energy);
}

// can_evolve_steps
// Input: The genes array and two pointers to genes in the Genes array, one pointed at the target gene
// and one pointed at the source gene.
// Description: Checks whether or not the given source gene can mutate into the 
// given target gene and returns either -1 or the number of steps if possible
// Output: Returns an int that represents whether the proposed evolution is 
// possible and in how many steps
int can_evolve_steps(Gene *dna, Gene *src, Gene *tgt){
    if (src->muts == nullptr){
        return -1;
    }
    if (src->muts->target->sq == tgt->sq){
        return 1;
    }
    if (src->seen == true){
        return -1;
    }
    src->seen = true;
    Gene *new_src = src->muts->target;
    int n = can_evolve_steps(dna, new_src, tgt);
    if(n != -1){
        n = n+1;
    }
    return n;
    
}

// evolutionary_path
// Input: The genes array and two pointers to genes in the Genes array, one pointed at the target gene
// and one pointed at the source gene
// Description: Checks whether or not the given source gene can mutate into the 
// given target gene and updates a string representing the path when necessary.
// Output: Returns a boolean that represents whether the proposed evolution is 
// possible. Modiefies a string to represent the path
bool evolutionary_path(Gene *dna, Gene *src, Gene *tgt, string *path_steps){
    if (src->muts == nullptr){
        return false;
    }
    if (src->muts->target->sq == tgt->sq){
        return true;
    }
    if (src->seen == true){
        return false;
    }
    src->seen = true;
    Gene *new_src = src->muts->target;
    *path_steps += (" -> " + src->muts->target->sq);
    return evolutionary_path(dna, new_src, tgt, path_steps);
}