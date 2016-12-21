#include "iostream"
#include "fstream"
#include "string"
#include "algorithm"
#include "stdlib.h"
#include "stdio.h"
#include "Input.h"

using namespace std;


/** @brief The Input Class transforms the sequences to be used in the Adapter Sequencer.      

	The Input Class is responsible to reverse complement Read 2 before dynamic programming is done. 
	@author Rayan Gan
	@date April 2015
 */

int Input::checkNucleotide(string line, string& seq)
{
     	bool onlynuc = true;
  	if(line[0]=='A'||line[0]=='C'||line[0]=='G'||line[0]=='T'||line[0]=='N')
    	{
		onlynuc = false;
	 	for(int a = 0; a < line.length(); a++)
	 	{	
			if(line[a]=='A'||line[a]=='C'||line[a]=='G'||line[a]=='T'||line[a]=='N')
				onlynuc = true;
			else
			{ 
				onlynuc = false;
				break;
				}
		}
		if(onlynuc == true)
		{
			seq = line;
		}
   	 }
   	 return 0;
}    

/**Reverse complement the input for sequences from the second input file.     
 * @param seq The sequence to be complemented
 */
int Input::complementInput(string& seq)
{
        for(int reverse = 0; reverse < seq.length(); reverse++)
        {
            switch(seq[reverse]){
                case 'A': seq[reverse]='T';break;
                case 'T': seq[reverse]='A';break;
                case 'G': seq[reverse]='C';break;
                case 'C': seq[reverse]='G';
            }
        }      
	return 0;
}

string Input::reform(string file, bool &fourline){
    bool startline=false;
    ifstream myfile;
    ofstream outfile;
    string line, str;
    int cmpvar=0;
    
    myfile.open(file.c_str(),ios::in);
    file="4line_"+file;
    outfile.open(file.c_str(),ios::out); 
    
    while (getline(myfile, line)){
        
                if (startline==true) {
                        if(line[0]=='@'||line[0]=='+'){
                            if (cmpvar==0) str+=line;
                            else str+="\n"+line;
                            cmpvar=0;
                        }      
                        else {
                            if (cmpvar==1) str+=line;
                            else str+="\n"+line;
                            cmpvar=1;
                        }
                    }
                    
                else {
                        if(line[0]=='@'){
                            startline=true;
                            str+=line;
                            cmpvar=0;
                        }
                    }
    }
            outfile<<str;
            outfile.close();
            fourline=true;
            return file;
}

