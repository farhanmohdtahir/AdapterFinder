//#include<gperftools/profiler.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <getopt.h>
#include "Input.h"
#include "NW.h"
#include "CS.h"
#include "time.h"


using namespace std;

/**	@mainpage Paired-End Adapter Finder v15.4.14

	The Paired-End Adapter Finder constructs two consensus adapter sequences from Reads 1 and Reads 2 of paired-end sequencing by using the Needleman-Wunsh algorithm to align the two reads and determining which region of the sequence is actually an adapter. The user has to input two fastq format files, corresponding to Read 1 and Read 2, and the result will be the consensus adapter sequences for both the Adapters in Read 1 and Read 2 respectively of the paired-end sequencing.
	
	@author Rayan Gan
	@date April 2015
 */
 
void print_usage();
void help();

int main(int argc, char *argv[]) 
{	
    //clock
   clock_t t1,t2;
    t1=clock();	 
    //ProfilerStart("result.prof");
	string file1 = "test.fastq", file2 = "test2.fastq", seq_1, seq_2, seq_1_al, seq_2_al;
  	int opt = 0, seqLength = 0, debugLevel = 0; //check=0;
  	double percentage = 0, confLevel = .0;
  	bool option = false;
  	int rowmax = 0, colmax = 0, confTrue, adapLenCount = 0, iteration = 0;
  	bool back = true;
  	    	    	int q=0,w=0;
 	string line, line2;
 	bool onlynuc = true;
 	bool onlynuc2 = true;
 	
          static struct option long_options[] = {
        {"help",      no_argument,       0,  'h' },
        {"f1",    required_argument,     0,  'a' },
        {"f2",    required_argument,     0,  'b' },
        {"seql",  optional_argument,     0,  'c' },
        {"perc",  optional_argument,     0,  'd' }, 
        {"conf",  optional_argument,     0,  'e' },  
        {"debug", optional_argument,     0,  'f' },                                
        {0,                0,            0,   0  }
        };

        int long_index =0;
        while ((opt = getopt_long_only(argc, argv,"", 
                   long_options, &long_index )) != -1) {
        switch (opt) {
             case 'h' : help();option = true;
                 break;
             case 'a' : file1 = optarg;option = true;
                 break;
             case 'b' : file2 = optarg;option = true;
                 break;
             case 'c' : seqLength = atoi(optarg);option = true;
                 break;
             case 'd' : percentage = atoi(optarg);option = true;
                 break;  
             case 'e' : confLevel = atoi(optarg);option = true;
                 break; 
             case 'f' : debugLevel = atoi(optarg);option = true;
                 break;                                                      
             default: print_usage(); 
                 exit(EXIT_FAILURE);
        }
        }
        //cout<<option;
    if(seqLength == 0) seqLength = 70;
    if(percentage == 0) percentage = 85;
    if(confLevel == 0) confLevel = 1;
    if(debugLevel == 0) debugLevel = 0;
       
    if(option == false) print_usage();
        
	if(debugLevel == 0 || debugLevel == 1 || debugLevel == 2){

    if(file1 != "" && file2 != "")
    {
       
    	Input ab;
	NW b;
 	CS c;
	CS d;
 	ifstream myfile (file1.c_str());
 	ifstream myfile2 (file2.c_str());
 	if (myfile.is_open() && myfile2.is_open())
  	{
  	 while (getline (myfile,line) && getline (myfile2,line2))
  	 {
  	 /*
  	    ab.checkNucleotide(line, seq_1);
  	    ab.checkNucleotide(line2, seq_2);	
  	 */    
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
			seq_1 = line;
                        q+=1;
		 }
    	    }
  	    if(line2[0]=='A'||line2[0]=='C'||line2[0]=='G'||line2[0]=='T'||line2[0]=='N')
   	    {
		 onlynuc2 = false;
		 for(int a = 0; a < line2.length(); a++)
		 {	
			if(line2[a]=='A'||line2[a]=='C'||line2[a]=='G'||line2[a]=='T'||line2[a]=='N')
				onlynuc2 = true;
 			else
			{ 
				onlynuc2 = false;
				break;
			}
		 }
		 if(onlynuc2 == true)
		 {
			seq_2 = line2;
                        w+=1;
		 }
		 
// Creating NW and CS objects
//           reverse( seq_1.begin(), seq_1.end() );     
//	    ab.complementInput(seq_1);     
//            cout<<endl;
          // cout<<"Seq2:"<<seq_2<<endl;
	   reverse( seq_2.begin(), seq_2.end() );
	    ab.complementInput(seq_2);    
            //cout<<seq_2<<endl ;
    	    double L1 = seq_1.length();
	    double L2 = seq_2.length();	   

	    b.nw(seq_1, seq_2, seq_1_al, seq_2_al, debugLevel);
	   // cout<<seq_1_al<<endl<<endl;
 
	    //cout<<b.rowmax<<endl;
	    if(b.percentage > percentage && (b.rowmax/L1*100) > seqLength && ((seq_2.length()-b.colmax)/L2*100) > seqLength && b.colmax != 0)
	    {
          	if(debugLevel == 1 || debugLevel == 2)
	        {
	    	    cout<<"\nSeq1:"<<seq_1<<endl;
	    	    cout<<"Seq2:"<<seq_2<<endl;
	        }
	    
		rowmax = b.rowmax;
		colmax = seq_2.length()-b.colmax;		
		c.cs(seq_1, rowmax);
		reverse( seq_2.begin(), seq_2.end() );
		d.cs(seq_2, colmax); 	

		if(debugLevel == 1 || debugLevel == 2)c.print_nucCount_phred();
		if(debugLevel == 1 || debugLevel == 2)d.print_nucCount_phred();

		confTrue = 0;
		c.checkConfidence(confLevel, confTrue, adapLenCount);
		d.checkConfidence(confLevel, confTrue, adapLenCount);
               // cout<<confTrue<<endl;
                
		adapLenCount++;
//                check++;
//                cout<<check;
		if(confTrue == 2)
		{
			cout <<"\n";
			c.print_cs(0);
			cout << "	";
			d.print_cs(1);
			cout << "\n\n";
                        t2=clock();
    float diff ((float)t2-(float)t1);
    float seconds = diff / CLOCKS_PER_SEC;
    cout<<seconds<<endl;
			exit(0);
		}	
	    }
	
// After NW and CS
               // cout<<debugLevel<<endl; 
                //cout<<b.percentage<<"   "<<percentage<<"   "<<b.rowmax/L1*100<<"   "<<seqLength<<"   "<<(seq_2.length()-b.colmax)/L2*100<<"   "<<b.colmax<<endl;	    
    	    }



  	  }
  	  myfile.close();
  	  myfile2.close();
                         // cout<<q<<endl<<w<<endl;
  	  cout << "Confidence level could not be achieved...\n";
          
 	}
        else cout << "Unable to open file"; 
    }}
    //ProfilerStop();  
           //clock
 t2=clock();
    float diff ((float)t2-(float)t1);
    float seconds = diff / CLOCKS_PER_SEC;
    cout<<seconds<<endl;
  return 0;
}

void print_usage() 
{
    cout<<"Usage: ./PEAdapterFinder -f1 filename -f2 filename -seql percentage length -perc percentage match -conf confidence level -debug debug level\n";
}

void help()
{
    cout<<"\nPaired-End Adapter Finder\n\n";
    cout<<"Please enter ./PEAdapterFinder -f1 file1 -f2 file2 -seql seqLength -perc percentage -conf confLevel -debug debugLevel\n";
    cout<<"\nfile1 - name of first fastq file \n";
    cout<<"file2 - name of second fastq file\n";
    cout<<"seqLength - minimum length percentage to get adapter sequence (default = 70, to change use '-seql=')\n";
    cout<<"percentage - minimum match percentage to get adapter sequence (default = 85, to change use '-perc=')\n";   
    cout<<"confLevel - minimum confidence level of nucleotides (default = 1, to change use '-conf=')\n";
    cout<<"debugLevel - debug level of programme (default = 0, to change use '-debug=' : 0 - only adapter sequences, 1 - nucleotide count and phred score, 3 - dynamic programming matrix and traceback matrix)\n";         
    cout<<"\nPlease refer to the documentation for more help...\n";
    cout<<"\nThank you. \n\n";
}
