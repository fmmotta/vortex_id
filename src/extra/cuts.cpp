#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <fstream>
using namespace std;

int main()
{	
  double Pi = 3.14159265359;
  ifstream myReadFile;
  ofstream myWriteFile, axisFile,myCheckFile;
  char nameout[67];
  double NFs,NFe,nSteps,nCurrent;
  
  //cout << "Start Cycle" << endl;
  //cin >> NFs;
  //cout << "End Cycle" << endl;
  //cin >> NFe;
  //cout << "Steps" << endl;
  //cin >> nSteps;

  NFs = 10.0075;
  //NFe = 10.0195;
  NFe = 10.0075;
  nSteps = 0.0120;
  //--	
	int NN,Nx,Nxmax,Ny,Nz, NzSkip;
  //cout << "Domain Definition [NN,Nx,Ny,Nz]" << endl;
  //cin >> NN >> Nx >> Ny >> Nz;
  //cout << "Skip how many points in z direction?" << endl;
  //cin >> NzSkip;

  NN = 1;
  Nx = 256;
  Nxmax = 256;
  Ny = 96;
  //Ny=96;
  //Nz = 384;
  Nz = 192;
  NzSkip = 32;
	
  //-- Grid distribution in Y direction
  //-- Grid istribution in X,Z are fixed and known: x-> 0,2pi
  
  double dx = 2.0*Pi/(double)Nxmax;
  double dz = Pi/(double)Nz;
  double Y[1000],Z[1000];
  double *mU,*mV,*mW;
  
  mU = new double [Nx*Ny*Nz];
  mV = new double [Nx*Ny*Nz];
  mW = new double [Nx*Ny*Nz];

  string line;

  //-----------------------------------------------------------  
  // Start getting the Y positions
  
  myReadFile.open("../../data/DNS_OPEN_FOAM/constant/polyMesh/points");
  for (int ii=0;ii<20;ii++) 
    getline(myReadFile,line);
	
  for (int j=0;j<Ny+1;j++){
    getline(myReadFile,line);
    int f1 = line.find(" ");
    int f2 = line.find(" ",f1+1);

    string s1 = line.substr(1,f1-1);
    string s2 = line.substr(f1+1,f2-(f1+1));
    string s3 = line.substr(f2+1,line.find(")")-(f2+1));
    Y[j] = atof(s2.c_str());

    for (int i=1;i<Nx+1;i++){
      getline(myReadFile,line);
      f1 = line.find(" ");
      f2 = line.find(" ",f1+1);
      s1 = line.substr(1,f1-1);
      s2 = line.substr(f1+1,f2-(f1+1));
      s3 = line.substr(f2+1,line.find(")")-(f2+1));			
      //-- stores the vertex position in y direction
    }
  }
  myReadFile.close();
  
  myWriteFile.open("Reference/yAxis.dat",ios::out);
  myWriteFile.precision(8);
  std::sort(Y, Y+Ny+1);
  for(int j=0;j<Ny;j++)
    myWriteFile << (Y[j]+Y[j+1])/2 << endl;
  myWriteFile.close();

  myWriteFile.open("Reference/zAxis.dat",ios::out);
  myWriteFile.precision(8);
  for(int k=0;k<Nz;k+=1)
    myWriteFile << (k+0.5)*dz << endl;
  myWriteFile.close();

  // End getting the Y positions
  //---------------------------------------------------------

  myReadFile.open("../../data/DNS_OPEN_FOAM/constant/polyMesh/points");
  myWriteFile.precision(8);
  for (int ii=0;ii<20;ii++) 
    getline(myReadFile,line);
	
  for (int j=0;j<Ny+1;j++){
    for (int i=0;i<Nx;i++){
      getline(myReadFile,line);
      int f1 = line.find(" ");
      int f2 = line.find(" ",f1+1);
      string s1 = line.substr(1,f1-1);
      string s2 = line.substr(f1+1,f2-(f1+1));
      string s3 = line.substr(f2+1,line.find(")")-(f2+1));			
      //-- stores the vertex position in y direction
      Y[j] = atof(s2.c_str());
    }
  }
  myReadFile.close();

  //-- Iteration through Hamid's files
  nCurrent = NFs;
  int nnn = 0;

  //cout << "got here\n";
  
  while (nCurrent<=NFe) {
    //-- Read first lines of the file	
    //cout << "nCurrent =" << nCurrent << endl;

    char FN[100],FN2[100];
    sprintf(FN,"../../data/DNS_OPEN_FOAM/%6.4f/U",nCurrent);
    myReadFile.open(FN, ios::in);
    
    if (myReadFile.is_open()) {
      //cout << "Opened File \n";
      //cout << FN << endl;
      nnn++;
      for (int ii=0;ii<22;ii++) 
      	getline(myReadFile,line);	
      
      double dU,dV,dW;		
      for (int n=0;n<NN;n++) {
        //--- Each file
		    /*
			  for (int k=0;k<Nz;k+=NzSkip) {

				  // create file to save X-Y axis values 
				  if((nCurrent == NFs && k==0))
			  		axisFile.open("planes/axis.dat");
		  			
		  		// create output file: x-y plane cuts 
		  		sprintf(nameout,"planes/t%.3d_z%.3d.dat", nCurrent,k);
		  		myWriteFile.open(nameout, ios::out);
			
			  	for (int j=0;j<Ny;j++) 
			  	{
			  		for (int i=0;i<Nx;i++)
			  		{
			  			//-- read the variable
			  			getline(myReadFile,line);						
			  			int f1 = line.find(" ");
			  			int f2 = line.find(" ",f1+1);
			  			string s1 = line.substr(1,f1-1);
			  			string s2 = line.substr(f1+1,f2-(f1+1));
			  			string s3 = line.substr(f2+1,line.find(")")-(f2+1));
			  			dU = atof(s1.c_str());
			  			dV = atof(s2.c_str());
			  			dW = atof(s3.c_str());					
			  			//-- dU is U(i,j,k)
			  			//-- dV is V(i,j,k)
			  			//-- dW is W(i,j,k)
			  			if (axisFile.is_open()) 
			  				axisFile << (double)i*dx << " " << Y[j] << "\n";						
                                                  myWriteFile.precision(8);
		  				myWriteFile << dU << " " << dV << "\n";
		  			}
			  	}
			  	// close file with X-Y axis 
			  	if((nCurrent == NFs && k==0))
			  		axisFile.close();
    			myWriteFile.close();
		  	}*/
			
			  for(int i=0;i<Nx;i++){
			    /* create file to save Z-Y axis values */
			    //cout << "begin processing plane " << i << "\n";
			    if((nCurrent == NFs && i==0))
                  axisFile.open("Reference/axis.dat");
					
			    /* create output file: x-y plane cuts */
			    			  
			    for(int k=0;k<Nz;k+=1){
			      for(int j=0;j<Ny;j+=1){
					
				      //-- read the variable
				      getline(myReadFile,line);						
     				  int f1 = line.find(" ");
		     		  int f2 = line.find(" ",f1+1);
			    	  string s1 = line.substr(1,f1-1);
				      string s2 = line.substr(f1+1,f2-(f1+1));
				      string s3 = line.substr(f2+1,line.find(")")-(f2+1));
				      dU = atof(s1.c_str());
				      dV = atof(s2.c_str());
				      dW = atof(s3.c_str());					
				      //-- dU is U(i,j,k)
				      //-- dV is V(i,j,k)
				      //-- dW is W(i,j,k)
               mU[Nx*Ny*k+Nx*j+i] = dU;
               mV[Nx*Ny*k+Nx*j+i] = dV;
               mW[Nx*Ny*k+Nx*j+i] = dW;

				      if (axisFile.is_open()) 
				        axisFile << (double)k*dz << " " << Y[j] << "\n";						
                                   myWriteFile.precision(8);
				      //myWriteFile << dW << " " << dV << "\n";
				    }
			    } 
			    //cout << "end processing plane " << i << "\n";
			  }
			  //--- Each file
        /*
        myCheckFile.open("planes/check.dat");
        for(int i=0;i<Nx;i++){
          for(int k=0;k<Nz;k+=1){
            for(int j=0;j<Ny;j+=1){
              dU = mU[Nx*Ny*k+Nx*j+i];
              dW = mW[Nx*Ny*k+Nx*j+i];
              dV = mV[Nx*Ny*k+Nx*j+i];
              myCheckFile << dU << " " << dV << " " << dW << "\n";
            }
          }
        }
        myCheckFile.close();
        */
        for(int i=0;i<Nx;i++){
         	sprintf(nameout,"Reference/t%.3f_x%d.dat", nCurrent,i);
			    myWriteFile.open(nameout, ios::out);

			    for(int k=0;k<Nz;k+=1){
			      for(int j=0;j<Ny;j+=1){
              dW = mW[Nx*Ny*k+Nx*j+i];
              dV = mV[Nx*Ny*k+Nx*j+i];
			      	myWriteFile << dW << " " << dV << "\n";
            }
			    }
          myWriteFile.close(); 
        }
		  }
		  //-- close the file
		  myReadFile.close();
		  nCurrent = nCurrent + nSteps;
		} 
		else {
			nCurrent = NFe + nSteps;
		}

	}	
	//--


   return 0;
}
