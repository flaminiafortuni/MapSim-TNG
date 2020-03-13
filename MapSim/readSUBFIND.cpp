#include "readSUBFIND.h"

// long replaced with int32_t

double DoubleSwap(double d)
{
  union
  {
    double d;
    unsigned char b[8];
  } dat1, dat2;
  dat1.d = d;
  dat2.b[0] = dat1.b[7];
  dat2.b[1] = dat1.b[6];
  dat2.b[2] = dat1.b[5];
  dat2.b[3] = dat1.b[4];
  dat2.b[4] = dat1.b[3];
  dat2.b[5] = dat1.b[2];
  dat2.b[6] = dat1.b[1];
  dat2.b[7] = dat1.b[0];
  return dat2.d;
}

int IntSwap (int i)
{
  unsigned char b1, b2, b3;
  b1 = i & 255;
  b2 = (i>>8) & 255;
  b3 = (i>>16) & 255;
  return ((int)b1<<16) + ((int)b2<<8) + b3;
}

int32_t LongSwap (int32_t i)
{
  unsigned char b1, b2, b3, b4;
  b1 = i & 255;
  b2 = (i>>8) & 255;
  b3 = (i>>16) & 255;
  b4 = (i>>24) & 255;
  return ((int)b1<<24) + ((int)b2<<16) + ((int)b3<<8) + b4;
}

long long LongLongSwap (long long i)
{
  unsigned char b1, b2, b3, b4, b5, b6, b7, b8;
  b1 = i & 255;
  b2 = (i>>8) & 255;
  b3 = (i>>16) & 255;
  b4 = (i>>24) & 255;
  b5 = (i>>32) & 255;
  b6 = (i>>40) & 255;
  b7 = (i>>48) & 255;
  b8 = (i>>56) & 255;
  return ((long long)b1<<56) + ((long long)b2<<48) + ((long long)b3<<40) + ((long long)b4<<32) + ((long long)b5<<24) + ((long long)b6<<16) + ((long long)b7<<8) + b8;
}

float FloatSwap(float f)
{
  union
  {
    float f;
    unsigned char b[4];
  } dat1, dat2;
  dat1.f = f;
  dat2.b[0] = dat1.b[3];
  dat2.b[1] = dat1.b[2];
  dat2.b[2] = dat1.b[1];
  dat2.b[3] = dat1.b[0];
  return dat2.f;
}

istream & operator>>(istream &input, DATA_subhalo &Data) { 
  input.read((char *)&Data, sizeof(Data));
#ifdef SWAP
  Data.Ngroups = LongSwap(Data.Ngroups);
  Data.TotNgroups = LongSwap(Data.TotNgroups);
  Data.Nids = LongSwap(Data.Nids);
  Data.TotNids = LongLongSwap(Data.TotNids);
  Data.NFiles = LongSwap(Data.NFiles);
  Data.Nsubhalos = LongSwap(Data.Nsubhalos);
  Data.TotNsubhalos = LongSwap(Data.TotNsubhalos);
#endif
  return input;
}

void read_SUBFIND (string &file,
		   vector<double> &mfof,
		   vector<double> &xfof, vector<double> &yfof, vector<double> &zfof,
		   vector<long> &firstsubinfof,
		   vector<long> &nsubinfof,
		   vector<double> &m200, vector<double> &r200,
		   vector<double> &msub,
		   vector<double> &xsub, vector<double> &ysub, vector<double> &zsub,
		   vector<long> &idPH, vector<double> &velDisp, vector<double> &vmax,
		   vector<double> &halfmassradius,vector<double> &rmax, vector<long> &groupN
		   ){
  

  // ckeck if it contains massive neutrinos particles, look for eV in the file name string
  string seV = "eV";
  bool neutrinos = false;
  if (file.find(seV) != std::string::npos) {
    cout << "found ... " << seV << " in the file name ... " << endl;
    neutrinos = true;
  }

  vector<double> mvir,rvir;
  vector<double> masstab0,masstab1,masstab2,masstab3,masstab4,masstab5;
  vector<double> velDisp200m,velDisp200c,velDispvir;
  string file_in = file+".0";
  ifstream fin (file_in.c_str());
  if (!fin) {cerr <<"Error in opening the input file: "<<file_in<<"!\n\a";exit(1);}
  
  DATA_subhalo data; fin>>data; 
  fin.clear(); fin.close();
  
  data.Ngroups = data.TotNgroups; 
  data.Nsubhalos = data.TotNsubhalos;

  cout << "  " << endl;
  cout << " total number of groups = " << data.TotNgroups << endl;
  cout << " total number of subhaloes = " << data.TotNsubhalos << endl;
  cout << " catalogues divided in " << data.NFiles << " files" << endl;
  cout << " " << endl;

  float num_float; int32_t num_long; int num_int;
  int64_t num_long_long;
  vector<double> vxsub, vysub, vzsub;
  // haloes
  vector<int> npartinfof;
  
  int n=0;
  int ns=0;
  if (data.Nsubhalos>0) { 
    int offset_fof=0;
    int offset_sub=0;
    int offset_ids=0;
    for (int ifile=0; ifile<data.NFiles; ifile++) { // loop on the input files 
      
      file_in = file+"."+conv(ifile,fINT);
      ifstream fin (file_in.c_str());
      if (!fin) {cerr <<"Error in opening the input file: "<<file_in<<"!\n\a";exit(1);}
      
      DATA_subhalo datai; fin>>datai; 
      
      // cout << " total number of groups in the file = " << datai.Ngroups << endl;
      n+=datai.Ngroups;
      if (datai.Ngroups>0) {
	// Group Len
	// cout << " haloes " << endl;
	for (int i=0; i<datai.Ngroups; i++){
	  fin.read((char *)&num_long, sizeof(num_long)); 
#ifdef SWAP
	  num_long = LongSwap(num_long);
#else
	  num_long = num_long;
#endif
	  npartinfof.push_back(num_long);
	}
	// Offset
	for (int i=0; i<datai.Ngroups; i++){
	  fin.read((char *)&num_long, sizeof(num_long)); 
#ifdef SWAP
	  num_long = LongSwap(num_long);
#else
	  num_long = num_long;
#endif
	}
	// FOF mass
	for (int i=0; i<datai.Ngroups; i++) {
	  fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	  num_float = FloatSwap(num_float);
#else
	  num_float = num_float;
#endif
	  mfof.push_back(num_float*1.0e+10);
	}
	// FOF position
	for (int i=0; i<datai.Ngroups; i++){
	  fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	  num_float = FloatSwap(num_float);
#else
	  num_float = num_float;
#endif
	  xfof.push_back(num_float);
	  fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	  num_float = FloatSwap(num_float);
#else
	  num_float = num_float;
#endif
	  yfof.push_back(num_float);
	  fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	  num_float = FloatSwap(num_float);
#else
	  num_float = num_float;
#endif
	  zfof.push_back(num_float);
	}
	// mass 200 the mean density
	for (int i=0; i<datai.Ngroups; i++){
	  fin.read((char *)&num_float, sizeof(num_float));
#ifdef SWAP
	  num_float = FloatSwap(num_float);
#else
	  num_float = num_float;
#endif
	}
	// radius 200 the mean density
	for (int i=0; i<datai.Ngroups; i++){
	  fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	  num_float = FloatSwap(num_float);
#else
	  num_float = num_float;
#endif
	}
	// mass 200 the critical density
	for (int i=0; i<datai.Ngroups; i++){
	  fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	  num_float = FloatSwap(num_float);
#else
	  num_float = num_float;
#endif
	  m200.push_back(num_float*1.0e+10);
	}
	// radius 200 the critical density
	for (int i=0; i<datai.Ngroups; i++){ 
	  fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	  num_float = FloatSwap(num_float);
#else
	  num_float = num_float;
#endif
	  r200.push_back(num_float);
	}
	// Mvir
	for (int i=0; i<datai.Ngroups; i++){
	  fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	  num_float = FloatSwap(num_float);
#else
	  num_float = num_float;
#endif
	  mvir.push_back(num_float*1.0e+10);
	}
	// Rvir
	for (int i=0; i<datai.Ngroups; i++){
	  fin.read((char *)&num_float, sizeof(num_float));
#ifdef SWAP
	  num_float = FloatSwap(num_float);
#else
	  num_float = num_float;
#endif
	  rvir.push_back(num_float);
	}

	// disp.vel 200m
        for (int i=0; i<datai.Ngroups; i++){
	  fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	  num_float = FloatSwap(num_float);
#else
	  num_float = num_float;
#endif
	  velDisp200m.push_back(num_float);
	}

	// disp.vel 200c
        for (int i=0; i<datai.Ngroups; i++){
	  fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	  num_float = FloatSwap(num_float);
#else
	  num_float = num_float;
#endif
	  velDisp200c.push_back(num_float);
	}

	// disp.vel vir
        for (int i=0; i<datai.Ngroups; i++){
	  fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	  num_float = FloatSwap(num_float);
#else
	  num_float = num_float;
#endif
	  velDispvir.push_back(num_float);
	}

	// contamination - float
	for (int i=0; i<datai.Ngroups; i++) {
	  fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	  num_float = FloatSwap(num_float);
#else
	  num_float = num_float;
#endif
	  //contd.push_back(num_float);
	}

	// contamination - long
	for (int i=0; i<datai.Ngroups; i++) {
	  fin.read((char *)&num_long, sizeof(num_long)); 
#ifdef SWAP
	  num_long = LongSwap(num_long);
#else
	  num_long = num_long;
#endif
	  //contl.push_back(num_long);
	}

	// N sub in FOF
	for (int i=0; i<datai.Ngroups; i++){   
	  fin.read((char *)&num_long, sizeof(num_long)); 
#ifdef SWAP
	  num_long = LongSwap(num_long);
#else
	  num_long = num_long;
#endif
	  nsubinfof.push_back(num_long);
	}
	// ID First Subhalo of the Halo
	for (int i=0; i<datai.Ngroups; i++){
	  fin.read((char *)&num_long, sizeof(num_long)); 
#ifdef SWAP
	  num_long = LongSwap(num_long);
#else
	  num_long = num_long;
#endif
	  firstsubinfof.push_back(num_long+offset_sub);
	}
      }
      //if(offset_fof==0){
      // for(int ii=0;ii<10;ii++){
      // int i = ii + offset_fof;
      // cout << i << "  " << nsubinfof[i] << "  " << mfof[i] << "  " << mvir[i] << "  " << rvir[i] << "  " << m200[i] << "  " << r200[i] 
      // << " | " << velDisp200m[i] << "  " << velDisp200c[i] << "  " << velDispvir[i] 
      // << " |  " << firstsubinfof[i] << "  " << xfof[i] << "  " << yfof[i] << "  " << zfof[i] << endl;
      // }
      // }
      ns+=datai.Nsubhalos;
      // cout << " subhaloes " << endl;
      if (datai.Nsubhalos) {
	// Subhalo Len
	for (int i=0; i<datai.Nsubhalos; i++){
	  fin.read((char *)&num_long, sizeof(num_long)); 
#ifdef SWAP
	  num_long = LongSwap(num_long);
#else
	  num_long = num_long;
#endif
	}
	// offset	
	for (int i=0; i<datai.Nsubhalos; i++){
	  fin.read((char *)&num_long, sizeof(num_long)); 
#ifdef SWAP
	  num_long = LongSwap(num_long);
#else
	  num_long = num_long;
#endif
	}
	// ID parent Halo
	for (int i=0; i<datai.Nsubhalos; i++){
	  fin.read((char *)&num_long, sizeof(num_long)); 
#ifdef SWAP
	  num_long = LongSwap(num_long);
#else
	  num_long = num_long;
#endif
	  idPH.push_back(num_long+offset_fof);
	}
	// total subhalo mass
	for (int i=0; i<datai.Nsubhalos; i++) {
	  fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	  num_float = FloatSwap(num_float);
#else
	  num_float = num_float;
#endif
	  msub.push_back(num_float*1.0e+10);  
	}
	// position in kcp/h 
	for (int i=0; i<datai.Nsubhalos; i++) {
	  fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	  num_float = FloatSwap(num_float); 
#else
	  num_float = num_float;
#endif
	  xsub.push_back(num_float); 
	  fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	  num_float = FloatSwap(num_float);
#else
	  num_float = num_float;
#endif
	  ysub.push_back(num_float);
	  fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	  num_float = FloatSwap(num_float);
#else
	  num_float = num_float;
#endif
	  zsub.push_back(num_float);
	}
	// velocity 
	for (int i=0; i<datai.Nsubhalos; i++) {
	  fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	  num_float = FloatSwap(num_float);
#else
	  num_float = num_float;
#endif
	  vxsub.push_back(num_float);
	  fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	  num_float = FloatSwap(num_float);
#else
	  num_float = num_float;
#endif
	  vysub.push_back(num_float);
	  fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	  num_float = FloatSwap(num_float);
#else
	  num_float = num_float;
#endif
	  vzsub.push_back(num_float);
	}
	// CM location index
	for (int i=0; i<datai.Nsubhalos; i++) 
	  for (int j=0; j<3; j++) 
	    fin.read((char *)&num_float, sizeof(num_float));
	// Spin
	for (int i=0; i<datai.Nsubhalos; i++) 
	  for (int j=0; j<3; j++) 
	    fin.read((char *)&num_float, sizeof(num_float));
	// Velocity Dispersion
	for (int i=0; i<datai.Nsubhalos; i++){
	  fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	  num_float = FloatSwap(num_float);
#else
	  num_float = num_float;
#endif
	  velDisp.push_back(num_float);
	}
	// Vmax
	for (int i=0; i<datai.Nsubhalos; i++){
	  fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	  num_float = FloatSwap(num_float);
#else
	  num_float = num_float;
#endif
	  vmax.push_back(num_float);
	}
	// Rmax of Vmax
	for (int i=0; i<datai.Nsubhalos; i++){
	  fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	  num_float = FloatSwap(num_float);
#else
	  num_float = num_float;
#endif
	  rmax.push_back(num_float);
	}
	// half mass radius
	for (int i=0; i<datai.Nsubhalos; i++){
	  fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	  num_float = FloatSwap(num_float);
#else
	  num_float = num_float;
#endif
	  halfmassradius.push_back(num_float);
	}
	// most bound ID
	for (int i=0; i<datai.Nsubhalos; i++)
	  //if(neutrinos) fin.read((char *)&num_long_long, sizeof(num_long_long)); 
	  //	  else 
	  fin.read((char *)&num_long, sizeof(num_long)); 
	// group number
	for (int i=0; i<datai.Nsubhalos; i++){
	  fin.read((char *)&num_long, sizeof(num_long)); 
#ifdef SWAP
	  num_long = LongSwap(num_long);
#else
	  num_long = num_long;
#endif
	  groupN.push_back(num_long);
	}
      }

      // mass tab 6 elements
      for (int i=0; i<datai.Nsubhalos; i++) {
	fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	num_float = FloatSwap(num_float); 
#else
	num_float = num_float;
#endif
	masstab0.push_back(num_float); 
	fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	num_float = FloatSwap(num_float);
#else
	num_float = num_float;
#endif
	masstab1.push_back(num_float);
	fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	num_float = FloatSwap(num_float);
#else
	num_float = num_float;
#endif
	masstab2.push_back(num_float);
	fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	num_float = FloatSwap(num_float); 
#else
	num_float = num_float;
#endif
	masstab3.push_back(num_float); 
	fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	num_float = FloatSwap(num_float);
#else
	num_float = num_float;
#endif
	masstab4.push_back(num_float);
	fin.read((char *)&num_float, sizeof(num_float)); 
#ifdef SWAP
	num_float = FloatSwap(num_float);
#else
	num_float = num_float;
#endif
	masstab5.push_back(num_float);
      }
      if(offset_sub==0){
	std:: cout << masstab0[0] << "   " <<  masstab1[0] << "   " <<  masstab2[0] << "   " << 
	  masstab3[0] << "   " <<  masstab4[0] << "   " <<  masstab5[0] << std:: endl;
	std:: cout << masstab0[1] << "   " <<  masstab1[1] << "   " <<  masstab2[1] << "   " << 
	  masstab3[1] << "   " <<  masstab4[1] << "   " <<  masstab5[1] << std:: endl;
	std:: cout << masstab0[2] << "   " <<  masstab1[2] << "   " <<  masstab2[2] << "   " << 
	  masstab3[2] << "   " <<  masstab4[2] << "   " <<  masstab5[2] << std:: endl;
	std:: cout << nsubinfof[0] << "    " << nsubinfof[1] << std:: endl;
	for(int ii=0;ii<128;ii++){
	  int i = ii + offset_sub;
	  cout << i << "  " << idPH[i] << "  " << msub[i] << "  " 
	       << xsub[i] << "  " << ysub[i] << "  " << zsub[i] << "  " << groupN[i] << "   " << halfmassradius[i] << endl;
	}
      }
      fin.clear(); fin.close();
      offset_fof+=datai.Ngroups;
      offset_sub+=datai.Nsubhalos;
      offset_ids+=datai.Nids;
    }
    cout << " TOTAL number of Groups read " << n << endl;
    cout << " TOTAL number of Subhaloes read " << ns << endl;
  }
}


