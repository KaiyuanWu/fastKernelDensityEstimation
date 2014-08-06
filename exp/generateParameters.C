void generateParameters(){
	double h=0.01	;
	double maxX=10	;
	int p=1		;
	int maxP=30	;
	int beta=0	;
	int maxBeta=20	;
	int maxRound=800;
	double v=0.	;
	double x	;
	double s0	;
	double s	;
	TFile* fout=TFile::Open("parameterDataBase.root","recreate");
	TTree* tree;
	for(beta=0;beta<=maxBeta;beta++){
		tree=new TTree(Form("pB%d",beta),Form("pB%d",beta));
		tree->Branch("v",&v,"v/D");
		tree->Branch("x",&x,"x/D");
		for(x=0.;x<=maxX;x+=h){
			v=1.;
			s0=1.;
			for(int i1=1;i1<maxRound;i1++){
				s=sqrt(i1+beta)/double(i1)*x;
				s0*=s;
				v+=s0;
			}
			tree->Fill();
		}
		cout<<"s0= "<<s0<<endl;
		tree->Write("",TObject::kOverwrite);
		for(int ip=p;ip<=maxP;ip++){
			tree=new TTree(Form("p%dB%d",ip,beta),Form("p%dB%d",ip,beta));
                	tree->Branch("v",&v,"v/D");	
			tree->Branch("x",&x,"x/D");
			for(x=0.;x<=maxX;x+=h){
				v=1.;
				s0=1.;
				for(int i1=1;i1<ip;i1++){
					s=sqrt(i1+beta)/double(i1)*x;
					s0*=s;
					v+=s0;
				}
				tree->Fill();
			}
			tree->Write("",TObject::kOverwrite);
			cout<<"ip= "<<ip<<endl;
		}
		cout<<"beta= "<<beta<<endl;
	}
	TNamed* name1=new TNamed("h",Form("%f",h));
	name1->Write("",TObject::kOverwrite);
	TNamed* name2=new TNamed("minP",Form("%d",p));
        name2->Write("",TObject::kOverwrite);
	TNamed* name3=new TNamed("maxP",Form("%d",maxP));
        name3->Write("",TObject::kOverwrite); 
	TNamed* name4=new TNamed("maxBeta",Form("%d",maxBeta));
        name4->Write("",TObject::kOverwrite);	
	TNamed* name5=new TNamed("maxX",Form("%f",maxX));
        name5->Write("",TObject::kOverwrite);
	fout->Close();			
}
