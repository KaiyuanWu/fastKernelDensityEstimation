void allow(){
	/*
	for(int m=20;m<30;m++){
	cout<<"+++++++++++++++ m= "<<m<<" +++++++++++++++++++++"<<endl;
	int total=0;
	for(int i1=1;i1<30;i1++){
		for(int i2=1;i2<30;i2++){
			for(int i3=1;i3<30;i3++){
				int k=i1*i2*i3;
				if(k<m){
					cout<<"("<<i1<<", "<<i2<<", "<<i3<<"), ";
					total++;
				}
			}
		}
	}	
	cout<<"\ntotal= "<<total<<endl;
	}
	*/
	TFile* fin=TFile::Open("parameterDataBase.root");
	int beta=0;
	TGraph* gr=new TGraph();
	double v;
	double maxv;
	TTree* tree;
	int c=1;
	for(int beta=0;beta<20;beta++){
	for(int nx=1;nx<400;nx++){
		double x=0.01*nx;
		tree=(TTree*) fin->Get(Form("p%dB%d",29,beta));
		tree->SetBranchAddress("v",&v);
		tree->GetEntry(nx);
		maxv=v;
		for(int ib=2;ib<=30;ib++){
                        tree=(TTree*) fin->Get(Form("p%dB%d",ib,beta));
			if(!tree){
				cout<<"Can not find "<<Form("p%dB%d",ib,beta)<<endl;
				return;
			}
                        tree->SetBranchAddress("v",&v);
			tree->GetEntry(nx);
			if(v>=maxv-0.001){
				//cout<<"("<<x<<", "<<ib<<", v="<<v<<", maxv="<<maxv<<") "<<endl;
				printf("%2d, ",ib);
				gr->SetPoint(nx-1,x,ib);
				break;
			}
		}
		if(c%40==0)
			printf("\n");
		c++;
		//cout<<maxv<<" "<<x<<endl;
		}
	}
	gr->Draw("AL");
}
