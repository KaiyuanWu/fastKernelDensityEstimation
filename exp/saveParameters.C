void saveParameters(){
	TFile* fin=TFile::Open("parameterDataBase.root");
	int n=400;
	double v;
	FILE* fout=fopen("initParameters.h","w");
	fprintf(fout,"///////////////////////// g(beta,p) //////////////////////\n");
	fprintf(fout,"const double paramBP[]={\n");
	for(int iP=2;iP<29;iP++){
		for(int iB=0;iB<20;iB++){
			TTree* tree=(TTree*) fin->Get(Form("p%dB%d",iP,iB));
			int entries=tree->GetEntries();
			tree->SetBranchAddress("v",&v);
			for(int iE=0;iE<n;iE++){
				tree->GetEntry(iE);
				if(iE%10==0&&iE!=0)
					fprintf(fout,"\n");
				if(iE!=n-1)
					fprintf(fout,"%.2e, ",v);
				else
					fprintf(fout,"%.2e",v);
			}
			if(iB==19&&iP==28)
				fprintf(fout,"};\n");
			else
				fprintf(fout,",\n");
		}
	}
	fprintf(fout,"///////////////////////// g(beta) //////////////////////\n");
	fprintf(fout,"const double paramB[]={\n");
	for(int iB=0;iB<20;iB++){
        	TTree* tree=(TTree*) fin->Get(Form("pB%d",iB));
                int entries=tree->GetEntries();
                tree->SetBranchAddress("v",&v);
                for(int iE=0;iE<n;iE++){
                tree->GetEntry(iE);
                if(iE%10==0&&iE!=0)
                	fprintf(fout,"\n");
                if(iE!=n-1)
                	fprintf(fout,"%.2e, ",v);
                else
                        fprintf(fout,"%.2e",v);
                }
		if(iB==19)
                	fprintf(fout,"};\n");
		else
			fprintf(fout,",\n");
	}
	fclose(fout);

	/*
	double maxv;
	int maxi;
	double v1,v2;
                for(int iB=0;iB<19;iB++){
			maxv=-1.0e100;
			maxi=0;
			TTree* tree1=(TTree*) fin->Get(Form("p%dB%d",20,iB));
			TTree* tree2=(TTree*) fin->Get(Form("pB%d",iB));
			tree1->SetBranchAddress("v",&v1);
			tree2->SetBranchAddress("v",&v2);
			int entries=tree1->GetEntries();
			for(int iE=0;iE<n;iE++){
				tree1->GetEntry(iE);
				tree2->GetEntry(iE);
				if(fabs(v1-v2)>maxv){
					maxv=fabs(v1-v2);
					maxi=iE;
				}
			}
			cout<<"iB= "<<iB<<", maxv= "<<maxv<<", maxi= "<<maxi<<endl;
		}
	*/
}
