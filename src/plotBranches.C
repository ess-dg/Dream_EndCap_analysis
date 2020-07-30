/// Comment about uppercase C 

void plotBranches() {
  // the are ROOT command lines to plot various branches

  //		cdt_new->AddFriend("cdt_new_cal","cdt_new_cal.root")
  //      cdt_new->Scan("ntime:nchopperTime:tdiff:ccycle:ntof","")
  //		cdt_new->Draw("tdiff/10000>>tt(800,0,80)","","")
  //		cdt_new->Draw("ntof/10000>>tt(900,-10,80)","","")
  //		cdt_new->Draw("ntof_wfm/10000>>tt(900,-10,80)","","")
  //		cdt_new->Draw("lambda_c>>tt(1000,0,10)","","")
  //		cdt_new->Draw("lambda_c_wfm_corr>>tt(1000,0,10)","","")
  //		cdt_new->Draw("nvoxel_y:nvoxel_z","nsumo<=4","")
  //		cdt_new->Draw("dspacing_c>>tt(600,0,6)","nsumo==6","")
  //		cdt_new->Draw("twotheta_c>>tt(350,135,170)","","")
  //		cdt_new->Draw("nvoxel_y:nvoxel_z>>tt(500,-1600,-1100,1000,200,1200)","","colz")
  //		cdt_new->Draw("nvoxel_y:nvoxel_x>>tt(300,-100,200,1000,200,1200)","","colz")

  // *************and some automatic plotting

	TCanvas *can=new TCanvas("can","can",100,100,700,700);

	gStyle->SetOptTitle(1);
	gStyle->SetOptStat(1);

	can->SetFillColor(0);
	can->SetGrid();
	Float_t small=1e-5;

	can->Divide(2,2,small,small);
	can->ToggleEventStatus();

	can->cd(1);
	cdt_new->Draw("ntof_wfm/10000>>tt(1000,-20,80)","","");
	can->cd(2);
	cdt_new->Draw("ntof_wfm_corr/10000>>tt1(1000,-20,80)","","");
	can->cd(3);
	cdt_new->Draw("tdiff/10000>>tt2(1200,-20,100)","","");
	can->cd(4);
	cdt_new->Draw("ntof/10000>>tt3(1200,-20,100)","","");

  //	newt->Draw("ncathode:nanode","nsumo==4","colz");
  //	can->cd(3);
  //	newt->Draw("ncathode:nanode","nsumo==5","colz");
  //	can->cd(4);
  //	newt->Draw("ncathode:nanode","nsumo==6","colz");

  /*	TCanvas *canT=new TCanvas("canT","canT",100,100,800,800);

          gStyle->SetOptTitle(1);
          gStyle->SetOptStat(1);

          canT->SetFillColor(0);
          canT->SetGrid();

          canT->Divide(2,2,small,small);

          canT->cd(1);
          newt->Draw("(ntime-nchopperTime)/10e6>>tt(2000,-0.1,0.1)","nboardID==1418045");
          canT->cd(2);
          newt->Draw("(ntime-nchopperTime)/10e6","nboardID==1416964");
          canT->cd(3);
          newt->Draw("(ntime-nchopperTime)/10e6","nboardID==1416799");
          canT->cd(4);
          newt->Draw("(ntime-nchopperTime)/10e6","nboardID==1416697");*/
}
