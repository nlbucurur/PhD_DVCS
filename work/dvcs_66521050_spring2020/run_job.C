void run_job() {
  gROOT->ProcessLine(".L analysis_DVCS_preskimmed_fiducials_NID.C");
  gROOT->ProcessLine(".L run_DVCS.C");
  gROOT->ProcessLine("run_DVCS(\"pDVCS\",\"spring2020\");");
}
run_job();
