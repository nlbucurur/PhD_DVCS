// clas12root -l -b -q run_all.C
#include <TROOT.h>

void run_all()
{
    gROOT->ProcessLine(".L analysis_DVCS_preskimmed_fiducials_NID.C+");
    gROOT->ProcessLine(".L run_DVCS.C+");
    gROOT->ProcessLine("run_DVCS(\"pDVCS\",\"spring2019\");");
    gROOT->ProcessLine("run_DVCS(\"pDVCS\",\"fall2019\");");
    gROOT->ProcessLine("run_DVCS(\"pDVCS\",\"spring2020\");");
}
