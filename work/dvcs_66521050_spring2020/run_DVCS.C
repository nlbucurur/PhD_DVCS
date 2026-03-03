#ifdef analysis_cxx
#undef analysis_cxx
#endif

#include "analysis_DVCS_preskimmed_fiducials_NID.h"
// #include "spring2019_preskimmed.C"

void run_DVCS(TString channel = "pDVCS", TString period = "spring2019")
{
    analysis a(channel.Data(), period.Data());
    a.Loop(false);
}