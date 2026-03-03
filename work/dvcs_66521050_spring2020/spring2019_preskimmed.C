#include <TChain.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TString.h>
#include <TObjArray.h>

#include <iostream>
#include <vector>
#include <algorithm>

static std::string sidisdvcs_dir_for_run(const std::string &runs)
{
    if (runs == "spring2019")
        return "/lustre24/expphy/cache/clas12/rg-b/production/recon/spring2019/torus-1/pass2/v0/dst/train/sidisdvcs";
    if (runs == "fall2019")
        return "/lustre24/expphy/cache/clas12/rg-b/production/recon/fall2019/torus+1/pass2/v1/dst/train/sidisdvcs";
    if (runs == "spring2020")
        return "/lustre24/expphy/cache/clas12/rg-b/production/recon/spring2020/torus-1/pass2/v1/dst/train/sidisdvcs";

    // Default: treat `runs` itself as a directory path
    return runs;
}

TChain *hipochain(std::string runs, std::string chains)
{
    (void)chains; // same train for pDVCS/nDVCS; keep argument for compatibility

    std::string dir = sidisdvcs_dir_for_run(runs);
    if (!dir.empty() && dir.back() != '/')
        dir += "/";

    auto *chain = new TChain("hipo");

    TSystemDirectory sysdir("sidisdvcs_dir", dir.c_str());
    TList *list = sysdir.GetListOfFiles();
    if (!list)
    {
        std::cerr << "ERROR: Could not list directory: " << dir << std::endl;
        return chain;
    }

    std::vector<TString> to_add;

    TIter next(list);
    while (TObject *obj = next())
    {
        auto *f = dynamic_cast<TSystemFile *>(obj);
        if (!f)
            continue;

        TString name = f->GetName();
        if (f->IsDirectory())
            continue;

        // keep only sidisdvcs_*.hipo
        if (!name.BeginsWith("sidisdvcs_"))
            continue;
        if (!name.EndsWith(".hipo"))
            continue;

        TString full = TString(dir.c_str()) + name;
        to_add.push_back(full);
    }

    std::sort(to_add.begin(), to_add.end());
    for (const auto &full : to_add)
        chain->Add(full);

    const int nfiles = chain->GetListOfFiles() ? chain->GetListOfFiles()->GetEntries() : 0;
    std::cout << "Added " << nfiles << " files from " << dir << std::endl;

    return chain;
}
