#include <TChain.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TString.h>
#include <TObjArray.h>

#include <iostream>
#include <vector>
#include <algorithm>

// Data
// To run: auto *ch_data = hipochain("spring2019", "pDVCS");   // uses sidisdvcs_*.hipo from production cache
// MC both energies
// To run: auto *ch_mc = hipochain("spring2019_mc", "", true); // uses all *.hipo from 10362+10363
// MC 10.6 only
// To run: auto *ch_mc_10_6 = hipochain("spring2019_mc_10_6", "", true); // uses all *.hipo from 10362 only
// MC 10.2 only
// To run: auto *ch_mc_10_2 = hipochain("spring2019_mc_10_2", "", true); // uses all *.hipo from 10363 only

// ------------------------
// Data directory resolver
// ------------------------
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

// ------------------------
// MC directory resolver
// ------------------------
static std::vector<std::string> mc_dirs_for_run(const std::string &runs)
{
    // Spring2019 MC:
    if (runs == "spring2019_mc" || runs == "spring2019-mc" || runs == "spring2019MC")
    {
        return {
            "/w/hallb-scshelf2102/clas12/nlbucuru/data_dvcs/mc_dvcs/10362", // 10.6 set 
            "/w/hallb-scshelf2102/clas12/nlbucuru/data_dvcs/mc_dvcs/10363"  // 10.2 set
        };
    }
    if (runs == "spring2019_mc_10_6" || runs == "spring2019_mc_10.6" || runs == "spring2019_10_6_mc")
    {
        return {"/w/hallb-scshelf2102/clas12/nlbucuru/data_dvcs/mc_dvcs/10362"};
    }
    if (runs == "spring2019_mc_10_2" || runs == "spring2019_mc_10.2" || runs == "spring2019_10_2_mc")
    {
        return {"/w/hallb-scshelf2102/clas12/nlbucuru/data_dvcs/mc_dvcs/10363"};
    }

    // Default: treat `runs` itself as a directory path
    return {runs};
}


// ------------------------
// Internal helper: add files from one directory
// ------------------------
static void add_files_from_dir(TChain *chain,
                               const std::string &dir_in,
                               const TString &begins_with,
                               const TString &ends_with)
{
    std::string dir = dir_in;
    if (!dir.empty() && dir.back() != '/')
        dir += "/";

    TSystemDirectory sysdir("input_dir", dir.c_str());
    TList *list = sysdir.GetListOfFiles();
    if (!list)
    {
        std::cerr << "ERROR: Could not list directory: " << dir << std::endl;
        return;
    }

    std::vector<TString> to_add;
    TIter next(list);
    while (TObject *obj = next())
    {
        auto *f = dynamic_cast<TSystemFile *>(obj);
        if (!f)
            continue;

        if (f->IsDirectory())
            continue;

        TString name = f->GetName();

        if (!begins_with.IsNull() && begins_with.Length() > 0 && !name.BeginsWith(begins_with))
            continue;
        if (!ends_with.IsNull() && ends_with.Length() > 0 && !name.EndsWith(ends_with))
            continue;

        TString full = TString(dir.c_str()) + name;
        to_add.push_back(full);
    }

    std::sort(to_add.begin(), to_add.end());
    for (const auto &full : to_add)
        chain->Add(full);
}

// ------------------------
// Main entry point
// ------------------------
// Backward-compatible: hipochain("spring2019","whatever") still works for DATA.
// New usage: hipochain("spring2019_mc","",true) for MC (chains all *.hipo in 10362+10363).
TChain *hipochain(std::string runs,
                  std::string chains,
                  bool is_mc = false,
                  std::string file_prefix = "")
{
    (void)chains; // keep argument for compatibility

    auto *chain = new TChain("hipo");

    if (!is_mc)
    {
        // DATA behavior (same as before): only sidisdvcs_*.hipo
        std::string dir = sidisdvcs_dir_for_run(runs);
        TString prefix = file_prefix.empty() ? "sidisdvcs_" : TString(file_prefix.c_str());
        add_files_from_dir(chain, dir, prefix, ".hipo");

        const int nfiles = chain->GetListOfFiles() ? chain->GetListOfFiles()->GetEntries() : 0;
        std::cout << "Added " << nfiles << " DATA files from " << dir << std::endl;
        return chain;
    }

    // MC behavior: chain all *.hipo (optionally restrict by prefix)
    const auto dirs = mc_dirs_for_run(runs);
    TString prefix = file_prefix.empty() ? "" : TString(file_prefix.c_str());

    for (const auto &d : dirs)
        add_files_from_dir(chain, d, prefix, ".hipo");

    const int nfiles = chain->GetListOfFiles() ? chain->GetListOfFiles()->GetEntries() : 0;
    std::cout << "Added " << nfiles << " MC files from ";
    for (size_t i = 0; i < dirs.size(); i++)
        std::cout << (i ? ", " : "") << dirs[i];
    std::cout << std::endl;

    return chain;
}
