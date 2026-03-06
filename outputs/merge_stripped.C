// merge_stripped.C
// Merge stripped ROOT outputs in outputs/ using TFileMerger.
// Usage examples (from outputs/):
//   clas12root -l -b -q 'merge_stripped.C("spring2019","pDVCS")'
//   clas12root -l -b -q 'merge_stripped.C("fall2019","pDVCS")'
//   clas12root -l -b -q 'merge_stripped.C("spring2020","pDVCS")'
//   clas12root -l -b -q 'merge_stripped.C("spring2019_mc_10_6","pDVCS")'
//
// Dry-run:
//   clas12root -l -b -q 'merge_stripped.C("spring2019","pDVCS",".","",true)'

#include <TFileMerger.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TString.h>
#include <TPRegexp.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TSystem.h>

#include <iostream>
#include <vector>
#include <algorithm>

static int ExtractTrailingIndex(const TString& name) {
  // Extract trailing number in "...__12.root" or "..._12.root"
  TPRegexp re("_+([0-9]+)\\.root$");
  TObjArray* m = re.MatchS(name);
  int idx = -1;
  if (m && m->GetLast() >= 1) {
    TString s = ((TObjString*)m->At(1))->GetString();
    idx = s.Atoi();
  }
  if (m) delete m;
  return idx;
}

static std::vector<TString> ListMatchingFiles(const char* dir, const TString& regexPattern) {
  std::vector<TString> matches;
  TPRegexp re(regexPattern);

  TSystemDirectory d("outputs_dir", dir);
  TList* list = d.GetListOfFiles();
  if (!list) {
    std::cerr << "ERROR: Could not list directory: " << dir << "\n";
    return matches;
  }

  TIter it(list);
  while (TObject* obj = it()) {
    auto* f = dynamic_cast<TSystemFile*>(obj);
    if (!f) continue;
    if (f->IsDirectory()) continue;

    TString name = f->GetName();
    if (!name.EndsWith(".root")) continue;
    if (!re.MatchB(name)) continue;

    matches.push_back(TString::Format("%s/%s", dir, name.Data()));
  }

  std::sort(matches.begin(), matches.end(),
    [](const TString& a, const TString& b){
      TString an = gSystem->BaseName(a);
      TString bn = gSystem->BaseName(b);
      int ia = ExtractTrailingIndex(an);
      int ib = ExtractTrailingIndex(bn);
      if (ia >= 0 && ib >= 0 && ia != ib) return ia < ib;
      return a < b;
    });

  return matches;
}

void merge_stripped(const char* run,
                    const char* channel = "pDVCS",
                    const char* dir = ".",
                    const char* outFile = "",
                    bool dryRun = false,
                    bool fastMethod = true)
{
  TString RUN(run), CH(channel);

  // Decide output filename
  TString outName;
  if (outFile && TString(outFile).Length() > 0) outName = outFile;
  else outName = TString::Format("stripped_data_%s_%s_merged.root", RUN.Data(), CH.Data());

  TString outPath = outName;
  if (!outPath.Contains("/")) outPath = TString::Format("%s/%s", dir, outName.Data());

  // Prefer "__N.root" chunks, then "_N.root" chunks, else single file.
  TString patt_double = TString::Format("^stripped_data_%s_%s__\\d+\\.root$", RUN.Data(), CH.Data());
  TString patt_single = TString::Format("^stripped_data_%s_%s_\\d+\\.root$",  RUN.Data(), CH.Data());
  TString patt_base   = TString::Format("^stripped_data_%s_%s\\.root$",        RUN.Data(), CH.Data());

  auto inputs = ListMatchingFiles(dir, patt_double);
  TString used = "__N chunks";
  if (inputs.empty()) {
    inputs = ListMatchingFiles(dir, patt_single);
    used = "_N chunks";
  }
  if (inputs.empty()) {
    inputs = ListMatchingFiles(dir, patt_base);
    used = "single file";
  }

  std::cout << "Run            : " << RUN.Data() << "\n";
  std::cout << "Channel        : " << CH.Data() << "\n";
  std::cout << "Directory      : " << dir << "\n";
  std::cout << "Selection      : " << used.Data() << "\n";
  std::cout << "Output         : " << outPath.Data() << "\n";
  std::cout << "Dry run        : " << (dryRun ? "YES" : "NO") << "\n\n";

  if (inputs.empty()) {
    std::cerr << "ERROR: No files found for this run/channel.\n";
    std::cerr << "Tried patterns:\n"
              << "  " << patt_double.Data() << "\n"
              << "  " << patt_single.Data() << "\n"
              << "  " << patt_base.Data()   << "\n";
    return;
  }

  // Never merge the output into itself
  std::vector<TString> filtered;
  filtered.reserve(inputs.size());
  for (auto& f : inputs) {
    if (f != outPath) filtered.push_back(f);
  }
  inputs.swap(filtered);

  std::cout << "Input files (" << inputs.size() << "):\n";
  for (auto& f : inputs) std::cout << "  " << f.Data() << "\n";
  std::cout << "\n";

  if (dryRun) return;

  TFileMerger merger;
  merger.SetFastMethod(fastMethod);
  merger.OutputFile(outPath.Data(), "RECREATE");

  for (auto& f : inputs) merger.AddFile(f.Data());

  if (!merger.Merge()) {
    std::cerr << "ERROR: merge failed.\n";
    return;
  }

  std::cout << "DONE: wrote " << outPath.Data() << "\n";
}