#define analysis_cxx
#include "analysis_DVCS_preskimmed_fiducials_NID.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void analysis::Loop(bool simu = false)
{

    Target_Vec.SetXYZT(0, 0, 0, Dmass);
    PTarget_Vec.SetXYZT(0, 0, 0, Pmass);
    NTarget_Vec.SetXYZT(0, 0, 0, Nmass);

    MC = simu;

    outfile = new TFile(std::string("stripped_data_" + run + ".root").c_str(), "RECREATE");

    AddBranches();

    TObjArray *files = chain->GetListOfFiles();

    cout << "size of files   " << files->GetLast() + 1 << endl;

    for (int noffiles = 0; noffiles < files->GetLast() + 1; noffiles++)
    {

        reader.open(files->At(noffiles)->GetTitle());

        reader.readDictionary(factory);

        iCandidate = -1;
        iMinCandidate = 0;
        MinCandidate = 999999;

        hipo::bank CONF(factory.getSchema("RUN::config"));
        hipo::bank HEL(factory.getSchema("REC::Event"));
        hipo::bank PART(factory.getSchema("REC::Particle"));
        hipo::bank CAL(factory.getSchema("REC::Calorimeter"));
        hipo::bank TRAK(factory.getSchema("REC::Track"));
        hipo::bank TRAJ(factory.getSchema("REC::Traj"));
        hipo::bank SCINT(factory.getSchema("REC::Scintillator"));
        hipo::bank SCINTX(factory.getSchema("REC::ScintExtras"));

        /*if (simu) {
            hipo::bank MCPART(factory.getSchema("MC::Particle"));
        }*/

        while (reader.next() == true)
        {

            ClearVectors(3);

            reader.read(event);
            event.getStructure(PART);
            event.getStructure(HEL);
            event.getStructure(CONF);
            event.getStructure(CAL);
            event.getStructure(TRAK);
            event.getStructure(TRAJ);
            event.getStructure(SCINT);
            event.getStructure(SCINTX);

            // event.getStructure(MCPART);
            /*
            if (simu) {
                RECPMCpid   =PART.getInt("pid",1);
            }

            if (RECPMCpid != 2212) {
                continue;
            }*/

            size = PART.getRows();

            RunNumber = CONF.getInt("run", 0);
            EventNumber = CONF.getInt("event", 0);
            Helicity = HEL.getInt("helicity", 0);

            Ebeam = 10.6;

            if (RunNumber >= 6420)
                Ebeam = 10.2;

            if (RunNumber > 10000)
                Ebeam = 10.4;

            ElectronBeam.SetXYZT(0, 0, Ebeam, Ebeam);

            Cal_Nentries = 0;
            Traj_Nentries = 0;
            TRK_Nentries = 0;
            Scint_Nentries = 0;
            ScintX_Nentries = 0;

            Cal_Nentries = CAL.getRows();
            TRK_Nentries = TRAK.getRows();
            Traj_Nentries = TRAJ.getRows();
            Scint_Nentries = SCINT.getRows();
            ScintX_Nentries = SCINTX.getRows();

            // detector type in Scint bank   CND       ( 3, "CND"),     CTOF      ( 4, "CTOF"),

            if (Scint_Nentries > 0 && ScintX_Nentries > 0)
            {
                for (int i = 0; i < Scint_Nentries; i++)
                {
                    if (SCINT.getInt("pindex", i) >= 0 && SCINT.getInt("pindex", i) < size && SCINT.getInt("detector", i) == 3)
                    {
                        part_Scint_CND_energy[SCINT.getInt("pindex", i)] = SCINT.getFloat("energy", i);
                        part_ScintX_CND_size[SCINT.getInt("pindex", i)] = SCINTX.getInt("size", i);
                        part_ScintX_CND_layermult[SCINT.getInt("pindex", i)] = SCINTX.getInt("layermult", i);
                        part_ScintX_CND_dedx[SCINT.getInt("pindex", i)] = SCINTX.getFloat("dedx", i);
                        part_Scint_CND_x[SCINT.getInt("pindex", i)] = SCINT.getFloat("x", i);
                        part_Scint_CND_y[SCINT.getInt("pindex", i)] = SCINT.getFloat("y", i);
                        part_Scint_CND_z[SCINT.getInt("pindex", i)] = SCINT.getFloat("z", i);
                        part_Scint_CND_t[SCINT.getInt("pindex", i)] = SCINT.getFloat("time", i);
                    }
                    if (SCINT.getInt("pindex", i) >= 0 && SCINT.getInt("pindex", i) < size && SCINT.getInt("detector", i) == 4)
                    {
                        part_Scint_CTOF_energy[SCINT.getInt("pindex", i)] = SCINT.getFloat("energy", i);
                        part_ScintX_CTOF_size[SCINT.getInt("pindex", i)] = SCINTX.getInt("size", i);
                        part_ScintX_CTOF_layermult[SCINT.getInt("pindex", i)] = SCINTX.getInt("layermult", i);
                        part_ScintX_CTOF_dedx[SCINT.getInt("pindex", i)] = SCINTX.getFloat("dedx", i);
                        part_Scint_CTOF_x[SCINT.getInt("pindex", i)] = SCINT.getFloat("x", i);
                        part_Scint_CTOF_y[SCINT.getInt("pindex", i)] = SCINT.getFloat("y", i);
                        part_Scint_CTOF_z[SCINT.getInt("pindex", i)] = SCINT.getFloat("z", i);
                        part_Scint_CTOF_t[SCINT.getInt("pindex", i)] = SCINT.getFloat("time", i);
                    }
                }
            }

            // Calorimeter bank (layer:  PCAL = 1, ECin = 4, ECout = 7)

            if (Cal_Nentries > 0)
            {
                for (int i = 0; i < Cal_Nentries; i++)
                {
                    if (CAL.getInt("pindex", i) >= 0 && CAL.getInt("pindex", i) < size && CAL.getInt("layer", i) == 1)
                    {
                        part_Cal_PCAL_sector[CAL.getInt("pindex", i)] = CAL.getInt("sector", i);
                        part_Cal_PCAL_energy[CAL.getInt("pindex", i)] = CAL.getFloat("energy", i);
                        part_Cal_PCAL_x[CAL.getInt("pindex", i)] = CAL.getFloat("x", i);
                        part_Cal_PCAL_y[CAL.getInt("pindex", i)] = CAL.getFloat("y", i);
                        part_Cal_PCAL_z[CAL.getInt("pindex", i)] = CAL.getFloat("z", i);
                        part_Cal_PCAL_lu[CAL.getInt("pindex", i)] = CAL.getFloat("lu", i);
                        part_Cal_PCAL_lv[CAL.getInt("pindex", i)] = CAL.getFloat("lv", i);
                        part_Cal_PCAL_lw[CAL.getInt("pindex", i)] = CAL.getFloat("lw", i);
                    }
                    if (CAL.getInt("pindex", i) >= 0 && CAL.getInt("pindex", i) < size && CAL.getInt("layer", i) == 4)
                    {
                        part_Cal_ECin_sector[CAL.getInt("pindex", i)] = CAL.getInt("sector", i);
                        part_Cal_ECin_energy[CAL.getInt("pindex", i)] = CAL.getFloat("energy", i);
                        part_Cal_ECin_x[CAL.getInt("pindex", i)] = CAL.getFloat("x", i);
                        part_Cal_ECin_y[CAL.getInt("pindex", i)] = CAL.getFloat("y", i);
                        part_Cal_ECin_z[CAL.getInt("pindex", i)] = CAL.getFloat("z", i);
                        part_Cal_ECin_lu[CAL.getInt("pindex", i)] = CAL.getFloat("lu", i);
                        part_Cal_ECin_lv[CAL.getInt("pindex", i)] = CAL.getFloat("lv", i);
                        part_Cal_ECin_lw[CAL.getInt("pindex", i)] = CAL.getFloat("lw", i);
                    }
                    if (CAL.getInt("pindex", i) >= 0 && CAL.getInt("pindex", i) < size && CAL.getInt("layer", i) == 7)
                    {
                        part_Cal_ECout_sector[CAL.getInt("pindex", i)] = CAL.getInt("sector", i);
                        part_Cal_ECout_energy[CAL.getInt("pindex", i)] = CAL.getFloat("energy", i);
                        part_Cal_ECout_x[CAL.getInt("pindex", i)] = CAL.getFloat("x", i);
                        part_Cal_ECout_y[CAL.getInt("pindex", i)] = CAL.getFloat("y", i);
                        part_Cal_ECout_z[CAL.getInt("pindex", i)] = CAL.getFloat("z", i);
                        part_Cal_ECout_lu[CAL.getInt("pindex", i)] = CAL.getFloat("lu", i);
                        part_Cal_ECout_lv[CAL.getInt("pindex", i)] = CAL.getFloat("lv", i);
                        part_Cal_ECout_lw[CAL.getInt("pindex", i)] = CAL.getFloat("lw", i);
                    }
                }
            }

            for (Int_t i = 0; i < size; i++)
            {
                part_Cal_energy_total[i] = part_Cal_PCAL_energy[i] + part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i];
            }

            // tracking banks  (detectors: DC = 6, BST = 2,  BMT = 1, FMT = 8)

            if (TRK_Nentries > 0)
            {
                for (int i = 0; i < TRK_Nentries; i++)
                {
                    if (TRAK.getInt("pindex", i) < size && TRAK.getInt("detector", i) == 6)
                    { // DC
                        part_DC_Track_chi2[TRAK.getInt("pindex", i)] = TRAK.getFloat("chi2", i);
                        part_DC_Track_NDF[TRAK.getInt("pindex", i)] = TRAK.getInt("NDF", i);
                        part_DC_Track_status[TRAK.getInt("pindex", i)] = TRAK.getInt("status", i);
                    }
                }
            }

            // trajectory crosses  (layers: 6 = DC region 1 start,  18 = DC region 2 center,  36 = DC region 3 end )

            if (Traj_Nentries > 0)
            {
                for (int i = 0; i < Traj_Nentries; i++)
                {
                    if (TRAJ.getInt("pindex", i) >= 0 && TRAJ.getInt("pindex", i) < size && TRAJ.getInt("detector", i) == 6 && TRAJ.getInt("layer", i) == 6)
                    {
                        part_DC_c1x[TRAJ.getInt("pindex", i)] = TRAJ.getFloat("x", i);
                        part_DC_c1y[TRAJ.getInt("pindex", i)] = TRAJ.getFloat("y", i);
                        part_DC_c1z[TRAJ.getInt("pindex", i)] = TRAJ.getFloat("z", i);
                        part_DC_region[TRAJ.getInt("pindex", i)] = 1;
                    }
                    if (TRAJ.getInt("pindex", i) >= 0 && TRAJ.getInt("pindex", i) < size && TRAJ.getInt("detector", i) == 6 && TRAJ.getInt("layer", i) == 18)
                    {
                        part_DC_c2x[TRAJ.getInt("pindex", i)] = TRAJ.getFloat("x", i);
                        part_DC_c2y[TRAJ.getInt("pindex", i)] = TRAJ.getFloat("y", i);
                        part_DC_c2z[TRAJ.getInt("pindex", i)] = TRAJ.getFloat("z", i);
                        part_DC_region[TRAJ.getInt("pindex", i)] = 2;
                    }
                    if (TRAJ.getInt("pindex", i) >= 0 && TRAJ.getInt("pindex", i) < size && TRAJ.getInt("detector", i) == 6 && TRAJ.getInt("layer", i) == 36)
                    {
                        part_DC_c3x[TRAJ.getInt("pindex", i)] = TRAJ.getFloat("x", i);
                        part_DC_c3y[TRAJ.getInt("pindex", i)] = TRAJ.getFloat("y", i);
                        part_DC_c3z[TRAJ.getInt("pindex", i)] = TRAJ.getFloat("z", i);
                        part_DC_region[TRAJ.getInt("pindex", i)] = 3;
                    }
                }
            }

            /// The sector of charged particles is extracted from the hit position in the local phi coordinate:

            for (int i = 0; i < size; ++i)
            {
                part_DC_sector[i] = determineSector(i);
            }

            iCandidate = -1;
            iMinCandidate = 0;
            MinCandidate = 999999;

            for (int j = 0; j < size; j++)
            {
                ClearVectors(4);
                RECPpx = PART.getFloat("px", j);
                RECPpy = PART.getFloat("py", j);
                RECPpz = PART.getFloat("pz", j);
                RECPcharge = PART.getInt("charge", j);
                RECPpid = PART.getInt("pid", j);
                RECPvx = PART.getFloat("vx", j);
                RECPvy = PART.getFloat("vy", j);
                RECPvz = PART.getFloat("vz", j);
                RECPstatus = PART.getInt("status", j);
                RECpidChi2 = PART.getFloat("chi2pid", j);
                RECPbeta = PART.getFloat("beta", j);

                // electron ID
                if (RECPcharge < 0 && RECPpid == 11) // && pmom.at(j)>1.)
                {
                    if (EC_hit_position_fiducial_cut_homogeneous(j) == true && DC_fiducial_cut_chi2(j, part_DC_region[j]) == true)
                    {

                        pmom = TMath::Sqrt(RECPpx * RECPpx + RECPpy * RECPpy + RECPpz * RECPpz);

                        El_Vec_temp.SetPxPyPzE(RECPpx, RECPpy, RECPpz, TMath::Sqrt(pmom * pmom + Elmass * Elmass));

                        El_Vec_.push_back(El_Vec_temp);

                        El_info_temp.push_back(RECPvx);
                        El_info_temp.push_back(RECPvy);
                        El_info_temp.push_back(RECPvz);
                        El_info_temp.push_back(RECPstatus);
                        El_info_temp.push_back(RECpidChi2);
                        El_info_temp.push_back(RECPbeta);

                        El_info_.push_back(El_info_temp);
                    }
                }

                // photon ID
                if (RECPcharge == 0 && RECPpid == 22)
                {
                    if (EC_hit_position_fiducial_cut(j) == true)
                    {
                        pmom = TMath::Sqrt(RECPpx * RECPpx + RECPpy * RECPpy + RECPpz * RECPpz);

                        Ph_Vec_temp.SetPxPyPzE(RECPpx, RECPpy, RECPpz, pmom);

                        Ph_Vec_.push_back(Ph_Vec_temp);

                        Ph_info_temp.push_back(RECPvx);
                        Ph_info_temp.push_back(RECPvy);
                        Ph_info_temp.push_back(RECPvz);
                        Ph_info_temp.push_back(RECPstatus);
                        Ph_info_temp.push_back(RECpidChi2);
                        Ph_info_temp.push_back(RECPbeta);

                        Ph_info_.push_back(Ph_info_temp);
                    }
                }

                // proton ID
                if (RECPpid == 2212 && RECPcharge > 0)
                {
                    if (DC_fiducial_cut_chi2(j, part_DC_region[j]) == true)
                    {
                        pmom = TMath::Sqrt(RECPpx * RECPpx + RECPpy * RECPpy + RECPpz * RECPpz);

                        Pr_Vec_temp.SetPxPyPzE(RECPpx, RECPpy, RECPpz, TMath::Sqrt(pmom * pmom + Pmass * Pmass));

                        Pr_Vec_.push_back(Pr_Vec_temp);

                        Pr_info_temp.push_back(RECPvx);
                        Pr_info_temp.push_back(RECPvy);
                        Pr_info_temp.push_back(RECPvz);
                        Pr_info_temp.push_back(RECPstatus);
                        Pr_info_temp.push_back(RECpidChi2);
                        Pr_info_temp.push_back(RECPbeta);
                        Pr_info_temp.push_back(part_Scint_CND_energy[j]);
                        Pr_info_temp.push_back(part_ScintX_CND_size[j]);
                        Pr_info_temp.push_back(part_ScintX_CND_layermult[j]);
                        Pr_info_temp.push_back(part_Scint_CTOF_energy[j]);
                        Pr_info_temp.push_back(part_ScintX_CTOF_size[j]);
                        Pr_info_temp.push_back(part_ScintX_CTOF_layermult[j]);
                        Pr_info_temp.push_back(part_ScintX_CND_dedx[j]);
                        Pr_info_temp.push_back(part_ScintX_CTOF_dedx[j]);

                        Pr_info_.push_back(Pr_info_temp);
                    }
                }

                // neutron ID
                if (RECPcharge == 0 && (RECPpid == 2112 || RECPpid == 0)) // && pmom.at(j)>0)
                {
                    if (RECPpid == 2112)
                    {
                        pmom = TMath::Sqrt(RECPpx * RECPpx + RECPpy * RECPpy + RECPpz * RECPpz);
                        N_Vec_temp.SetPxPyPzE(RECPpx, RECPpy, RECPpz, TMath::Sqrt(pmom * pmom + Nmass * Nmass));
                    }
                    else
                    {
                        px = 0;
                        py = 0;
                        pz = 0;
                        th = 0;
                        ph = 0;

                        pmom = TMath::Sqrt(TMath::Power(Nmass * RECPbeta, 2) / (1 - RECPbeta * RECPbeta));
                        if (part_Scint_CND_x[j] == 0 && part_Scint_CND_y[j] == 0 && part_Scint_CND_z[j] == 0)
                        {
                            direction.SetXYZ(part_Scint_CTOF_x[j], part_Scint_CTOF_y[j], part_Scint_CTOF_z[j]);
                        }
                        else
                        {
                            direction.SetXYZ(part_Scint_CND_x[j], part_Scint_CND_y[j], part_Scint_CND_z[j]);
                        }
                        th = direction.Theta();
                        ph = direction.Phi();
                        px = pmom * TMath::Sin(th) * TMath::Cos(ph);
                        py = pmom * TMath::Sin(th) * TMath::Sin(ph);
                        pz = pmom * TMath::Cos(th);
                        N_Vec_temp.SetPxPyPzE(px, py, pz, TMath::Sqrt(pmom * pmom + Nmass * Nmass));
                    }

                    N_Vec_.push_back(N_Vec_temp);

                    N_info_temp.push_back(RECPvx);
                    N_info_temp.push_back(RECPvy);
                    N_info_temp.push_back(RECPvz);
                    N_info_temp.push_back(RECPstatus);
                    N_info_temp.push_back(RECpidChi2);
                    N_info_temp.push_back(RECPbeta);
                    N_info_temp.push_back(part_Scint_CND_energy[j]);
                    N_info_temp.push_back(part_ScintX_CND_size[j]);
                    N_info_temp.push_back(part_ScintX_CND_layermult[j]);
                    N_info_temp.push_back(part_Scint_CTOF_energy[j]);
                    N_info_temp.push_back(part_ScintX_CTOF_size[j]);
                    N_info_temp.push_back(part_ScintX_CTOF_layermult[j]);
                    N_info_temp.push_back(part_ScintX_CND_dedx[j]);
                    N_info_temp.push_back(part_ScintX_CTOF_dedx[j]);

                    N_info_temp.push_back(part_Scint_CND_x[j]);
                    N_info_temp.push_back(part_Scint_CND_y[j]);
                    N_info_temp.push_back(part_Scint_CND_z[j]);
                    N_info_temp.push_back(part_Scint_CND_t[j]);
                    N_info_temp.push_back(part_Scint_CTOF_x[j]);
                    N_info_temp.push_back(part_Scint_CTOF_y[j]);
                    N_info_temp.push_back(part_Scint_CTOF_z[j]);
                    N_info_temp.push_back(part_Scint_CTOF_t[j]);

                    N_info_.push_back(N_info_temp);
                }

                // pion ID
                if (RECPpid == 211)
                {

                    pmom = TMath::Sqrt(RECPpx * RECPpx + RECPpy * RECPpy + RECPpz * RECPpz);

                    Pip_Vec_temp.SetPxPyPzE(RECPpx, RECPpy, RECPpz, TMath::Sqrt(pmom * pmom + Pipmass * Pipmass));

                    Pip_Vec_.push_back(Pip_Vec_temp);

                    Pip_info_temp.push_back(RECPvx);
                    Pip_info_temp.push_back(RECPvy);
                    Pip_info_temp.push_back(RECPvz);
                    Pip_info_temp.push_back(RECPstatus);
                    Pip_info_temp.push_back(RECpidChi2);
                    Pip_info_temp.push_back(RECPbeta);

                    Pip_info_.push_back(Pip_info_temp);
                }

                if (RECPpid == -211)
                {

                    pmom = TMath::Sqrt(RECPpx * RECPpx + RECPpy * RECPpy + RECPpz * RECPpz);

                    Pim_Vec_temp.SetPxPyPzE(RECPpx, RECPpy, RECPpz, TMath::Sqrt(pmom * pmom + Pimmass * Pimmass));

                    Pim_Vec_.push_back(Pim_Vec_temp);

                    Pim_info_temp.push_back(RECPvx);
                    Pim_info_temp.push_back(RECPvy);
                    Pim_info_temp.push_back(RECPvz);
                    Pim_info_temp.push_back(RECPstatus);
                    Pim_info_temp.push_back(RECpidChi2);
                    Pim_info_temp.push_back(RECPbeta);

                    Pim_info_.push_back(Pim_info_temp);
                }

                // check for other tracks matching a fake neutron deposit in the cnd
                if (RECPcharge > 0 && part_Scint_CND_energy[j] != 0)
                {

                    fN_info_temp.push_back(part_Scint_CND_energy[j]);
                    fN_info_temp.push_back(part_ScintX_CND_size[j]);
                    fN_info_temp.push_back(part_ScintX_CND_layermult[j]);
                    fN_info_temp.push_back(part_Scint_CTOF_energy[j]);
                    fN_info_temp.push_back(part_ScintX_CTOF_size[j]);
                    fN_info_temp.push_back(part_ScintX_CTOF_layermult[j]);
                    fN_info_temp.push_back(part_ScintX_CND_dedx[j]);
                    fN_info_temp.push_back(part_ScintX_CTOF_dedx[j]);

                    fN_info_temp.push_back(part_Scint_CND_x[j]);
                    fN_info_temp.push_back(part_Scint_CND_y[j]);
                    fN_info_temp.push_back(part_Scint_CND_z[j]);
                    fN_info_temp.push_back(part_Scint_CND_t[j]);
                    fN_info_temp.push_back(part_Scint_CTOF_x[j]);
                    fN_info_temp.push_back(part_Scint_CTOF_y[j]);
                    fN_info_temp.push_back(part_Scint_CTOF_z[j]);
                    fN_info_temp.push_back(part_Scint_CTOF_t[j]);

                    fN_info_.push_back(fN_info_temp);
                }
            }
            // pDVCS tree
            if (channels == "pDVCS")
            {
                build_tree(PTarget_Vec, El_Vec_, El_info_, Pr_Vec_, Pr_info_, Ph_Vec_, Ph_info_);
                if (strip_Q2.size() != 0)
                    pDVCS_tree->Fill();
            }
            // nDVCS tree
            if (channels == "nDVCS") // && Pim_Vec_.size()==0 && Pip_Vec_.size()==0 && Pr_Vec_.size()==0)
            {
                build_tree(NTarget_Vec, El_Vec_, El_info_, N_Vec_, N_info_, Ph_Vec_, Ph_info_, fN_info_);
                if (strip_Q2.size() != 0)
                    nDVCS_tree->Fill();
            }
        }
    }
    cout << endl;

    if (channels == "pDVCS")
        pDVCS_tree->Write();
    if (channels == "nDVCS")
        nDVCS_tree->Write();

    outfile->Close();
}

void analysis::build_tree(TLorentzVector NucTarget_Vec, vector<TLorentzVector> electron, vector<vector<double>> electron_inf, vector<TLorentzVector> nucleon, vector<vector<double>> nucleon_inf, vector<TLorentzVector> photon, vector<vector<double>> photon_inf)
{
    vector<vector<double>> dummyVector;
    build_tree(NucTarget_Vec, electron, electron_inf, nucleon, nucleon_inf, photon, photon_inf, dummyVector);
}

void analysis::build_tree(TLorentzVector NucTarget_Vec, vector<TLorentzVector> electron, vector<vector<double>> electron_inf, vector<TLorentzVector> nucleon, vector<vector<double>> nucleon_inf, vector<TLorentzVector> photon, vector<vector<double>> photon_inf, vector<vector<double>> fnucleon_inf)
{
    iCandidate = -1;
    iMinCandidate = 0;
    MinCandidate = 999999;
    ClearVectors(2);

    if (electron.size() >= 1 && nucleon.size() >= 1 && photon.size() >= 1)
    {
        for (int j = 0; j < electron.size(); j++)
        {

            El_Vec = electron.at(j);

            for (int k = 0; k < nucleon.size(); k++)
            {

                Nuc_Vec = nucleon.at(k);

                for (int l = 0; l < photon.size(); l++)
                {

                    Ph_Vec = photon.at(l);

                    mass_temp = NucTarget_Vec.T();
                    Q2 = 4 * Ebeam * El_Vec.P() * TMath::Power(TMath::Sin(El_Vec.Theta() / 2), 2.);
                    W = TMath::Sqrt(mass_temp * mass_temp + 2 * mass_temp * (Ebeam - El_Vec.P()) - Q2);
                    Xbj = Q2 / (2 * mass_temp * (Ebeam - El_Vec.P()));

                    if (Q2 < 1)
                        continue;

                    if (W < 2)
                        continue;

                    El_px = El_Vec.Px();
                    El_py = El_Vec.Py();
                    El_pz = El_Vec.Pz();
                    El_E = El_Vec.E();
                    El_P = El_Vec.P();
                    El_Theta = El_Vec.Theta() * 180. / TMath::Pi();
                    El_Phi = El_Vec.Phi() * 180. / TMath::Pi();
                    El_vx = electron_inf.at(j).at(0);
                    El_vy = electron_inf.at(j).at(1);
                    El_vz = electron_inf.at(j).at(2);
                    El_status = electron_inf.at(j).at(3);
                    El_chi2pid = electron_inf.at(j).at(4);
                    El_beta = electron_inf.at(j).at(5);

                    Ph_px = Ph_Vec.Px();
                    Ph_py = Ph_Vec.Py();
                    Ph_pz = Ph_Vec.Pz();
                    Ph_E = Ph_Vec.E();
                    Ph_P = Ph_Vec.P();
                    Ph_Theta = Ph_Vec.Theta() * 180. / TMath::Pi();
                    Ph_Phi = Ph_Vec.Phi() * 180. / TMath::Pi();
                    Ph_vx = photon_inf.at(l).at(0);
                    Ph_vy = photon_inf.at(l).at(1);
                    Ph_vz = photon_inf.at(l).at(2);
                    Ph_status = photon_inf.at(l).at(3);
                    Ph_chi2pid = photon_inf.at(l).at(4);
                    Ph_beta = photon_inf.at(l).at(5);

                    Nuc_px = Nuc_Vec.Px();
                    Nuc_py = Nuc_Vec.Py();
                    Nuc_pz = Nuc_Vec.Pz();
                    Nuc_E = Nuc_Vec.E();
                    Nuc_P = Nuc_Vec.P();
                    Nuc_Theta = Nuc_Vec.Theta() * 180. / TMath::Pi();
                    Nuc_Phi = Nuc_Vec.Phi() * 180. / TMath::Pi();
                    Nuc_vx = nucleon_inf.at(k).at(0);
                    Nuc_vy = nucleon_inf.at(k).at(1);
                    Nuc_vz = nucleon_inf.at(k).at(2);
                    Nuc_status = nucleon_inf.at(k).at(3);
                    Nuc_chi2pid = nucleon_inf.at(k).at(4);
                    Nuc_beta = nucleon_inf.at(k).at(5);
                    Nuc_CND_energy = nucleon_inf.at(k).at(6);
                    Nuc_CND_size = nucleon_inf.at(k).at(7);
                    Nuc_CND_layermult = nucleon_inf.at(k).at(8);
                    Nuc_CTOF_energy = nucleon_inf.at(k).at(9);
                    Nuc_CTOF_size = nucleon_inf.at(k).at(10);
                    Nuc_CTOF_layermult = nucleon_inf.at(k).at(11);
                    Nuc_CND_dedx = nucleon_inf.at(k).at(12);
                    Nuc_CTOF_dedx = nucleon_inf.at(k).at(13);

                    if (Ph_P < 2)
                        continue;

                    if (El_P < 1)
                        continue;

                    if (Nuc_P < 0.3)
                        continue;

                    iCandidate++;

                    // checking if neutron had a matched track
                    matched_track_to_fN = false;
                    for (int m = 0; m < fnucleon_inf.size(); m++)
                    {
                        Dx = abs(nucleon_inf.at(k).at(14) - fnucleon_inf.at(m).at(8));
                        Dy = abs(nucleon_inf.at(k).at(15) - fnucleon_inf.at(m).at(9));
                        Dz = abs(nucleon_inf.at(k).at(16) - fnucleon_inf.at(m).at(10));
                        Dt = abs(nucleon_inf.at(k).at(17) - fnucleon_inf.at(m).at(11));
                        if (Dx < 20 && Dy < 20 && Dz < 10 && Dt < 1.2)
                        {
                            matched_track_to_fN = true;

                            Nuc_CND_energy += fnucleon_inf.at(m).at(0);
                            Nuc_CND_size += fnucleon_inf.at(m).at(1);
                            Nuc_CTOF_energy += fnucleon_inf.at(m).at(3);
                            Nuc_CTOF_size += fnucleon_inf.at(m).at(4);
                            if (fnucleon_inf.at(m).at(2) > Nuc_CND_layermult)
                            {
                                Nuc_CND_layermult = fnucleon_inf.at(m).at(2);
                            }
                        }
                    }

                    if (matched_track_to_fN)
                    {
                        matched_fN.push_back(1);
                    }
                    else
                    {
                        matched_fN.push_back(0);
                    }

                    strip_Q2.push_back(Q2);
                    strip_W.push_back(W);
                    strip_Xbj.push_back(Xbj);

                    strip_El_px.push_back(El_px);
                    strip_El_py.push_back(El_py);
                    strip_El_pz.push_back(El_pz);
                    strip_El_E.push_back(El_E);
                    strip_El_P.push_back(El_P);
                    strip_El_Theta.push_back(El_Theta);
                    strip_El_Phi.push_back(El_Phi);
                    strip_El_vx.push_back(El_vx);
                    strip_El_vy.push_back(El_vy);
                    strip_El_vz.push_back(El_vz);
                    strip_El_status.push_back(El_status);
                    strip_El_chi2pid.push_back(El_chi2pid);
                    strip_El_beta.push_back(El_beta);

                    strip_Ph_px.push_back(Ph_px);
                    strip_Ph_py.push_back(Ph_py);
                    strip_Ph_pz.push_back(Ph_pz);
                    strip_Ph_E.push_back(Ph_E);
                    strip_Ph_P.push_back(Ph_P);
                    strip_Ph_Theta.push_back(Ph_Theta);
                    strip_Ph_Phi.push_back(Ph_Phi);
                    strip_Ph_vx.push_back(Ph_vx);
                    strip_Ph_vy.push_back(Ph_vy);
                    strip_Ph_vz.push_back(Ph_vz);
                    strip_Ph_status.push_back(Ph_status);
                    strip_Ph_chi2pid.push_back(Ph_chi2pid);
                    strip_Ph_beta.push_back(Ph_beta);

                    strip_Nuc_px.push_back(Nuc_px);
                    strip_Nuc_py.push_back(Nuc_py);
                    strip_Nuc_pz.push_back(Nuc_pz);
                    strip_Nuc_E.push_back(Nuc_E);
                    strip_Nuc_P.push_back(Nuc_P);
                    strip_Nuc_Theta.push_back(Nuc_Theta);
                    strip_Nuc_Phi.push_back(Nuc_Phi);
                    strip_Nuc_vx.push_back(Nuc_vx);
                    strip_Nuc_vy.push_back(Nuc_vy);
                    strip_Nuc_vz.push_back(Nuc_vz);
                    strip_Nuc_status.push_back(Nuc_status);
                    strip_Nuc_chi2pid.push_back(Nuc_chi2pid);
                    strip_Nuc_beta.push_back(Nuc_beta);
                    strip_Nuc_CND_energy.push_back(Nuc_CND_energy);
                    strip_Nuc_CND_size.push_back(Nuc_CND_size);
                    strip_Nuc_CND_layermult.push_back(Nuc_CND_layermult);
                    strip_Nuc_CND_dedx.push_back(Nuc_CND_dedx);
                    strip_Nuc_CTOF_energy.push_back(Nuc_CTOF_energy);
                    strip_Nuc_CTOF_size.push_back(Nuc_CTOF_size);
                    strip_Nuc_CTOF_layermult.push_back(Nuc_CTOF_layermult);
                    strip_Nuc_CTOF_dedx.push_back(Nuc_CTOF_dedx);

                    VelectronIn = ElectronBeam.Vect();
                    VelectronOut = El_Vec.Vect();
                    VnucleonOut = Nuc_Vec.Vect();
                    VphotonOut = Ph_Vec.Vect();
                    Vvirtualphoton = (ElectronBeam - El_Vec).Vect();

                    Vlepto = VelectronIn.Cross(VelectronOut);
                    Vhadro = VnucleonOut.Cross(Vvirtualphoton);
                    VhadroPP = VnucleonOut.Cross(VphotonOut);

                    Phi_Nuc_temp = 180. / TMath::Pi() * Vlepto.Angle(Vhadro);
                    Phi_Ph_temp = 180. / TMath::Pi() * Vlepto.Angle(VhadroPP);

                    if (Vlepto.Dot(VnucleonOut) > 0.)
                        Phi_Nuc_temp = 360. - Phi_Nuc_temp;
                    if (Vlepto.Dot(VphotonOut) < 0.)
                        Phi_Ph_temp = 360. - Phi_Ph_temp;

                    Phi_Nuc.push_back(Phi_Nuc_temp);
                    Phi_Ph.push_back(Phi_Ph_temp);

                    t_Nuc.push_back((Nuc_Vec - NucTarget_Vec).M2());
                    t_Ph.push_back((ElectronBeam - El_Vec - Ph_Vec).M2());

                    BalV = ElectronBeam + NucTarget_Vec - Ph_Vec - El_Vec - Nuc_Vec;

                    mm2_eNg.push_back((ElectronBeam + Target_Vec - Nuc_Vec - El_Vec - Ph_Vec).M2());
                    mm2_eNg_N.push_back((ElectronBeam + NucTarget_Vec - Nuc_Vec - El_Vec - Ph_Vec).M2());
                    mm2_eg.push_back((ElectronBeam + NucTarget_Vec - El_Vec - Ph_Vec).M2());
                    mm2_eNX_N.push_back((ElectronBeam + NucTarget_Vec - Nuc_Vec - El_Vec).M2());

                    bestCandidateFlag.push_back(0);
                    // add criteria here: thetagammax delta phi delta t
                    if (abs((ElectronBeam + NucTarget_Vec - Nuc_Vec - El_Vec - Ph_Vec).M2()) < MinCandidate)
                    {
                        MinCandidate = abs((ElectronBeam + NucTarget_Vec - Nuc_Vec - El_Vec - Ph_Vec).M2());

                        bestCandidateFlag.at(iMinCandidate) = 0;
                        bestCandidateFlag.at(iCandidate) = 1;

                        iMinCandidate = iCandidate;
                    }

                    Xbal.push_back(BalV.X());
                    Ybal.push_back(BalV.Y());
                    Zbal.push_back(BalV.Z());
                    Ebal.push_back(BalV.E());

                    delta_Phi.push_back(Phi_Nuc_temp - Phi_Ph_temp);
                    delta_t.push_back((Nuc_Vec - NucTarget_Vec).M2() - (ElectronBeam - El_Vec - Ph_Vec).M2());
                    miss_mom_eNg.push_back(TMath::Sqrt(pow(BalV.X(), 2) + pow(BalV.Y(), 2) + pow(BalV.Z(), 2)));
                    p_perp.push_back(TMath::Sqrt(pow(BalV.X(), 2) + pow(BalV.Y(), 2)));

                    missing_gamma_N = (ElectronBeam + NucTarget_Vec - El_Vec - Nuc_Vec).Vect();

                    theta_gamma_e.push_back(180. / TMath::Pi() * TMath::ACos((VphotonOut.Dot(VelectronOut)) / (VphotonOut.Mag() * VelectronOut.Mag())));
                    theta_gamma_X.push_back(180. / TMath::Pi() * TMath::ACos((missing_gamma_N.Dot(VphotonOut)) / (missing_gamma_N.Mag() * VphotonOut.Mag())));
                    theta_N_e.push_back(180. / TMath::Pi() * TMath::ACos((VnucleonOut.Dot(VelectronOut)) / (VnucleonOut.Mag() * VelectronOut.Mag())));

                    missing_N_beta.push_back((ElectronBeam + NucTarget_Vec - El_Vec - Ph_Vec).Beta());
                    recon_N_beta.push_back(Nuc_Vec.Beta());

                    exclusivity_chi2 = 0;
                    exclusivity_chi2 += pow((ElectronBeam + NucTarget_Vec - Nuc_Vec - El_Vec).M2(), 2);
                    exclusivity_chi2 += pow(180. / TMath::Pi() * TMath::ACos((missing_gamma_N.Dot(VphotonOut)) / (missing_gamma_N.Mag() * VphotonOut.Mag())), 2);
                    exclusivity_chi2 += pow(Phi_Nuc_temp - Phi_Ph_temp, 2);
                    exclusivity_chi2 += pow((Nuc_Vec - NucTarget_Vec).M2() - (ElectronBeam - El_Vec - Ph_Vec).M2(), 2);
                    Exclusive_4DChi2.push_back(exclusivity_chi2);
                }
            }
        }
    }
    Minchi2 = 999999;

    for (int m = 0; m < Exclusive_4DChi2.size(); m++)
    {
        if (Exclusive_4DChi2.at(m) < Minchi2)
        {
            Minchi2 = Exclusive_4DChi2.at(m);
        }
    }

    for (int m = 0; m < Exclusive_4DChi2.size(); m++)
    {
        if (Exclusive_4DChi2.at(m) == Minchi2)
        {
            best4DChi2Flag.push_back(1);
        }
        else
        {
            best4DChi2Flag.push_back(0);
        }
    }
}