/**
 * DumpWaveformInfo_PDVD.C
 *
 * Author: Angelo Ralaikoto
 *
 * Functionality:
 *  - Loads a channel map JSON file (e.g., PDVD_PDS_Mapping_v07082025.json)
 *  - Opens a ROOT file containing waveform histograms from the analyzer: "OpDetDigiAnaDUNE" (e.g., detsimHist.root)
 *  - Looks for histograms with names like: 
 *        event_<N>_opchannel_<OfflineChannel>_waveform_<M>
 *  - Extracts OfflineChannel from the name and retrieves full metadata from the JSON:
 *        channel, pd_type, wls, eff_Ar, eff_Xe, name, Slot, Link, DaphneChannel
 *  - Plots each waveform and saves it as a PNG in waveform_plots/
 *  - Creates waveform_channels_info_PDVD.txt summarizing waveform metadata
 *
 * Usage inside ROOT:
 *   root -l
 *   .x DumpWaveformInfo_PDVD.C
 */

#include <TFile.h>
#include <TKey.h>
#include <TClass.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TDirectoryFile.h>
#include <TROOT.h>
#include <TStyle.h>

#include <nlohmann/json.hpp>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <regex>
#include <map>
#include <sstream>

using json = nlohmann::json;
using namespace std;

struct ChannelInfo {
    int channel;
    std::string pd_type;
    std::string wls;
    double eff_Ar;
    double eff_Xe;
    std::string name;
    int Slot;
    int Link;
    int DaphneChannel;
    int OfflineChannel;
};

std::map<int, ChannelInfo> LoadChannelMap(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: cannot open JSON file " << filename << "\n";
        return {};
    }

    json j;
    file >> j;

    std::map<int, ChannelInfo> offlineMap;

    for (const auto& entry : j) {
        for (const auto& hw : entry["HardwareChannel"]) {
            ChannelInfo info;
            info.channel = entry["channel"];
            info.pd_type = entry["pd_type"];
            info.wls = entry["wls"];
            info.eff_Ar = entry["eff_Ar"];
            info.eff_Xe = entry["eff_Xe"];
            info.name = entry["name"];
            info.Slot = hw["Slot"];
            info.Link = hw["Link"];
            info.DaphneChannel = hw["DaphneChannel"];
            info.OfflineChannel = hw["OfflineChannel"];
            offlineMap[info.OfflineChannel] = info;
        }
    }

    return offlineMap;
}

void DumpWaveformInfo(const char* jsonFile = "../../../../dunecore/dunecore/ChannelMap/PDVD_PDS_Mapping_v07082025.json",
                      const char* rootFile = "../../../../detsim_single_protoDUNE_hist.root") {
    gROOT->SetBatch(kTRUE); // no GUI popups
    gStyle->SetOptStat(0);  // no stats box
    gSystem->mkdir("waveform_plots", kTRUE);

    std::map<int, ChannelInfo> channelMap = LoadChannelMap(jsonFile);

    std::ofstream outFile("waveform_channels_info_PDVD.txt");
    if (!outFile.is_open()) {
        std::cerr << "Error: cannot create waveform_channels_info_PDVD.txt\n";
        return;
    }

    TFile* file = TFile::Open(rootFile);
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening ROOT file\n";
        return;
    }

    TDirectoryFile* ana = (TDirectoryFile*)file->Get("ana");
    if (!ana) {
        std::cerr << "Directory 'ana' not found in ROOT file\n";
        return;
    }

    TIter next(ana->GetListOfKeys());
    TKey* key;
    std::regex pattern("event_\\d+_opchannel_(\\d+)_waveform_\\d+");
    int count = 0;

    while ((key = (TKey*)next())) {
        std::string name = key->GetName();
        std::smatch match;

        if (std::regex_search(name, match, pattern)) {
            int offlineChannel = std::stoi(match[1]);

            if (channelMap.find(offlineChannel) != channelMap.end()) {
                ChannelInfo info = channelMap[offlineChannel];

                TH1* hist = dynamic_cast<TH1*>(key->ReadObj());
                if (!hist) continue;

                TCanvas* c = new TCanvas();
                hist->Draw();

                std::stringstream pngNameStream;
                pngNameStream << "waveform_plots/waveform_"
                              << "CH" << info.channel << "_"
                              << "Offline" << info.OfflineChannel << "_"
                              << info.pd_type << "_"
                              << info.name << "_"
                              << "Slot" << info.Slot << "_"
                              << "Link" << info.Link << "_"
                              << "Daphne" << info.DaphneChannel << ".png";

                std::string pngName = pngNameStream.str();
                c->SaveAs(pngName.c_str());
                delete c;

                outFile << "Waveform: " << name << "\n"
                        << " -> Channel: " << info.channel << "\n"
                        << " -> pd_type: " << info.pd_type << "\n"
                        << " -> wls: " << info.wls << "\n"
                        << " -> eff_Ar: " << info.eff_Ar << "\n"
                        << " -> eff_Xe: " << info.eff_Xe << "\n"
                        << " -> name: " << info.name << "\n"
                        << " -> Slot: " << info.Slot << "\n"
                        << " -> Link: " << info.Link << "\n"
                        << " -> DaphneChannel: " << info.DaphneChannel << "\n"
                        << " -> OfflineChannel: " << info.OfflineChannel << "\n\n";

                count++;
            }
        }
    }

    outFile.close();
    file->Close();

    std::cout << "Done. Saved " << count << " plots and created waveform_channels_info_PDVD.txt\n";
}

