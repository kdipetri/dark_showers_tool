/*
    This File is a simplified version fo the old higgsPortal.cc from Anthony.
    The point of this file is to run a pythia card and output to a hepmc file

    Matias Mantiñan
*/


#include <iostream>
#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TTree.h"
#include "Pythia8Plugins/HepMC3.h"

// declare pythia8 namespace
using namespace Pythia8;
using namespace std;



// Progress bar
void printProgressBar(int progress, int total, int barWidth) {
    float ratio = progress / static_cast<float>(total);
    int pos = static_cast<int>(barWidth * ratio);
    
    cout << "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) cout << "=";
        else if (i == pos) cout << ">";
        else cout << " ";
    }
    cout << "] " << int(ratio * 100.0) << " %\r";
    cout.flush();
}


int main(int argc, char* argv[]){
    // check user is inputting correct parameters
    if(argc != 4){
        std::cout << "Usage: ./card_runner.exe <pythiaCard> <outFileName> <maxEvents>" << std::endl;
        return 1;
    }

    // read in user parameters
    string pythiaCard = argv[1];
    string outFileName = argv[2];
    int maxEvents = stoi(argv[3]);

    // Generator.
    Pythia pythia;

    // Shorthand for the event record in pythia.
    Event& event = pythia.event;
    const Info& info = pythia.info;

    // Read in commands from external file.
    pythia.readFile(pythiaCard.c_str());

    // Initialize.
    pythia.init();

    // Interface for conversion from Pythia8::Event to HepMC event.
    HepMC3::Pythia8ToHepMC3 toHepMC;
    // Specify file where HepMC events will be stored.
    HepMC3::WriterAscii ascii_io((outFileName + ".hepmc").c_str());


    int barWidth = 50;

    for (int iE = 0; iE < maxEvents; ++iE) {
        printProgressBar(iE,maxEvents,barWidth); // print progress bar

        if(!pythia.next()) continue;
        
        // Construct new empty HepMC event and fill it.
        // Default units are ( HepMC3::Units::GEV, HepMC3::Units::MM)
        // but can be changed in the GenEvent constructor.
        HepMC3::GenEvent hepmcevt;
        toHepMC.fill_next_event( pythia, &hepmcevt );

        // Write the HepMC event to file.
        ascii_io.write_event(hepmcevt);

    }

    cout << endl; // last step for the progressbar

    return 0;
}