/*******************************************************************************
 * Dual-head FlexPET Full Energy Calibration Program
 * 
 * Author: Tran Cong Thien / PraanaTech
 * 
 * 
 * =============================================================================
 * DESCRIPTION
 * =============================================================================
 * This program performs energy calibration for the dual-head FlexPET
 * detector system. It combines two previously separate workflows:
 *   1. Peak analysis from binary measurement files (BinaryRead9)
 *   2. Energy calibration curve fitting (Read_BinOut_Calibrate)
 * 
 * The FlexPET detector consists of:
 *   - 4 GMSL (Gigabit Multimedia Serial Link) boards (gmslID: 0-3)
 *   - Each GMSL has 4 STiC (Silicon Timing readout Chip) ASICs (sticID: 0-3)
 *   - Each STiC reads 64 channels (channel: 0-63)
 *   - Total: 4 × 4 × 64 = 1024 channels
 *   
 * The detector measures gamma photon energy using Time-over-Threshold (ToT),
 * which has a non-linear relationship with deposited energy:
 *   ToT = P0 * exp(E/P1) + P2
 * 
 * =============================================================================
 * WORKFLOW
 * =============================================================================
 * 1. Scan input folder for .dat/.bin files named by gamma energy (e.g., 511.0.dat)
 * 2. For each file:
 *    a. Auto-select optimal peak finding parameters based on gamma energy
 *    b. Read binary data and decode using LFSR (Linear Feedback Shift Register)
 *    c. Fill energy histograms for each channel
 *    d. Find peaks using ROOT's TSpectrum
 *    e. Fit Gaussian to dominant peak after background subtraction
 *    f. Store peak parameters (ToT position, resolution, area)
 * 3. Combine peak data from all sources
 * 4. Fit exponential calibration curve for each channel
 * 5. Generate output files and quality plots
 * 
 * =============================================================================
 * INPUT
 * =============================================================================
 *   argv[1]: Folder containing .dat/.bin files named by gamma energy
 *   argv[2]: Minimum ToT ADC value (default: 100)
 *   argv[3]: Maximum ToT ADC value (default: 800)
 *   argv[4]: Minimum peak area threshold (default: 100)
 *   argv[5]: Maximum peak area threshold (default: 1000000)
 * 
 * Note: Peak finding parameters (sigP, sigB) are automatically determined
 *       based on the gamma energy extracted from each filename:
 *         > 1000 keV : sigP=0.1, sigB=3 (e.g., Co-60)
 *         600-1000   : sigP=0.6, sigB=7 (e.g., Cs-137)
 *         400-600    : sigP=0.5, sigB=7 (e.g., Na-22 511 keV)
 *         < 400 keV  : sigP=0.6, sigB=7 (e.g., Ba-133)
 * 
 * =============================================================================
 * OUTPUT
 * =============================================================================
 *   - [InputFolder]/<energy>_Analysis/  : Per-source spectrum analysis
 *   - [InputFolder]_ECal/               : Combined calibration results
 * 
 * =============================================================================
 * COMPILE
 * =============================================================================
 *   g++ FlexPET_EnergyCal_v2.cpp -o FlexPET_EnergyCal_v2 `root-config --cflags --glibs` -lSpectrum
 *       
 * =============================================================================
 * USAGE
 * =============================================================================
 *   ./FlexPET_EnergyCal_v2 <InputFolder> [xmin] [xmax]
 ******************************************************************************/

//=============================================================================
// ROOT FRAMEWORK HEADERS
//=============================================================================
#include "TCanvas.h"      // For creating drawing canvases
#include "TStyle.h"       // For global style settings
#include "TH1F.h"         // 1D histogram with float bins
#include "TH1D.h"         // 1D histogram with double bins
#include "TH2F.h"         // 2D histogram with float bins
#include "TH2D.h"         // 2D histogram with double bins
#include "TF1.h"          // 1D function for fitting
#include "TGraphErrors.h" // Graph with error bars for calibration curves
#include "TAxis.h"        // Axis manipulation
#include "TSpectrum.h"    // Peak finding and background estimation
#include "TFile.h"        // ROOT file I/O
#include "TROOT.h"        // ROOT system interface

//=============================================================================
// STANDARD C++ HEADERS
//=============================================================================
#include <iostream>       // Console I/O (cout, cerr)
#include <fstream>        // File I/O (ifstream, ofstream)
#include <cstdlib>        // Standard library (strtod, atoi)
#include <cstring>        // String manipulation (sprintf)
#include <cmath>          // Math functions (abs, exp, pow)
#include <vector>         // Dynamic arrays
#include <string>         // String class
#include <algorithm>      // Sorting algorithms

//=============================================================================
// POSIX HEADERS FOR FILE SYSTEM OPERATIONS
//=============================================================================
#include <sys/stat.h>     // Directory creation (mkdir)
#include <dirent.h>       // Directory scanning (opendir, readdir)

using namespace std;

//=============================================================================
// GLOBAL CONSTANTS
//=============================================================================
// Detector geometry constants
const int N_GMSL   = 4;   // Number of GMSL boards
const int N_STIC   = 4;   // Number of STiC ASICs per GMSL
const int N_CHANNEL = 64; // Number of channels per STiC

// Histogram binning constants
const int NBINS_ENERGY = 8192;    // Number of bins for energy histograms
const int MAX_ENERGY_ADC = 32768; // Maximum ADC value

// LFSR lookup table size (15-bit LFSR)
const int LFSR_SIZE = 1 << 15;    // 2^15 = 32768

//=============================================================================
// DATA STRUCTURE: Peak analysis results for one channel
//=============================================================================
/**
 * @struct PeakData
 * @brief Stores the peak fitting results for a single channel from one source
 * 
 * This structure holds the results of Gaussian peak fitting for a single
 * detector channel. Multiple PeakData entries (from different gamma sources)
 * are used to fit the energy calibration curve.
 */
struct PeakData {
    int gmslID;        // GMSL board ID (0-3)
    int sticID;        // STiC ASIC ID (0-3)
    int channel;       // Channel number (0-63)
    double peakArea;   // Integrated peak area (counts)
    double ToT_ADC;    // Peak centroid in ToT ADC units
    double resolution; // Peak width (sigma) in ADC units
    double realEnergy; // Known gamma energy in keV
};

//=============================================================================
// FUNCTION: Peak position comparator for sorting
//=============================================================================
/**
 * @brief Comparison function for qsort to sort peaks by X position (energy)
 * @param arg1 Pointer to first peak data array [X, Y]
 * @param arg2 Pointer to second peak data array [X, Y]
 * @return -1 if arg1 < arg2, +1 if arg1 > arg2, 0 if equal
 * 
 * Used to sort peaks found by TSpectrum in ascending order of position,
 * allowing identification of the highest-energy peak.
 */
int comparePeaksByPosition(const void *arg1, const void *arg2)
{
    const int *lhs = static_cast<const int*>(arg1);
    const int *rhs = static_cast<const int*>(arg2);
    
    if (lhs[0] < rhs[0]) return -1;
    if (lhs[0] > rhs[0]) return +1;
    return 0;
}

//=============================================================================
// FUNCTION: Channel to detector X-Y coordinate mapping
//=============================================================================
/**
 * @brief Maps (channel, sticID, gmslID) to detector X coordinate
 * @param channel Channel number within STiC (0-63)
 * @param sid STiC ASIC ID (0-3)
 * @param gid GMSL board ID (0-3)
 * @return X coordinate in detector coordinate system (0-31)
 * 
 * The AS-PET detector uses a complex channel mapping due to the physical
 * arrangement of SiPM arrays and readout ASICs. This function converts
 * the electronic channel address to physical detector position.
 * 
 * The mapping accounts for:
 *   - 8×8 channel layout within each STiC
 *   - STiC orientation (some are flipped)
 *   - GMSL board position (left/right detector head)
 */
int mapChannelToX(int channel, int sid, int gid)
{
    // First, map channel to local X within 8×8 STiC array
    int x = 0;
    
    // Column 0 channels
    if ((channel==58)||(channel==57)||(channel==54)||(channel==52)||
        (channel==46)||(channel==44)||(channel==38)||(channel==37)) { x = 0; }
    // Column 1 channels
    if ((channel==60)||(channel==59)||(channel==56)||(channel==50)||
        (channel==48)||(channel==42)||(channel==40)||(channel==36)) { x = 1; }
    // Column 2 channels
    if ((channel==61)||(channel==55)||(channel==53)||(channel==47)||
        (channel==45)||(channel==39)||(channel==32)||(channel==35)) { x = 2; }
    // Column 3 channels
    if ((channel==63)||(channel==62)||(channel==51)||(channel==49)||
        (channel==43)||(channel==41)||(channel==34)||(channel==33)) { x = 3; }
    // Column 4 channels
    if ((channel==0)||(channel==1)||(channel==12)||(channel==14)||
        (channel==20)||(channel==22)||(channel==29)||(channel==30)) { x = 4; }
    // Column 5 channels
    if ((channel==2)||(channel==8)||(channel==10)||(channel==16)||
        (channel==18)||(channel==24)||(channel==31)||(channel==28)) { x = 5; }
    // Column 6 channels
    if ((channel==3)||(channel==4)||(channel==7)||(channel==13)||
        (channel==15)||(channel==21)||(channel==23)||(channel==27)) { x = 6; }
    // Column 7 channels
    if ((channel==5)||(channel==6)||(channel==9)||(channel==11)||
        (channel==17)||(channel==19)||(channel==25)||(channel==26)) { x = 7; }
    
    // Apply STiC orientation correction (some STiCs are mounted flipped)
    if (sid / 2 == 1) { 
        x = 7 - x;  // Mirror X for STiC 2 and 3
    }
    
    // Calculate global X position combining GMSL and STiC offsets
    // Formula: GMSL contributes 16-column offset, STiC contributes 8-column offset
    x = (gid % 2) * 16 + (sid % 2) * 8 + x;
    
    return x;
}

/**
 * @brief Maps (channel, sticID, gmslID) to detector Y coordinate
 * @param channel Channel number within STiC (0-63)
 * @param sid STiC ASIC ID (0-3)
 * @param gid GMSL board ID (0-3)
 * @return Y coordinate in detector coordinate system (0-15)
 * 
 * Similar to mapChannelToX but for Y coordinate. The Y mapping also
 * depends on the XOR of STiC and GMSL IDs due to the detector geometry.
 */
int mapChannelToY(int channel, int sid, int gid)
{
    // First, map channel to local Y within 8×8 STiC array
    int y = 0;
    
    // Row 0 channels
    if ((channel==58)||(channel==60)||(channel==61)||(channel==63)||
        (channel==0)||(channel==2)||(channel==3)||(channel==5)) { y = 0; }
    // Row 1 channels
    if ((channel==57)||(channel==59)||(channel==55)||(channel==62)||
        (channel==1)||(channel==8)||(channel==4)||(channel==6)) { y = 1; }
    // Row 2 channels
    if ((channel==54)||(channel==56)||(channel==53)||(channel==51)||
        (channel==12)||(channel==10)||(channel==7)||(channel==9)) { y = 2; }
    // Row 3 channels
    if ((channel==52)||(channel==50)||(channel==47)||(channel==49)||
        (channel==14)||(channel==16)||(channel==13)||(channel==11)) { y = 3; }
    // Row 4 channels
    if ((channel==46)||(channel==48)||(channel==45)||(channel==43)||
        (channel==20)||(channel==18)||(channel==15)||(channel==17)) { y = 4; }
    // Row 5 channels
    if ((channel==44)||(channel==42)||(channel==39)||(channel==41)||
        (channel==22)||(channel==24)||(channel==21)||(channel==19)) { y = 5; }
    // Row 6 channels
    if ((channel==38)||(channel==40)||(channel==32)||(channel==34)||
        (channel==29)||(channel==31)||(channel==23)||(channel==25)) { y = 6; }
    // Row 7 channels
    if ((channel==26)||(channel==36)||(channel==35)||(channel==33)||
        (channel==30)||(channel==28)||(channel==27)||(channel==37)) { y = 7; }
    
    // Apply orientation correction based on XOR of STiC and GMSL quadrant
    // This accounts for the physical mounting orientation
    int quadrantXOR = (sid / 2) ^ (gid / 2);
    if (quadrantXOR == 1) { 
        y = 7 - y;  // Mirror Y for certain quadrant combinations
    }
    
    // Calculate global Y position
    y = quadrantXOR * 8 + y;
    
    return y;
}

//=============================================================================
// FUNCTION: Initialize LFSR lookup table
//=============================================================================
/**
 * @brief Initializes the Linear Feedback Shift Register (LFSR) lookup table
 * @param m_lut Pointer to lookup table array of size 2^15
 * 
 * The STiC ASIC encodes timing values using a 15-bit LFSR (pseudo-random
 * counter) instead of a linear counter. This provides better noise immunity.
 * This function creates a lookup table to decode LFSR values back to linear.
 * 
 * LFSR polynomial: new_bit = NOT(bit13 XOR bit14)
 * This creates a maximal-length sequence of 2^15 - 1 states.
 */
void initLFSR(int16_t* m_lut)
{
    // Mark invalid state (all 1s should never occur in valid data)
    m_lut[0x7FFF] = -1;
    
    uint16_t lfsr = 0x0000;  // Starting state
    
    // Generate all 2^15 - 1 states and store their linear equivalents
    for (int16_t n = 0; n < (1 << 15) - 1; ++n)
    {
        m_lut[lfsr] = n;  // Map LFSR value to linear counter value
        
        // Calculate next LFSR state using feedback polynomial
        const uint8_t bits13_14 = lfsr >> 13;  // Extract bits 13 and 14
        uint8_t new_bit;
        
        // new_bit = NOT(bit13 XOR bit14)
        switch (bits13_14)
        {
            case 0x00:  // 00 -> NOT(0 XOR 0) = NOT(0) = 1
            case 0x03:  // 11 -> NOT(1 XOR 1) = NOT(0) = 1
                new_bit = 0x01;
                break;
            case 0x01:  // 01 -> NOT(0 XOR 1) = NOT(1) = 0
            case 0x02:  // 10 -> NOT(1 XOR 0) = NOT(1) = 0
                new_bit = 0x00;
                break;
        }
        
        // Shift register left and insert new bit
        lfsr = (lfsr << 1) | new_bit;
        lfsr &= 0x7FFF;  // Keep only 15 bits
    }
}

//=============================================================================
// FUNCTION: Determine optimal peak finding parameters based on gamma energy
//=============================================================================
/**
 * @brief Determines optimal sigP and sigB parameters based on gamma energy
 * @param energy Gamma energy in keV
 * @param sigP Output: Peak finding sigma parameter for TSpectrum
 * @param sigB Output: Background estimation sigma parameter
 * 
 * Different gamma energies require different peak finding sensitivities:
 * 
 * | Energy Range      | sigP | sigB | Typical Sources                    |
 * |-------------------|------|------|------------------------------------|
 * | > 1000 keV        | 0.1  | 3    | Co-60 (1173, 1332 keV)            |
 * | 600 - 1000 keV    | 0.6  | 7    | Cs-137 (662 keV)                  |
 * | 400 - 600 keV     | 0.5  | 7    | Na-22 (511 keV)                   |
 * | < 400 keV         | 0.6  | 7    | Ba-133 (356, 303, 276 keV)        |
 * 
 * Lower sigP = more sensitive peak finding (finds more peaks)
 * Higher sigB = smoother background estimation
 * 
 * For high-energy peaks (>1 MeV), the spectrum is cleaner with fewer
 * background features, so we use lower sigP to find the peak precisely
 * and lower sigB since the background is simpler.
 */
void getOptimalPeakParameters(double energy, double& sigP, double& sigB)
{
    if (energy > 1000.0) {
        // High energy gamma (e.g., Co-60: 1173 keV, 1332 keV)
        // Cleaner spectrum, need precise peak finding
        sigP = 0.1;
        sigB = 3.0;
    } 
    else if (energy > 600.0) {
        // Medium-high energy (e.g., Cs-137: 662 keV)
        sigP = 0.6;
        sigB = 7.0;
    }
    else if (energy > 400.0) {
        // Medium energy (e.g., Na-22: 511 keV)
        sigP = 0.5;
        sigB = 7.0;
    }
    else {
        // Low energy (e.g., Ba-133: 356, 303, 276 keV)
        // More Compton background, need careful peak finding
        sigP = 0.6;
        sigB = 7.0;
    }
}

//=============================================================================
// FUNCTION: Extract gamma energy from filename
//=============================================================================
/**
 * @brief Extracts the gamma energy value from a filename
 * @param filename Filename string (e.g., "511.0.dat")
 * @return Gamma energy in keV, or -1.0 if parsing fails
 * 
 * The program expects data files to be named by their dominant gamma energy.
 * This function parses the filename to extract that energy value.
 * 
 * Examples:
 *   "511.0.dat"    -> 511.0
 *   "306.82.dat"   -> 306.82
 *   "1332.492.bin" -> 1332.492
 *   "invalid.dat"  -> -1.0
 */
double extractEnergyFromFilename(const char* filename)
{
    string fname(filename);
    
    // Find and remove the file extension (.dat or .bin)
    size_t dotPos = fname.rfind(".dat");
    if (dotPos == string::npos) {
        dotPos = fname.rfind(".bin");
    }
    if (dotPos == string::npos) {
        return -1.0;  // No valid extension found
    }
    
    // Extract the part before the extension
    string energyStr = fname.substr(0, dotPos);
    
    // Attempt to convert to double
    try {
        return stod(energyStr);
    } catch (...) {
        return -1.0;  // Conversion failed
    }
}

//=============================================================================
// FUNCTION: Create directory if it doesn't exist
//=============================================================================
/**
 * @brief Creates a directory if it doesn't already exist
 * @param path Directory path to create
 * 
 * Uses POSIX mkdir() with permissions 0755 (rwxr-xr-x).
 * Silently does nothing if directory already exists.
 */
void createDirectory(const char* path)
{
    struct stat st = {0};
    if (stat(path, &st) == -1) {
        mkdir(path, 0755);
    }
}

//=============================================================================
// FUNCTION: Process single binary data file
//=============================================================================
/**
 * @brief Processes a single binary measurement file and extracts peak data
 * @param filePath Full path to the binary data file
 * @param realEnergy Known gamma energy in keV for this source
 * @param xmin Minimum ToT ADC value to include
 * @param xmax Maximum ToT ADC value to include
 * @param sigP Peak finding sigma parameter for TSpectrum
 * @param sigB Background estimation sigma parameter
 * @param peakLT Minimum peak area threshold
 * @param peakHT Maximum peak area threshold
 * @param peakResults Output vector to store peak data for all channels
 * @param outputDir Directory to save spectrum plots
 * @param m_lut LFSR lookup table for decoding
 * @return true if processing succeeded, false otherwise
 * 
 * This function performs the complete peak analysis workflow:
 * 1. Read binary data and decode events
 * 2. Fill energy (ToT) histograms for each channel
 * 3. Find peaks using TSpectrum
 * 4. Subtract background and fit Gaussian to highest peak
 * 5. Store results and generate diagnostic plots
 */
bool processDatFile(const char* filePath, double realEnergy, 
                    double xmin, double xmax, 
                    double sigP, double sigB, 
                    double peakLT, double peakHT,
                    vector<PeakData>& peakResults, 
                    const char* outputDir,
                    int16_t* m_lut)
{
    cout << "\n======================================================" << endl;
    cout << "Processing: " << filePath << endl;
    cout << "Real Energy: " << realEnergy << " keV" << endl;
    cout << "======================================================" << endl;
    
    //-------------------------------------------------------------------------
    // STEP 1: Open file and count total events
    //-------------------------------------------------------------------------
    FILE *ptr = fopen(filePath, "rb");
    if (!ptr) {
        cerr << "Error: Cannot open file: " << filePath << endl;
        return false;
    }
    
    // Each event is 8 bytes in the binary file
    const int BYTES_PER_EVENT = 8;
    unsigned char buffer[BYTES_PER_EVENT];
    
    // Count total events by reading through the file
    int eventsTOTAL = 0;
    while (fread(buffer, BYTES_PER_EVENT, 1, ptr) == 1) {
        eventsTOTAL++;
    }
    cout << "Total events: " << eventsTOTAL << endl;
    fclose(ptr);
    
    // Reopen file for actual processing
    ptr = fopen(filePath, "rb");
    if (!ptr) {
        cerr << "Error: Cannot reopen file: " << filePath << endl;
        return false;
    }
    
    //-------------------------------------------------------------------------
    // STEP 2: Create histograms for each channel
    //-------------------------------------------------------------------------
    // Energy histograms for raw data
    TH1F *hE[N_GMSL][N_STIC][N_CHANNEL];
    // Background estimation histograms
    TH1F *hBkg[N_GMSL][N_STIC][N_CHANNEL];
    // Signal histograms (after background subtraction)
    TH1F *hSig[N_GMSL][N_STIC][N_CHANNEL];
    
    char hname[100];
    for (int g = 0; g < N_GMSL; g++) {
        for (int s = 0; s < N_STIC; s++) {
            for (int ch = 0; ch < N_CHANNEL; ch++) {
                sprintf(hname, "hE_G%dS%dCh%d_E%.1f", g, s, ch, realEnergy);
                hE[g][s][ch] = new TH1F(hname, "Energy", NBINS_ENERGY, 0, MAX_ENERGY_ADC);
                
                sprintf(hname, "hBkg_G%dS%dCh%d_E%.1f", g, s, ch, realEnergy);
                hBkg[g][s][ch] = new TH1F(hname, "Background", NBINS_ENERGY, 0, MAX_ENERGY_ADC);
                
                sprintf(hname, "hSig_G%dS%dCh%d_E%.1f", g, s, ch, realEnergy);
                hSig[g][s][ch] = new TH1F(hname, "Signal", NBINS_ENERGY, 0, MAX_ENERGY_ADC);
            }
        }
    }
    
    // 2D maps for visualization
    TH2F *hSingle0 = new TH2F("hSingle0", "Singles Board 0;X index;Y index", 32, 0, 32, 16, 0, 16);
    TH2F *hSingle1 = new TH2F("hSingle1", "Singles Board 1;X index;Y index", 32, 0, 32, 16, 0, 16);
    TH2F *hPeak0 = new TH2F("hPeak0", "Peak Area Board 0;X index;Y index", 32, 0, 32, 16, 0, 16);
    TH2F *hPeak1 = new TH2F("hPeak1", "Peak Area Board 1;X index;Y index", 32, 0, 32, 16, 0, 16);
    TH2F *hRes0 = new TH2F("hRes0", "Resolution Board 0;X index;Y index", 32, 0, 32, 16, 0, 16);
    TH2F *hRes1 = new TH2F("hRes1", "Resolution Board 1;X index;Y index", 32, 0, 32, 16, 0, 16);
    TH1F *hResAll = new TH1F("hResAll", "Resolution Distribution;Resolution (FWHM/mean);Channels", 100, 0.0, 0.5);
    
    //-------------------------------------------------------------------------
    // STEP 3: Read binary data and fill histograms
    //-------------------------------------------------------------------------
    // Binary data format (8 bytes per event):
    // Bytes 0-2: Energy coarse counter (LFSR encoded)
    // Bytes 3-5: Time coarse counter (LFSR encoded)
    // Byte 5: Also contains channel number (upper 6 bits)
    // Byte 6: Contains STiC ID (upper 2 bits)
    // Byte 7: Contains GMSL ID (lower 2 bits) and board ID
    
    unsigned int ecoarse, tcoarse, energy;
    unsigned int channel, sticID, gmslID;
    
    while (fread(buffer, BYTES_PER_EVENT, 1, ptr) == 1)
    {
        // Check for sync packet (not data) - skip if detected
        // Sync packets have a specific pattern in the data
        if ((buffer[0]/4 + buffer[7]*256/4 + (buffer[6]%2)*256*256/4) == 32767) {
            continue;  // Skip sync packets
        }
        
        // Decode energy coarse counter using LFSR lookup
        // Extract 15 bits from bytes 0-2: shift and combine
        ecoarse = m_lut[buffer[0]/32 + buffer[1]*256/32 + (buffer[2]%16)*256*256/32];
        
        // Decode time coarse counter using LFSR lookup
        // Extract 15 bits from bytes 3-5: shift and combine
        tcoarse = m_lut[buffer[3]/4 + buffer[4]*256/4 + (buffer[5]%2)*256*256/4];
        
        // Extract channel number (6 bits from byte 5, shifted right by 2)
        channel = buffer[5] / 4;  // Equivalent to buffer[5] >> 2, but keeping lower 6 bits
        
        // Extract STiC ID (2 bits from byte 6, upper bits)
        sticID = buffer[6] / 64;  // Equivalent to buffer[6] >> 6
        
        // Extract GMSL ID (2 bits from byte 7, lower bits)
        gmslID = buffer[7] % 4;   // Equivalent to buffer[7] & 0x03
        
        // Calculate energy as Time-over-Threshold (ToT)
        // ToT = Ecoarse - Tcoarse, with wrap-around handling
        if (ecoarse >= tcoarse) {
            energy = ecoarse - tcoarse;
        } else {
            // Handle counter wrap-around (15-bit counter max = 32766)
            energy = ecoarse + 32766 - tcoarse;
        }
        
        // Apply energy window cut
        if (energy < xmin || energy > xmax) {
            continue;  // Skip events outside energy window
        }
        
        // Validate channel indices
        if (gmslID >= N_GMSL || sticID >= N_STIC || channel >= N_CHANNEL) {
            continue;  // Skip invalid channel addresses
        }
        
        // Fill energy histogram for this channel
        hE[gmslID][sticID][channel]->Fill(energy);
        
        // Fill 2D singles map (board 0 = GMSL 0,1; board 1 = GMSL 2,3)
        int px = mapChannelToX(channel, sticID, gmslID);
        int py = mapChannelToY(channel, sticID, gmslID);
        
        if (gmslID < 2) {
            hSingle0->Fill(px, py);
        } else {
            hSingle1->Fill(px, py);
        }
    }
    fclose(ptr);
    
    //-------------------------------------------------------------------------
    // STEP 4: Peak finding and fitting for each channel
    //-------------------------------------------------------------------------
    // Peak amplitude and position arrays for intermediate storage
    double peakPosition[N_GMSL][N_STIC][N_CHANNEL];
    double peakAmplitude[N_GMSL][N_STIC][N_CHANNEL];
    
    // Factor to distinguish between main peak and secondary peaks (e.g., 1274 keV from Na-22)
    const double SECONDARY_PEAK_FACTOR = 1.7;
    
    for (int g = 0; g < N_GMSL; g++) {
        for (int s = 0; s < N_STIC; s++) {
            cout << "Processing G" << g << "S" << s << "..." << endl;
            
            for (int ch = 0; ch < N_CHANNEL; ch++) {
                // Skip empty channels
                if (hE[g][s][ch]->GetEntries() < 10) {
                    continue;
                }
                
                //-------------------------------------------------------------
                // 4a. Find peaks using TSpectrum
                //-------------------------------------------------------------
                TSpectrum *spectrum = new TSpectrum();
                
                // Search for peaks with threshold controlled by sigP
                // Lower sigP = more sensitive, finds more peaks
                int nPeaks = spectrum->Search(hE[g][s][ch], 1, "", sigP);
                
                if (nPeaks == 0) {
                    delete spectrum;
                    continue;
                }
                
                // Get peak positions and amplitudes
                Double_t *peakX = spectrum->GetPositionX();
                Double_t *peakY = spectrum->GetPositionY();
                
                // Store peaks in array for sorting
                int parray[nPeaks][2];
                for (int p = 0; p < nPeaks; p++) {
                    parray[p][0] = static_cast<int>(peakX[p]);  // Position
                    parray[p][1] = static_cast<int>(peakY[p]);  // Amplitude
                }
                
                // Sort peaks by position (ascending energy)
                qsort(parray, nPeaks, 2 * sizeof(int), comparePeaksByPosition);
                
                // Highest position peak (typically the one we want for calibration)
                int highestPeakPos = parray[nPeaks - 1][0];
                int highestPeakAmp = parray[nPeaks - 1][1];
                
                // For sources like Na-22, the 511 keV peak may be higher than 1274 keV
                // Find the dominant peak that's not just the highest position
                peakPosition[g][s][ch] = highestPeakPos;
                peakAmplitude[g][s][ch] = highestPeakAmp;
                
                for (int p = nPeaks - 1; p >= 0; p--) {
                    if (parray[p][1] > SECONDARY_PEAK_FACTOR * highestPeakAmp) {
                        peakPosition[g][s][ch] = parray[p][0];
                        peakAmplitude[g][s][ch] = parray[p][1];
                        break;
                    }
                }
                
                //-------------------------------------------------------------
                // 4b. Background subtraction using TSpectrum
                //-------------------------------------------------------------
                Double_t source[NBINS_ENERGY];
                
                // Copy histogram contents to array
                for (int bin = 0; bin < NBINS_ENERGY; bin++) {
                    source[bin] = hE[g][s][ch]->GetBinContent(bin + 1);
                }
                
                // Estimate background using TSpectrum algorithm
                // sigB controls smoothing; larger = smoother background
                spectrum->Background(source, NBINS_ENERGY, static_cast<int>(sigB),
                                    TSpectrum::kBackDecreasingWindow,
                                    TSpectrum::kBackOrder2,
                                    kFALSE,
                                    TSpectrum::kBackSmoothing3,
                                    kFALSE);
                
                // Store background in histogram
                for (int bin = 0; bin < NBINS_ENERGY; bin++) {
                    hBkg[g][s][ch]->SetBinContent(bin + 1, source[bin]);
                }
                
                // Create signal histogram (data - background)
                hSig[g][s][ch]->Add(hE[g][s][ch]);
                hSig[g][s][ch]->Add(hBkg[g][s][ch], -1);
                
                //-------------------------------------------------------------
                // 4c. Fit Gaussian to the highest peak
                //-------------------------------------------------------------
                // Define fit range around peak position
                double fitLeft = peakPosition[g][s][ch] * (1.0 - sigP / 8.0);
                double fitRight = peakPosition[g][s][ch] * (1.0 + sigP / 8.0);
                
                // Create Gaussian function: A * exp(-0.5 * ((x-mean)/sigma)^2)
                TF1 *gaussFit = new TF1("gaussFit", "[0]*exp(-0.5*((x-[1])/[2])^2)", 
                                        fitLeft, fitRight);
                
                // Initial parameter estimates:
                // [0] = amplitude (use peak height)
                // [1] = mean (use peak position)
                // [2] = sigma (estimate from typical resolution ~12% FWHM)
                double sigmaEstimate = peakPosition[g][s][ch] * 0.12 / 2.355;
                gaussFit->SetParameters(peakAmplitude[g][s][ch], 
                                        peakPosition[g][s][ch], 
                                        sigmaEstimate);
                
                // Perform fit on background-subtracted signal
                // Options: R = use function range, Q = quiet (no printout)
                hSig[g][s][ch]->Fit(gaussFit, "RQ");
                
                //-------------------------------------------------------------
                // 4d. Calculate peak area and quality metrics
                //-------------------------------------------------------------
                double fittedMean = gaussFit->GetParameter(1);
                double fittedSigma = abs(gaussFit->GetParameter(2));
                double fittedAmplitude = gaussFit->GetParameter(0);
                
                // Calculate peak area by numerical integration (±3 sigma)
                double integrateStart = fittedMean - 3.0 * fittedSigma;
                double integrateEnd = fittedMean + 3.0 * fittedSigma;
                double stepSize = 6.0 * fittedSigma / 1000.0;
                
                double peakArea = 0.0;
                for (int step = 0; step < 1000; step++) {
                    double x = integrateStart + step * stepSize;
                    double dx = (x - fittedMean) / fittedSigma;
                    peakArea += fittedAmplitude * exp(-0.5 * dx * dx) * stepSize;
                }
                // Convert to bin units (histogram bin width = 4 ADC)
                peakArea = peakArea / 4.0;
                
                // Calculate resolution as FWHM/mean
                double resolution = fittedSigma * 2.355 / fittedMean;
                
                //-------------------------------------------------------------
                // 4e. Fill 2D maps and store results
                //-------------------------------------------------------------
                int px = mapChannelToX(ch, s, g);
                int py = mapChannelToY(ch, s, g);
                
                if (g < 2) {
                    hPeak0->SetBinContent(px + 1, py + 1, peakArea);
                    hRes0->SetBinContent(px + 1, py + 1, resolution);
                } else {
                    hPeak1->SetBinContent(px + 1, py + 1, peakArea);
                    hRes1->SetBinContent(px + 1, py + 1, resolution);
                }
                hResAll->Fill(resolution);
                
                //-------------------------------------------------------------
                // 4f. Store peak data if fit quality is acceptable
                //-------------------------------------------------------------
                // Quality criteria:
                // - Valid number of peaks found (1-7)
                // - Fitted sigma not too large (< 2 × sigB)
                // - Peak area within specified thresholds
                
                bool goodFit = (nPeaks >= 1 && nPeaks <= 7) &&
                              (fittedSigma < sigB * 2) &&
                              (peakArea >= peakLT && peakArea <= peakHT);
                
                if (goodFit) {
                    PeakData pd;
                    pd.gmslID = g;
                    pd.sticID = s;
                    pd.channel = ch;
                    pd.peakArea = peakArea;
                    pd.ToT_ADC = fittedMean;
                    pd.resolution = fittedSigma;
                    pd.realEnergy = realEnergy;
                    peakResults.push_back(pd);
                }
                
                // Clean up
                delete gaussFit;
                delete spectrum;
            }
        }
    }
    
    //-------------------------------------------------------------------------
    // STEP 5: Save diagnostic plots
    //-------------------------------------------------------------------------
    // Create spectrum output directory
    char specDir[500];
    sprintf(specDir, "%s/Spectrum", outputDir);
    createDirectory(specDir);
    
    // Suppress ROOT info messages during plotting
    gErrorIgnoreLevel = kWarning;
    
    // Save spectrum plots for each GMSL/STiC combination
    for (int g = 0; g < N_GMSL; g++) {
        for (int s = 0; s < N_STIC; s++) {
            for (int setNum = 0; setNum < 4; setNum++) {
                // Create canvas with 4×4 subplots (16 channels per page)
                sprintf(hname, "%s/Spectrum/G%dS%d_ChSet%d.png", outputDir, g, s, setNum);
                TCanvas *canvas = new TCanvas("canvas", "Energy Spectra", 1000, 1000);
                canvas->Divide(4, 4, 0, 0);
                
                for (int i = 0; i < 16; i++) {
                    int ch = 16 * setNum + i;
                    canvas->cd(i + 1);
                    
                    int px = mapChannelToX(ch, s, g);
                    int py = mapChannelToY(ch, s, g);
                    
                    char title[100];
                    sprintf(title, "G%dS%d-Ch%d (X%d,Y%d)", g, s, ch, px, py);
                    hE[g][s][ch]->SetTitle(title);
                    hE[g][s][ch]->SetAxisRange(xmin, xmax, "X");
                    hE[g][s][ch]->GetXaxis()->SetLabelSize(0.06);
                    hE[g][s][ch]->GetYaxis()->SetLabelSize(0.06);
                    hE[g][s][ch]->GetXaxis()->SetNdivisions(5);
                    hE[g][s][ch]->GetYaxis()->SetNdivisions(5);
                    
                    // Draw raw data (blue), background (black), signal (green)
                    hE[g][s][ch]->SetLineColor(kBlue);
                    hE[g][s][ch]->Draw();
                    hBkg[g][s][ch]->SetLineColor(kBlack);
                    hBkg[g][s][ch]->Draw("same");
                    hSig[g][s][ch]->SetLineColor(kGreen + 2);
                    hSig[g][s][ch]->Draw("same");
                }
                
                canvas->SaveAs(hname);
                delete canvas;
            }
        }
    }
    
    // Save 2D singles map
    TCanvas *cSingles = new TCanvas("cSingles", "Singles Map", 1000, 1000);
    gStyle->SetPalette(55);
    cSingles->Divide(1, 2);
    cSingles->cd(1);
    hSingle0->Draw("colz");
    cSingles->cd(2);
    hSingle1->Draw("colz");
    sprintf(hname, "%s/Singles_2D.png", outputDir);
    cSingles->SaveAs(hname);
    delete cSingles;
    
    // Save peak area map
    TCanvas *cPeak = new TCanvas("cPeak", "Peak Area Map", 1000, 1000);
    cPeak->Divide(1, 2);
    cPeak->cd(1);
    hPeak0->Draw("colz");
    cPeak->cd(2);
    hPeak1->Draw("colz");
    sprintf(hname, "%s/Peak_2D.png", outputDir);
    cPeak->SaveAs(hname);
    delete cPeak;
    
    // Save resolution map
    TCanvas *cRes = new TCanvas("cRes", "Resolution Map", 1000, 1000);
    cRes->Divide(1, 2);
    cRes->cd(1);
    hRes0->SetMinimum(0);
    hRes0->SetMaximum(0.15);
    hRes0->Draw("colz");
    cRes->cd(2);
    hRes1->SetMinimum(0);
    hRes1->SetMaximum(0.15);
    hRes1->Draw("colz");
    sprintf(hname, "%s/Resolution_2D.png", outputDir);
    cRes->SaveAs(hname);
    delete cRes;
    
    // Save resolution distribution
    TCanvas *cResHist = new TCanvas("cResHist", "Resolution Distribution", 800, 600);
    hResAll->Draw();
    sprintf(hname, "%s/Resolution_Distribution.png", outputDir);
    cResHist->SaveAs(hname);
    delete cResHist;
    
    //-------------------------------------------------------------------------
    // STEP 6: Clean up memory
    //-------------------------------------------------------------------------
    for (int g = 0; g < N_GMSL; g++) {
        for (int s = 0; s < N_STIC; s++) {
            for (int ch = 0; ch < N_CHANNEL; ch++) {
                delete hE[g][s][ch];
                delete hBkg[g][s][ch];
                delete hSig[g][s][ch];
            }
        }
    }
    delete hSingle0;
    delete hSingle1;
    delete hPeak0;
    delete hPeak1;
    delete hRes0;
    delete hRes1;
    delete hResAll;
    
    cout << "Valid peak data extracted: " << peakResults.size() << " channels" << endl;
    return true;
}

//=============================================================================
// FUNCTION: Perform energy calibration fitting
//=============================================================================
/**
 * @brief Combines peak data from multiple sources and fits calibration curves
 * @param allPeakData Vector of peak data vectors, one per source file
 * @param outputDir Directory to save calibration results
 * 
 * This function takes the peak data from multiple gamma sources and fits
 * the exponential calibration function for each detector channel:
 *   ToT(E) = P0 × exp(E/P1) + P2
 * 
 * The fitted parameters allow conversion from measured ToT to energy:
 *   E = P1 × ln((ToT - P2) / P0)
 */
void performEnergyCalibration(vector<vector<PeakData>>& allPeakData, 
                              const char* outputDir)
{
    cout << "\n======================================================" << endl;
    cout << "Performing Energy Calibration Fitting" << endl;
    cout << "======================================================" << endl;
    
    int numSources = allPeakData.size();
    cout << "Number of calibration sources: " << numSources << endl;
    
    if (numSources < 3) {
        cerr << "Warning: At least 3 sources recommended for exponential fit." << endl;
        cerr << "         With " << numSources << " sources, fit quality may be poor." << endl;
    }
    
    //-------------------------------------------------------------------------
    // Create output directories
    //-------------------------------------------------------------------------
    char fittingDir[500];
    sprintf(fittingDir, "%s/fitting", outputDir);
    createDirectory(fittingDir);
    
    //-------------------------------------------------------------------------
    // Organize peak data into multi-dimensional array
    //-------------------------------------------------------------------------
    // Array structure: data[source][gmsl][stic][channel][feature]
    // Features: 0=peakArea, 1=ToT_ADC, 2=resolution, 3=realEnergy
    const int NUM_FEATURES = 4;
    double data[numSources][N_GMSL][N_STIC][N_CHANNEL][NUM_FEATURES];
    
    // Initialize all values to zero
    memset(data, 0, sizeof(data));
    
    // Fill array from peak data vectors
    for (int src = 0; src < numSources; src++) {
        for (const auto& pd : allPeakData[src]) {
            data[src][pd.gmslID][pd.sticID][pd.channel][0] = pd.peakArea;
            data[src][pd.gmslID][pd.sticID][pd.channel][1] = pd.ToT_ADC;
            data[src][pd.gmslID][pd.sticID][pd.channel][2] = pd.resolution;
            data[src][pd.gmslID][pd.sticID][pd.channel][3] = pd.realEnergy;
        }
    }
    
    //-------------------------------------------------------------------------
    // Open output file for calibration parameters
    //-------------------------------------------------------------------------
    char calFileName[500];
    sprintf(calFileName, "%s/Energy_Calibration.txt", outputDir);
    ofstream outFile(calFileName);
    
    outFile << "# ============================================================" << endl;
    outFile << "# FlexPET Energy Calibration Parameters" << endl;
    outFile << "# ============================================================" << endl;
    outFile << "# Calibration function: ToT = P0 × exp(E/P1) + P2" << endl;
    outFile << "# Inverse (for energy reconstruction):" << endl;
    outFile << "#   E = P1 × ln((ToT - P2) / P0)" << endl;
    outFile << "# ============================================================" << endl;
    outFile << "# Columns:" << endl;
    outFile << "#   1: gmslID    - GMSL board ID (0-3)" << endl;
    outFile << "#   2: sticID    - STiC ASIC ID (0-3)" << endl;
    outFile << "#   3: channel   - Channel number (0-63)" << endl;
    outFile << "#   4: P0        - Amplitude parameter" << endl;
    outFile << "#   5: P0_err    - P0 uncertainty" << endl;
    outFile << "#   6: P1        - Decay constant (keV)" << endl;
    outFile << "#   7: P1_err    - P1 uncertainty" << endl;
    outFile << "#   8: P2        - Offset parameter" << endl;
    outFile << "#   9: P2_err    - P2 uncertainty" << endl;
    outFile << "#  10: Chi2      - Fit chi-squared" << endl;
    outFile << "# ============================================================" << endl;
    outFile << "# gmslID\tsticID\tchannel\tP0\tP0_err\tP1\tP1_err\tP2\tP2_err\tChi2" << endl;
    
    //-------------------------------------------------------------------------
    // Create quality assessment histograms
    //-------------------------------------------------------------------------
    TH1D *hChi2Dist = new TH1D("hChi2Dist", "Chi2 Distribution;Chi2;Channels", 
                               100, 0, 2.0);
    TH2D *hChi2Map0 = new TH2D("hChi2Map0", "Chi2 Map Board 0;X index;Y index", 
                               32, 0, 32, 16, 0, 16);
    TH2D *hChi2Map1 = new TH2D("hChi2Map1", "Chi2 Map Board 1;X index;Y index", 
                               32, 0, 32, 16, 0, 16);
    
    // Graphs for storing calibration data points
    TGraphErrors *grCal[N_GMSL][N_STIC][N_CHANNEL];
    
    int goodFitCount = 0;
    char hname[200];
    
    //-------------------------------------------------------------------------
    // Fit calibration curve for each channel
    //-------------------------------------------------------------------------
    for (int g = 0; g < N_GMSL; g++) {
        for (int s = 0; s < N_STIC; s++) {
            for (int ch = 0; ch < N_CHANNEL; ch++) {
                // Initialize fit parameters
                double p0 = 0, p1 = 0, p2 = 0;
                double p0_err = 0, p1_err = 0, p2_err = 0;
                double chi2 = 1000;  // Default to high chi2 (indicates bad fit)
                
                // Create graph for this channel
                grCal[g][s][ch] = new TGraphErrors();
                
                // Collect data points from all sources
                int nPoints = 0;
                for (int src = 0; src < numSources; src++) {
                    // Check if this channel has valid data for this source
                    if (data[src][g][s][ch][0] > 0) {  // peakArea > 0 indicates valid data
                        double energy = data[src][g][s][ch][3];     // Real energy (keV)
                        double tot = data[src][g][s][ch][1];        // ToT (ADC)
                        double tot_err = data[src][g][s][ch][2];    // ToT uncertainty (sigma)
                        
                        grCal[g][s][ch]->SetPoint(nPoints, energy, tot);
                        grCal[g][s][ch]->SetPointError(nPoints, 1.0, tot_err);
                        nPoints++;
                    }
                }
                
                //-------------------------------------------------------------
                // Perform fitting based on number of data points
                //-------------------------------------------------------------
                if (nPoints >= 3) {
                    // Fit exponential function: P0 × exp(E/P1) + P2
                    // Typical initial values for LYSO-based detector
                    TF1 *expFit = new TF1("expFit", "[0]*(exp(x/[1]))+[2]", 190, 1400);
                    expFit->SetParameters(-380.0, -432.0, 565.0);
                    
                    // Perform fit (R=use range, Q=quiet)
                    grCal[g][s][ch]->Fit(expFit, "RQ");
                    
                    // Extract fit parameters and errors
                    p0 = expFit->GetParameter(0);
                    p1 = expFit->GetParameter(1);
                    p2 = expFit->GetParameter(2);
                    p0_err = expFit->GetParError(0);
                    p1_err = expFit->GetParError(1);
                    p2_err = expFit->GetParError(2);
                    chi2 = expFit->GetChisquare();
                    
                    delete expFit;
                    
                    // Update chi2 distribution
                    hChi2Dist->Fill(chi2);
                    
                } else if (nPoints == 2) {
                    // With only 2 points, can only do linear fit
                    // Mark chi2 as invalid (999) to indicate insufficient data
                    TF1 *linFit = new TF1("linFit", "[0]*x+[1]", 190, 1400);
                    grCal[g][s][ch]->Fit(linFit, "RQ");
                    
                    // Store linear fit in P2 (offset), set others to 0
                    p0 = 0;
                    p1 = 0;
                    p2 = linFit->GetParameter(1);
                    chi2 = 999;  // Special value indicating linear fit only
                    
                    delete linFit;
                }
                
                //-------------------------------------------------------------
                // Write results to file
                //-------------------------------------------------------------
                outFile << g << "\t" << s << "\t" << ch << "\t"
                        << p0 << "\t" << p0_err << "\t"
                        << p1 << "\t" << p1_err << "\t"
                        << p2 << "\t" << p2_err << "\t"
                        << chi2 << endl;
                
                // Count good fits
                if (chi2 < 1.0) {
                    goodFitCount++;
                }
                
                // Fill chi2 map
                int px = mapChannelToX(ch, s, g);
                int py = mapChannelToY(ch, s, g);
                
                if (g < 2) {
                    hChi2Map0->SetBinContent(px + 1, py + 1, chi2);
                } else {
                    hChi2Map1->SetBinContent(px + 1, py + 1, chi2);
                }
            }
        }
    }
    
    outFile.close();
    cout << "Calibration file saved: " << calFileName << endl;
    cout << "Good fits (Chi2 < 1.0): " << goodFitCount << " channels" << endl;
    
    //-------------------------------------------------------------------------
    // Save fitting plots
    //-------------------------------------------------------------------------
    gErrorIgnoreLevel = kWarning;  // Suppress ROOT info messages
    
    for (int g = 0; g < N_GMSL; g++) {
        for (int s = 0; s < N_STIC; s++) {
            for (int setNum = 0; setNum < 4; setNum++) {
                sprintf(hname, "%s/fitting/G%dS%d_ChSet%d.png", outputDir, g, s, setNum);
                TCanvas *canvas = new TCanvas("canvas", "Calibration Curves", 1000, 1000);
                canvas->Divide(4, 4, 0, 0);
                
                for (int i = 0; i < 16; i++) {
                    int ch = 16 * setNum + i;
                    canvas->cd(i + 1);
                    canvas->cd(i + 1)->SetGrid(1);
                    
                    int px = mapChannelToX(ch, s, g);
                    int py = mapChannelToY(ch, s, g);
                    
                    char title[100];
                    sprintf(title, "G%dS%d-Ch%d (X%d,Y%d)", g, s, ch, px, py);
                    grCal[g][s][ch]->SetTitle(title);
                    grCal[g][s][ch]->GetXaxis()->SetTitle("Energy (keV)");
                    grCal[g][s][ch]->GetYaxis()->SetTitle("ToT (ADC)");
                    grCal[g][s][ch]->GetXaxis()->SetLabelSize(0.06);
                    grCal[g][s][ch]->GetYaxis()->SetLabelSize(0.06);
                    grCal[g][s][ch]->GetXaxis()->SetRangeUser(0, 1500);
                    grCal[g][s][ch]->GetYaxis()->SetRangeUser(200, 700);
                    grCal[g][s][ch]->SetMarkerStyle(22);
                    grCal[g][s][ch]->SetMarkerSize(1.1);
                    grCal[g][s][ch]->SetMarkerColor(kBlue);
                    grCal[g][s][ch]->Draw("AP");
                }
                
                canvas->SaveAs(hname);
                delete canvas;
            }
        }
    }
    
    //-------------------------------------------------------------------------
    // Save chi2 distribution histogram
    //-------------------------------------------------------------------------
    double meanChi2 = hChi2Dist->GetMean();
    
    sprintf(hname, "Chi2 Distribution (Good fits: %d channels)", goodFitCount);
    TCanvas *cChi2Hist = new TCanvas("cChi2Hist", "Chi2 Distribution", 800, 600);
    cChi2Hist->SetGrid(1);
    hChi2Dist->SetTitle(hname);
    hChi2Dist->Draw();
    sprintf(hname, "%s/Chi2_distribution.png", outputDir);
    cChi2Hist->SaveAs(hname);
    sprintf(hname, "%s/Chi2_distribution.pdf", outputDir);
    cChi2Hist->SaveAs(hname);
    delete cChi2Hist;
    
    //-------------------------------------------------------------------------
    // Save chi2 2D maps
    //-------------------------------------------------------------------------
    TCanvas *cChi2Map = new TCanvas("cChi2Map", "Chi2 2D Maps", 1000, 1000);
    gStyle->SetPalette(55);
    cChi2Map->Divide(1, 2);
    
    cChi2Map->cd(1);
    hChi2Map0->SetMinimum(0);
    hChi2Map0->SetMaximum(meanChi2 * 5);
    hChi2Map0->Draw("colz");
    
    cChi2Map->cd(2);
    hChi2Map1->SetMinimum(0);
    hChi2Map1->SetMaximum(meanChi2 * 5);
    hChi2Map1->Draw("colz");
    
    sprintf(hname, "%s/Chi2_2D.png", outputDir);
    cChi2Map->SaveAs(hname);
    sprintf(hname, "%s/Chi2_2D.pdf", outputDir);
    cChi2Map->SaveAs(hname);
    delete cChi2Map;
    
    //-------------------------------------------------------------------------
    // Clean up memory
    //-------------------------------------------------------------------------
    delete hChi2Dist;
    delete hChi2Map0;
    delete hChi2Map1;
    
    for (int g = 0; g < N_GMSL; g++) {
        for (int s = 0; s < N_STIC; s++) {
            for (int ch = 0; ch < N_CHANNEL; ch++) {
                delete grCal[g][s][ch];
            }
        }
    }
}

//=============================================================================
// MAIN PROGRAM
//=============================================================================
/**
 * @brief Main entry point for the combined energy calibration program
 * @param argc Number of command-line arguments
 * @param argv Array of command-line argument strings
 * @return 0 on success, -1 on error
 */
int main(int argc, char *argv[])
{
    //-------------------------------------------------------------------------
    // Print program header
    //-------------------------------------------------------------------------
    cout << "\n================================================================" << endl;
    cout << "  Flex Combined Energy Calibration Program - Version 2" << endl;
    cout << "================================================================\n" << endl;
    
    //-------------------------------------------------------------------------
    // Check command-line arguments
    //-------------------------------------------------------------------------
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " <InputFolder> [options]" << endl;
        cout << "\nRequired:" << endl;
        cout << "  InputFolder : Folder containing .dat/.bin files named by gamma energy" << endl;
        cout << "                Example: 511.0.dat, 306.82.dat, 661.657.dat" << endl;
        cout << "\nOptional parameters:" << endl;
        cout << "  [2] xmin    : Minimum ToT ADC value (default: 100)" << endl;
        cout << "  [3] xmax    : Maximum ToT ADC value (default: 800)" << endl;
        cout << "  [4] peakLT  : Minimum peak area threshold (default: 100)" << endl;
        cout << "  [5] peakHT  : Maximum peak area threshold (default: 1000000)" << endl;
        cout << "\nNote: Peak finding parameters (sigP, sigB) are automatically set" << endl;
        cout << "      based on the gamma energy in each filename:" << endl;
        cout << "        > 1000 keV : sigP=0.1, sigB=3 (e.g., Co-60)" << endl;
        cout << "        600-1000   : sigP=0.6, sigB=7 (e.g., Cs-137)" << endl;
        cout << "        400-600    : sigP=0.5, sigB=7 (e.g., Na-22 511keV)" << endl;
        cout << "        < 400 keV  : sigP=0.6, sigB=7 (e.g., Ba-133)" << endl;
        cout << "\nExample:" << endl;
        cout << "  ./FlexPET_EnergyCal_v2 CalibrationData 100 800" << endl;
        return -1;
    }
    
    //-------------------------------------------------------------------------
    // Parse command-line arguments
    //-------------------------------------------------------------------------
    char inputFolder[500];
    sprintf(inputFolder, "%s", argv[1]);
    
    double xmin   = (argc >= 3) ? strtod(argv[2], NULL) : 100;
    double xmax   = (argc >= 4) ? strtod(argv[3], NULL) : 800;
    double peakLT = (argc >= 5) ? strtod(argv[4], NULL) : 100;
    double peakHT = (argc >= 6) ? strtod(argv[5], NULL) : 1000000;
    
    cout << "Configuration:" << endl;
    cout << "  Input folder        : " << inputFolder << endl;
    cout << "  ToT ADC range       : " << xmin << " - " << xmax << endl;
    cout << "  Peak area threshold : " << peakLT << " - " << peakHT << endl;
    cout << "  Peak parameters     : Auto-selected based on gamma energy" << endl;
    
    //-------------------------------------------------------------------------
    // Initialize LFSR lookup table for binary data decoding
    //-------------------------------------------------------------------------
    int16_t lfsrLUT[LFSR_SIZE];
    initLFSR(lfsrLUT);
    cout << "\nLFSR lookup table initialized (" << LFSR_SIZE << " entries)" << endl;
    
    //-------------------------------------------------------------------------
    // Create output directory for calibration results
    //-------------------------------------------------------------------------
    char calOutputDir[500];
    sprintf(calOutputDir, "%s_ECal", inputFolder);
    createDirectory(calOutputDir);
    cout << "Output directory: " << calOutputDir << endl;
    
    //-------------------------------------------------------------------------
    // Scan input folder for data files
    //-------------------------------------------------------------------------
    vector<string> dataFiles;
    vector<double> gammaEnergies;
    
    cout << "\nScanning for data files..." << endl;
    
    DIR *dir = opendir(inputFolder);
    if (dir == NULL) {
        cerr << "Error: Cannot open directory: " << inputFolder << endl;
        return -1;
    }
    
    struct dirent *entry;
    while ((entry = readdir(dir)) != NULL) {
        string filename(entry->d_name);
        
        // Check for .dat or .bin extension
        if (filename.length() > 4) {
            string ext = filename.substr(filename.length() - 4);
            if (ext == ".dat" || ext == ".bin") {
                // Try to extract gamma energy from filename
                double energy = extractEnergyFromFilename(entry->d_name);
                if (energy > 0) {
                    dataFiles.push_back(filename);
                    gammaEnergies.push_back(energy);
                    cout << "  Found: " << filename << " -> " << energy << " keV" << endl;
                }
            }
        }
    }
    closedir(dir);
    
    //-------------------------------------------------------------------------
    // Validate input files
    //-------------------------------------------------------------------------
    if (dataFiles.empty()) {
        cerr << "\nError: No valid data files found in " << inputFolder << endl;
        cerr << "Files must be named by their gamma energy (e.g., 511.0.dat)" << endl;
        return -1;
    }
    
    cout << "\nFound " << dataFiles.size() << " data file(s)" << endl;
    
    if (dataFiles.size() < 3) {
        cout << "Warning: For optimal calibration, use at least 3 different gamma sources." << endl;
    }
    
    //-------------------------------------------------------------------------
    // Process each data file
    //-------------------------------------------------------------------------
    vector<vector<PeakData>> allPeakData;
    
    for (size_t i = 0; i < dataFiles.size(); i++) {
        // Build full file path
        char filePath[500];
        sprintf(filePath, "%s/%s", inputFolder, dataFiles[i].c_str());
        
        // Create output directory for this file's analysis
        string baseName = dataFiles[i].substr(0, dataFiles[i].length() - 4);
        char analysisDir[500];
        sprintf(analysisDir, "%s/%s_Analysis", inputFolder, baseName.c_str());
        createDirectory(analysisDir);
        
        // Get optimal peak finding parameters based on gamma energy
        double sigP, sigB;
        getOptimalPeakParameters(gammaEnergies[i], sigP, sigB);
        cout << "\nUsing parameters for " << gammaEnergies[i] << " keV: "
             << "sigP=" << sigP << ", sigB=" << sigB << endl;
        
        // Process the file
        vector<PeakData> peakResults;
        bool success = processDatFile(filePath, gammaEnergies[i], 
                                      xmin, xmax, sigP, sigB, 
                                      peakLT, peakHT,
                                      peakResults, analysisDir, lfsrLUT);
        
        if (success && !peakResults.empty()) {
            allPeakData.push_back(peakResults);
            
            // Save peak results to text file
            char peakFileName[500];
            sprintf(peakFileName, "%s/HighestPeak.txt", analysisDir);
            ofstream peakOut(peakFileName);
            
            peakOut << "# Peak analysis results for " << dataFiles[i] << endl;
            peakOut << "# Gamma energy: " << gammaEnergies[i] << " keV" << endl;
            peakOut << "# gmslID\tsticID\tchannel\tpeakArea\tToT_ADC\tresolution\trealEnergy" << endl;
            
            for (const auto& pd : peakResults) {
                peakOut << pd.gmslID << "\t" 
                       << pd.sticID << "\t" 
                       << pd.channel << "\t"
                       << pd.peakArea << "\t" 
                       << pd.ToT_ADC << "\t" 
                       << pd.resolution << "\t" 
                       << pd.realEnergy << endl;
            }
            peakOut.close();
        }
    }
    
    //-------------------------------------------------------------------------
    // Perform energy calibration fitting
    //-------------------------------------------------------------------------
    if (allPeakData.size() >= 2) {
        performEnergyCalibration(allPeakData, calOutputDir);
    } else {
        cerr << "\nError: Need at least 2 valid data files for calibration." << endl;
        cerr << "Only " << allPeakData.size() << " file(s) produced valid peak data." << endl;
    }
    
    //-------------------------------------------------------------------------
    // Print summary
    //-------------------------------------------------------------------------
    cout << "\n================================================================" << endl;
    cout << "  Energy Calibration Complete!" << endl;
    cout << "================================================================" << endl;
    cout << "\nOutput locations:" << endl;
    cout << "  Spectrum analysis : " << inputFolder << "/<energy>_Analysis/" << endl;
    cout << "  Calibration data  : " << calOutputDir << "/" << endl;
    cout << "\nMain output file:" << endl;
    cout << "  " << calOutputDir << "/Energy_Calibration.txt" << endl;
    
    return 0;
}
