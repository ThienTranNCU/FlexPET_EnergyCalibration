# FlexPET Full Energy Calibration Program

## Author

Tran Cong Thien / PraanaTech

## Overview

This program performs energy calibration for the dual-head FlexPET detector system. It combines two previously separate workflows into a single automated program:
1. Peak analysis from binary measurement files (BinaryRead9)
2. Energy calibration curve fitting (Read_BinOut_Calibrate)

**Key Feature:** The program automatically selects optimal peak finding parameters based on the gamma energy in each filename, eliminating the need for manual parameter tuning.

## Features

1. **Automatic file discovery**: Scans input folder for `.dat`/`.bin` files named by their gamma energy
2. **Auto-parameter selection**: Determines optimal sigP and sigB based on gamma energy
3. **Peak analysis**: Performs TSpectrum-based peak finding and Gaussian fitting for each channel
4. **Background subtraction**: Applies TSpectrum background estimation for clean peak fitting
5. **Energy calibration**: Fits exponential calibration curve: `ToT = P0 × exp(E/P1) + P2`
6. **Quality metrics**: Generates Chi² maps and distributions for fit quality assessment
7. **Comprehensive outputs**: Saves spectrum plots, 2D maps, and calibration parameters

## File Naming Convention

Data files must be named by their **dominant gamma energy** in keV:
- `511.0.dat` - Na-22 annihilation peak (511 keV)
- `306.82.dat` - Ba-133 peak (306.82 keV)
- `661.657.dat` - Cs-137 peak (661.657 keV)
- `1332.492.dat` - Co-60 peak (1332.492 keV)
- `201.83.dat` - Additional Ba-133 peak (optional)

## Compilation

```bash
g++ FlexPET_EnergyCal_v2.cpp -o FlexPET_EnergyCal_v2 `root-config --cflags --glibs` -lSpectrum
```

## Usage

```bash
./FlexPET_EnergyCal_v2 <InputFolder> [xmin] [xmax] [peakLT] [peakHT]
```

### Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| InputFolder | Folder containing .dat/.bin files | Required |
| xmin | Minimum ToT ADC value | 100 | optional |
| xmax | Maximum ToT ADC value | 800 | optional |
| peakLT | Minimum peak area threshold | 100 | optional |
| peakHT | Maximum peak area threshold | 1000000 | optional |

### Automatic Parameter Selection

The program automatically determines the optimal peak finding parameters (sigP, sigB) based on the gamma energy extracted from each filename:

| Energy Range | sigP | sigB | Typical Sources |
|--------------|------|------|-----------------|
| > 1000 keV | 0.1 | 3 | Co-60 (1173, 1332 keV) |
| 600 - 1000 keV | 0.6 | 7 | Cs-137 (662 keV) |
| 400 - 600 keV | 0.5 | 7 | Na-22 (511 keV) |
| < 400 keV | 0.6 | 7 | Ba-133 (356, 303, 276 keV) |

**Why different parameters?**
- **sigP (peak finding sensitivity)**: Lower values find peaks more precisely but may pick up noise. High-energy peaks (>1 MeV) have cleaner spectra, so lower sigP works well.
- **sigB (background smoothing)**: Higher values create smoother background estimates. Lower energy spectra have more Compton scattering background, requiring higher sigB.

## Example

```bash
# Create folder structure with properly named files
mkdir CalibrationData
cp measurement_ba133.bin CalibrationData/306.82.dat
cp measurement_na22.bin CalibrationData/511.0.dat
cp measurement_cs137.bin CalibrationData/661.657.dat
cp measurement_co60.bin CalibrationData/1332.492.dat

# Run calibration (simplest form - uses all defaults)
./FlexPET_EnergyCal_v2 CalibrationData

# Run with custom ToT window
./FlexPET_EnergyCal_v2 CalibrationData 100 800

# Run with custom ToT window and peak area thresholds
./FlexPET_EnergyCal_v2 CalibrationData 100 800 100 1000000
```

## Output Structure

```
CalibrationData/
├── 306.82_Analysis/
│   ├── Spectrum/
│   │   └── G*S*_ChSet*.png    (16 files per GMSL/STiC combination)
│   ├── Singles_2D.png
│   ├── Peak_2D.png
│   ├── Resolution_2D.png
│   ├── Resolution_Distribution.png
│   └── HighestPeak.txt
├── 511.0_Analysis/
│   └── ... (same structure)
├── 661.657_Analysis/
│   └── ...
└── 1332.492_Analysis/
    └── ...

CalibrationData_ECal/
├── fitting/
│   └── G*S*_ChSet*.png        (calibration curve plots)
├── Energy_Calibration.txt     (main calibration parameters)
├── Chi2_distribution.png
├── Chi2_distribution.pdf
├── Chi2_2D.png
└── Chi2_2D.pdf
```

## Output File Formats

### HighestPeak.txt (per source file)
```
# Peak analysis results for 511.0.dat
# Gamma energy: 511.0 keV
# gmslID  sticID  channel  peakArea  ToT_ADC  resolution  realEnergy
0         0       0        1234.5    456.7    12.3        511.0
...
```

### Energy_Calibration.txt
```
# ============================================================
# FlexPET Energy Calibration Parameters
# ============================================================
# Calibration function: ToT = P0 × exp(E/P1) + P2
# Inverse (for energy reconstruction):
#   E = P1 × ln((ToT - P2) / P0)
# ============================================================
# gmslID  sticID  channel  P0      P0_err  P1      P1_err  P2      P2_err  Chi2
0         0       0        -380.8  12.3    -432.6  15.7    565.6   8.9     0.45
...
```

## Calibration Function

The program fits the ToT-Energy relationship using:

```
ToT(E) = P0 × exp(E/P1) + P2
```

Where:
- **ToT**: Time-over-Threshold in ADC counts
- **E**: Gamma energy in keV
- **P0, P1, P2**: Fitting parameters (channel-dependent)

To convert measured ToT to energy, invert this function:

```
E = P1 × ln((ToT - P2) / P0)
```

## Quality Assessment

- **Chi² < 1.0**: Good calibration fit
- **Chi² 1.0 - 2.0**: Acceptable fit
- **Chi² > 2.0**: Poor fit, check spectrum quality
- **Chi² = 999**: Insufficient data points (only 2 sources available)

The program reports the number of channels with Chi² < 1.0 and generates 2D Chi² maps for visual inspection of problematic channels.

## Requirements

- ROOT 6.x with TSpectrum
- C++11 or later
- Linux/Unix environment

## Detector Geometry

The FlexPET detector uses a hierarchical addressing scheme:
- **GMSL ID (0-3)**: Identifies the GMSL serializer board
- **STiC ID (0-3)**: Identifies the STiC ASIC on each GMSL
- **Channel (0-63)**: Identifies the individual SiPM channel

The program includes channel mapping functions to convert from (GMSL, STiC, Channel) to physical (X, Y) detector coordinates for 2D visualization.

## Comparison with Original Programs

| Feature | BinaryRead9 + Read_BinOut_Calibrate | Combined Program v2 |
|---------|-------------------------------------|---------------------|
| Steps required | 2 separate programs | 1 program |
| File naming | Manual energy input | Auto-extract from filename |
| Parameter selection | Manual per-file | Automatic based on energy |
| Directory creation | Manual | Automatic |
| Parameter passing | Per-file basis | Single set for all files |
| Output organization | Scattered | Structured by analysis |

## Troubleshooting

### No peaks found
- Check if ToT window (xmin, xmax) covers the expected peak positions
- Verify that the file contains valid data from the expected source

### Poor calibration fits (high Chi²)
- Ensure at least 3 different gamma sources are used
- Check individual spectrum plots for detector issues
- Verify file naming matches actual source energies

### Missing channels in output
- Some channels may have insufficient statistics
- Check Singles_2D.png to identify dead or noisy channels


