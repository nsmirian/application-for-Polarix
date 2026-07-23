# Application for PolariX

MATLAB application for the operation and analysis of the **PolariX X-band Transverse Deflecting Structure (TDS)** at **FLASH**, DESY.

The application provides a graphical user interface for configuring measurements, acquiring beam diagnostics, and analyzing the longitudinal properties of relativistic electron bunches using the PolariX transverse deflecting cavity.

---

# Overview

The PolariX transverse deflecting cavity (TDS) is one of the most important longitudinal diagnostics at FLASH. It converts the temporal distribution of an electron bunch into a transverse spatial distribution on an observation screen, allowing direct measurements of

* longitudinal bunch profile,
* bunch length,
* longitudinal phase space,
* slice energy spread,
* beam arrival time,
* longitudinal beam dynamics.

This application was developed to simplify machine operation and data acquisition through an intuitive MATLAB graphical interface.

---

# Features

The application provides tools for

* operating the PolariX TDS;
* configuring measurement parameters;
* controlling measurement scans;
* visualizing beam images in real time;
* reconstructing longitudinal bunch profiles;
* measuring bunch length;
* analyzing longitudinal phase-space images;
* saving measurement results;
* exporting processed data.

---

# Application

The repository contains the MATLAB App Designer application

```text
Application-for-Polarix.mlapp
```

which provides an integrated graphical interface for PolariX measurements.

---

# Requirements

* MATLAB
* MATLAB App Designer
* Access to the FLASH control system
* Access to PolariX diagnostics
* Appropriate permissions for machine operation

Depending on the local installation, additional MATLAB toolboxes or control-system interfaces may be required.

---

# Installation

Clone the repository

```bash
git clone https://github.com/nsmirian/application-for-Polarix.git
cd application-for-Polarix
```

or download it as a ZIP archive.

Add the repository to the MATLAB path

```matlab
addpath('/path/to/application-for-Polarix')
```

Open

```text
Application-for-Polarix.mlapp
```

from MATLAB App Designer and run the application.

---

# Typical Workflow

## 1. Prepare the machine

Before performing a measurement

* establish stable beam operation;
* verify that the PolariX cavity is operational;
* confirm the observation screen is available;
* verify beam transport through the diagnostic section.

---

## 2. Configure the measurement

Typical parameters include

* TDS voltage;
* RF phase;
* beam energy;
* observation screen;
* camera settings;
* image acquisition;
* averaging parameters.

---

## 3. Acquire beam images

The application records beam images while monitoring

* beam centroid;
* beam size;
* image intensity;
* machine status.

Real-time visualization allows immediate evaluation of the measurement quality.

---

## 4. Longitudinal reconstruction

Using the calibrated TDS streak, the transverse beam image is converted into the longitudinal time coordinate.

The application can be used to reconstruct

* bunch profile,
* bunch duration,
* current profile,
* longitudinal phase space.

---

## 5. Data analysis

The measured data can be analyzed to determine

* RMS bunch length;
* FWHM bunch length;
* peak current;
* longitudinal beam structure;
* slice properties;
* beam stability.

---

# Measurement Principle

A transverse deflecting cavity produces a time-dependent transverse kick

[
y' \propto V_{\mathrm{TDS}}\sin(\omega t+\phi),
]

where

* (V_{\mathrm{TDS}}) is the cavity voltage,
* (\omega) is the RF frequency,
* (\phi) is the RF phase.

Near the zero-crossing phase,

[
y \propto t,
]

allowing the temporal coordinate to be mapped onto the transverse beam profile observed on the screen.

By combining the TDS with a downstream dipole magnet, the full longitudinal phase space can be measured.

---

# Output

Typical output includes

* raw beam images;
* processed beam images;
* bunch profiles;
* longitudinal phase-space distributions;
* bunch-length measurements;
* calibration results;
* exported figures;
* saved measurement data.

---

# Typical Applications

The software can be used for

* injector optimization;
* bunch compression studies;
* RF phase optimization;
* FEL optimization;
* beam arrival-time studies;
* longitudinal phase-space characterization;
* accelerator physics experiments.

---

# Repository Structure

```text
application-for-Polarix/
│
├── Application-for-Polarix.mlapp
├── README.md
└── LICENSE
```

---

# Citation

If this software contributes to published work, please cite

```bibtex
@software{mirian_polarix,
  author = {Najmeh S. Mirian},
  title = {Application for PolariX: MATLAB GUI for X-Band Transverse Deflecting Cavity Diagnostics at FLASH},
  url = {https://github.com/nsmirian/application-for-Polarix},
  note = {MATLAB application for PolariX operation and longitudinal beam diagnostics}
}
```

A DOI should be added when a software release is archived on Zenodo.

---

# Author

**Dr. Najmeh S. Mirian**

Developed for beam diagnostics with the PolariX X-band Transverse Deflecting Structure at FLASH, DESY.

---

# License

This project is distributed under the GNU General Public License v3.0.

See the LICENSE file for details.

---

# Disclaimer

This software is intended for accelerator research and machine studies. Users are responsible for validating the measurements and operating the accelerator in accordance with local machine-protection and safety procedures.
