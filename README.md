# Security Proof for Variable-Length Quantum Key Distribution

This is a public version of the code used in *[Security Proof for Variable-Length Quantum Key Distribution](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.6.023002)* \[[arxiv link](https://arxiv.org/abs/2311.01600)]. This was built for [v2.0.1](https://github.com/Optical-Quantum-Communication-Theory/openQKDsecurity/releases/tag/v2.0.1) of the Open QKD Security package.

The code computes the key rates for fixed-length and variable-length qubit BB84 protocols studied in the paper and generates the plots in the paper. To compute expected key rates, we sample from the honest behaviour of the channel a large number of times. 


## Install instructions
> [!CAUTION]
> This repository is for archival and transparency purposes; we do not guarantee compatibility with other versions of the Open QKD Security package beyond the ones listed above.

### as zip
1. Download the linked version of the code from above and follow all [install instructions](https://github.com/Optical-Quantum-Communication-Theory/openQKDsecurity/commit/bb1c6490c6bffb0661cef52f6b48de41b5e78027).
2. Download the latest release on the side bar and unzip in your preferred directory and add this folder to the Matlab path.
3.  Run `mainQubitAdaptive.m` for Fig.1 of the paper. Run `mainFullyAdaptiveQubitBB84.m` for Fig. 2 of the paper.


### with git
1. Clone this repository and its exact submodules navigate to your desired directory and run,
```
git clone --recurse-submodules https://github.com/Optical-Quantum-Communication-Theory/SecurityProofForVariableLengthQuantumKeyDistribution
```
2. Follow all further [install instructions](https://github.com/Optical-Quantum-Communication-Theory/openQKDsecurity/commit/bb1c6490c6bffb0661cef52f6b48de41b5e78027).
3. Add this repository's folder to the Matlab path.
4. Run `mainQubitAdaptive.m` for Fig.1 of the paper. Run `mainFullyAdaptiveQubitBB84.m` for Fig. 2 of the paper.
