# HMcode2020 Emulator

Emulated [HMcode2020](https://arxiv.org/abs/2009.01858) non-linear power
spectra for fast weak lensing analysis.

## Installation

To install it, just clone the repository, go to the folder and do 

```bash
pip install . [--user]
```

## Requirements
Required python packages:
* numpy
* scipy
* cosmopower

For tutorials:
* matplotlib


## Usage

```python
import HMcode2020Emu as hmcodeemu


params = {
    'omega_cdm'     :  [0.315],
    'As'            :  [np.exp(3.07)*1.e-10],
    'omega_baryon'  :  [0.05],
    'ns'            :  [0.96],
    'hubble'        :  [0.67],
    'neutrino_mass' :  [0.0],
    'w0'            :  [-1.0],
    'wa'            :  [0.0],
    'log10TAGN'     :  [7.8],
    'z'             :  [0.]
}
emulator = hmcodeemu.Matter_powerspectrum()

k_lin, pk_lin_total = emulator.get_linear_pk(nonu=False, **params)
k_lin, pk_lin_nonu = emulator.get_linear_pk(nonu=True, **params)

k_nonlin, pk_nonlin_total = emulator.get_nonlinear_pk(baryonic_boost=True, **params)
```
Note that for neutrino calculations we assume 2 massless neutrinos and 1 massive neutrino with the mass equal to 'neutrino_mass'.
You can also see an example in the tutorials-folder.

## Parameter ranges
| parameter     | limits                |
| :---:         | :---:                 |
| omega_cdm     | [0.15, 0.6]           |
| omega_baryon  | [0.03, 0.07]          |
| hubble        | [0.5, 0.9]            |
| As            | [0.547e-9, 5.459e-9]  |
| ns            | [0.6, 1.2]            |
| neutrino_mass | [0., 0.5]             |
| w0            | [-1.5, -0.3]          |
| wa            | [-2., 2.]             |
| log10TAGN     | [7.6, 8.3]            |
| z             | [0.0, 3.]             |
| k_lin         | [1e-4, 50]            |
| k_nonlin      | [0.01, 5]             |


## License
[MIT](https://choosealicense.com/licenses/mit/)
