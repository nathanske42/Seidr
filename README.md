#### SEIDR REPOSITORY
----------------------

A smorgasbord of code for the Seidr concept - an instrument to perform photonic latnern fed kernel nulling at the VLTI.

Kernel nuller code courtesy of Frantz Martinache (https://github.com/fmartinache/knuller-sim), and photonic lantern code courtesy of Barnaby Norris (unpublished).


## Installation

Code has been tested in a clean install in python 3.10
To install all the required packages, run the following command in the terminal:

```bash
pip install -r requirements.txt
```

## Citation

If you found this code useful, please cite the following paper:

```bibtex
@inproceedings{taras2024kernel,
  title={Kernel nulling at VLTI with photonic lanterns for optimal fibre injection},
  author={Taras, Adam K and Norris, Barnaby and Chhabra, Sorabh and Cvetojevic, Nick and Foriel, Vincent and Ireland, Michael and Kraus, Stefan and Leon-Saval, Sergio and Martinache, Frantz and Paul, Jyotirmay and Spaldin, Eckhart and Sweeney, David and Tuthill, Peter},
  booktitle={Optical and Infrared Interferometry and Imaging IX},
  volume={13095},
  pages={242--250},
  year={2024},
  organization={SPIE}
}
```

## Structure

```
seidr
│   README.md
│   requirements.txt    
│
└───black_hole_followup
│   │   reflight_planet_snr_seidr_02.py
│   │   x_simple_contrast_and_sep.py
│   
│   
└───device_sim
│   │   glintcalc.py
│   │   run_glintcalc.py
│   │   run_glintcalc_vlti.py
│
│
└───end_to_end_simulator
│   │   lanternfibre.py
│   │   x_phase_screens.py
│
│
└───null_depth_sims
│   │   main.py
│   │   SeidrSim.py
│   │   correlatednoise.py
│   │   x_zernike_to_LP.py
│   │   ...
│   │
│   └───knuller-sim
│   │   │   ...
│   │
│   └───glint_example
│   │   │   ...
│   │
│   │
└───planets
│   │
│   └───exoplanet_catalogue
│   │   │   ...
│   │
│   └───gaia
│   │   │   ...
│   │
│   └───old_star_catalogue
│   │   │   ...
│   │
│   └───planet_2_electric_boogaloo
│   │   │   planet_calcs.py
│   │   │   run_over_catalogue.py
│   │   │   ...
│   │

```

## Script Functions

Script				 | Description											  |
---------------|---------------------------------------------------------------																
main.py				 | Simulates examples of a nuller and a kernel nuller,    |
					     | using planetary and telescope geometric information,   |
					     | as well random amplitude/phase errors, as input.		  |
...............................................................................
planet_calcs.py 	 | Defines the functions needed to calculate explanet     |
					 | separations and contrasts.							  |
...............................................................................
run_over_catalgue.py | Loads the Exoplanet Catalogue from PADC and calculates |
					 | separations and contrasts, as well as error bars, for  |
					 | for each exoplanet. Produces example of Figs. 1 and 5  |
					 | from Taras et. al. 2024.  							  |
