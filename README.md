# Using Machine Learning to Enhance Hyperparameter Optimization in Pandemic Modeling: Case study of COVID-19 Dynamics in Ghana

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20298147.svg)](https://zenodo.org/records/20298148)

This repository contains information and code to reproduce the results presented
in the article
```bibtex
@online{AIM2026,
      title={Using Machine Learning to Enhance Hyperparameter Optimization in Pandemic Modeling: Case study of COVID-19 Dynamics in Ghana}, 
      author={Thomas Izgin and Andreas Meister and Isaac Azure},
      year={2026},
      eprint={TODO},
      archivePrefix={arXiv},
      primaryClass={math.NA},
      url={https://arxiv.org/abs/TODO}, 
}
```

If you find these results useful, please cite the article mentioned above. If you
use the implementations provided here, please **also** cite this repository as
```bibtex
@misc{AIM202&repository,
  title={Reproducibility repository for
         "Using Machine Learning to Enhance Hyperparameter Optimization in Pandemic Modeling: Case study of COVID-19 Dynamics in Ghana"},
  author={Izgin, Thomas and Meister, Andreas and Azure, Isaac},
  year={2026},
  howpublished={\url{https://github.com/IzginThomas/COVID-19}},
  doi={10.5281/zenodo.20298147}
}
```

## Abstract

  In this study, five distinct COVID-19 models developed in different countries, each designed to reflect the prevailing epidemiological condition at the time of formulation, are examined. The models are reformulated while still maintaining their original structure, using their common transmissions from one compartment to the other. Modified Patankar--Runge--Kutta (MPRK) methods are then applied to approximate the solutions of the resulting system of nonlinear ordinary differential equations (ODEs) representing each model to produce unconditionally positive approximations and to preserve the conservative part of the ODEs. In particular, we incorporate the numerical solution into a cost function to improve the estimates for the non-autonomous model hyperparameters. In a first step we obtain piecewise constant parameters that fit real data. Later we perform a WENO reconstruction in a post-process to approximate the true time-dependent coefficients inside the ODEs. As a proof-of-concept, we apply our approach to improve the parameters of a paper concerned with modeling COVID-19 in Ghana, where we can make 5-day predictions within a 10% error range.


## Numerical experiments

To reproduce the numerical experiments presented in this article, you need
to install [Julia](https://julialang.org/downloads/).
The numerical experiments presented in this article were performed using
Julia v1.12.1.

First, you need to download this repository, e.g., by cloning it with `git`
or by downloading an archive via the GitHub interface. Then, you need to start
julia in the `code` directory of this repository and follow the instructions
described in the `README.md` file therein.


## Authors
- [Thomas Izgin](https://uni-kassel.de/go/izgin) (University of Kassel, Germany)
- [Andreas Meister](https://www.uni-kassel.de/fb10/institute/mathematik/arbeitsgruppen/analysis-und-angewandte-mathematik/prof-dr-andreas-meister/team/detailseite.html?tx_ukpersons_personfunctiondetail%5Bcfpid%5D=41557&tx_ukpersons_personfunctiondetail%5BpersonFunction%5D=129&cHash=fb750321314a2c0920ef8cbe5db3c57d) (University of Kassel, Germany)
- [Isaac Azure](https://webapps.knust.edu.gh/staff/dirsearch/profile/summary/4d1bc4eb13e6.html) (Kwame Nkrumah University of Science and Technology - KNUST, Ghana)


## License

The code in this repository is published under the MIT license, see the
`LICENSE` file.


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!



# COVID-19



