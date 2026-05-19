# Using Machine Learning to Design Time Step Size Controllers for Stable Time Integrators

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17361026.svg)](https://zenodo.org/records/17988179)

This repository contains information and code to reproduce the results presented
in the article
```bibtex
@online{IR2025,
      title={Using Bayesian Optimization to Design Time Step Size Controllers with Application to Modified Patankar--Runge--Kutta Methods}, 
      author={Thomas Izgin and Hendrik Ranocha},
      year={2023},
      eprint={2312.01796},
      archivePrefix={arXiv},
      primaryClass={math.NA},
      url={https://arxiv.org/abs/2312.01796}, 
}
```

If you find these results useful, please cite the article mentioned above. If you
use the implementations provided here, please **also** cite this repository as
```bibtex
@misc{IR2025repository,
  title={Reproducibility repository for
         "Using Machine Learning to Design Time Step Size Controllers for Stable Time Integrators"},
  author={Izgin, Thomas and Ranocha, Hendrik},
  year={2025},
  howpublished={\url{https://github.com/IzginThomas/ML-TimeStepControl}},
  doi={10.5281/zenodo.17988178}
}
```

## Abstract

We present a new method for developing time step controllers based on a technique from the field of machine learning. This method is applicable to stable time integrators that have an embedded scheme, i.e., that have local error estimation similar to Runge-Kutta pairs.
To design good time step size controllers using these	error estimates, we propose to use Bayesian optimization. In particular, we design a novel objective function that captures important properties	such as tolerance convergence and computational stability. We apply our new approach to several modified Patankar--Runge--Kutta (MPRK) schemes and a Rosenbrock-type scheme, equipping them with controllers based	on digital signal processing which extend classical PI and PID controllers.
We demonstrate that the optimization process yields controllers that are at least as good as the best controllers chosen from a wide range of
suggestions available for classical explicit and implicit time integration	methods by providing work-precision diagrams on a variety of ordinary and partial differential equations.


## Numerical experiments

To reproduce the numerical experiments presented in this article, you need
to install [Matlab](https://de.mathworks.com/products/matlab.html).
The numerical experiments presented in this article were performed using
Matlab R2025b.

First, you need to download this repository, e.g., by cloning it with `git`
or by downloading an archive via the GitHub interface. Then, you need to start
Matlab in the `code` directory of this repository and follow the instructions
described in the `README.md` file therein.


## Authors
- [Thomas Izgin](https://uni-kassel.de/go/izgin) (University of Kassel, Germany)
- [Hendrik Ranocha](https://ranocha.de) (Johannes Gutenberg University Mainz, Germany)


## License

The code in this repository is published under the MIT license, see the
`LICENSE` file.


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!



# ML-TimeStepControl



