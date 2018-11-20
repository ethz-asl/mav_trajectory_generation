=========================
mav_trajectory_generation
=========================

This repository contains tools for polynomial trajectory generation and optimization based on methods described in [1].
These techniques are especially suitable for rotary-wing micro aerial vehicles (MAVs).
This README provides a brief overview of our trajectory generation utilities with some examples.

**Authors**: Markus Achtelik, Michael Burri, Helen Oleynikova, Rik Bähnemann, Marija Popović  
**Maintainer**: Rik Bähnemann, brik@ethz.ch  
**Affiliation**: Autonomous Systems Lab, ETH Zurich  

.. toctree::
   :maxdepth: 3
   :caption: Table of Contents
   :glob:

   pages/*

   api/library_root

Bibliography
============

This implementation is largely based on the work of C. Richter *et al*, who should be cited if this is used in a scientific publication (or the preceding conference papers):  
[1] C. Richter, A. Bry, and N. Roy, “**Polynomial trajectory planning for aggressive quadrotor flight in dense indoor environments,**” in *International Journal of Robotics Research*, Springer, 2016.

.. code-block:: latex

	@incollection{richter2016polynomial,
	  title={Polynomial trajectory planning for aggressive quadrotor flight in dense indoor environments},
	  author={Richter, Charles and Bry, Adam and Roy, Nicholas},
	  booktitle={Robotics Research},
	  pages={649--666},
	  year={2016},
	  publisher={Springer}
	}

Furthermore, the nonlinear optimization features our own extensions, described in:  

Michael Burri, Helen Oleynikova, Markus Achtelik, and Roland Siegwart, “**Real-Time Visual-Inertial Mapping, Re-localization and Planning Onboard MAVs in Previously Unknown Environments**”. In *IEEE Int. Conf. on Intelligent Robots and Systems* (IROS), September 2015.

.. code-block:: latex

	@inproceedings{burri2015real-time,
	  author={Burri, Michael and Oleynikova, Helen and  and Achtelik, Markus W. and Siegwart, Roland},
	  booktitle={Intelligent Robots and Systems (IROS 2015), 2015 IEEE/RSJ International Conference on},
	  title={Real-Time Visual-Inertial Mapping, Re-localization and Planning Onboard MAVs in Unknown Environments},
	  year={2015},
	  month={Sept}
	}


