# VLPCalibration

Programs can calibrate VLP systems that use PD-based RSS.

One can use either the command `cmake .` or vs code to build these programs, and then use the command  `./run.sh` or vs code to run.

To this date, there are two programs: `calibration.cc` and `calibration_lamber.cc`.

`calibration.cc` uses the model ![](http://latex.codecogs.com/svg.latex?d=aRSS^b), which supposes the PD is in a horizontal position and the height difference is constant. Please refer to [the paper](https://ieeexplore.ieee.org/document/9728724/).

`calibration_lamber.cc` uses the model ![](http://latex.codecogs.com/svg.latex?RSS=a(m+1)cos^m(\theta)cos^M(\phi)/d^2), which is the complete Lambertian formulation.
