# Joint Magnetometer-IMU Calibration
This is a MATLAB implementation of the paper [Joint Magnetometer-IMU Calibration via Maximum A Posteriori Estimation](https://www.arxiv.org/abs/2505.16662).
## Install dependencies
- [YALMIP](https://github.com/yalmip/YALMIP)
- [Manopt](https://github.com/NicolasBoumal/manopt)

## Install
```
git clone --recurse-submodules https://github.com/Huang-Chuan/Mag-IMU-JointCalibration.git
```

## Usage
The configuration files in the folder `config.json` can be used to change the parameters of the algorithm.

Run `main.m` to see the results of the experiments. 


## Citation

If you are using MAINS for academic work, please use the following citation:

```bibtex
@ARTICLE{11505718,
  author={Huang, Chuan and Hendeby, Gustaf and Skog, Isaac},
  journal={IEEE Sensors Journal}, 
  title={Joint Magnetometer-IMU Calibration via Maximum A Posteriori Estimation}, 
  year={2026},
  volume={},
  number={},
  pages={1-1},
  keywords={Aerospace and electronic systems;Sensor systems;Feeds;Kalman filters;Filters;Filtering;Central Processing Unit;Nonlinear filters;Circuits and systems;Electronic circuits;inertial sensors;magnetometers;in-situ calibration;MAP estimation;IMU preintegration},
  doi={10.1109/JSEN.2026.3688575}}
```
