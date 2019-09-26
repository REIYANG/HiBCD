# Hybrid inexact Block Coordinate Descent (HiBCD)
Matlab code for "Hybrid Inexact BCD for Coupled Structured Matrix Factorization in Hyperspectral Super-Resolution", submitted to IEEE Transaction on Signal Processing, 2019.

## Usage
1. Semi-real dataset experiment: run "Fig2_SAMmap.m", "Fig3_PSNRcurve.m", "Table12_one_shot.m" and "Table3_step_size.m" in folder "semireal experiment"
   * Note: please download the real HS image from the link provided below, and run script "data_generation" to get the data matrix "Chikusei_1080.mat"

2. Synthetic dataset experiment: run "Table45_monte_carlo" in folder "synthetic experiment"

## Reference
Ruiyuan Wu, Hoi-To Wai, and Wing-Kin Ma. "Hybrid Inexact BCD for Coupled Structured Matrix Factorization in Hyperspectral Super-Resolution" [[pdf]](https://arxiv.org/pdf/1909.09183.pdf)

### Abstract
This paper develops a first-order optimization method for coupled structured matrix factorization (CoSMF) problems that arise in the context of hyperspectral super-resolution (HSR) in remote sensing. To best leverage the problem structures for computational efficiency, we introduce a hybrid inexact block coordinate descent (HiBCD) scheme wherein one coordinate is updated via the fast proximal gradient (FPG) method, while another via the Frank-Wolfe (FW) method. The FPG-type methods are known to take less number of iterations to converge, by numerical experience, while the FW-type methods can offer lower per-iteration complexity in certain cases; and we wish to take the best of both. We show that the limit points of this HiBCD scheme are stationary. Our proof treats HiBCD as an optimization framework for a class of multi-block structured optimization problems, and our stationarity claim is applicable not only to CoSMF but also to many other problems. Previous optimization research showed the same stationarity result for inexact block coordinate descent with either FPG or FW updates only. Numerical results indicate that the proposed HiBCD scheme is computationally much more efficient than the state-of-the-art CoSMF schemes in HSR.

## Sources of the tested algorithms
1. Hyperspectral and multispectral data fusion toolbox/CNMF: http://naotoyokoya.com/Download.html

2. FUMI: https://github.com/qw245/FUMI

3. SuperRes-PALM: https://github.com/lanha/SupResPALM

## Source of the real HS image:
1. Airborne hyperspectral data taken over Chikusei: http://naotoyokoya.com/Download.html

### Please advise to remove immediately if any infringement caused.
