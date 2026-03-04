# VELMAP

**VELMAP:** Referencing multi-frame InSAR velocity results to a GNSS reference frame using a meshing approach. The original implementation was developed by Wang and Wright (2012) and further expanded in Wang et al. (2019). 

## Updates in this repository

- Support for **SBOI (along-track) InSAR inputs**, allowing joint inversion with standard **LoS InSAR** for referencing to a GNSS dataset.
- **Iterative processing over multiple time windows** to estimate cumulative deformation fields, in addition to scalar velocity estimation.
- Additional plotting utilities, including:
  - deformation standard deviation
  - full-resolution InSAR outputs
- Improved execution of the main routine (`main(config.inv)`), allowing **queue-based job submission** for large-scale processing.

## Contact

For questions or collaboration, contact:

**Muhammet Nergizci**  
University of Leeds  
Email: eemne@leeds.ac.uk

## Related Work

Nergizci et al., 2025.  
*Interseismic and Postseismic Deformation of the 2023 Kahramanmaraş Earthquakes from Burst Overlap Interferometry (BOI).*  
AGU Fall Meeting 2025. (Manuscript in preparation)

## References

Wang, H., & Wright, T. J., 2012.  
Satellite geodetic imaging reveals internal deformation of western Tibet.  
*Geophysical Research Letters*, 39(7).  
https://doi.org/10.1029/2012GL051222

Wang, H., Wright, T. J., Liu-Zeng, J., & Peng, L., 2019.  
Strain rate distribution in South-Central Tibet from two decades of InSAR and GPS.  
*Geophysical Research Letters*, 46(10), 5170–5179.  
https://doi.org/10.1029/2019GL081916

---

> **Note**
>
> This repository is a temporary development repository. After the code is finalized, all updates will be transferred to the official COMET repository:  
> https://github.com/nerc-comet/velmap
<!-- TODO: push the updates to COMET-VELMAP -->
