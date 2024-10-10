<!-- PROJECT SHIELDS -->
[![arXiv][arxiv-shield]][arxiv-url]
[![MIT License][license-shield]][license-url]
[![Webpage][webpage-shield-NC]][webpage-url-NC]
[![Webpage][webpage-shield-RS]][webpage-url-RS]

# End-to-end guarantees for indirect data-driven control of bilinear systems with finite stochastic data
This repository contains the code from our paper "End-to-end guarantees for indirect data-driven control of bilinear systems with finite stochastic data" which can be acessed [here](https://arxiv.org/abs/2409.18010). 

If you use this project for academic work, please consider citing our publication 

    Chatzikiriakos, N., Strässer, R., Allgöwer, F., & Iannelli, A. (2024). 
    End-to-end guarantees for indirect data-driven control of bilinear systems with finite stochastic data. 
    submitted, Preprint:arXiv:2409.18010.

## Structure
- AnalysisBounds: Code for the analysis of the finite sample indentification error bounds (Section 5.1)
- ControllerDesign: Code for the controller design (Sections 5.2 & 5.3) 

### Analysis of Indentification error bounds
While the main code is provided in python the generated data can be plotted using the corresponding matlab scirpts 
#### Requirements
To install all relevant packages execute 
```bash 
pip install -r requirements.txt
```
#### LSE and error bound estimates
- To recreate Fig.1 a)&b) run main.py once and plot data using plotError_T.m.
- To recreate Fig.1 c)&d) run main.py with different system dimensions and plot data using plotError_nx.m.

### Controller design  
The controller design is carried out in matlab. Each of the exmaples has its own main.m, indicated by the suffix.
#### Requirements 
The code uses the YALMIB toolbox with the solver MOSEK.
## Contact
🧑‍💻 Nicolas Chatzikiriakos - [nicolas.chatzikiriakos@ist.uni-stuttgart.de](mailto:nicolas.chatzikiriakos@ist.uni-stuttgart.de)

🧑‍💻 Robin Strässer - [robin.straesser@ist.uni-stuttgart.de](mailto:robin.straesser@ist.uni-stuttgart.de)


[license-shield]: https://img.shields.io/badge/License-MIT-T?style=flat&color=blue
[license-url]: https://github.com/col-tasas/2024-bilinear-end-to-end/blob/main/LICENSE
[webpage-shield-NC]: https://img.shields.io/badge/Webpage-Nicolas%20Chatzikiriakos-T?style=flat&logo=codementor&color=green
[webpage-url-NC]: https://www.ist.uni-stuttgart.de/institute/team/Chatzikiriakos/
[webpage-shield-RS]: https://img.shields.io/badge/Webpage-Robin%20Strässer-T?style=flat&logo=codementor&color=green
[webpage-url-RS]: https://www.ist.uni-stuttgart.de/institute/team/Straesser/
[arxiv-shield]: https://img.shields.io/badge/arXiv-2409.18010-t?style=flat&logo=arxiv&logoColor=white&color=red
[arxiv-url]: https://arxiv.org/abs/2409.18010
[researchgate-shield-NC]: https://img.shields.io/badge/ResearchGate-Nicolas%20Chatzikiriakos-T?style=flat&logo=researchgate&color=darkgreen
[researchgate-url-NC]: https://www.researchgate.net/profile/Nicolas-Chatzikiriakos
[researchgate-shield-RS]: https://img.shields.io/badge/ResearchGate-Robin%20Strässer-T?style=flat&logo=researchgate&color=darkgreen
[researchgate-url-RS]: https://www.researchgate.net/profile/Robin-Straesser


