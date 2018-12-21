#  The winning pear
After harvest, the respiration metabolism of pome fruit, apple and pear, still remains active. In order to maintain fruit quality for a long period of time, consumption of oxygen and production of carbon dioxide need to be controlled. In practice, this is done by cold and controlled atmosphere storage: low temperature in combination with a reduced oxygen concentration and a slightly increased carbon dioxide concentration slow down the respiration metabolism. However, suboptimal or extreme storage conditions can cause physiological disorders in fruit. For example, if the oxygen concentration is too low and the carbon dioxide concentration is too high, this can lead to core breakdown in Conference pears (tissue browning around the core, development of cavities). It is believed that this phenomenon is due to altered respiration (a switch from aerobic to anaerobic respiration or fermentation) and gas exchange properties of the tissue (diffusivity of metabolic gasses).

As the exchange of metabolic gasses, such as oxygen and carbon dioxide, is crucial for maintaining normal metabolic/physiological functioning, it is important to study/understand how these gasses are transported and distributed within the fruit structure.

At present, no good methods are available to measure internal gas concentrations in fruit. Therefore, in recent years, a scientific computing approach has been adopted to simulate and predict internal gas concentrations/distributions. Furthermore, this approach allows to study the effect of fruit geometry (shape and size) or controlled storage conditions on local oxygen and carbon diox- ide concentrations, while reducing experimental costs.

[(Source)](/doc/STATEMENT.pdf)


## Results
Plot information: left O2 and right CO2 concentration (last image is the mesh used to compute FEM)

<img src="/matlab/results/disorder-inducing.png?raw=true" width="280"> <img src="/matlab/results/optimal-ca.png?raw=true" width="280"> <img src="/matlab/results/pre-cooling.png?raw=true" width="280"> <img src="/matlab/results/refrigerator.png?raw=true" width="280"> <img src="/matlab/results/shelf-life.png?raw=true" width="280"> <img src="/matlab/results/_mesh_generation.png?raw=true" width="280">

## C++ version
The program is coded in C++. Tu execute it, just run the Makefile (one may want to change the compiler or the flags to suit it to its system). The Eigen library is only constituted of headers - which are included in this git - and needs no installation.
The program will then be compiled and put into the bin folder. Just run it.

To change the parameters and the mesh, one may want to rerun the matlab with another mesh and put it into the imports folder. Furthermore, the problem values can be changed in the perameters.hpp file and need recompilation.

Last detail. The paths are absolute and one may change them in main.cpp. They are suited to my computer as for now.
