**Effectiveness and efficiency of pre-season administration of long-acting monoclonal antibodies for infants born to RSV vaccinated mothers: a modelling study**

Julia Mayer, Ayaka Monoi, Fabienne Krauer, Kevin van Zandvoort, Beate Kampmann, Matthieu Domenech de Cellès, Stefan Flasche

Code accompanying the manuscript "Effectiveness and efficiency of pre-season administration of long-acting monoclonal antibodies for infants born to RSV vaccinated mothers: a modelling study".

# Structure of the code base:
The folder **src** contains the code used for the main analysis. It is made of 7 files.

0. **Odin2 model** and **Simulation model** contain the catalytic MSR model used to estimate the proportion of seroconverted children <5 years, stratified by season of birth and age in months. They cannot be run on their own.
1. **Main** fits the *odin2* MSR model above to seroprevalence data using a Bayesian framework and runs diagnostic checks. The *simulation model* is then used to output the proportion of proportion of seroconverted children <5 years at a given age, stratified by season of birth.
2. **VE estimates** estimates the efficacy of the maternal RSV vaccine at a given age using a Bayesian framework. It is based on work by Monoi et al.
3. **mAB IE estimates** estimates the efficacy of nirvsevimab at a given time after administration using a Bayesian framework. It is based on work by Monoi et al.
4. **Outputs** is the main file. It runs the *main* file and uses it to calculate the number of seroconverted children <1 year, stratified by season of birth and age in months, in Germany. Using estimates of disease progression by Mahmud et al., the number of RSV-related hospitalisations and ICU admission are then estimated, stratified by season of birth and age in months. The impact of maternal vaccination and of maternal vaccination + nirsevimab is then modelled using estimates from the *VE estimates* and *mAB IE estimates* files.
5. **Numerical outputs** calculates the numerical estimates presented in the paper from the *outputs*, as well as the numbers used to generate the graphical figures.
6. **Plots** uses the the *numerical outputs* to produce the graphical figures presented in the paper.

The **sensitivity analysis** folder contains the files used for the sensitivity analysis. We evaluated the following:
1. Assume an exponential waning of naturally-derived maternal immunity: files **Main 1M** and **Odin2 model 1M**.
2. Assume that the waning of naturally-derived maternal immunity follows an Erlang-3 distribution: files **Main 3M** and **Odin2 model 3M**.
3. Assume that the proportion of children born with naturally-derived maternal immunity varies by season of birth: files **Main 4pi** and **Odin2 model 4 pi**.
