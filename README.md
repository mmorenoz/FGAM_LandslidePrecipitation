# Functional regression for space-time prediction of precipitation-induced shallow landslides in South Tyrol, Italy
This is the R code to reproduce the analysis in "Functional regression for space-time prediction of precipitation-induced shallow landslides in South Tyrol, Italy."
The provided script implements functional GAMs for integrating hourly precipitation time series (INCA dataset) and static factors in landslide prediction. Cross-validation routines (random, spatial, and factor) and visualization are also included.


### Abstract
Shallow landslides are geomorphic hazards in mountainous terrains across the globe. Their occurrence can be attributed to the interplay of static and dynamic landslide controls. In previous studies, data-driven approaches have been employed to model shallow landslides on a regional scale, focusing on analyzing the spatial aspects and time-varying conditions separately. Still, the joint assessment of shallow landslides in space and time using data-driven methods remains challenging. This study aims to predict the occurrence of precipitation-induced shallow landslides in space and time (i.e., the ‘where’ and the ‘when’) within the Italian province of South Tyrol (7,400 km²). In this context, we investigate the benefits of considering precipitation leading to landslide events as a functional predictor, in contrast to conventional approaches that treat precipitation (or summary measures thereof) as a scalar predictor. We built upon hourly precipitation analysis data and past landslide occurrences from 2012 to 2021. We implemented a novel functional generalized additive model to establish statistical relationships between the spatiotemporal occurrence of shallow landslides, various static factors included as scalar predictors, and the hourly precipitation pattern preceding a potential landslide event used as a functional predictor. We evaluated the resulting predictions through several stratified cross-validation routines, achieving high model performance scores. To showcase the model capabilities, we performed a hindcast for the storm event in the Passeier Valley on August 4th and 5th, 2016. This novel approach enables the prediction of landslides in space and time for large areas by accounting for static and dynamic functional landslide controls, seasonal effects, statistical uncertainty, and underlying data limitations such as the spatially inconsistent mapping of landslide data.


### Acknowledgements
The research that led to these results is related to the PROSLIDE project (https://www.mountainresearch.at/proslide/), which received funding from the research program Research Südtirol/Alto Adige 2019 of the autonomous province of Bolzano. We thank the Faculty of Geo-information Science and Earth Observation (ITC) – University of Twente for covering the open-access publication fees We thank. We thank Dr. Serkan Girgin for his support with using the CRIB platform. We thank the Office for Meteorology and Avalanche Prevention, especially Mauro Tollardo, for supporting and providing data. Finally, we thank the provincial Office for Geology and Building Materials Testing for assisting in preparing landslide data.


### Further dissemination
Previous work was presented at the EGU General Assembly 2023. The abstract is available at:

> Moreno, M., Steger, S., Lombardo, L., Opitz, T., Crespi, A., Marra, F., Vugt, L. de, Zieher & Westen, C. (2023). Functional regression for space-time prediction of precipitation-induced shallow landslides in South Tyrol, Italy (Nos. EGU23-9538). EGU23. Copernicus Meetings. https://doi.org/10.5194/egusphere-egu23-9538


### Financial support
This research has been supported by the autonomous province of Bolzano (grant no. 9/34), the Faculty of Geo-information Science and Earth Observation (ITC-UTWENTE), and the Biostatistics and Spatial Processes unit (INRAE).
