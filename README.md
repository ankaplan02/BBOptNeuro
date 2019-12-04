# Batch Bayesian Optimization for Optimizing a Neurostimulator

This repository consists of functions that conduct analysis and simulations for the Bayesian adaptive clinical trial presented in the article "Batch Bayesian Optimization for Optimizing a Neurostimulator" by Adam Kaplan and Thomas A. Murray. The motivation for this trial is to efficiently optimize a medical device that has one or two dimensional configuration space (i.e., frequency and the pair of frequency and pulse width, respectively) for one patient. In our specific case, optimization was for a patient using an implanted spinal cord stimulation device. There are three components to the device calibration: (1) Efficient exploration of the device configuration space through Batch Bayesian Optimization, (2) a probability model for the latent preferences for the configurations, and (3) two stopping rules to stop the trial early for one of two scenarios: whether the patient's outcomes suggest **that there is preference for a particular configuration at all**, or that **further calibration does not provide improvement** over the currently estimated-to-be most preferred configuration. We refer to these two situations as *preference neutrality* and *calibration convergence*. Further information is provided in the article.   

```

```
