# LB1_project
#### TASK:
Build a Profile Hidden Markov Model for the Kunitz-type protease inhibitor domain.

## Building an Hidden Markov Model for the detection of the Kunitz domain starting from structural data
### Abstract
Kunitz domains are conserved domains largely found in nature, involved in protease inhibition and multiple biological functions, which are highly studied for pharmaceutical purposes.
In this work is proposed a workflow for the construction of an Hidden Markov Model profile to predict the presence of the Kunitz domain inside unseen protein sequences. The model has been generated thanks to the HMMER software and trained on highly resolved structures. It has been optimized using a 2-fold cross-validation procedure and tested on SwissProt sequences, considering as main statistical parameter for the evaluation the Matthews Correlation Coefficient. The model ultimately proved itself to be promising for the classification of the domains in new sequences.
