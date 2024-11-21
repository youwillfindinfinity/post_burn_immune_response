# An in silico modeling approach to understanding the dynamics of the post-burn immune response
Paper published [here](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1303776/full)

**Authors:**  
H. Ibrahim Korkmaz<sup>1,2,3,4*</sup>†, Vivek M. Sheraton<sup>5,6,7*</sup>†, Roland V. Bumbuc<sup>1,2,5,6,7</sup>, Meifang Li<sup>5</sup>, Anouk Pijpe<sup>1,3</sup>, Patrick P.G. Mulder<sup>4,8</sup>, Bouke K.H.L. Boekema<sup>1,4</sup>, Evelien de Jong<sup>9</sup>, Stephan G.F. Papendorp<sup>9</sup>, Ruud Brands<sup>10,11</sup>, Esther Middelkoop<sup>1,3,4</sup>, Peter M.A. Sloot<sup>5</sup>, Paul P.M. van Zuijlen<sup>1,3,4,12</sup>

† These authors contributed equally to this work and share first authorship.

**Affiliations:**
1. Department of Plastic, Reconstructive and Hand Surgery, Amsterdam Movement Sciences (AMS) Institute, Amsterdam UMC, Location VUmc, Amsterdam, The Netherlands.
2. Department of Molecular Cell Biology and Immunology, Amsterdam Infection and Immunity (AII) Institute, Amsterdam UMC, Location VUmc, Amsterdam, The Netherlands.
3. Burn Center and Department of Plastic and Reconstructive Surgery, Red Cross Hospital, Beverwijk, The Netherlands.
4. Preclinical Research, Association of Dutch Burn Centres (ADBC), Beverwijk, The Netherlands.
5. Computational Science Lab, Informatics Institute, University of Amsterdam, UvA - LAB42, Amsterdam, The Netherlands.
6. Center for Experimental and Molecular Medicine (CEMM), Amsterdam UMC, Amsterdam, The Netherlands.
7. Laboratory for Experimental Oncology and Radiobiology, ONCODE, Amsterdam UMC, Location AMC, Amsterdam, The Netherlands.
8. Laboratory of Medical Immunology, Department of Laboratory Medicine, Radboud University Medical Center, Nijmegen, The Netherlands.
9. Department of Intensive Care, Red Cross Hospital, Beverwijk, The Netherlands.
10. Complexity Institute, Nanyang Technological University, Singapore, Singapore.
11. Alloksys Life Sciences BV, Wageningen, Netherlands.
12. Paediatric Surgical Centre, Emma Children’s Hospital, Amsterdam UMC, location AMC, Amsterdam, The Netherlands.

**Correspondence:**
H. Ibrahim Korkmaz  
Email: h.korkmaz@amsterdamumc.nl

**Keywords:**  
burns, wound healing, inflammation, immune response, computational modeling.

![Image](https://github.com/youwillfindinfinity/post_burn_immune_response/blob/main/Paper_figures/Figure%201/flow%20chart%20of%20conceptual%20model.png)

This paper was submitted, acceptesd, peer-reviewed and publushed in *Frontiers in Immunology: New Approach Methods in Immunology*.

![Image](https://github.com/youwillfindinfinity/post_burn_immune_response/blob/8f2af0e6082c1acdeabf84115769948ae6816db4/endothelial_experiment_fixed_modulus/E4/2D%20data/iter_4_il8_700kmcs.png)

## Abstract

Introduction: Burns are characterized by a massive and prolonged acute inflammation, which persists for up to months after the initial trauma. Due to the complexity of the inflammatory process, Predicting the dynamics of wound healing process can be challenging for burn injuries. The aim of this study was to develop simulation models for the post-burn immune response based on (pre)clinical data.

Methods: The simulation domain was separated into blood and tissue compartments. Each of these compartments contained solutes and cell agents. Solutes comprise pro-inflammatory cytokines, anti-inflammatory cytokines and inflammation triggering factors. The solutes diffuse around the domain based on their concentration profiles. The cells include mast cells, neutrophils, and macrophages, and were modeled as independent agents. The cells are motile and exhibit chemotaxis based on concentrations gradients of the solutes. In addition, the cells secrete various solutes that in turn alter the dynamics and responses of the burn wound system.

Results: We developed an Glazier-Graner-Hogeweg method-based model (GGH) to capture the complexities associated with the dynamics of inflammation after burn injuries, including changes in cell counts and cytokine levels. Through simulations from day 0 – 4 post-burn, we successfully identified key factors influencing the acute inflammatory response, i.e., the initial number of endothelial cells, the chemotaxis threshold, and the level of chemoattractants.

Conclusion: Our findings highlight the pivotal role of the initial endothelial cell count as a key parameter for intensity of inflammation and progression of acute inflammation, 0 – 4 days post-burn.

![Image](https://github.com/youwillfindinfinity/post_burn_immune_response/blob/main/Simulation%20results/FIgures%20and%20comparisons/cell_count_data_E1.png)


## Software Implementation

All [source code](https://github.com/youwillfindinfinity/post_burn_immune_response/tree/8724f978fbc22f7aa8acb01098d3450ec607c277/Code) used to generate the [results](https://github.com/youwillfindinfinity/post_burn_immune_response/tree/65fa78e7835fa46f80284881e243a8567832f673/Simulation%20results), [raw figures](https://github.com/youwillfindinfinity/post_burn_immune_response/tree/8724f978fbc22f7aa8acb01098d3450ec607c277/Simulation%20results/FIgures%20and%20comparisons) and [some data(depending on size)](https://github.com/youwillfindinfinity/post_burn_immune_response/tree/8724f978fbc22f7aa8acb01098d3450ec607c277/Simulation%20results) in the paper can be accessed. The calculations and figure generation are all ran in Python 3.10 and all simulations are carried jointly through CompuCell3D. Figures used in the paper are saved in [figures](https://github.com/youwillfindinfinity/post_burn_immune_response/tree/8724f978fbc22f7aa8acb01098d3450ec607c277/Paper_figures).

## Getting the Figures and Videos of the Simulations

You can download a copy of all the files in this repository by cloning the repository:

git clone https://github.com/youwillfindinfinity/post_burn_immune_response


or by clicking on `Code` > `Download ZIP`.

The authors reserve the rights to the article content, which is currently being submitted for publication in the *Frontiers in Immunology: New Approach Methods in Immunology* special issue journal.

Any use of the code, figures or data, must be cited as follows:

@article{korkmaz2024silico,
  title={An in silico modeling approach to understanding the dynamics of the post-burn immune response},
  author={Korkmaz, H Ibrahim and Sheraton, Vivek M and Bumbuc, Roland V and Li, Meifang and Pijpe, Anouk and Mulder, Patrick PG and Boekema, Bouke KHL and de Jong, Evelien and Papendorp, Stephan GF and Brands, Ruud and others},
  journal={Frontiers in immunology},
  volume={15},
  pages={1303776},
  year={2024},
  publisher={Frontiers Media SA}
}
