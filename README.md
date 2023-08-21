# Tryptic-Cleavage

By Ahmad Arnaout
Supervisor: Professor Conrad Bessant

This research dissertation is submitted for the MSc in Bioinformatics at Queen Mary, University of London
MSc Bioinformatics| Ahmad Arnaout


## Abstract ##

The most widely used approach in quantitative proteomics research involves utilizing proteolytic peptides as substitutes for parent proteins. These peptides play a pivotal role by serving as indicators for both protein identity and quantity. The accuracy of this technique hinges on the breakdown of the parent protein into peptides, preserving their accurate representation. This breakdown is predominantly achieved through the use of trypsin, an enzyme known for its exceptional specificity in hydrolysing peptide bonds exclusively at the carboxyl-termini of Lysine and Arginine residues. This high degree of specificity renders trypsin the preferred enzyme for the bottom-up proteomic approach. Nevertheless, a challenge arises during the process of proteolysis, as proteins occasionally undergo incomplete digestion, resulting in the production of what are termed missed-cleaved peptides. In response, developing computational methods to predict missed cleavages is challenging, however, in this research project machine learning models including Random Forest and XGBoost were developed to predict missed cleavage sites. Existing computational methodology only predicts missed cleavages using amino acid sequences and overlooks 3D structures of the protein. This study emphasizes structural aspects and their properties to predict missed cleavages, by including features such as secondary motifs and solvent accessibility to train machine learning models for improved prediction accuracy. As a result, the research project yielded a random forest model achieving 0.81 accuracy in predicting missed cleavage, while the XGBoost model showed poorer performance. In addition, structural features including relative ASA, residue depth, phi, and psi were discovered to be the more predictive features that influence the likelihood of a successful tryptic cleavage. 


