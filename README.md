# RTP_Co-adapt_Model

Source Code File 1: RTP_Coadapt_ModelV1p2.m
Computational model of retinotopic mapping. MATLAB source code of our computational model (base version: Gebhardt et al., Development,139,335 (2012)) including co-adaptation.

Source Code File 2: wkeitpd.m
Function called by RTP_Coadapt_ModelV1p2.m while calculating probabilistic step decisions of simulated fiber terminals. wkeitpd.m calculates a Gaussian shaped distribution with standard deviation, sigma, describing the sojourn probability of a terminal experiencing a given local guidance potential, D.

Source Code File 3: wkeit01.m
Function called by RTP_Coadapt_ModelV1p2.m to calculate probabilistic step decisions of simulated fiber terminals. wkeit01.m calculates the probability, p, of changing the current position, based on the sojourn probabilities for the old and the new target position given by wkeitpd.m.

Source Code File 4: Substrate_Tectum.m
Target matrix called by RTP_Coadapt_ModelV1p2.m simulating the a-p graded distributions of EphAs and ephrin-As on the tectum. 

Source Code File 5: Substrate_Tectal_Innervation.m
Target matrix called by RTP_Coadapt_ModelV1p2.m simulating the a-p graded distributions of EphAs and ephrin-As on the tectum with an additional, cue-free stretch in front of the tectums’ anterior pole.

Source Code File 6: Substrate_Gap_Assay.m
Target matrix called by RTP_Coadapt_ModelV1p2.m simulating the various ephrin-A/EphA single and double-cue in-vitro gap substrates.


Source Code File 7: Gcanalyze_v2_0.m
MATLAB GUI source code of a tool developed to count axon fibers/bundles that grow (in parallel) from retinal explant strips for the quantification of in vitro gap-assays.
