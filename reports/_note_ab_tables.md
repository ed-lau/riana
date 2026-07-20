# ¹⁸O note re-defaulting — A/B ladder tables

## Peptide level (curation yield, within-protein CV)

| run | arm | gate | n_fitted | n_admit | pct | med_R2_fitted | geomCV | n_prot_ge3 | med_k |
|---|---|---|---|---|---|---|---|---|---|
| boomi_ipsc_o18 | railoff | r2 | 12204 | 5704 | 46.7 | 0.778 | 0.142 | 615 | 0.0478 |
| boomi_ipsc_o18 | railoff | kcv | 12204 | 6249 | 51.2 | 0.778 | 0.152 | 659 | 0.0477 |
| boomi_ipsc_o18 | railon | r2 | 10944 | 6131 | 56.0 | 0.837 | 0.151 | 650 | 0.0471 |
| boomi_ipsc_o18 | railon | kcv | 10944 | 6671 | 61.0 | 0.837 | 0.168 | 691 | 0.0472 |
| boomi_ipsc_d2o | railoff | r2 | 15462 | 5708 | 36.9 | 0.702 | 0.14 | 576 | 0.0355 |
| boomi_ipsc_d2o | railoff | kcv | 15462 | 6598 | 42.7 | 0.702 | 0.154 | 655 | 0.0352 |
| boomi_ipsc_d2o | railon | r2 | 13719 | 6646 | 48.4 | 0.79 | 0.154 | 653 | 0.0343 |
| boomi_ipsc_d2o | railon | kcv | 13719 | 7395 | 53.9 | 0.79 | 0.168 | 715 | 0.0345 |
| juber_ac16_o18 | railoff | r2 | 5175 | 468 | 9.0 | 0.071 | 0.159 | 36 | 0.0305 |
| juber_ac16_o18 | railoff | kcv | 5175 | 514 | 9.9 | 0.071 | 0.172 | 40 | 0.0311 |
| juber_ac16_o18 | railon | r2 | 3009 | 594 | 19.7 | 0.314 | 0.16 | 51 | 0.0297 |
| juber_ac16_o18 | railon | kcv | 3009 | 654 | 21.7 | 0.314 | 0.172 | 56 | 0.0302 |
| juber_ac16_d2o | railoff | r2 | 11281 | 959 | 8.5 | 0.174 | 0.136 | 76 | 0.0335 |
| juber_ac16_d2o | railoff | kcv | 11281 | 1018 | 9.0 | 0.174 | 0.143 | 80 | 0.0338 |
| juber_ac16_d2o | railon | r2 | 6254 | 1190 | 19.0 | 0.384 | 0.163 | 96 | 0.0319 |
| juber_ac16_d2o | railon | kcv | 6254 | 1271 | 20.3 | 0.384 | 0.179 | 103 | 0.0324 |
| timeseries_lauren5_7_ipsc_mesoderm_o18 | railoff | r2 | 13101 | 1400 | 10.7 | 0.312 | 0.173 | 127 | 0.035 |
| timeseries_lauren5_7_ipsc_mesoderm_o18 | railoff | kcv | 13101 | 2138 | 16.3 | 0.312 | 0.202 | 211 | 0.0342 |
| timeseries_lauren5_7_ipsc_mesoderm_o18 | railon | r2 | 11677 | 2359 | 20.2 | 0.365 | 0.184 | 254 | 0.0322 |
| timeseries_lauren5_7_ipsc_mesoderm_o18 | railon | kcv | 11677 | 2958 | 25.3 | 0.365 | 0.201 | 330 | 0.0328 |
| timeseries_lauren9 | railoff | r2 | 28233 | 4328 | 15.3 | 0.399 | 0.139 | 486 | 0.0371 |
| timeseries_lauren9 | railoff | kcv | 28233 | 6063 | 21.5 | 0.399 | 0.176 | 685 | 0.0364 |
| timeseries_lauren9 | railon | r2 | 25372 | 6449 | 25.4 | 0.457 | 0.171 | 723 | 0.0349 |
| timeseries_lauren9 | railon | kcv | 25372 | 7756 | 30.6 | 0.457 | 0.187 | 855 | 0.0356 |

## Protein level (rollup)

| run | arm | gate | model | n_prot | med_k |
|---|---|---|---|---|---|
| boomi_ipsc_o18 | railoff | r2 | ols | 974 | 0.048 |
| boomi_ipsc_o18 | railoff | r2 | wls | 974 | 0.0488 |
| boomi_ipsc_o18 | railoff | kcv | ols | 1033 | 0.0478 |
| boomi_ipsc_o18 | railoff | kcv | wls | 1033 | 0.0486 |
| boomi_ipsc_o18 | railon | r2 | ols | 1024 | 0.0465 |
| boomi_ipsc_o18 | railon | r2 | wls | 1024 | 0.0475 |
| boomi_ipsc_o18 | railon | kcv | ols | 1090 | 0.0465 |
| boomi_ipsc_o18 | railon | kcv | wls | 1090 | 0.0476 |
| boomi_ipsc_d2o | railoff | r2 | ols | 874 | 0.0364 |
| boomi_ipsc_d2o | railoff | r2 | wls | 874 | 0.0359 |
| boomi_ipsc_d2o | railoff | kcv | ols | 951 | 0.0359 |
| boomi_ipsc_d2o | railoff | kcv | wls | 951 | 0.0356 |
| boomi_ipsc_d2o | railon | r2 | ols | 969 | 0.0352 |
| boomi_ipsc_d2o | railon | r2 | wls | 969 | 0.0349 |
| boomi_ipsc_d2o | railon | kcv | ols | 1036 | 0.0354 |
| boomi_ipsc_d2o | railon | kcv | wls | 1036 | 0.0351 |
| juber_ac16_o18 | railoff | r2 | ols | 72 | 0.0296 |
| juber_ac16_o18 | railoff | r2 | wls | 72 | 0.0294 |
| juber_ac16_o18 | railoff | kcv | ols | 85 | 0.0309 |
| juber_ac16_o18 | railoff | kcv | wls | 85 | 0.0308 |
| juber_ac16_o18 | railon | r2 | ols | 99 | 0.0295 |
| juber_ac16_o18 | railon | r2 | wls | 99 | 0.029 |
| juber_ac16_o18 | railon | kcv | ols | 110 | 0.03 |
| juber_ac16_o18 | railon | kcv | wls | 110 | 0.0295 |
| juber_ac16_d2o | railoff | r2 | ols | 145 | 0.0312 |
| juber_ac16_d2o | railoff | r2 | wls | 145 | 0.0307 |
| juber_ac16_d2o | railoff | kcv | ols | 156 | 0.0316 |
| juber_ac16_d2o | railoff | kcv | wls | 156 | 0.0312 |
| juber_ac16_d2o | railon | r2 | ols | 189 | 0.0315 |
| juber_ac16_d2o | railon | r2 | wls | 189 | 0.0313 |
| juber_ac16_d2o | railon | kcv | ols | 203 | 0.0318 |
| juber_ac16_d2o | railon | kcv | wls | 203 | 0.0318 |
| timeseries_lauren5_7_ipsc_mesoderm_o18 | railoff | r2 | ols | 441 | 0.0342 |
| timeseries_lauren5_7_ipsc_mesoderm_o18 | railoff | r2 | wls | 441 | 0.0339 |
| timeseries_lauren5_7_ipsc_mesoderm_o18 | railoff | kcv | ols | 632 | 0.0326 |
| timeseries_lauren5_7_ipsc_mesoderm_o18 | railoff | kcv | wls | 632 | 0.0327 |
| timeseries_lauren5_7_ipsc_mesoderm_o18 | railon | r2 | ols | 711 | 0.0315 |
| timeseries_lauren5_7_ipsc_mesoderm_o18 | railon | r2 | wls | 711 | 0.0313 |
| timeseries_lauren5_7_ipsc_mesoderm_o18 | railon | kcv | ols | 820 | 0.0318 |
| timeseries_lauren5_7_ipsc_mesoderm_o18 | railon | kcv | wls | 820 | 0.0318 |
| timeseries_lauren9 | railoff | r2 | ols | 825 | 0.0372 |
| timeseries_lauren9 | railoff | r2 | wls | 825 | 0.0371 |
| timeseries_lauren9 | railoff | kcv | ols | 1051 | 0.0366 |
| timeseries_lauren9 | railoff | kcv | wls | 1051 | 0.0362 |
| timeseries_lauren9 | railon | r2 | ols | 1163 | 0.0353 |
| timeseries_lauren9 | railon | r2 | wls | 1163 | 0.0347 |
| timeseries_lauren9 | railon | kcv | ols | 1305 | 0.0356 |
| timeseries_lauren9 | railon | kcv | wls | 1305 | 0.0352 |

## Cross-label head-to-head (§5.2)

| pair | arm | gate | model | n_pep | rho_pep | n_prot | rho_prot | med_log2_d2o_over_o18 |
|---|---|---|---|---|---|---|---|---|
| AICS52 (boomi) | railoff | r2 | ols | 3002 | 0.679 | 724 | 0.74 | -0.464 |
| AICS52 (boomi) | railoff | r2 | wls | 3002 | 0.679 | 724 | 0.711 | -0.464 |
| AICS52 (boomi) | railoff | kcv | ols | 3558 | 0.667 | 795 | 0.724 | -0.466 |
| AICS52 (boomi) | railoff | kcv | wls | 3558 | 0.667 | 795 | 0.691 | -0.466 |
| AICS52 (boomi) | railon | r2 | ols | 3442 | 0.677 | 789 | 0.718 | -0.458 |
| AICS52 (boomi) | railon | r2 | wls | 3442 | 0.677 | 789 | 0.691 | -0.458 |
| AICS52 (boomi) | railon | kcv | ols | 3950 | 0.662 | 851 | 0.71 | -0.458 |
| AICS52 (boomi) | railon | kcv | wls | 3950 | 0.662 | 851 | 0.685 | -0.458 |
| AC16 (juber) | railoff | r2 | ols | 260 | 0.673 | 62 | 0.667 | 0.169 |
| AC16 (juber) | railoff | r2 | wls | 260 | 0.673 | 62 | 0.654 | 0.169 |
| AC16 (juber) | railoff | kcv | ols | 295 | 0.664 | 70 | 0.697 | 0.161 |
| AC16 (juber) | railoff | kcv | wls | 295 | 0.664 | 70 | 0.692 | 0.161 |
| AC16 (juber) | railon | r2 | ols | 330 | 0.641 | 79 | 0.62 | 0.189 |
| AC16 (juber) | railon | r2 | wls | 330 | 0.641 | 79 | 0.609 | 0.189 |
| AC16 (juber) | railon | kcv | ols | 373 | 0.649 | 89 | 0.636 | 0.167 |
| AC16 (juber) | railon | kcv | wls | 373 | 0.649 | 89 | 0.625 | 0.167 |
| SCVI480 (lauren) | railoff | r2 | ols | 749 | 0.668 | 371 | 0.678 | 0.175 |
| SCVI480 (lauren) | railoff | r2 | wls | 749 | 0.668 | 371 | 0.673 | 0.175 |
| SCVI480 (lauren) | railoff | kcv | ols | 1264 | 0.649 | 535 | 0.676 | 0.181 |
| SCVI480 (lauren) | railoff | kcv | wls | 1264 | 0.649 | 535 | 0.644 | 0.181 |
| SCVI480 (lauren) | railon | r2 | ols | 1206 | 0.657 | 591 | 0.642 | 0.206 |
| SCVI480 (lauren) | railon | r2 | wls | 1206 | 0.657 | 591 | 0.618 | 0.206 |
| SCVI480 (lauren) | railon | kcv | ols | 1708 | 0.646 | 690 | 0.631 | 0.199 |
| SCVI480 (lauren) | railon | kcv | wls | 1708 | 0.646 | 690 | 0.61 | 0.199 |

## §6 mesoderm Δk — significance counts

| arm | gate | model | n_tested | n_sig | pct_sig | faster | slower | med_abs_dk |
|---|---|---|---|---|---|---|---|---|
| railoff | legacy | ols | 869 | 519 | 59.7 | 253 | 266 | 0.0116 |
| railon | legacy | ols | 750 | 419 | 55.9 | 208 | 211 | 0.0112 |
| railoff | r2 | ols | 378 | 241 | 63.8 | 125 | 116 | 0.0106 |
| railoff | r2 | wls | 378 | 178 | 47.1 | 65 | 113 | 0.0108 |
| railoff | kcv | ols | 552 | 349 | 63.2 | 182 | 167 | 0.0102 |
| railoff | kcv | wls | 552 | 230 | 41.7 | 87 | 143 | 0.0111 |
| railon | r2 | ols | 561 | 333 | 59.4 | 173 | 160 | 0.0099 |
| railon | r2 | wls | 561 | 213 | 38.0 | 80 | 133 | 0.0107 |
| railon | kcv | ols | 652 | 369 | 56.6 | 196 | 173 | 0.0103 |
| railon | kcv | wls | 652 | 225 | 34.5 | 84 | 141 | 0.0111 |
