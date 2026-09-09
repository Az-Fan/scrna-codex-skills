# Fixed figure design

- FeaturePlot: `#ECECEC` to `#B2182B`, independent per-gene scales, up to six genes per page.
- Dot plot: dot area is detected percentage; colour is scaled mean expression with `#2166AC`–white–`#B2182B`.
- Violin: condition colours are deterministic, with no inferential annotation.
- Sample expression: jittered biological-sample points plus a black descriptive median; pseudobulk values are preferred when available.
- Heatmap: gene-wise z scores across samples, diverging blue-white-red palette.
- DE effect: horizontal log2 fold change plot centred at zero; adjusted P value controls point size only.
- Volcano: all finite DE rows in grey and requested genes in orange-red with labels.

Do not add pie charts, significance stars, cell-level P values, rainbow gradients, or condition colours that change between runs.

