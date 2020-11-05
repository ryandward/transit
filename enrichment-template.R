require('pacman')
pacman::p_load(data.table, R.utils, scales, ggplot2, pracma, ggrepel, devtools, ggbiplot)

#Differentially essential
differentials <- fread('/home/ryandward/Class/Genetics885/transit/Escherichia_coli_BW25113_Acinetobacter_baumannii_ATCC_17978_diff_ess.tsv_GO_all_ecocyc_functional_classification.tsv')
differentials[, c("GO", "Process") := tstrsplit(`Process~name`, "~")]
differentials <- differentials[,.(GO, Process, num_of_Genes, gene_group, `Benjamini and Hochberg (FDR)`)]
differentials[, significance := ifelse(`Benjamini and Hochberg (FDR)` < 0.01, "FDR<=0.01", "Not significant")]
setorder(differentials,`Benjamini and Hochberg (FDR)`)

this_plot <- ggplot(differentials, aes(x = log10(num_of_Genes), y = -log10(`Benjamini and Hochberg (FDR)`))) +
	geom_point(aes(color = significance)) +
	scale_color_manual(values = c("red", "grey")) +
	theme_bw(base_size = 12) + theme(legend.position = "bottom") +
	geom_text_repel(
		data = subset(differentials, `Benjamini and Hochberg (FDR)` < 0.01),
		aes(label = Process),
		size = 4,
		box.padding = unit(0.35, "lines"),
		point.padding = unit(0.3, "lines"),
		max.iter = 1000
	) +
	ggtitle("GO Pathway Enrichment of Differentially Essential Genes in\nAcinetobacter baumannii 17978 vs. Escherichia coli BW25113") +
	theme(plot.title = element_text(hjust = 0.5))

plot(this_plot)

#Differentially non-essential
differentials <- fread('/home/ryandward/Class/Genetics885/transit/Escherichia_coli_BW25113_Acinetobacter_baumannii_ATCC_17978_diff_noness.tsv_GO_all_ecocyc_functional_classification.tsv')
differentials[, c("GO", "Process") := tstrsplit(`Process~name`, "~")]
differentials <- differentials[,.(GO, Process, num_of_Genes, gene_group, `Benjamini and Hochberg (FDR)`)]
differentials[, significance := ifelse(`Benjamini and Hochberg (FDR)` < 0.01, "FDR<=0.01", "Not significant")]
setorder(differentials,`Benjamini and Hochberg (FDR)`)

this_plot <- ggplot(differentials, aes(x = log10(num_of_Genes), y = -log10(`Benjamini and Hochberg (FDR)`))) +
	geom_point(aes(color = significance)) +
	scale_color_manual(values = c("red", "grey")) +
	theme_bw(base_size = 12) + theme(legend.position = "bottom") +
	geom_text_repel(
		data = subset(differentials, `Benjamini and Hochberg (FDR)` < 0.01),
		aes(label = Process),
		size = 4,
		box.padding = unit(0.35, "lines"),
		point.padding = unit(0.3, "lines"),
		max.iter = 1000
	) +
	ggtitle("GO Pathway Enrichment of Differentially Essential Genes in\nAcinetobacter baumannii 17978 vs. Escherichia coli BW25113") +
	theme(plot.title = element_text(hjust = 0.5))

plot(this_plot)

#Both essential
differentials <- fread('/home/ryandward/Class/Genetics885/transit/Escherichia_coli_BW25113_Acinetobacter_baumannii_ATCC_17978_both_ess.tsv_GO_all_ecocyc_functional_classification.tsv')
differentials[, c("GO", "Process") := tstrsplit(`Process~name`, "~")]
differentials <- differentials[,.(GO, Process, num_of_Genes, gene_group, `Benjamini and Hochberg (FDR)`)]
differentials[, significance := ifelse(`Benjamini and Hochberg (FDR)` < 0.001, "FDR<=0.001", "Not significant")]
setorder(differentials,`Benjamini and Hochberg (FDR)`)

this_plot <- ggplot(differentials, aes(x = log10(num_of_Genes), y = -log10(`Benjamini and Hochberg (FDR)`))) +
	geom_point(aes(color = significance)) +
	scale_color_manual(values = c("red", "grey")) +
	theme_bw(base_size = 12) + theme(legend.position = "bottom") +
	geom_text_repel(
		data = subset(differentials, `Benjamini and Hochberg (FDR)` < 0.001),
		aes(label = Process),
		size = 4,
		box.padding = unit(0.35, "lines"),
		point.padding = unit(0.3, "lines"),
		max.iter = 1000
	) +
	ggtitle("GO Pathway Enrichment of Differentially Essential Genes in\nAcinetobacter baumannii 17978 vs. Escherichia coli BW25113") +
	theme(plot.title = element_text(hjust = 0.5))

plot(this_plot)






