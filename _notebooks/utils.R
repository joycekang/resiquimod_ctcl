suppressPackageStartupMessages({
    suppressWarnings({
    library(Seurat)
    library(ggplot2)
    library(ggrepel)
    library(data.table)
    library(tidyverse)
    library(Matrix)
    library(harmony)
    library(presto)
    library(viridis)
    library(MAST)
    library(SingleCellExperiment)
    library(singlecellmethods)
    library(lme4)
    library(ggpubr)
    library(ggrastr)
    library(cowplot)
    library(limma)
    library(mixtools)
    library(rstatix)
    library(RColorBrewer)
    library(patchwork)
    library(dplyr)
    library(scales)
    library(rlang)

    #GSEA
    library(clusterProfiler)
    library(msigdbr)
    library(org.Hs.eg.db)})
})

# Helper functions
fig.size = function (height, width) {
    options(repr.plot.height = height, repr.plot.width = width, repr.plot.res = 200)
}

lower_str = function(original_string) {
    gsub("[^A-Za-z0-9]", "", tolower(original_string))
}

# Color palettes
custom.colors = c(  'Keratinocyte undiff.'='#3256a8',
                    'Keratinocyte diff.'='#3b9edb',
                    'Fibroblast'='#6533ed',
                    'T cell'='#e39309',
                    'Vascular EC'='#e64a3c',
                    'Myeloid'='#c75656',
                    'Pericyte'='#9698dc',
                    'Hair follicle'='#6da15f',
                    'Melanocyte'='#f283e3',
                    'Lymphatic EC'='#6ef0c9',
                    'Mast'='#dbb81a',
                    'Dermal adipocyte' = '#ffb2be')

color_scheme <- c("Non-lesional" = "#66c2a5", "Treated\nlesional" = "#3288bd", "Untreated\nlesional" = "#d53e4f")

mycolors = colorRampPalette(brewer.pal(8, "Paired"))(14)

# Function to find valley between two normal distributions
find_valley <- function(mu1, mu2, sd1, sd2, lambda1, lambda2, range_vals) {
  # Define mixture density function
  mix_density <- function(x) {
    lambda1 * dnorm(x, mu1, sd1) + lambda2 * dnorm(x, mu2, sd2)
  }
  
  # Find minimum point between the means
  valley_range <- seq(min(mu1, mu2), max(mu1, mu2), length.out = 100)
  densities <- sapply(valley_range, mix_density)
  valley <- valley_range[which.min(densities)]
  return(valley)
}

plot_genes_hex <- function(
    seurat_obj,
    features = NULL,       # List of features (genes or module score names)
    feature_type = "gene", # Either "gene" or "module"
    save_path = NULL,      # Path to save PDF
    height = 3.5,
    width = 10,
    hex_bins = 50,
    color_palette = c("mediumpurple4", "maroon", "orange", "lightgoldenrod1"),
    quantile_limit = 0.99,
    ncol = NULL            # Parameter for controlling layout
) {
    # Get UMAP coordinates
    umap_coords <- seurat_obj@meta.data[, c('merged_UMAP_1', 'merged_UMAP_2')]
    colnames(umap_coords) <- c("UMAP1", "UMAP2")
    
    # Get expression data based on feature type
    if (feature_type == "gene") {
        # Get data from assay slot
        exp_data <- GetAssayData(seurat_obj, layer = "data")[features, ]
        exp_data <- as.data.frame(t(as.matrix(exp_data)))
    } else if (feature_type == "module") {
        # Get data from meta.data slot
        exp_data <- seurat_obj@meta.data[, features, drop = FALSE]
    } else {
        stop("feature_type must be either 'gene' or 'module'")
    }
    
    # Combine expression data with UMAP coordinates
    plot_data <- cbind(umap_coords, exp_data)
    
    # Convert to long format for plotting
    long_data <- gather(plot_data, 
                       key = "feature", 
                       value = "expression", 
                       -UMAP1, -UMAP2)
    
    # Determine number of columns if not specified
    if (is.null(ncol)) {
        n_plots <- length(features)
        ncol <- if (n_plots == 2) 2 else 1
    }
    
    # Create plot with hexagon-based normalization
    p <- long_data %>%
        group_split(feature) %>%
        map(function(df) {
            # First create the hex bins to get their statistics
            hex_plot <- ggplot(df, aes(x = UMAP1, y = UMAP2, z = expression)) +
                stat_summary_hex(bins = hex_bins)
            
            # Extract the computed statistics from the hex layer
            hex_data <- layer_data(hex_plot, 1)
            
            # For module scores, make the limits symmetric around zero if there are negative values
            if (feature_type == "module" && min(hex_data$value) < 0) {
                max_value <- max(hex_data$value)
                min_value = min(hex_data$value)
                value_limits <- c(min_value, max_value)
                values <- seq(0, 1, length.out = length(color_palette))
            } else {
                value_limits <- c(0, max(hex_data$value))
                values <- seq(0, 1, length.out = length(color_palette))
            }
            
            # Create the final plot with limits based on hex summary statistics
            ggplot(df, aes(x = UMAP1, y = UMAP2, z = expression)) +
                stat_summary_hex(bins = hex_bins) +
                scale_fill_gradientn(
                    colours = color_palette,
                    values = values,
                    limits = value_limits,
                    na.value = color_palette[length(color_palette)]
                ) +
                facet_wrap(~feature) +
                theme_void() +
                labs(fill = ifelse(feature_type == "gene", 
                                 'Mean normalized\nexpression',
                                 'Module score')) +
                theme(
                    strip.background = element_blank(),
                    strip.text = element_text(size = 12, face = "bold"),
                    legend.text = element_text(size = 10),
                    legend.title = element_text(size = 10)
                )
        }) %>%
        plot_grid(plotlist = ., align = 'hv', ncol = ncol)
    
    # Save if path provided
    if (!is.null(save_path)) {
        pdf(save_path, height = height, width = width, useDingbats = FALSE)
        print(p)
        dev.off()
    }
    return(p)
}

plot_split_hex <- function(
    seurat_obj,
    features = NULL,       # List of features (genes or module score names)
    feature_type = "gene", # Either "gene" or "module"
    split_by = NULL,       # Metadata column to split by
    facet_by = NULL,       # New parameter for faceting
    facet_ncol = NULL,     # New parameter for number of facet columns
    save_path = NULL,      # Path to save PDF
    height = 3.5,
    width = 10,
    hex_bins = 50,
    color_palette = c("mediumpurple4", "maroon", "orange", "lightgoldenrod1"),
    quantile_limit = 0.99,
    plot_ncol = 2
) {
    # Get UMAP coordinates
    umap_coords <- seurat_obj@meta.data[, c('merged_UMAP_1', 'merged_UMAP_2')]
    colnames(umap_coords) <- c("UMAP1", "UMAP2")
    
    # Get expression data based on feature type
    if (feature_type == "gene") {
        exp_data <- GetAssayData(seurat_obj, layer = "data")[features, ]
        exp_data <- as.data.frame(t(as.matrix(exp_data)))
    } else if (feature_type == "module") {
        exp_data <- seurat_obj@meta.data[, features, drop = FALSE]
    } else {
        stop("feature_type must be either 'gene' or 'module'")
    }
    
    # Get split variable
    if (is.null(split_by)) {
        stop("split_by parameter must be specified")
    }
    split_data <- seurat_obj@meta.data[[split_by]]
    
    # Get facet variable if specified
    if (!is.null(facet_by)) {
        facet_data <- seurat_obj@meta.data[[facet_by]]
    }
    
    # Combine all data
    plot_data <- cbind(umap_coords, exp_data, split = split_data)
    if (!is.null(facet_by)) {
        plot_data$facet <- facet_data
    }
    
    # Convert to long format for plotting
    if (!is.null(facet_by)) {
        long_data <- gather(plot_data, 
                          key = "feature", 
                          value = "expression", 
                          -UMAP1, -UMAP2, -split, -facet)
    } else {
        long_data <- gather(plot_data, 
                          key = "feature", 
                          value = "expression", 
                          -UMAP1, -UMAP2, -split)
    }
    
    # First pass to get consistent color scales per feature
    feature_limits <- long_data %>%
        group_split(feature) %>%
        map(function(df) {
            # Create hex bins for all split categories for this feature
            all_hexes <- df %>%
                group_split(split) %>%
                map(function(split_df) {
                    hex_plot <- ggplot(split_df, aes(x = UMAP1, y = UMAP2, z = expression)) +
                        stat_summary_hex(bins = hex_bins)
                    hex_data <- layer_data(hex_plot, 1)
                    return(hex_data$value)
                }) %>%
                unlist()
            
            if (feature_type == "module" && min(all_hexes) < 0) {
                return(c(min(all_hexes), max(all_hexes)))
            } else {
                return(c(0, max(all_hexes)))
            }
        })
    names(feature_limits) <- unique(long_data$feature)
    
    # Create plots
    plot_list <- long_data %>%
        group_split(feature) %>%
        map(function(feature_df) {
            current_feature <- unique(feature_df$feature)
            value_limits <- feature_limits[[current_feature]]
            
            # Create base plot
            base_plot <- ggplot(feature_df, aes(x = UMAP1, y = UMAP2, z = expression)) +
                stat_summary_hex(bins = hex_bins) +
                scale_fill_gradientn(
                    colours = color_palette,
                    values = seq(0, 1, length.out = length(color_palette)),
                    limits = value_limits,
                    na.value = color_palette[length(color_palette)]
                ) +
                theme_void(base_size = 11) +
                labs(fill = ifelse(feature_type == "gene", 
                                 'Mean normalized\nexpression',
                                 'Module score')) +
                theme(
                    plot.title = element_text(hjust = 0.5, size = 12),
                    legend.text = element_text(size = 10),
                    legend.title = element_text(size = 10),
                    strip.text = element_text(size = 12)
                )
            
            # Add faceting if specified
            if (!is.null(facet_by)) {
                base_plot <- base_plot + 
                    facet_wrap(~facet + split, ncol = facet_ncol)
            } else {
                base_plot <- base_plot +
                    facet_wrap(~split, ncol = facet_ncol)
            }
            
            # Add feature label
            title <- ggdraw() + 
                draw_label(current_feature, 
                          fontface = "bold",
                          x = 0,
                          hjust = 0)
            
            plot_grid(title, base_plot,
                     ncol = 1,
                     rel_heights = c(0.1, 1))
        })
    
    # Combine all features into final plot
    p <- plot_grid(plotlist = plot_list,
                  ncol = plot_ncol,
                  align = "v")
    
    # Save if path provided
    if (!is.null(save_path)) {
        pdf(save_path, height = height, width = width, useDingbats = FALSE)
        print(p)
        dev.off()
    }
    
    return(p)
}

plot_gene_violins <- function(
    seurat_obj,
    genes,              # List/vector of genes to plot
    split_by,           # Metadata column to split violins by
    group_by,
    ncol = 2,           # Number of columns in plot grid
    pt_size = 0.1       # Size of individual points
) {
    # Verify genes exist in the dataset
    missing_genes <- genes[!genes %in% rownames(seurat_obj)]
    if(length(missing_genes) > 0) {
        warning(sprintf("Following genes not found: %s", 
                paste(missing_genes, collapse = ", ")))
        genes <- genes[genes %in% rownames(seurat_obj)]
    }
    
    if(length(genes) == 0) {
        stop("No valid genes provided")
    }
    
    # Verify split_by column exists
    if(!split_by %in% colnames(seurat_obj@meta.data)) {
        stop(sprintf("Column '%s' not found in metadata", split_by))
    }
    
    # Calculate number of rows needed
    nrows <- ceiling(length(genes)/ncol)
    
    # Create plots for each gene
    plots <- lapply(genes, function(gene) {
        # Extract data for plotting
        plot_data <- data.frame(
            expression = GetAssayData(seurat_obj, layer = "data")[gene,],
            split_var = seurat_obj@meta.data[[split_by]],
            group_var = seurat_obj@meta.data[[group_by]]
        )
        
        # Create violin plot
        p <- ggplot(plot_data, 
                   aes(x = group_var, y = expression, 
                       fill = split_var)) +
            geom_violin(alpha = 1, 
                       scale = "width",  linewidth = 0.4,
                       position = position_dodge(width = 1)) +
            #geom_jitter_rast(aes(),
            #           color = 'black', 
            #           size = pt_size,
            #           alpha = 0.25,
            #           position = position_jitterdodge(dodge.width = 1)) +
            theme_bw(base_size = 11) + scale_fill_manual(values = c("Wk0" = "#e5f59f", "Wk8" = "#4b95db", "Wk24" = "#de5e47" )) +
            theme(
                axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
                axis.title.x = element_blank(),
                plot.title = element_text(face = "bold", hjust = 0.5),
                legend.position = "bottom",
                legend.title = element_blank()) +
            labs(title = gene, y = "Expression")
        
        return(p)
    })
    
    # Combine plots
    combined_plot <- wrap_plots(plots, ncol = ncol)
    return(combined_plot)
}

cell_type_levels = c('Keratinocyte undiff.', 'Keratinocyte diff.', 'Melanocyte', 'Hair follicle',
                     'Dermal adipocyte', 'Fibroblast', 'Pericyte', 'Vascular EC', 'Lymphatic EC',
                     'Mast', 'Myeloid', 'T cell')

# Load snRNA-seq object with updated cell type labels from the merged dataset.
# Returns a Seurat object with harmonized cell_type labels and UMAP coordinates.
load_snRNA_obj <- function(
    rds_path = '../_data/snRNAseq/obj_43854cells_labeled_v2.rds',
    labels_path = '../_data/obj_merged_snRNA_labels.csv'
) {
    obj <- readRDS(rds_path)
    labels <- read.csv(labels_path, row.names = 1)
    obj$cell_type <- NULL
    stopifnot(all(obj$Cell == labels$Cell))
    rownames(labels) <- labels$Cell
    labels$Cell <- NULL
    obj <- AddMetaData(obj, labels)
    obj@meta.data$cell_type <- factor(obj@meta.data$cell_type, levels = cell_type_levels)
    return(obj)
}

# Load FFPE scRNA-seq object with updated cell type labels from the merged dataset.
load_ffpe_obj <- function(
    rds_path = '../_data/obj_27448cells_labeled_v2.rds',
    labels_path = '../_data/obj_merged_scRNA_labels.csv'
) {
    obj <- readRDS(rds_path)
    labels <- read.csv(labels_path, row.names = 1)
    obj$cell_type <- NULL
    stopifnot(all(obj$Cell == labels$Cell))
    rownames(labels) <- labels$Cell
    labels$Cell <- NULL
    obj <- AddMetaData(obj, labels)
    obj@meta.data$cell_type <- factor(obj@meta.data$cell_type, levels = cell_type_levels)
    return(obj)
}

# Load MSigDB gene sets (hallmark + KEGG canonical + GO:BP).
# Returns a list with $hallmark and $all (combined) data frames.
load_gene_sets <- function() {
    hallmark  <- msigdbr(species = "Homo sapiens", category = "H")
    canonical <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
    go        <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")

    all_gene_sets <- rbind(
        hallmark  %>% mutate(collection = "Hallmark"),
        canonical %>% mutate(collection = "Canonical-Kegg"),
        go        %>% mutate(collection = "GO-BP")
    )
    hallmark_gene_sets <- hallmark %>% mutate(collection = "Hallmark")

    return(list(hallmark = hallmark_gene_sets, all = all_gene_sets))
}

# Volcano plot for one cell type's DE results.
# markers_ct: data frame with avg_log2FC, p_val_adj, gene columns
# tolabel:    subset of markers_ct rows to highlight and label
# title:      plot title (typically the cell type name)
# nx, ny:     nudge for label repulsion
volcano_plot <- function(markers_ct, tolabel, title = '', nx = 1.5, ny = 5) {
    ggplot(markers_ct, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
        geom_point(size = 0.2) +
        geom_point(data = tolabel, size = 0.2, col = 'red') +
        geom_text_repel(
            data = tolabel,
            aes(label = gene),
            size = 2,
            segment.size = 0.1, force = 6, min.segment.length = 0.1, max.overlaps = 20,
            nudge_x = nx * sign(tolabel$avg_log2FC),
            nudge_y = ny, segment.color = 'grey48', seed = 1
        ) +
        geom_vline(xintercept = 0, col = 'grey42', linewidth = 0.2) +
        theme_bw(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5)) +
        labs(title = title, x = 'Average log2FC', y = '-log10(Padj)')
}

# Run GSEA for all cell types in a named list of DE results and save plots.
# markers_list: named list of DE result data frames (one per cell type)
# out_dir:      output directory (e.g. "../_results/GSEA_Tx/hallmark_gs")
# label:        filename prefix (e.g. "allW8vsW0")
# gene_set:     msigdbr gene set data frame (use load_gene_sets()$hallmark or $all)
# plot_width:   width of saved PDF in inches
run_GSEA_on_list <- function(markers_list, out_dir, label = '', gene_set, plot_width = 9) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    results_list <- list()
    for (cell_type in names(markers_list)) {
        message(paste("Processing", cell_type))
        results_list[[cell_type]] <- run_gsea_analysis(markers_list[[cell_type]], cell_type, gene_set)
        ggsave(
            filename = file.path(out_dir, paste0(label, "_gsea_plot_", cell_type, ".pdf")),
            plot = results_list[[cell_type]]$plot,
            width = plot_width,
            height = 5
        )
    }
    invisible(results_list)
}

# Save a named list of DE result data frames to CSVs for IPA pathway analysis.
# out_dir: directory path (created if needed), one CSV per cell type.
save_to_csvs <- function(reslist, out_dir) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    for (name in names(reslist)) {
        DEres <- reslist[[name]]
        DEres$cell_type <- name
        DEres$gene <- rownames(DEres)
        write.csv(DEres, file.path(out_dir, paste0(name, '.csv')), quote = FALSE)
    }
}

# Create a function to run GSEA and create plot for one cell type
run_gsea_analysis <- function(de_results, cell_type, gene_set) {
    
    # Extract and prepare gene list
    log2fc <- de_results$avg_log2FC
    names(log2fc) <- rownames(de_results)
    log2fc <- sort(log2fc, decreasing = TRUE)
    
    # Run GSEA
    gsea_result <- GSEA(
        geneList = log2fc,
        TERM2GENE = gene_set[, c("gs_name", "gene_symbol")],
        pvalueCutoff = 0.05,
        minGSSize = 10,
        maxGSSize = 500
    )
    
    # Process results
    gsea_df <- as.data.frame(gsea_result)
    gsea_df$padj <- -log10(gsea_df$p.adjust)
    gsea_df$Direction <- ifelse(gsea_df$NES > 0, "Up", "Down")
    
    # Get top pathways
    top_pathways <- rbind(
        gsea_df %>% 
            filter(Direction == "Up") %>% 
            top_n(8, wt = abs(NES)),
        gsea_df %>% 
            filter(Direction == "Down") %>% 
            top_n(8, wt = abs(NES))
    )
    
    # Create plot
    p <- ggplot(top_pathways, aes(x = NES, y = reorder(Description, NES), fill = NES >0)) +
        geom_bar(stat = "identity") +
        geom_vline(xintercept = 0) +
        theme_minimal(base_size = 11) +
        theme(
            axis.text.y = element_text(size = 10),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, hjust = 0.5)
        ) +
        labs(title = paste("GSEA -", cell_type),
            x = "Normalized Enrichment Score (NES)",
            y = "Pathway") +
          scale_fill_manual(
              values = c(`TRUE` = "orange", `FALSE` = "navy"),
              labels = c(`TRUE` = "NES > 0", `FALSE` = "NES < 0"),
              name = NULL) + guides(fill = "none")
    
    # Return both results and plot
    return(list(
        gsea_result = gsea_result,
        plot = p,
        top_pathways = top_pathways
    ))
}