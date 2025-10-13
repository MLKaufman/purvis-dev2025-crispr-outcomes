### Create the target location and reference sequence

```{r target-location}
# Define the target location and reference sequence as granges
target_location <- GRanges(
  seqnames = params$target_chr,
  ranges = IRanges(start = params$target_start, end = params$target_end),
  strand = params$target_strand
)

# Extract reference sequence
reference_raw <- system(sprintf("samtools faidx %s %s:%s-%s",
                               params$ref_genome,
                               seqnames(target_location)[1], 
                               start(target_location)[1], 
                               end(target_location)[1]),
                       intern = TRUE)[[2]]

# Convert to DNAString and reverse complement
reference <- Biostrings::reverseComplement(Biostrings::DNAString(reference_raw))

# Check lengths match
target_width <- width(target_location)
ref_length <- nchar(reference_raw)

message("Target width: ", target_width)
message("Reference length: ", ref_length)
message("Reference sequence: ", reference_raw)

# If they don't match, adjust the target to match the reference
if(target_width != ref_length) {
  message("Adjusting target location to match reference sequence length")
  target_location <- GRanges(
    seqnames = params$target_chr,
    ranges = IRanges(start = params$target_start, end = params$target_start + ref_length - 1),
    strand = params$target_strand
  )
  message("New target width: ", width(target_location))
}
```

### Process guide sequences for analysis

```{r guide-processing}
# Parse guide sequences from params
guides_raw <- params$guides
message("Raw guide sequences from params:")
print(guides_raw)

# Extract guide names and sequences
guide_info <- data.frame(
  name = character(),
  sequence = character(),
  stringsAsFactors = FALSE
)

for(guide_entry in guides_raw) {
  # Split on '=' and clean up whitespace
  parts <- trimws(strsplit(guide_entry, "=")[[1]])
  if(length(parts) == 2) {
    guide_info <- rbind(guide_info, data.frame(
      name = parts[1],
      sequence = parts[2],
      stringsAsFactors = FALSE
    ))
  }
}

message("Processed guide information:")
print(guide_info)

# Find guide locations within reference sequence
reference_seq <- as.character(reference)

guide_locations <- data.frame(
  name = character(),
  sequence = character(),
  sequence_searched = character(),
  start_pos = integer(),
  end_pos = integer(),
  strand = character(),
  match_type = character(),
  found = logical(),
  stringsAsFactors = FALSE
)

for(i in 1:nrow(guide_info)) {
  guide_name <- guide_info$name[i]
  guide_seq <- guide_info$sequence[i]
  guide_seq_rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(guide_seq)))
  
  message("Searching for guide: ", guide_name)
  message("  Original sequence: ", guide_seq)
  message("  Reverse complement: ", guide_seq_rc)
  
  found_match <- FALSE
  
  # Search for original sequence on forward strand
  forward_match <- regexpr(guide_seq, reference_seq, fixed = TRUE)
  if(forward_match[1] != -1) {
    start_pos <- forward_match[1]
    end_pos <- start_pos + attr(forward_match, "match.length") - 1
    guide_locations <- rbind(guide_locations, data.frame(
      name = guide_name,
      sequence = guide_seq,
      sequence_searched = guide_seq,
      start_pos = start_pos,
      end_pos = end_pos,
      strand = "+",
      match_type = "exact_forward",
      found = TRUE,
      stringsAsFactors = FALSE
    ))
    message("  Found exact sequence on forward strand at position ", start_pos, "-", end_pos)
    found_match <- TRUE
  }
  
  # Search for reverse complement on forward strand
  if(!found_match) {
    rc_forward_match <- regexpr(guide_seq_rc, reference_seq, fixed = TRUE)
    if(rc_forward_match[1] != -1) {
      start_pos <- rc_forward_match[1]
      end_pos <- start_pos + attr(rc_forward_match, "match.length") - 1
      guide_locations <- rbind(guide_locations, data.frame(
        name = guide_name,
        sequence = guide_seq,
        sequence_searched = guide_seq_rc,
        start_pos = start_pos,
        end_pos = end_pos,
        strand = "+",
        match_type = "reverse_complement_forward",
        found = TRUE,
        stringsAsFactors = FALSE
      ))
      message("  Found reverse complement on forward strand at position ", start_pos, "-", end_pos)
      found_match <- TRUE
    }
  }
  
  # If still not found, search on reverse strand
  if(!found_match) {
    reference_seq_rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(reference_seq)))
    
    # Search for original sequence on reverse strand
    reverse_match <- regexpr(guide_seq, reference_seq_rc, fixed = TRUE)
    if(reverse_match[1] != -1) {
      # Convert position back to forward strand coordinates
      start_pos_rc <- reverse_match[1]
      end_pos_rc <- start_pos_rc + attr(reverse_match, "match.length") - 1
      # Convert to forward strand positions
      start_pos <- nchar(reference_seq) - end_pos_rc + 1
      end_pos <- nchar(reference_seq) - start_pos_rc + 1
      
      guide_locations <- rbind(guide_locations, data.frame(
        name = guide_name,
        sequence = guide_seq,
        sequence_searched = guide_seq,
        start_pos = start_pos,
        end_pos = end_pos,
        strand = "-",
        match_type = "exact_reverse",
        found = TRUE,
        stringsAsFactors = FALSE
      ))
      message("  Found exact sequence on reverse strand at position ", start_pos, "-", end_pos, " (forward coordinates)")
      found_match <- TRUE
    }
    
    # Search for reverse complement on reverse strand
    if(!found_match) {
      rc_reverse_match <- regexpr(guide_seq_rc, reference_seq_rc, fixed = TRUE)
      if(rc_reverse_match[1] != -1) {
        # Convert position back to forward strand coordinates
        start_pos_rc <- rc_reverse_match[1]
        end_pos_rc <- start_pos_rc + attr(rc_reverse_match, "match.length") - 1
        # Convert to forward strand positions
        start_pos <- nchar(reference_seq) - end_pos_rc + 1
        end_pos <- nchar(reference_seq) - start_pos_rc + 1
        
        guide_locations <- rbind(guide_locations, data.frame(
          name = guide_name,
          sequence = guide_seq,
          sequence_searched = guide_seq_rc,
          start_pos = start_pos,
          end_pos = end_pos,
          strand = "-",
          match_type = "reverse_complement_reverse",
          found = TRUE,
          stringsAsFactors = FALSE
        ))
        message("  Found reverse complement on reverse strand at position ", start_pos, "-", end_pos, " (forward coordinates)")
        found_match <- TRUE
      }
    }
  }
  
  # If no match found anywhere
  if(!found_match) {
    guide_locations <- rbind(guide_locations, data.frame(
      name = guide_name,
      sequence = guide_seq,
      sequence_searched = "none",
      start_pos = NA,
      end_pos = NA,
      strand = NA,
      match_type = "not_found",
      found = FALSE,
      stringsAsFactors = FALSE
    ))
    message("  WARNING: Neither original nor reverse complement found in reference")
  }
}

DT::datatable(guide_locations, caption = "Guide sequence locations within target region")
```

### CrisprSet object

This object is created from the sorted BAM files and the target location. It will contain the read counts and variant information for each sample at the specified target region.

```{r crispr-set}
# Grab sorted bams from processed_data using the exact same pattern
sorted_bam_files <- list.files(params$processed_data_path, pattern = "_sorted\\.bam$", full.names = TRUE)

if(length(sorted_bam_files) == 0) {
  stop("No sorted BAM files found in ", params$processed_data_path, ". Check alignment step.")
}

message("BAM files found for CrisprSet:")
for(bam_file in sorted_bam_files) {
  message("  - ", basename(bam_file))
}

bam_sample_names <- gsub("_sorted\\.bam$", "", basename(sorted_bam_files))
message("Sample names: ", paste(bam_sample_names, collapse = ", "))

# Check all files exist and are indexed
for(bam_file in sorted_bam_files) {
  if(!file.exists(bam_file)) {
    stop("BAM file not found: ", bam_file)
  }
  if(!file.exists(paste0(bam_file, ".bai"))) {
    message("Indexing missing for: ", basename(bam_file))
    indexBam(bam_file)
  }
}

crispr_set <- readsToTarget(sorted_bam_files, 
                           target = target_location, 
                           reference = reference,
                           names = bam_sample_names, 
                           target.loc = params$target_loc_forward)

DT::datatable(variantCounts(crispr_set))
DT::datatable(as.data.frame(consensusSeqs(crispr_set)))
```

### Guide-specific analysis

This section performs guide-specific analysis using the guide sequences defined in the parameters. For each guide found in the reference sequence, we create targeted CrisprSet objects and analyze cutting efficiency.

```{r guide-specific-analysis}
# Create guide-specific CrisprSet objects
guide_crispr_sets <- list()

# Filter to only guides that were found in the reference
found_guides <- guide_locations[guide_locations$found == TRUE, ]

if(nrow(found_guides) > 0) {
  for(i in 1:nrow(found_guides)) {
    guide_name <- found_guides$name[i]
    guide_start <- found_guides$start_pos[i]
    guide_end <- found_guides$end_pos[i]
    guide_strand <- found_guides$strand[i]
    match_type <- found_guides$match_type[i]
    
    # Calculate the cut site based on strand and match type
    # For Cas9: cut site is typically 3bp upstream of PAM
    # PAM is usually at 3' end of spacer sequence
    
    if(match_type %in% c("exact_forward", "reverse_complement_forward")) {
      # Guide found on forward strand (+ orientation)
      if(guide_strand == "+") {
        cut_site <- guide_end - 3  # 3bp upstream of 3' end
      } else {
        cut_site <- guide_start + 3  # Adjust for reverse orientation
      }
    } else if(match_type %in% c("exact_reverse", "reverse_complement_reverse")) {
      # Guide found on reverse strand (- orientation)
      if(guide_strand == "-") {
        cut_site <- guide_start + 3  # 3bp downstream of 5' end on reverse
      } else {
        cut_site <- guide_end - 3   # Adjust for forward orientation
      }
    } else {
      # Default calculation
      if(guide_strand == "+") {
        cut_site <- guide_end - 3
      } else {
        cut_site <- guide_start + 3
      }
    }
    
    message("Analyzing guide: ", guide_name)
    message("  Guide position: ", guide_start, "-", guide_end)
    message("  Strand: ", guide_strand)
    message("  Match type: ", match_type)
    message("  Predicted cut site: ", cut_site)
    
    # Create CrisprSet with target.loc at the predicted cut site
    tryCatch({
      guide_crispr_set <- readsToTarget(sorted_bam_files, 
                                       target = target_location, 
                                       reference = reference,
                                       names = bam_sample_names, 
                                       target.loc = cut_site)
      
      guide_crispr_sets[[guide_name]] <- guide_crispr_set
      message("  Successfully created CrisprSet for ", guide_name)
      
    }, error = function(e) {
      message("  Error creating CrisprSet for ", guide_name, ": ", e$message)
    })
  }
} else {
  message("No guides found in reference sequence - skipping guide-specific analysis")
}

message("Created ", length(guide_crispr_sets), " guide-specific CrisprSet objects")
```

### Guide-specific mutation efficiency

```{r guide-efficiency}
if(length(guide_crispr_sets) > 0) {
  guide_efficiency <- data.frame(
    guide = character(),
    efficiency = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(guide_name in names(guide_crispr_sets)) {
    tryCatch({
      eff <- mutationEfficiency(guide_crispr_sets[[guide_name]], 
                               filter.cols = "ControlDHS15primers", 
                               exclude.cols = "ControlDHS15primers")
      
      guide_efficiency <- rbind(guide_efficiency, data.frame(
        guide = guide_name,
        efficiency = eff,
        stringsAsFactors = FALSE
      ))
      
      message("Mutation efficiency for ", guide_name, ": ", round(eff, 3))
      
    }, error = function(e) {
      message("Error calculating efficiency for ", guide_name, ": ", e$message)
    })
  }
  
  if(nrow(guide_efficiency) > 0) {
    DT::datatable(guide_efficiency, caption = "Guide-specific mutation efficiency")
  }
} else {
  message("No guide-specific CrisprSet objects available for efficiency analysis")
}
```

### Guide-specific variant plots

```{r guide-variant-plots}
#| fig.width: 16
#| fig.height: 12

if(length(guide_crispr_sets) > 0) {
  # Create plots for each guide
  for(guide_name in names(guide_crispr_sets)) {
    message("Creating variant plot for guide: ", guide_name)
    
    tryCatch({
      cat("\n\n#### ", guide_name, "\n\n")
      
      # Create the plot
      plotVariants(guide_crispr_sets[[guide_name]], 
                   txdb = txdb, 
                   gene.text.size = 8,
                   main = paste("Variants for guide:", guide_name),
                   row.ht.ratio = c(1,8), 
                   col.wdth.ratio = c(4,2),
                   plotAlignments.args = list(line.weight = 0.5, ins.size = 2, 
                                             legend.symbol.size = 4),
                   plotFreqHeatmap.args = list(plot.text.size = 3, x.size = 8, 
                                               group = sample_metadata$treatment, 
                                               legend.text.size = 8, 
                                               legend.key.height = grid::unit(0.5, "lines")))
      
      # Save individual plot
      pdf(file.path(params$outs_path, paste0("variant_plot_", guide_name, ".pdf")), 
          width = 12, height = 8)
      plotVariants(guide_crispr_sets[[guide_name]], main = paste("Variants for guide:", guide_name))
      dev.off()
      
    }, error = function(e) {
      message("Error creating plot for ", guide_name, ": ", e$message)
    })
  }
} else {
  message("No guide-specific CrisprSet objects available for plotting")
}
```

### Compare guides efficiency summary

```{r guide-comparison}
if(length(guide_crispr_sets) > 1) {
  # Create a comparison summary
  comparison_data <- data.frame(
    Guide = character(),
    Total_Reads = integer(),
    Variant_Reads = integer(),
    Efficiency = numeric(),
    Top_Variant_Type = character(),
    stringsAsFactors = FALSE
  )
  
  for(guide_name in names(guide_crispr_sets)) {
    crispr_set <- guide_crispr_sets[[guide_name]]
    
    tryCatch({
      # Get basic statistics
      total_reads <- sum(totalCounts(crispr_set))
      variant_counts <- variantCounts(crispr_set)
      
      # Calculate total variant reads (exclude "no_variant" if present)
      if("no_variant" %in% rownames(variant_counts)) {
        variant_reads <- total_reads - sum(variant_counts["no_variant", ])
      } else {
        variant_reads <- sum(variant_counts) - sum(variant_counts[1, ]) # Assume first row is reference
      }
      
      efficiency <- variant_reads / total_reads * 100
      
      # Find most common variant type
      if(nrow(variant_counts) > 1) {
        variant_sums <- rowSums(variant_counts)
        top_variant <- names(variant_sums)[which.max(variant_sums[-1]) + 1] # Exclude reference
      } else {
        top_variant <- "No variants"
      }
      
      comparison_data <- rbind(comparison_data, data.frame(
        Guide = guide_name,
        Total_Reads = total_reads,
        Variant_Reads = variant_reads,
        Efficiency = round(efficiency, 2),
        Top_Variant_Type = top_variant,
        stringsAsFactors = FALSE
      ))
      
    }, error = function(e) {
      message("Error in comparison analysis for ", guide_name, ": ", e$message)
    })
  }
  
  if(nrow(comparison_data) > 0) {
    DT::datatable(comparison_data, 
                  caption = "Guide comparison summary",
                  options = list(pageLength = 10))
  }
} else {
  message("Need at least 2 guides for comparison analysis")
}
```

### Variant plots

This section creates variant plots for the CRISPR editing outcomes. It uses the `CrispRVariants` package to visualize the variants detected in the target region across samples. The main analysis uses the center of the target region, while guide-specific analyses above use predicted cut sites for each guide.

```{r variant-plots}
gtf_fname <- params$gtf_file

# Check if GTF file exists
if(file.exists(gtf_fname)) {
  message("Creating TxDb from GTF file...")
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_fname, format = "gtf")
  
  # Check if there are any genes in our target region
  target_genes <- genes(txdb, filter = list(tx_chrom = seqnames(target_location)[1]))
  overlapping_genes <- subsetByOverlaps(target_genes, target_location)
  
  if(length(overlapping_genes) == 0) {
    message("No genes found overlapping target region. Proceeding without gene annotation.")
    txdb <- NULL
  } else {
    message("Found ", length(overlapping_genes), " genes overlapping target region")
    print(overlapping_genes)
  }
} else {
  message("GTF file not found at: ", gtf_fname)
  message("Proceeding without gene annotation")
  txdb <- NULL
}
```


```{r variant-plots2}
#| fig.width: 16
#| fig.height: 16

group <- sample_metadata$treatment

plotVariants(crispr_set, txdb = txdb, gene.text.size = 8, 
    row.ht.ratio = c(1,8), col.wdth.ratio = c(4,2),
    plotAlignments.args = list(line.weight = 0.5, ins.size = 2, 
                               legend.symbol.size = 4),
    plotFreqHeatmap.args = list(plot.text.size = 3, x.size = 8, group = group, 
                                legend.text.size = 8, 
                                legend.key.height = grid::unit(0.5, "lines"))) 

pdf("variant_plot.pdf", width = 8, height = 6)
plotVariants(crispr_set)
dev.off()
```

### Mutation efficiency

```{r mutation-efficiency}
eff <- mutationEfficiency(crispr_set, filter.cols = "ControlDHS15primers", exclude.cols = "ControlDHS15primers")
eff
```

```{r crispr-set-reverse}
#| fig.width: 16
#| fig.height: 16
crispr_set_rev <-  readsToTarget(sorted_bam_files, target = target_location, reference = reference,
                            names = bam_sample_names, target.loc = params$target_loc_reverse, 
                                orientation = "opposite")
plotVariants(crispr_set_rev)
```

### Calculate frequency of each variant and total variants

```{r variant-frequency}
#| eval: false

# Calculate frequency of each variant
variant_freq <- variantFrequency(crispr_set, filter.cols = "ControlDHS15primers", exclude.cols = "ControlDHS15primers")
datatable(variant_freq)

```
