#!/usr/bin/Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Description:
# 	merge output of dotter2gpseq.py and add dataset and cell_type information.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

suppressMessages(library(argparser))

# INPUT ========================================================================

# Create arguent parser
parser = arg_parser('Description: merge dotter2gpseq.py output, add dataset and
cell type information. The software looks for dotter2gpseq.py output in
subfolders of the specified input directories. These subfolders must be named as
dataset_series, with series being in XXX format with leading zeros.

Example 1: output in current directory.
./merge_data.R -m meta_HAP1.tsv meta_IMR90.tsv -i HAP1/dots_auto IMR90/dots_auto

Example 2: output to "/home/user/out" directory.
./merge_data.R -m meta_HAP1.tsv meta_IMR90.tsv -i HAP1/dots_auto IMR90/dots_auto
	-o /home/user/out',
	name = 'dotter2gpseq_merge.R'
)

# Define mandatory arguments
parser = add_argument(parser, arg = '--meta', short = '-m', nargs = Inf,
	help = paste0('List of metadata tables.',
		' Needed columns: dataset, series, cell_line, label.'))
parser = add_argument(parser, arg = '--indir', short = '-i', nargs = Inf,
	help = 'List of input folders, same order as metadata.')
parser = add_argument(parser, arg = '--outdir', short = '-o', nargs = 1,
	help = 'Output folder, created if missing. Default to current one.',
	default = ".")
parser = add_argument(parser, arg = '--aspect', short = '-a', nargs = 3,
	help = paste0('Physical size of Z, Y and X voxel sides.',
		' Default: 300.0 130.0 130.0'), default = c(300.0, 130.0, 130.0))

# Parse arguments
p = parse_args(parser)

# Attach argument values to variables
attach(p['' != names(p)])

# Additional checks
if( 1 == length(meta) ) { if( is.na(meta) ) {
	stop("at least one matadata table must be provided.")
}}
if( 1 == length(indir) ) { if( is.na(indir) ) {
	stop("at least one input directory must be provided.")
}}
if( is.na(outdir) ) {
	stop("at least one output folder must be provided.")
}
if( length(meta) != length(indir)) {
	stop("provided metadata and input directories do not match.")
}

# RUN ==========================================================================

# Map metadata to input directory
data = cbind(meta, indir)

# Iterate over cell types
l1 = lapply(1:nrow(data), FUN = function(i) {
	meta = data[i, 1]
	indir = data[i, 2]

	# Read metadata table
	md = read.delim(meta, as.is = T, header = T)

	# Check that all columns are present
	l0 = lapply(c("dataset", "series", "cell_line", "label"),
		FUN = function(x) {
		if( !x %in% colnames(md) ) {
			stop(paste0("metadata missing the '", x, "' column.",
				"\nFile at: ", meta))
		}
	})
	
	# Iterate through metadata
	l2 = lapply(1:nrow(md), FUN = function(j) {
		# Extract dataset and series information
		dataset = md$dataset[j]
		series = sprintf("%03d", md$series[j])
		flag = paste0(dataset, "_", series)
		cell_type = md$cell_line[j]
		label = md$label[j]

		# Log current status
		cat(paste0("Working on ", flag, ".\n"))

		# Identify input folder and skip if missing
		ipath = paste0(indir, "/", flag)
		if( !dir.exists(ipath) ) {
			cat(paste0("Warning: cannot find folder for ", flag,
				". Skipped.\nFolder at: ", ipath, "\n"))
			return(NULL)
		}

		# Identify input files
		flist = list.files(ipath)
		nuclei = flist[grepl("nuclei.out", flist)]
		dots = flist[grepl("wCentr.out", flist) & ! grepl("noAllele", flist)]
		
		# Read input files
		if( 0 == length(nuclei) ) {
			cat(paste0("Warning: cannot find nuclei information in ",
				flag, ".\n"))
		} else {
			nuclei = read.delim(paste0(ipath, "/", nuclei),
				as.is = T, header = T)
		}
		if( 0 == length(dots) ) {
			cat(paste0("Warning: cannot find dot information in ",
				flag, ".\n"))
		} else {
			dots = read.delim(paste0(ipath, "/", dots),
				as.is = T, header = T)
		}

		# Skip if missing file
		if( 0 == length(dots) | 0 == length(nuclei) ) {
			return(NULL)
		}

		# Add dataset, series and cell_type information
		dots$dataset = rep(dataset, nrow(dots))
		dots$label = rep(label, nrow(dots))
		dots$cell_type = rep(cell_type, nrow(dots))
		nuclei$dataset = rep(dataset, nrow(nuclei))
		nuclei$cell_type = rep(cell_type, nrow(nuclei))

		# Prepare
		aldata = as.numeric(dots$Allele)
		aldata = dots[0 < aldata & !is.na(aldata),]
		if( 0 != nrow(aldata) ) {
			aldata$universal = paste(
				aldata$File, aldata$Channel, aldata$cell_ID, sep = "~")
			alleles = do.call(rbind, by(aldata, aldata$universal,
				FUN = function(subt) {
				d = subt[1, c("File", "Channel", "cell_ID", "G1")]
				d_3d = subt[1, c("x", "y", "z")] - subt[2, c("x", "y", "z")]
				d$d_3d = sqrt(sum(((d_3d) * aspect)^2))
				d$d_lamin = abs(diff(subt$lamin_dist))
				d$d_lamin_norm = abs(diff(subt$lamin_dist_norm))
				d$d_centr = abs(diff(subt$centr_dist))
				d$d_centr_norm = abs(diff(subt$centr_dist_norm))
				d$angle = subt$angle[1]
				d$dataset = dataset
				d$label = label
				d$cell_type = cell_type
				return(d)
			}))
			rownames(alleles) = c()
		} else {
			cat(paste0("Warning: no allele couples found in ", flag, ".\n"))
			alleles = NULL
		}

		# Output
		return(list(dots = dots, nuclei = nuclei, alleles = alleles))
	})

	# Remove skipped
	l2 = l2[!is.null(l2)]

	# Output
	alleles = lapply(l2, FUN = function(x) x[[3]])
	return(list(
		dots = do.call(rbind, lapply(l2, FUN = function(x) x[[1]])),
		nuclei = do.call(rbind, lapply(l2, FUN = function(x) x[[2]])),
		alleles = do.call(rbind, alleles[!is.null(alleles)])
	))
})

# Merge
dots = do.call(rbind, lapply(l1, FUN = function(x) x[[1]]))
nuclei = do.call(rbind, lapply(l1, FUN = function(x) x[[2]]))
alleles = do.call(rbind, lapply(l1, FUN = function(x) x[[3]]))

# Write output
if( !dir.exists(outdir) ) dir.create(outdir)
write.table(dots, paste0(outdir, "/", Sys.Date(), "_dots.merged.tsv"),
	col.names = T, row.names = F, quote = F, sep = "\t")
write.table(alleles, paste0(outdir, "/", Sys.Date(), "_alleles.merged.tsv"),
	col.names = T, row.names = F, quote = F, sep = "\t")
write.table(nuclei, paste0(outdir, "/", Sys.Date(), "_nuclei.merged.tsv"),
	col.names = T, row.names = F, quote = F, sep = "\t")

# END --------------------------------------------------------------------------

################################################################################
