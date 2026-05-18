#!/usr/bin/awk -f
# =============================================================================
# fix_bed_chr_names.awk
# =============================================================================
# Rename chromosomes in a BED file from NCBI RefSeq accession format
# (e.g. NC_000001.11) to UCSC-style (e.g. chr1), using a two-column
# mapping file.
#
# Usage:
#   awk -f scripts/fix_bed_chr_names.awk \
#       resources/beds/chr_map.txt \
#       resources/beds/Agilent_SureSelect_V6.bed \
#     > resources/beds/Agilent_SureSelect_V6_chr_fixed.bed
#
# chr_map.txt format (tab-separated, no header):
#   NC_000001.11    chr1
#   NC_000002.12    chr2
#   ...
#   NC_000023.11    chrX
#   NC_000024.10    chrY
#   NC_012920.1     chrM
#
# Notes:
#   - Lines starting with '#' in the BED file are passed through unchanged.
#   - Chromosomes not found in the map are passed through unchanged (with a
#     warning to stderr) so the file is never silently truncated.
# =============================================================================

NR == FNR {
    # First file: build the chromosome name mapping table
    map[$1] = $2
    next
}

# Second file: the BED
/^#/ { print; next }   # pass comment/track/browser lines through unchanged

{
    if ($1 in map) {
        $1 = map[$1]
    } else {
        print "WARNING: chromosome not in map: " $1 > "/dev/stderr"
    }
    print
}
