

## Function for reading in GFF files.
read.gff = function(fname, type = 'CDS')
{
    # seqid - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seq ID must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
    # source - name of the program that generated this feature, or the data source (database or project name)
    # type - type of feature. Must be a term or accession from the SOFA sequence ontology
    # start - Start position of the feature, with sequence numbering starting at 1.
    # end - End position of the feature, with sequence numbering starting at 1.
    # score - A floating point value.
    # strand - defined as + (forward) or - (reverse).
    # phase - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
    
    rawgff = fread(fname, sep='\t', header=FALSE, blank.lines.skip=TRUE)
    attrs = rawgff[V3 == type, V9] %>% strsplit(., ';')
    attrs = attrs %>% lapply(function(x){ y = strsplit(x, '='); sapply(y, function(z) setNames(z[2], z[1]))})
    
    loci = rawgff[V3 == type,-9]
    colnames(loci) =  c('seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase')
    
    ns = attrs %>% lapply(names) %>% unlist() %>% unique()
    attribute_dt = attrs %>% lapply(function(x) x[ns] %>% setNames(ns)) %>% lapply(as.list) %>% lapply(as.data.table) %>% do.call('rbind',.)
    
    return(cbind(loci, attribute_dt))
}