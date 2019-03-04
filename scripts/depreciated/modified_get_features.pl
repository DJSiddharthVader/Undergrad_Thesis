#!/usr/bin/perl
#
# Extract subsequences from a Genbank file to a fasta formatted file.
# Complement sequences are written correctly 5' to 3'.
#
# CDS or gene can be selected by -c -g
# or any features can be "grep'ed"
#

$cds=1;
$dna=0;
$rna=0;
$aa=1;
$terse=0;
$search = "";
$terseTitle_str = "";
$Title = "";

while(@ARGV) {
    $arg = shift @ARGV;
    if ($arg =~ /^\-/)  {
    FOO: {
        if($arg eq '-a') { $cds=1; $dna=0; $aa=1; last FOO; }
        if($arg eq '-b') { $cds=1; $dna=1; $aa=1; last FOO; }
        if($arg eq '-c') { $cds=1; $dna=1; $aa=0; last FOO; }
        if($arg eq '-d') { $cds=0; $dna=1; $aa=0; last FOO; }
        if($arg eq '-r') { $cds=0; $dna=1; $aa=0; $rna=1; last FOO; }
        if($arg eq '-s') { $search=shift @ARGV; last FOO; }
        if($arg eq '-t') { $terse=1; last FOO; }
        if($arg eq '-h') {
		print << "help_text";

Program extracts "feature" subsequences from genbank files.

Currently only the featureTags CDS (default), gene, repeat_region, RNA are
stored (but easily changed).
CURRENTLY ONLY WORKS FOR PROKARYOTES - NO JOINING
The options are
    -a selects CDS (protein)   [default]
    -b selects CDS (protein and then nucleotides)
    -c selects CDS nucleotides of protein coding sequences
    -d selects nucleotides
    -h print this message.
    -r selects RNA genes
    -s xxxx  selects search pattern (perl expressions accepted)
    -t print a terse title

    Syntax get_features [options] genbank.file
    output is sent to stdout.

    Note flags must be individually entered.

    Note there is sufficient variability in the format of GENBANK
    files even within the prokaryotes that some tweaking might be
    required.

    Without a search term all sequences are returned but if the search
    term is non-empty then only a match is returned.
    If preferred, the grep can be done after a run.  e.g.
    to select a particular CDS/gene pipe the output to "agrep -d "\>"
    search_term"; e.g.
        foreach f (*.gbk)
            ~/src/get_features.pl -r $f | agrep -d "\>" rrnA > `basename $f .gbk`.rrnA
        end


help_text
		exit(1); last FOO; }
	}
    }
}

$i=0;
$compliment[0] = 0;

open(FILE,"<$arg") || die "Can't open $arg \n\n";
while(<FILE>) {
    ROUND:
    if($_ =~ /^\s{5}(\S+)\s+[><]*(\d+)..[><]*(\d+)/) {
    	$featureTag[$i]=$1;  $x[$i] = $2;  $y[$i] = $3;
	$complement[$i] = 0;
	$title[$i] = "";
	$prot = "";
        $pseudo[$i] = 0;
	while(<FILE>) {
	    chomp;
	    if($_ =~ /^\s{5}(\S+)/ || $_ =~ /^\S+/) { $i++; goto ROUND; } # new tagged feature starting
#	    if($_ =~ /\/membrane protien/i) { goto ROUND; }
#	    if($_ =~ /\/phage/i) { goto ROUND; }
	    if($_ =~ /\/pseudo$/i) { $pseudo[$i]=1; }
	    if($_ =~ /\/locus_tag=\"(\S+)\"/) { $title[$i] .= " " . $1; }
	    if($_ =~ /\/gene=\"(\S+)\"/) { $title[$i] .= " " . $1; }
	    if($_ =~ /\/product=\"(.*)\"/) { $title[$i] .= " " . $1; }
	    if($_ =~ /\/protein_id=\"(.*)\"/) { $title[$i] .= " " . $1; }
	    if($_ =~ /\/protein_id=\"(.*)\"/) { $terseTitle[$i] .= $1; }
	    if($title[$i] =~ /\bphage\b/i) { goto ROUND; }
	    if($title[$i] =~ /\bprophage\b/i) { goto ROUND; }
	    if($title[$i] =~ /\bmobile.element\b/i) { goto ROUND; }
	    if($title[$i] =~ /\binsertion.sequence\b/i) { goto ROUND; }
	    if($title[$i] =~ /\bviral.element\b/i) { goto ROUND; }
	    if($title[$i] =~ /\btransposase\b/i) { goto ROUND; }
            #Added by Sid
            if($title[$i] =~ /phage(?!.*(shock|adsorption|receptor))/i) { goto ROUND; }
            if($title[$i] =~ /[^sA-Z]IS\d\d/i) { goto ROUND; }
            if($title[$i] =~ /viral.([^A-Z]|enhancin)/i) { goto ROUND; }
            if($title[$i] =~ /\bholin\b/i) { goto ROUND; }
            #/Added by Sid

	    if($_ =~ /\/translation=\"(\S+)\"/) { $prot[$i]=$1; last; }
	    if($_ =~ /\/translation=\"(\S+)/) {
		$prot[$i]=$1;
		$_ = "";  # reset
		until($_ =~ /(\S+)\"/) {
		    $_ = <FILE>;
		    chomp;
		    if($_ =~ /\s+\"/) { last; }
		    $_ =~ s/\s+//g;
		    $prot[$i] .= $_;
		    $prot[$i] =~ s/\"//g;
		}
		last;
	    }
	}
	$i++;
    }
    if($_ =~ /^\s{5}(\S+)\s+complement\([><]*(\d+)..[><]*(\d+)/) {
    	$featureTag[$i]=$1;  $x[$i] = $2;  $y[$i] = $3;
	$complement[$i] = 1;
	$title[$i] = "";
	$prot = "";
        $pseudo[$i] = 0;
	while(<FILE>) {
	    chomp;
	    if($_ =~ /^\s{5}(\S+)/ || $_ =~ /^\S+/) { $i++; goto ROUND; } # new tagged feature starting
#	    if($_ =~ /\/membrane protien/i) { goto ROUND; }
#	    if($_ =~ /\/phage/i) { goto ROUND; }
	    if($_ =~ /\/pseudo$/i) { $pseudo[$i]=1; }
	    if($_ =~ /\/locus_tag=\"(\S+)\"/) { $title[$i] .= " " . $1; }
	    if($_ =~ /\/gene=\"(\S+)\"/) { $title[$i] .= " " . $1; }
	    if($_ =~ /\/product=\"(.*)\"/) { $title[$i] .= " " . $1; }
	    if($_ =~ /\/protein_id=\"(.*)\"/) { $title[$i] .= " " . $1; }
	    if($_ =~ /\/protein_id=\"(.*)\"/) { $terseTitle[$i] .= $1; }
	    if($title[$i] =~ /\bphage\b/i) { goto ROUND; }
	    if($title[$i] =~ /\bprophage\b/i) { goto ROUND; }
	    if($title[$i] =~ /\bmobile.element\b/i) { goto ROUND; }
	    if($title[$i] =~ /\binsertion.sequence\b/i) { goto ROUND; }
	    if($title[$i] =~ /\bviral.element\b/i) { goto ROUND; }
	    if($title[$i] =~ /\btransposase\b/i) { goto ROUND; }
            #Added by Sid
            if($title[$i] =~ /\btransposase\b/i) { goto ROUND; }
            if($title[$i] =~ /phage(?!.*(shock|adsorption|receptor))/i) { goto ROUND; }
            if($title[$i] =~ /[^sA-Z]IS\d\d/i) { goto ROUND; }
            if($title[$i] =~ /viral.([^A-Z]|enhancin)/i) { goto ROUND; }
            if($title[$i] =~ /\bholin\b/i) { goto ROUND; }
            #/Added by Sid
#
    	    if($_ =~ /\/translation=\"(\S+)\"/) { $prot[$i]=$1; last; }
	    elsif($_ =~ /\/translation=\"(\S+)/) {
		$prot[$i]=$1;
		$_ = "";  # reset
		until($_ =~ /(\S+)\"/) {
		    $_ = <FILE>;
		    chomp;
		    if($_ =~ /\s+\"/) { last; }
		    $_ =~ s/\s+//g;
		    $prot[$i] .= $_;
		    $prot[$i] =~ s/\"//g;
		}
		last;
	    }
	}
	$i++;
    }
    if($_ =~ /^LOCUS/) {
        ($a, $b, $c) = split;
        $title_str = "> $b";
        $terseTitle_str = ">$b";
        $_ = <FILE> until($_ =~ /^DEFINITION\s+(.*)\s+$/);
        $def = $1;
        $_ = <FILE> until($_ =~ /^SOURCE\s+(.*)/);
	$source = $1;
        $def =~ s/[()]//;
        $source =~ s/[()]//;
	if($source =~ $def) { $def = ""; }
	else {
	    if($def =~ $source) { $source = ""; }
	}
	# $title_str = $title_str . " $source $def";
        $title_str = $title_str . " $source; ";
    }
    if($_ =~ /^ORIGIN\b/) {
        if($dna==1) {
	    $line = "";
	    $seq = "";
	    while (!(($line = <FILE>) =~ /\/\//)) {
		$line =~ s/[0-9 ]//g;
		chomp $line;
		$seq .= $line;
	    }
	    #########
	    for($j=0; $j<$i; $j++) {
		if($compliment[$j] == 0) { $title[$j] .= ";  gene = from $x[$j] to $y[$j]"; }
		else { $title[$j] .= ";  gene = from $y[$j] to $x[$j]"; }
		$sub_seq = substr($seq,$x[$j]-1,$y[$j]-$x[$j]+1);
		if($compliment[$j] == 1) {
		       $sub_seq = reverse $sub_seq;
		       $sub_seq =~ s/a/%/ig;
		       $sub_seq =~ s/t/A/ig;
		       $sub_seq =~ s/%/T/ig;
		       $sub_seq =~ s/c/%/ig;
		       $sub_seq =~ s/g/C/ig;
		       $sub_seq =~ s/%/G/ig;
		}
		$featureSeq[$j]=$sub_seq;
	    }
	}
    }
}
close(FILE);
#########
if($search ne "") { ### carry out a search
    for($j=0; $j<$i; $j++) {
        if($terse) { $Title = "$terseTitle_str:$terseTitle[$j]\n"; }
        else { $Title = "$title_str $title[$j]\n"; }
#	if(grep $search $tag[$j] || grep $search $prot[$j]
#	    || grep $search $featureSeq[$j] || grep $search $title[$j]) {
	if($featureTag[$j]=~ /$search/ || $prot[$j]=~ /$search/
	    || $featureSeq[$j]=~ /$search/ || $title[$j]=~ /$search/) {
	    if($dna) {
		print "$Title";
		for($k=0; $k<length($featureSeq[$j]); $k=$k+50) {
		    if($y[$j] > $k+50) { $line = substr($featureSeq[$j],$k,50); }
		    else { $line = substr($featureSeq[$j],$k,length($featureSeq[$j])-$k+1);}
		    print "$line\n";
		}
	    }
	    if($aa) {
		if(length($prot[$j]) <= 0) {
                    if($pseudo[$j]) {
                        print STDERR "Warning entry $j has no length -- possible pseudo. \n $title[$j] \n\n";
                    } else {
                        print STDERR "Warning entry $j has no length. \n $title[$j] \n\n";
                    }
                    next;
                }
		print "$Title";
		for($k=0; $k<length($prot[$j]); $k=$k+50) {
		    if($y[$j] > $k+50) { $line = substr($prot[$j],$k,50); }
		    else { $line = substr($prot[$j],$k,length($prot[$j])-$k+1);}
		    print "$line\n";
		}
	    }
	}
    }
    exit(0);
} else { ### no search -> print out all
    for($j=0; $j<$i; $j++) {
        if($terse) { $Title = "$terseTitle_str:$terseTitle[$j]\n"; }
        else { $Title = "$title_str $title[$j]\n"; }
	if($rna && $featureTag[$j] =~ /RNA/i) {
	    print "$Title";
	    for($k=0; $k<length($featureSeq[$j]); $k=$k+50) {
		if($y[$j] > $k+50) { $line = substr($featureSeq[$j],$k,50); }
		else { $line = substr($featureSeq[$j],$k,length($featureSeq[$j])-$k+1);}
		print "$line\n";
	    }
	}
	if($cds && $featureTag[$j] =~ /CDS/i) {
	    if($aa) {
		if(length($prot[$j]) <= 0) {
                    if($pseudo[$j]) {
                        print STDERR "Warning entry $j has no length -- possible pseudo. \n $title[$j] \n\n";
                    } else {
                        print STDERR "Warning entry $j has no length. \n $title[$j] \n\n";
                    }
                    next;
                }
		print "$Title";
		for($k=0; $k<length($prot[$j]); $k=$k+50) {
		    if($y[$j] > $k+50) { $line = substr($prot[$j],$k,50); }
		    else { $line = substr($prot[$j],$k,length($prot[$j])-$k+1);}
		    print "$line\n";
		}
	    }
	    if($dna) {
		print "$Title";
		for($k=0; $k<length($featureSeq[$j]); $k=$k+50) {
		    if($y[$j] > $k+50) { $line = substr($featureSeq[$j],$k,50); }
		    else { $line = substr($featureSeq[$j],$k,length($featureSeq[$j])-$k+1);}
		    print "$line\n";
		}
	    }
	}
	if($cds==0 && $featureTag[$j] =~ /gene/i && $rna==0) {
	    if($aa) { ### this shouldn't be possible; not all genes code proteins
		if(length($prot[$j]) <= 0) { print STDERR "Warning entry $j has no length. \n $title[$j] \n\n"; next; }
		print "$Title";
		for($k=0; $k<length($prot[$j]); $k=$k+50) {
		    if($y[$j] > $k+50) { $line = substr($prot[$j],$k,50); }
		    else { $line = substr($prot[$j],$k,length($prot[$j])-$k+1);}
		    print "$line\n";
		}
	    }
	    if($dna) {
		print "$Title";
		for($k=0; $k<length($featureSeq[$j]); $k=$k+50) {
		    if($y[$j] > $k+50) { $line = substr($featureSeq[$j],$k,50); }
		    else { $line = substr($featureSeq[$j],$k,length($featureSeq[$j])-$k+1);}
		    print "$line\n";
		}
	    }
	}
    }
}
