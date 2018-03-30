# RNA-Seq: Functional and strutural annotation tutorial
This repository is a usable, publicly available differential expression and functional annotation tutorial.
All steps have been provided for the UConn CBC Xanadu cluster here with appropriate headers for the Slurm scheduler that can be modified simply to run.  Commands should never be executed on the submit nodes of any HPC machine.  If working on the Xanadu cluster, you should use sbatch scriptname after modifying the script for each stage.  Basic editing of all scripts can be performed on the server with tools, such as nano, vim, or emacs.  If you are new to Linux, please use <a href="https://bioinformatics.uconn.edu/unix-basics/">this</a> handy guide for the operating system commands.  In this guide, you will be working with common bioinformatic file formats, such as <a href="https://en.wikipedia.org/wiki/FASTA_format">FASTA</a>, <a href="https://en.wikipedia.org/wiki/FASTQ_format">FASTQ</a>, <a href="https://en.wikipedia.org/wiki/SAM_(file_format)">SAM/BAM</a>, and <a href="https://en.wikipedia.org/wiki/General_feature_format">GFF3/GTF</a>. You can learn even more about each file format <a href="https://bioinformatics.uconn.edu/resources-and-events/tutorials/file-formats-tutorial/">here</a>. If you do not have a Xanadu account and are an affiliate of UConn/UCHC, please apply for one <a href="https://bioinformatics.uconn.edu/contact-us/">here</a>.
	
<div id="toc_container">
<p class="toc_title">Contents</p>
<ul class="toc_list">
<li><a href="#First_Point_Header">1 Overview<a/></li>
<li><a href="#Second_Point_Header">2 Downloading the data</a></li>
<li><a href="#Third_Point_Header">3 Identifying Regions of Genomic Repetition with RepeatModeler</a></li>
<li><a href="#Fourth_Point_Header">4 Masking Regions of Genomic Repetition with RepeatMasker</a></li>
<li><a href="#Fifth_Point_Header">5 Mapping RNA-Seq reads with HISAT2</a></li>
<li><a href="#Sixth_Point_Header">6 BRAKER2: Identifying Genes with RNA-Seq data</a></li>
<li><a href="#EnTAP">7 EnTAP: Functional Annotation for Genomes</a></li>
 <li><a href="#Integration">8 Integrating the DE Results with the Annotation Results</a></li>
<li><a href="#Citation">Citations</a></li>
</ul>
</div>

<h2 id="First_Point_Header">Overview</h2>
In this tutorial we will be performing functional and structural annotation of <i>Arabidopsis thaliana</i> using <a href="https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP135923transcriptomic"></a> data of its leaves at 4 weeks of age and software on the Xanadu cluster. The Xanadu cluster may be accessed via SSH with the following terminal command:

<pre style="color: silver; background: black;">ssh your.user.name@xanadu-submit-ext-cam.uchc.edu</pre>

It is important that after connecting via SSH the directory is set to

<pre style="color: silver; background: black;">cd /home/CAM/$your.user.name</pre> 

before proceeding. Your home directory contains 10TB of storage and will not pollute the capacities of other users on the cluster. 

The workflow may be cloned into the appropriate directory using the terminal command:
<pre style="color: silver; background: black;">$git clone https://github.com/wolf-adam-eily/annotation.git
$cd annotation
$ls  </pre>

<h2 id="Second_Point_Header">Downloading the data</h2>
Besides our RNA-Seq reads, the only other data required for our annotation is the <i>Arabidopsis thaliana</i> reference genome (this does not include the databases installed with the various software used in this tutorial). While the reference genome may be found on the NCBI website, we will be using <a href="ftp://ftp.ensemblgenomes.org/pub/plants/release-38/fasta/arabidopsis_thaliana/dna/">Ensembl</a>. After clicking on the link, you will see a variety of file types, including "rm", "sm", and "toplevel". The "rm" and "sm" are the hard-masked and soft-masked fastas, respectively. Hard-masked fastas have replaced low complexity and genomic repetition region nucleotides with 'N', preventing these regions from aligning and mapping during the various stages of analysis. Soft-masked fastas have replaced low complexity and genomic repetition region nucleotides with the lower-case correspondents, such as "aatgcgt" rather than "AATGCGT". Lastly, "toplevel" files contain the haplotype information of the sampled organisms. Because we will be masking the genome ourselves, we are only interested in the raw fastas, which we will download with the "wget" command and the "&#42;" operator:

<pre style="color: silver; background: black;">wget ftp://ftp.ensemblgenomes.org/pub/plants/release-38/fasta/arabidopsis_thaliana/dna/&#42;.dna.chromosome."[0-9]".fa.gz</pre>

The "&#42;" operator effectively instructs "wget" to accept any file with the succeeding file-name, "dna.chromosome", while "[0-9]" instructs "wget" to retrieve only those files whose file-name succeeding "&#42;.dna.chromosome" is a number 0-9 followed by ".fa.gz". Let's look at our directory now:

