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
Besides our RNA-Seq reads, the only other data required for our annotation is the <i>Arabidopsis thaliana</i> reference genome (this does not include the databases installed with the various software used in this tutorial). While the reference genome may be found on the NCBI website, we will be using <a href="ftp://ftp.ensemblgenomes.org/pub/plants/release-38/fasta/arabidopsis_thaliana/dna/">Ensembl</a>. After clicking on the link, you will see a variety of file types, including "rm", "sm", and "toplevel". The "rm" and "sm" are the hard-masked and soft-masked fastas, respectively. Hard-masked fastas have replaced low complexity and genomic repetition region nucleotides with 'N', preventing these regions from aligning and mapping during the various stages of analysis. Soft-masked fastas have replaced low complexity and genomic repetition region nucleotides with the lower-case correspondents, such as "aatgcgt" rather than "AATGCGT". Lastly, "toplevel" files contain the haplotype information of the sampled organisms. Because we will be masking the genome ourselves, we are only interested in the raw fastas, which we will download with the "wget" command and the "&#42;" operator.

<pre style="color: silver; background: black;">wget ftp://ftp.ensemblgenomes.org/pub/plants/release-38/fasta/arabidopsis_thaliana/dna/&#42;.dna.chromosome."[0-9]".fa.gz>
</pre>

The "&#42;" operator effectively instructs "wget" to accept any file with the succeeding file-name, "dna.chromosome", while "[0-9]" instructs "wget" to retrieve only those files whose file-name succeeding "&#42;.dna.chromosome" is a number 0-9 followed by ".fa.gz".

<pre style="color: silver; background: black;">ls
Arabidopsis_thaliana.TAIR10.dna.chromosome.1.fa.gz
Arabidopsis_thaliana.TAIR10.dna.chromosome.2.fa.gz
Arabidopsis_thaliana.TAIR10.dna.chromosome.3.fa.gz
Arabidopsis_thaliana.TAIR10.dna.chromosome.4.fa.gz
Arabidopsis_thaliana.TAIR10.dna.chromosome.5.fa.gz</pre>

Now, we must download our RNA-Seq using the <a href="https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc">SRA-toolkit</a>. We will be running this command as a <a href="https://bioinformatics.uconn.edu/resources-and-events/tutorials/xanadu/#Xanadu_6">Slurm scheduler</a> script. For more information, please visit the link provided. To initialize our Slurm script, we use the "nano" command as following:

<pre style="color: silver; background: black;">nano sra_download.sh
 GNU nano 2.3.1            File: sra_download.sh                               




















                                  [ New File ]
^G Get Help  ^O WriteOut  ^R Read File ^Y Prev Page ^K Cut Text  ^C Cur Pos
^X Exit      ^J Justify   ^W Where Is  ^V Next Page ^U UnCut Text^T To Spell</pre>

We may now begin writing our "&#35;SBATCH" arguments followed by our actual script commands (edit the email argument or you will not receive your process updates!):

<pre style="color: silver; background: black;">  GNU nano 2.3.1            File: sra_download.sh                     Modified  

#!/bin/bash
#SBATCH --job-name=sra_download.sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --partition=general
#SBATCH --mail-type=END
#SBATCH --mail-user=your.email@uconn.edu
#SBATCH --mem=50G
#SBATCH -o sra_download_%j.out
#SBATCH -e sra_download_%j.err
module load sratoolkit
fastq-dump SRR6852085
fastq-dump SRR6852086






^G Get Help  ^O WriteOut  ^R Read File ^Y Prev Page ^K Cut Text  ^C Cur Pos
^X Exit      ^J Justify   ^W Where Is  ^V Next Page ^U UnCut Text^T To Spell</pre>

Press CTRL+X, "y", "enter" to save your script. You may only use software previously loaded onto the Xanadu server, and cannot load your own. To view which modules you may use on Xanadu enter:
<pre style="color: silver; background: black;">module avail</pre>
Because of this, we must load the modules we wish to use in each script. In this case, we are going to utilize the sratookit with the command:
<pre style="color: silver; background: black;">module load sratoolkit</pre>
Let's submit our script and view the results:
<pre style="color: silver; background: black;">sbatch sra_download.sh</pre>

<h2 id="Third_Point_Header">Identifying Regions of Genomic Repetition with RepeatModeler</h2>
The largest proportion of genomes are low complexity regions, often consisting of <a href="https://en.wikipedia.org/wiki/Repeated_sequence_(DNA)">repetitive elements</a>. While these regions play crucial roles in safe-guarding the genome from deleterious mutations, novel protein synthesis, reproduction, and other processes, by virtue of their low-complexity they are quite common across organisms, even somewhat distantly unrelated organisms. Because of this, it can be hazardous to include these regions in alignment processes, as there is run a risk of false positives in the alignment profile. However, discarding low-complexity regions may also run the risk of removing high quality gene models alongside them. It is important to bear in mind your research goals and ambitions, choosing your modifications wisely. We will first be identifying our regions of low complexity using the <a href="http://www.repeatmasker.org/RepeatModeler/">RepeatModeler</a>. Before we identify our repeat regions, we must first compile our database using the "BuildDatabase" command of RepeatModeler. We may see our options with the following code in the bash:

<pre style="color: silver; background: black;">module load RepeatModeler
BuildDatabase
No query sequence file indicated

NAME
    BuildDatabase - Format FASTA files for use with RepeatModeler

SYNOPSIS
      BuildDatabase [-options] -name "mydb" <seqfile(s) in fasta format>
     or 
      BuildDatabase [-options] -name "mydb" 
                                  -dir <dir containing fasta files &#42;.fa, &#42;.fasta,
                                         &#42;.fast, &#42;.FA, &#42;.FASTA, &#42;.FAST, &#42;.dna,
                                         and  &#42;.DNA > 
     or
      BuildDatabase [-options] -name "mydb" 
                                  -batch <file containing a list of fasta files></pre>

Let's build our database using nano to create our Slurm script:
<pre style="color: silver; background: black;">nano repeat_modeler_db.sh
  GNU nano 2.3.1          File: repeat_modeler_db.sh                            

#!/bin/bash
#SBATCH --job-name=repeat_modeler_db
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --partition=general
#SBATCH --mail-type=END
#SBATCH --mail-user=your.email@uconn.edu
#SBATCH --mem=50G
#SBATCH -o repeat_modeler_%j.out
#SBATCH -e repeat_modeler_%j.err
module load RepeatModeler
gunzip &#42;.fa.gz
cat &#42;.fa > athaliana.txt
mv athaliana.txt athaliana.fa
BuildDatabase -name "athaliana_db" -engine ncbi athaliana.fa






^G Get Help  ^O WriteOut  ^R Read File ^Y Prev Page ^K Cut Text  ^C Cur Pos
^X Exit      ^J Justify   ^W Where Is  ^V Next Page ^U UnCut Text^T To Spell</pre>

<pre style="color: silver; background: black;">sbatch repeat_modeler_db.sh
ls
Arabidopsis_thaliana.TAIR10.dna.chromosome.1.fa  Arabidopsis_thaliana.TAIR10.dna.chromosome.4.fa  athaliana_db.nin  athaliana_db.nog          athaliana.fa          SRR6852085.fastq
Arabidopsis_thaliana.TAIR10.dna.chromosome.2.fa  Arabidopsis_thaliana.TAIR10.dna.chromosome.5.fa  athaliana_db.nnd  athaliana_db.nsq          repeat_modeler_db.sh  SRR6852086.fastq
Arabidopsis_thaliana.TAIR10.dna.chromosome.3.fa  athaliana_db.nhr                                 athaliana_db.nni  athaliana_db.translation  sra_download.sh</pre>

We are now ready to run the RepeatModeler. But first, let's have a look at our options:

<pre style="color: silver; background: black;">RepeatModeler
No database indicated

NAME
    RepeatModeler - Model repetitive DNA

SYNOPSIS
      RepeatModeler [-options] -database <XDF Database>

DESCRIPTION
    The options are:

    -h(elp)
        Detailed help

    -database
        The prefix name of a XDF formatted sequence database containing the
        genomic sequence to use when building repeat models. The database
        may be created with the WUBlast "xdformat" utility or with the
        RepeatModeler wrapper script "BuildXDFDatabase".

    -engine <abblast|wublast|ncbi>
        The name of the search engine we are using. I.e abblast/wublast or
        ncbi (rmblast version).

    -pa #
        Specify the number of shared-memory processors available to this
        program. RepeatModeler will use the processors to run BLAST searches
        in parallel. i.e on a machine with 10 cores one might use 1 core for
        the script and 9 cores for the BLAST searches by running with "-pa
        9".

    -recoverDir <Previous Output Directory>
        If a run fails in the middle of processing, it may be possible
        recover some results and continue where the previous run left off.
        Simply supply the output directory where the results of the failed
        run were saved and the program will attempt to recover and continue
        the run.</pre>
	
The options are quite straightforward. Let's go ahead and run our RepeatModeler process (the "nano" initialization of the Slurm script will be excluded from this point forward. <b>DO NOT SUBMIT THESE VIA THE BASH TERMINAL</b>. You must still initialize your Slurm script, or, alternatively, submit the batch files provided in this repository after addending your own email).

<pre style="color: silver; background: black;">
module load RepeatModeler
RepeatModeler -engine ncbi -pa 8 -database athaliana_db</pre>
