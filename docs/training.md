# bpertussisciwgs: Training

The following steps demonstrate the basic use of in this workflow on the CDC HPC (rosalind).  


1. Clone the `Bpertussis-CIWGS` repo from GitLab into your home directory or a sub-directory:  
    ```consol
    cd ~/
    git clone https://github.com/CDCgov/bpertussis-ciwgs.git
    ```
    **OR** *If you have previously cloned the repo,* confirm that your local copy is up-to-date and on the master branch:
    ```consol
    cd ~/path-to-your-local/Bpertussis-CIWGS/
    git checkout master
    git pull
    git status
    ```
      

    Should return message:  
    
    ```consol
    On branch master
    Your branch is up-to-date with 'origin/master'.
    ```

1. Start an interactive job:
    ```consol
    qlogin
    ```

1. Load the Nextflow module (or activate your Nextflow conda environment):  
    ```consol
    module load Nextflow
    ```

1. Run the bpertussisciwgs workflow with a test samplesheet provided with this repo:  
    ```consol
    nextflow run ./main.nf \
    --input ./assets/samplesheet.test.csv \
    --outdir ./mytestrun \
    --mlst ./assets/wgMLST-v2-allele1-clean-20190103.fasta \
    --skesa true \
    -profile sge,singularity;
    ```

1. Monitor pipeline progress with the on-screen messages from Nextflow and ensure no error messages are produced.  

1. Review the summary metrics table:  
    ```consol
    cat ./mytestrun/sample_summary.tsv | column -t | less

    ```
1. Stop your interactive HPC job:
    ```consol
    exit
    ```
