
//Globals
DATA_DIR="/home/sid/thesis_SidReed/data/"
SCRIPT_DIR="/home/sid/thesis_SidReed/scripts/scrapingData"
CONDA_ENV="thesis"
PROCESSES=8

//Input Files
CRISPR_ONE_FILE="${DATA_DIR}/inputScrapingData/CRISPRone.html"
GENOME_LIST="${DATA_DIR}/inputScrapingData/prokaryotes.csv"

//Output Files/Dirs
GENOMIC_DIR="${DATA_DIR}/genomic_data"
CRISPR_ONE_DATA="${DATA_DIR}/CRISPROne_data.json"
NCBI_DATA="${DATA_DIR}/NCBI_data.json"
ALL_DATA="${DATA_DIR}/all_data.tsv"
URL_LIST_DIR="${DATA_DIR}/url_lists


scrapeCRISPROne = {
    exec """
        python "$SCRIPT_DIR"/crispr_one_scraper.py \
                "$CRISPR_ONE_FILE" \
                "$CRISPR_ONE_DATA" \
                "$PROCESSES"
        """
}
scrapeNCBI = {
    exec """
        python "$SCRIPT_DIR"/NCBINameScraper.py \
                "$CRISPR_ONE_DATA" \
                "$NCBI_DATA" \
                "$PROCESSES"
        """
}
combineData = {
    exec """
        bash "$SCRIPT_DIR"/makeAnnotationDF.py \
                "$CRISPR_ONE_DATA" \
                "$NCBI_DATA" \
                "$GENOME_LIST" \
                "$URL_LIST_DIR" \
                "$ALL_DATA"
        """
}
downloadFastas = {
    exec """
        bash "$SCRIPT_DIR"/get_gbffs.sh \
                "$URL_LIST_DIR" \
                "$GENOMIC_DIR"
        """
}

run {
    scrapeCRISPROne + scrapeNCBI + combineData + downloadFastas
}

