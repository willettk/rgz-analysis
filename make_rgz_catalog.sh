#!/bin/bash

# Script for running the data reduction pipeline for Radio Galaxy Zoo. 

# Written by Kyle Willett (University of Minnesota). It has been tested 
# to run on the UMN physics servers (specifically tabernacle), but may 
# require tweaking on other systems.

# Working directory. Must be changed if not running on UMN servers.

RGZ_PATH='/data/tabernacle/larry/RGZdata'

# Start MongoDB in the background

echo "Starting MongoDB"

# The first two arguments manage memory on the UMN system. If this doesn't apply, just delete them and run:
#   mongod --fork --logpath $RGZ_PATH"/mongodb/log/mongodb.log" --dbpath $RGZ_PATH"/data/db"

numactl --interleave=all mongod --fork --logpath $RGZ_PATH"/mongodb/log/mongodb.log" --dbpath $RGZ_PATH"/data/db"

# Restore the raw files from the mongoexport

echo "Restoring MongoDB files"

# To get the most recent set of RGZ classifications:
#
# - Get email with link to file on AWS
# - Click on link in email and download zipped file. Format is sanitized_radio_YYYY-MM-DD.tar.gz
#       (alternatively, files are on the zooniverse-code AWS bucket if you have access)
# - Put file in $RGZ_PATH/mongodb/exports/
# - Unzip and untar

# If there are multiple data sets in this folder, use the most recent

arr=($(find $RGZ_PATH"/mongodb/exports/sanitized_radio"* -type d))
BACKUP_PATH=${arr[${#arr[@]}-1]}
echo ${#arr[@]}" backups of catalog found"
echo "Using "$BACKUP_PATH

# Import the raw RGZ data, overwriting any sets previously in Mongo 

mongoimport --db radio --drop --collection radio_subjects $BACKUP_PATH'/radio_subjects.json'
mongoimport --db radio --drop --collection radio_classifications $BACKUP_PATH'/radio_classifications.json'
mongoimport --db radio --drop --collection radio_groups $BACKUP_PATH'/radio_groups.json'

# Activate Python 2.7 via the modules system (this is UMN specific; comment it out if you already have Python 2.7)

echo "Activating Python environments"
module load python27

# Load a Python virtual environment that has the necessary packages and their dependencies:
#   astropy, astroquery, matplotlib, mechanize, numpy, pandas, PIL, pymongo, requests, scipy
# Comment out if these packages are already installed in your local version of Python

source $RGZ_PATH'/veastropy/bin/activate'

# Run the consensus algorithm in Python

echo "Running consensus algorithm on new subjects"
python2.7 $RGZ_PATH"/rgz-analysis/python/consensus.py"

# Run the cross-matching catalog in Python

echo "Updating catalog"
python2.7 $RGZ_PATH"/rgz-analysis/python/RGZcatalog.py" --consensus $RGZ_PATH"/rgz-analysis/csv/consensus_rgz_first.csv"

# Create a flat static version of the catalog

echo "Outputting static catalog"
python2.7 $RGZ_PATH"/rgz-analysis/python/static_catalog.py"

# Export the Mongo files (raw data + consensus + matched catalog) so that a hard copy is saved to disk

echo "Dumping Mongo files"
mongodump --db radio --out $RGZ_PATH'/rgz_mongo'

# Deactivate the Python virtual environment and unload the version of Python (only necessary if activated above)

echo "Closing Python environments"
deactivate
module unload python27

# Close the background MongoDB instance

echo "Killing MongoDB"
kill $(ps aux | grep '[m]ongod --fork' | awk '{print $2}')

# Finished!
