#!/bin/bash

# Script for running the data reduction pipeline for Radio Galaxy Zoo. 
# Written by Kyle Willett (willett@physics.umn.edu) 

# Sets where both input and output files are stored. Must be changed if not running on UMN servers.

RGZ_PATH='/data/tabernacle/larry/RGZdata'

# Start MongoDB in the background

echo "Starting MongoDB"
numactl --interleave=all mongod --fork --logpath $RGZ_PATH"/mongodb/log/mongodb.log" --dbpath $RGZ_PATH"/data2/db"

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
#BACKUP_PATH=$RGZ_PATH"/mongodb/exports/sanitized_radio_2016-08-02"
echo ${#arr[@]}" backups of catalog found"
echo "Using "$BACKUP_PATH

# Import the raw RGZ data, overwriting any sets previously in Mongo 

mongoimport --db radio --drop --collection radio_subjects $BACKUP_PATH'/radio_subjects.json'
mongoimport --db radio --drop --collection radio_classifications $BACKUP_PATH'/radio_classifications.json'
mongoimport --db radio --drop --collection radio_groups $BACKUP_PATH'/radio_groups.json'

# Activate Python 2.7 via the modules system

echo "Activating Python environments"
module load python27

# Load a Python virtual environment that has the necessary packages and their dependencies:
#   astropy, astroquery, matplotlib, mechanize, numpy, pandas, PIL, pymongo, requests, scipy

source $RGZ_PATH'/veastropy/bin/activate'

# Run the consensus algorithm in Python

echo "Running consensus algorithm on new subjects"
python2.7 $RGZ_PATH"/rgz-analysis/python/consensus.py"

# Run the cross-matching catalog in Python

echo "Updating catalog"
python2.7 $RGZ_PATH"/rgz-analysis/python/RGZcatalog.py"

# Create a flat static version of the catalog

echo "Outputting static catalog"
python2.7 $RGZ_PATH"/rgz-analysis/python/static_catalog.py"

# Export the Mongo files (raw data + consensus + matched catalog) so that we have a hard copy saved to disk

echo "Dumping Mongo files"
mongodump --db radio --out $RGZ_PATH'/rgz_mongo'

# Deactivate the Python virtual environment and unload the version of Python

echo "Closing Python environments"
deactivate
module unload python27

# Close the background MongoDB instance

echo "Killing MongoDB"
#kill $(ps aux | grep '[m]ongod --fork' | awk '{print $2}')
mongod --shutdown --dbpath $RGZ_PATH"/data2/db"

# Finished!
