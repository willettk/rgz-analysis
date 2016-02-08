#!/bin/bash

# Start MongoDB
echo "Starting MongoDB"
numactl --interleave=all mongod --fork --logpath /data/tabernacle/larry/RGZdata/mongodb/log/mongodb.log --dbpath /data/tabernacle/larry/RGZdata/data/db

# Restore the raw files
echo "Restoring MongoDB files"

arr=($(find mongodb/exports/radio* -type d))
backup_path=${arr[${#arr[@]}-1]}
echo ${#arr[@]}" backups of catalog found"
echo "Using "$backup_path

#mongorestore --db radio --drop --collection radio_users $backup_path'/radio_users.bson'
#mongorestore --db radio --drop --collection radio_subjects $backup_path'/radio_subjects.bson'
#mongorestore --db radio --drop --collection radio_classifications $backup_path'/radio_classifications.bson'

# Restore the latest version of the catalog
#mongorestore --db radio --drop --collection catalog '/data/tabernacle/larry/RGZdata/rgz_mongo/radio/catalog.bson'

# Activate necessary Python environments

echo "Activating Python environments"
module load python27
source /home/phys/willett/veastropy/bin/activate

# Run consensus algorithm

echo "Running consensus algorithm on new subjects"
python2.7 rgz-analysis/python/consensus.py

# Update catalog

echo "Updating catalog"
python2.7 rgz-analysis/RGZcatalogCode/RGZcatalog.py --consensus rgz-analysis/csv/consensus_rgz_first.csv

# Create static version (full flat catalog, removing peak/component data)

echo "Outputting static catalog"
python2.7 rgz-analysis/python/static_catalog.py

# Close Python environment

echo "Closing Python environments"
deactivate
module unload python27

# Export the Mongo files

echo "Dumping Mongo files"
mongodump --db radio --out rgz_mongo

# Kill the MongoDB process

echo "Killing MongoDB"
kill $(ps aux | grep '[m]ongod --fork' | awk '{print $2}')
